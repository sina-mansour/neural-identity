#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-14:00:00
#SBATCH --mem=10G
#SBATCH --partition=physical
#SBATCH -o slurm/volumetric_tractography/%x_%j.out

# This script uses mrtrix to create the streamlines of the white matter fibers.
# A probablistic CSD algorithm is used.
# 
# Usage: ./single_subject_tractography.sh <subject-id> <dataset> <streamlines> <nthreads>
# 
# make sure to replace <subject-id> with a subject that is already downloaded
# <dataset> should be either HCP_1200 or HCP_Retest and <streamlines> is the
# desired number of streamlines. <nthreads> is the number of threads used for
# multithreading. Here's an example:
# 
# ./single_subject_tractography.sh 877168 HCP_1200 5M 6
# 
# These softwares should be installed for the script to run properly:
# mrtrix, FSL
# 
# Detailed explanation of the method:
# - 5tt used to create a white matter mask to end the streamlines as they
#   exit white matter (subcortical regions are treated as white matter)
# - response function estimated using dhollander method
# - msmt CSD is used to get white matter fod for tractography
# - a surface mask (from .gii white matter surfaces) is generated and
#   used as tractography seeds
# - ACT is not done
# - the final streamlines should then be rejected if the endpoints are
#   far from surface vertices in order to avoid false tracks

module load FSL/5.0.11-intel-2017.u2-GCC-6.2.0-CUDA9
module load MRtrix/20190207-GCC-6.2.0
module load Python/3.6.4-intel-2017.u2
source venv/bin/activate

scratch_path="/scratch/punim0695"

# some colors for fancy logging :D
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# read the arguments
subject=$1
dataset=$2
streamlines=$3
# nthreads=$4
workingdir=`pwd`

# check if all files needed are already downloaded
# store_path="${workingdir}/DataStore/HCP_data/${dataset}/${subject}"
store_path="/data/gpfs/projects/punim0695/fMRI/DataStore/HCP_data/${dataset}/${subject}"

files=(
	# Diffusion files
	"T1w/T1w_acpc_dc_restore_1.25.nii.gz"
	"T1w/Diffusion/data.nii.gz"
	"T1w/Diffusion/bvals"
	"T1w/Diffusion/bvecs"
	"T1w/Diffusion/nodif_brain_mask.nii.gz"
	# Structural files
	"T1w/T1w_acpc_dc_restore.nii.gz"
	"T1w/T1w_acpc_dc_restore_brain.nii.gz"
	# MNI individual surfaces
	# "MNINonLinear/fsaverage_LR32k/${subject}.L.white.32k_fs_LR.surf.gii"
	# "MNINonLinear/fsaverage_LR32k/${subject}.R.white.32k_fs_LR.surf.gii"
	# native individual surfaces
	# "T1w/fsaverage_LR32k/${subject}.L.white.32k_fs_LR.surf.gii"
	# "T1w/fsaverage_LR32k/${subject}.R.white.32k_fs_LR.surf.gii"
	# Warp files
	"MNINonLinear/xfms/acpc_dc2standard.nii.gz"
	"MNINonLinear/xfms/standard2acpc_dc.nii.gz"
)

echo -e "${GREEN}[INFO]${NC} `date`: Checking if initial files are present"
for file in ${files[@]};
do
	if [ ! -f "${store_path}/${file}" ];
	then
		echo -e "${RED}[ERROR]${NC} File not found: ${store_path}/${file}"
		exit 1
	fi
done
echo -e "${GREEN}[INFO]${NC} `date`: File checks passed succesfuly"

# Create a temporary directory to store files
# tmp_path="${workingdir}/DataStore/tmp/${dataset}/${subject}"
# tmp_path="${scratch_path}/DataStore/tmp/${dataset}/${subject}"
tmp_path="/data/scratch/projects/punim0695/volumetric_tractography/${dataset}/${subject}"
out_path="/data/gpfs/projects/punim0695/fMRI/DataStore/tmp/volumetric_tractography/${dataset}/${subject}"

mkdir -p ${tmp_path}
mkdir -p ${out_path}

cd "${tmp_path}"

# First convert the initial diffusion image to .mif
echo -e "${GREEN}[INFO]${NC} `date`: Converting images to .mif"
mrconvert "${store_path}/T1w/Diffusion/data.nii.gz" "${tmp_path}/dwi.mif" \
          -fslgrad "${store_path}/T1w/Diffusion/bvecs" "${store_path}/T1w/Diffusion/bvals" \
          -datatype float32 -strides 0,0,0,1

# Five tissue type segmentation
echo -e "${GREEN}[INFO]${NC} `date`: Running five tissue type segmentation"
mrconvert "${store_path}/T1w/T1w_acpc_dc_restore_brain.nii.gz" "${tmp_path}/T1.mif"
5ttgen fsl "${tmp_path}/T1.mif" "${tmp_path}/5tt.mif" -premasked

# # Generate a meanb0.mif image
echo -e "${GREEN}[INFO]${NC} Generating mean bzero image"
dwiextract "${tmp_path}/dwi.mif" - -bzero | mrmath - mean "${tmp_path}/meanb0.mif" -axis 3

# Note: there are two alternative for response function estimation: msmt_5tt and dhollander
# This thread by dhollander explains why dhollander can be better:
# http://community.mrtrix.org/t/dwi2response-dhollander-questions/1095/5

# Estimate the response function using the dhollander method
echo -e "${GREEN}[INFO]${NC} `date`: Estimation of response function using dhollander"
dwi2response dhollander "${tmp_path}/dwi.mif" \
                        "${tmp_path}/wm.txt" \
                        "${tmp_path}/gm.txt" \
                        "${tmp_path}/csf.txt" \
                        -voxels "${tmp_path}/voxels.mif"

# # Estimate the response function using the msmt_5tt method
# echo -e "${GREEN}[INFO]${NC} `date`: Estimation of response function using msmt_5tt"
# dwi2response msmt_5tt "${tmp_path}/dwi.mif" \
#                       "${tmp_path}/5tt.mif" \
#                       "${tmp_path}/wm.txt" \
#                       "${tmp_path}/gm.txt" \
#                       "${tmp_path}/csf.txt" \
#                       -voxels "${tmp_path}/voxels.mif"

# Multi-Shell, Multi-Tissue Constrained Spherical Deconvolution
echo -e "${GREEN}[INFO]${NC} `date`: Running Multi-Shell, Multi-Tissue Constrained Spherical Deconvolution"
# First, creating a dilated brain mask (https://github.com/sina-mansour/UKB-connectomics/issues/4)
maskfilter -npass 2 "${store_path}/T1w/Diffusion/nodif_brain_mask.nii.gz" dilate "${tmp_path}/nodif_brain_mask_dilated_2.nii.gz"
# Now, perfoming CSD with the dilated mask
dwi2fod msmt_csd "${tmp_path}/dwi.mif" \
                 -mask "${tmp_path}/nodif_brain_mask_dilated_2.nii.gz" \
                 "${tmp_path}/wm.txt" "${tmp_path}/wmfod.mif" \
                 "${tmp_path}/gm.txt" "${tmp_path}/gmfod.mif" \
                 "${tmp_path}/csf.txt" "${tmp_path}/csffod.mif" \

# mtnormalise to perform multi-tissue log-domain intensity normalisation (~5sec)
# First, creating an eroded brain mask (https://github.com/sina-mansour/UKB-connectomics/issues/5)
maskfilter -npass 2 "${store_path}/T1w/Diffusion/nodif_brain_mask.nii.gz" erode "${tmp_path}/nodif_brain_mask_eroded_2.nii.gz"
mtnormalise "${tmp_path}/wmfod.mif" "${tmp_path}/wmfod_norm.mif" \
            "${tmp_path}/gmfod.mif" "${tmp_path}/gmfod_norm.mif" \
            "${tmp_path}/csffod.mif" "${tmp_path}/csffod_norm.mif" \
            -mask "${tmp_path}/nodif_brain_mask_eroded_2.nii.gz"


# # Create a mask of white matter gray matter interface
echo -e "${GREEN}[INFO]${NC} `date`: Creating gray matter white matter interface mask"
5tt2gmwmi "${tmp_path}/5tt.mif" "${tmp_path}/gmwmiSeed.mif"


# Create streamlines
# - 10 million streams are generated
# - with ACT
# - maxlength is set to 250 as the original 100 would result in streamlines
#   being at most 100 * voxelsize which will be 125mm and will result in loss
#   of long streamlines as the HCP images have high resolution
# - FOD amplitude cutoff is set to 0.06 as was recommended in ISMRM tutorial (2015)
# - ACT would appropriately crop streamline endpoints more precisely as they
#   cross white matter - cortical gray matter interface
# - using multiple threads to run faster
echo -e "${GREEN}[INFO]${NC} `date`: Creating streamlines"
tckgen -seed_gmwmi "${tmp_path}/gmwmiSeed.mif" \
	-act "${tmp_path}/5tt.mif" \
       -select "${streamlines}" \
       -maxlength 250 \
       -cutoff 0.06 \
       "${tmp_path}/wmfod_norm.mif" \
       "${out_path}/volumetric_probabilistic_tracks_${streamlines}.tck"
       # -nthreads 6 \

# Compute SIFT2 weights
tcksift2 "${out_path}/volumetric_probabilistic_tracks_${streamlines}.tck" \
	  "${tmp_path}/wmfod_norm.mif" \
	  "${out_path}/volumetric_probabilistic_sift_weights_${streamlines}.txt"

# resample endpoints
tckresample -endpoints "${out_path}/volumetric_probabilistic_tracks_${streamlines}.tck" \
 	     "${out_path}/volumetric_probabilistic_track_endpoints_${streamlines}.tck"


echo -e "${GREEN}[INFO]${NC} `date`: Removing unnecessary files."
rm -r "${tmp_path}/"

echo -e "${GREEN}[INFO]${NC} `date`: Script finished!"

deactivate


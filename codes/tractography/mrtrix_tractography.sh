#!/bin/bash

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


# some colors for fancy logging :D
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# read the arguments
subject=$1
dataset=$2
streamlines=$3
nthreads=$4
workingdir=`pwd`

# check if all files needed are already downloaded
store_path="./DataStore/HCP_data/${dataset}/${subject}"

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
	"MNINonLinear/fsaverage_LR32k/${subject}.L.white.32k_fs_LR.surf.gii"
	"MNINonLinear/fsaverage_LR32k/${subject}.R.white.32k_fs_LR.surf.gii"
	# native individual surfaces
	"T1w/fsaverage_LR32k/${subject}.L.white.32k_fs_LR.surf.gii"
	"T1w/fsaverage_LR32k/${subject}.R.white.32k_fs_LR.surf.gii"
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
tmp_path="./DataStore/tmp/tractography/${dataset}/${subject}"

mkdir -p ${tmp_path}

# First convert the initial diffusion image to .mif
echo -e "${GREEN}[INFO]${NC} `date`: Converting images to .mif"
mrconvert "${store_path}/T1w/Diffusion/data.nii.gz" "${tmp_path}/dwi.mif" \
          -fslgrad "${store_path}/T1w/Diffusion/bvecs" "${store_path}/T1w/Diffusion/bvals" \
          -datatype float32 -strides 0,0,0,1

# Five tissue type segmentation
echo -e "${GREEN}[INFO]${NC} `date`: Running five tissue type segmentation"
mrconvert "${store_path}/T1w/T1w_acpc_dc_restore_brain.nii.gz" "${tmp_path}/T1.mif"
5ttgen fsl "${tmp_path}/T1.mif" "${tmp_path}/5tt.mif" -premasked

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

# Multi-Shell, Multi-Tissue Constrained Spherical Deconvolution
echo -e "${GREEN}[INFO]${NC} `date`: Running Multi-Shell, Multi-Tissue Constrained Spherical Deconvolution"
dwi2fod msmt_csd "${tmp_path}/dwi.mif" \
                 -mask "${store_path}/T1w/Diffusion/nodif_brain_mask.nii.gz" \
                 "${tmp_path}/wm.txt" "${tmp_path}/wmfod.mif" \
                 "${tmp_path}/gm.txt" "${tmp_path}/gmfod.mif" \
                 "${tmp_path}/csf.txt" "${tmp_path}/csffod.mif" \

# Create a mask for surface seeds
echo -e "${GREEN}[INFO]${NC} `date`: Creating cortical surface mask"
python3 -c "\
import surface_mask as sm; \
sm.binary_volume_mask_from_surface_vertices(\
    '${store_path}/T1w/fsaverage_LR32k/${subject}.L.white.32k_fs_LR.surf.gii',\
    '${store_path}/T1w/fsaverage_LR32k/${subject}.R.white.32k_fs_LR.surf.gii',\
    '${store_path}/T1w/T1w_acpc_dc_restore.nii.gz',\
    '${tmp_path}/Cortical_surface_mask.nii.gz'\
)"

# create white matter + subcortical binary mask to trim streamlines with
# first extract the white matter and subcortical tissues from 5tt image
echo -e "${GREEN}[INFO]${NC} `date`: Creating white matter plus subcortical mask"
mrconvert --coord 3 2 -axes 0,1,2 "${tmp_path}/5tt.mif" "${tmp_path}/5tt-white_matter.mif"
mrconvert --coord 3 1 -axes 0,1,2 "${tmp_path}/5tt.mif" "${tmp_path}/5tt-subcortical.mif"
# add both tissues together
mrmath "${tmp_path}/5tt-white_matter.mif" \
       "${tmp_path}/5tt-subcortical.mif" \
       sum \
       "${tmp_path}/5tt-wm+sc.mif"
# binarise to create the trim mask
mrcalc "${tmp_path}/5tt-wm+sc.mif" 0 -gt 1 0 -if "${tmp_path}/5tt-trim.mif"

# Create streamlines
# - 10 million streams are generated
# - no ACT
# - maxlength is set to 250 as the original 100 would result in streamlines
#   being at most 100 * voxelsize which will be 125mm and will result in loss
#   of long streamlines as the HCP images have high resolution
# - FOD amplitude cutoff is set to 0.06 as was recommended in ISMRM tutorial (2015)
# - trim mask is used to crop streamline endpoints more precisely as they
#   cross white matter - cortical gray matter interface
# - using multiple threads to run it faster
# - the final streamlines should be filtered to remove spurious data that's not close
#   to any cortical brainordinate
echo -e "${GREEN}[INFO]${NC} `date`: Creating streamlines"
tckgen -seed_image "${tmp_path}/Cortical_surface_mask.nii.gz" \
       -mask "${tmp_path}/5tt-trim.mif" \
       -select "${streamlines}" \
       -maxlength 250 \
       -cutoff 0.06 \
       -nthreads 6 \
       "${tmp_path}/wmfod.mif" \
       "${tmp_path}/surface_tracks_${streamlines}.tck"

echo -e "${GREEN}[INFO]${NC} `date`: Removing unnecessary files."
rm "${tmp_path}/dwi.mif"

echo -e "${GREEN}[INFO]${NC} `date`: Script finished!"

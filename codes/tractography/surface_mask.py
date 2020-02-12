"""
this module contains a single function that only generates a volume mask of a surface from a cifti.
"""

import os
import numpy as np
import nibabel as nib
import nibabel.processing as nibprocessing


_main_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))


def _fwhm2sigma(fwhm):
    return fwhm / np.sqrt(8 * np.log(2))


def _gaussian_smoothing_cut_value(sigma, distance, dim):
    """
    return the value of the gaussian function (with variance sigma) at the wanted distance from the mean
    """
    return (np.exp(-(distance ** 2) / (2 * (sigma ** 2))) / ((sigma * np.sqrt(2 * np.pi)) ** dim))


def binary_volume_mask_from_surface_vertices(
        left_surface_file,
        right_surface_file,
        nifti_file,
        outfile,
        thickness=5,
        fwhm=2,
        cifti_file=os.path.join(_main_dir, 'templates/cifti/ones.dscalar.nii')):
    """
    This function generates a binary nifti volume mask (in the space of the
    nifti file) from a surface mesh file (.gii) by smoothing and thresholding
    the smoothed nifti (the threshold aims to make it nearly as thick as
    thickness value in mm)
    """
    # load the sample nifti file
    nifti = nib.load(nifti_file)

    # load the brain models from a sample cifti file
    cifti = nib.load(cifti_file)
    brain_models = [x for x in cifti.header.get_index_map(1).brain_models]

    # load surface files
    left_surface = nib.load(left_surface_file)
    right_surface = nib.load(right_surface_file)

    # get the list of surface vertices that have value in the cifti file (valid cortical surfaces)
    surfacexyz = np.vstack([
        left_surface.darrays[0].data[brain_models[0].vertex_indices],
        right_surface.darrays[0].data[brain_models[1].vertex_indices]
    ])

    # find the closest voxel to the vertex
    voxels_ijk = np.round(nib.affines.apply_affine(np.linalg.inv(nifti.affine), surfacexyz))

    # create an integer weighted mask from the closest voxel list
    data = nifti.get_data() * 0
    for (i, j, k) in voxels_ijk.astype(int):
        try:
            data[i, j, k] += 1
        except IndexError:
            print((i, j, k))

    # save the projections of vertices to the volume as a nifti image
    vertex_projection_mask = nib.nifti1.Nifti1Image(data, nifti.affine, nifti.header, nifti.extra)

    # smooth the image in volume space
    smoothed_projection_mask = nibprocessing.smooth_image(vertex_projection_mask, fwhm)

    # threshold the smoothed data and binarise
    threshold = _gaussian_smoothing_cut_value(_fwhm2sigma(fwhm), thickness / 2, 2)
    data = smoothed_projection_mask.get_data()
    # data = (data>threshold)*data
    data = (data > threshold) * 1.

    # save the final thresholded mask
    binarised_thresholded_projection_mask = nib.nifti1.Nifti1Image(data, nifti.affine, nifti.header, nifti.extra)

    nib.save(binarised_thresholded_projection_mask, outfile)

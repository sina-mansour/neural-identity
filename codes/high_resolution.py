"""
This module provides the methods used in our paper to conduct high-resolution
connectomics. The computations are based on creation of sparse representations
for the high resolution connectomes.

The module contains two main functionalities:

    - Codes for creating high-resolution sparse functional connectivity
      from raw fMRI time-series

    - Codes for creating high-resolution sparse structural connectivity
      from streamlines generated by tractography of diffusion MRI

Python implementation of high-resolution connectomic analyses
Author: Sina Mansour L.
Contact: sina.mansour.lakouraj@gmail.com
"""

import os
import numpy as np
import scipy.sparse as sparse
import scipy.spatial as spatial
from scipy.interpolate import RegularGridInterpolator
import sklearn.preprocessing
import nibabel as nib
import gdist


_main_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def _join_path(*args):
    return os.path.join(*args)


def _write_sparse(sp_obj, file_path):
    sparse.save_npz(file_path, sp_obj)


def _load_sparse(file_path):
    return sparse.load_npz(file_path)


def _print_log(message, mode='info'):
    if mode == 'info':
        print ('{}: \033[0;32m[INFO]\033[0m {}'.format(time_str(), message))
    if mode == 'err':
        print ('{}: \033[0;31m[ERROR]\033[0m {}'.format(time_str(), message))
        quit()
    # if mode == 'progress':
    #     print (' ' * 79, end="\r")
    #     print ('{}: \033[0;32m[PROGRESS]\033[0m {}'.format(time_str(), message), end="\r")
    # if mode == 'end_progress':
    #     print (' ' * 79, end="\r")
    #     print ('{}: \033[0;32m[PROGRESS RESULT]\033[0m {}'.format(time_str(), message))
    sys.stdout.flush()


def _handle_process_with_que(que, func, args, kwds):
    que.put(func(*args, **kwds))


def _run_in_separate_process(func, *args, **kwds):
    que = mp.Queue()
    p = mp.Process(target=_handle_process_with_que, args=(que, func, args, kwds))
    p.start()
    # join removed as it caused the parent to sleep for no reason
    # p.join()
    return que.get()


def _normalize_time_series(ts):
    return (ts - np.mean(ts, axis=0)) / np.std(ts, axis=0)


def _fwhm2sigma(fwhm):
    return fwhm / np.sqrt(8 * np.log(2))


def _max_smoothing_distance(sigma, epsilon, dim):
    """
    return the distance of the smoothing kernel that will miss a epsilon proportion of the
    smoothed signal energy
    """
    return sigma * (-stats.norm.ppf((1 - (1 - epsilon) ** (1 / dim)) / 2))


def _diagonal_stack_sparse_matrices(m1, m2):
    """
    Inputs are expected to be CSR matrices

    this is what the output looks like:

    | M1  0 |
    | 0  M2 |

    """
    return sparse.vstack((
        sparse.hstack((
            m1,
            sparse.csr_matrix((m1.shape[0], m2.shape[1]), dtype=m1.dtype)
        )).tocsr(),
        sparse.hstack((
            sparse.csr_matrix((m2.shape[0], m1.shape[1]), dtype=m1.dtype),
            m2
        )).tocsr()
    ))


# High-resolution functional connectivity:

def compute_sparse_functional_connectivity_from_timeseries(
        time_series,
        sparse_mask,
        chunk_size=1000):
    """
    This function computes a sparse vertex based connectome from the functional time-series.
    The sparse connectome gets created over the sparse mask provided.
    Note: This function can be used to compute sparse functional connectomes on both volumes
    (voxels), and surfaces (vertices). Howecer the sparse mask provided must be of the same
    dimensions as the input time-series.

    Args:
        time_series (ndarray):
            Nvert x Ntime array storing the timeseries of all vertices in the high-resolution.
        sparse_mask (sparse matrix of type numpy bool):
            Nvert x Nvert sparse mask to sparsify the dense connectome
        chunk_size (int):
            The chunk of the dense connectome rows to be loaded to sparsify at a time. This
            chunk is used to limit the memory footprint of loading the complete dense connectome
            at once. Smaller chunck size reduces the memory footprint at the cost of increase of
            time required to generate the sparsified connectome. The chunk size will not impact
            the value of the computed sparse connectivities.


    Returns:
        sparse_functional_connectome:
            Nvert x Nvert sparse matrix storing the functional connectivity values (pearson
            correlation) of the sparsified connectome.
    """
    nts = _normalize_time_series(time_series)

    sparse_chunck_list = []

    for i in range(0, int(nts.shape[1]), chunk_size):
        # portion of connectivity
        pcon = (np.matmul(nts[:, i:i + chunk_size].T, nts) / nts.shape[0])

        # sparsified connectivity portion
        spcon = sparse_mask[i:i + chunk_size, :].multiply(pcon)

        sparse_chunck_list.append(spcon)

    scon = sparse.vstack(sparse_chunck_list)

    return scon


def compute_sparse_functional_connectivity_from_dtseries(
        dtseries_file,
        sparse_mask=_load_sparse(_join_path(_main_dir, 'data/templates/sparse_mask/functional_sparse_mask_1%_density.npz'))):
    """
    Compute the high-resolution sparse functional connectivity of a dtseries file using a
    sparsification mask.
    Note: The default mask is the 1% group average mask used in our paper which can be used
          on the cifti dtseries functional MRI data provided in the human connectome project.
    Note: If you wish to extract high-resolution connecomes from other sources, it is
          recommended to use the <compute_sparse_functional_connectivity_from_timeseries> method
          instead and provide your own sparse mask.

    Args:
        dtseries_file (string):
            The path to dtseries file containing Nvert x Ntime array storing the functional
            activity timeseries of all vertices in the high-resolution.
        sparse_mask (sparse matrix of type numpy bool):
            Nvert x Nvert sparse mask to sparsify the dense connectome (default=1% group average
            computed from HCP S1200 group average)


    Returns:
        sparse_functional_connectome:
            Nvert x Nvert sparse matrix storing the functional connectivity values (pearson
            correlation) of the sparsified connectome.
    """
    time_series = nib.load(dtseries_file)

    return compute_sparse_functional_connectivity_from_timeseries(time_series, sparse_mask)


# High-resolution structural connectivity

def _get_xyz_hem_surface(hem_surface_file,
                         brain_model_index,
                         cifti_file):
    """
    returns the xyz mm coordinates of all brainordinates in that hemisphere's surface mesh
    """
    hem_surface = nib.load(hem_surface_file)
    if cifti_file is not None:
        img = nib.load(cifti_file)

        brain_models = [x for x in img.header.get_index_map(1).brain_models]

        return hem_surface.darrays[0].data[brain_models[brain_model_index].vertex_indices]
    else:
        return hem_surface.darrays[0].data


def _get_xyz_surface(left_surface_file, right_surface_file, cifti_file):
    """
    returns the xyz mm coordinates of all brainordinates in the surface mesh
    """
    # left cortex
    leftxyz = _get_xyz_hem_surface(left_surface_file, 0, cifti_file=cifti_file)

    # right cortex
    rightxyz = _get_xyz_hem_surface(right_surface_file, 1, cifti_file=cifti_file)

    return np.vstack([leftxyz, rightxyz])


def _apply_warp_to_points_mm_native_to_mm_MNI(native_mms, warpfile):
    """
    This function is used to warp a list of points from native mm space to MNI space
    using a warpfield file. make sure to put the reverse warp file (standard2acpc)
    Note: points are given as a m*3 array.
    """
    warp = nib.load(warpfile)

    x = np.linspace(0, warp.shape[0] - 1, warp.shape[0])
    y = np.linspace(0, warp.shape[1] - 1, warp.shape[1])
    z = np.linspace(0, warp.shape[2] - 1, warp.shape[2])

    xinterpolate = RegularGridInterpolator((x, y, z), warp.get_data()[:, :, :, 0])
    yinterpolate = RegularGridInterpolator((x, y, z), warp.get_data()[:, :, :, 1])
    zinterpolate = RegularGridInterpolator((x, y, z), warp.get_data()[:, :, :, 2])

    native_voxs = nib.affines.apply_affine(np.linalg.inv(warp.affine), native_mms)

    dx_mm, dy_mm, dz_mm = (-xinterpolate(native_voxs), yinterpolate(native_voxs), zinterpolate(native_voxs))

    return native_mms + np.array([dx_mm, dy_mm, dz_mm]).T


def _get_streamline_warped_endpoints_and_dists(track_file,
                                               rewarp_file,
                                               left_surface_file,
                                               right_surface_file,
                                               sample_cifti_file):
    """
    return the warped streamline endpoint distances from closest vertex on cortical
    surface mesh and the closest vertex index.
    """
    # load the track file streamlines
    _print_log('loading track file.')
    tracks = nib.streamlines.load(track_file)
    _print_log('track file loaded: {}'.format(track_file))

    # extract streamline endpoints
    starts = np.array([stream[0] for stream in tracks.streamlines])
    ends = np.array([stream[-1] for stream in tracks.streamlines])
    _print_log('endpoints extracted: #{}'.format(len(starts)))

    if rewarp_file is not None:
        # calculate endpoint coordinates in the MNI space
        warped_starts = _apply_warp_to_points_mm_native_to_mm_MNI(starts, rewarp_file)
        warped_ends = _apply_warp_to_points_mm_native_to_mm_MNI(ends, rewarp_file)
        _print_log('endpoints warped: #{}'.format(len(starts)))
    else:
        warped_starts = starts
        warped_ends = ends

    # extract cortical surface coordinates
    surface_xyz = _get_xyz_surface(left_surface_file, right_surface_file, sample_cifti_file)

    # store the coordinates in a kd-tree data structure to locate closest point faster
    kdtree = spatial.cKDTree(surface_xyz)

    # locate closest surface points to every endpoint
    start_dists, start_indices = kdtree.query(warped_starts)
    end_dists, end_indices = kdtree.query(warped_ends)
    _print_log('closest brainordinates located')

    return (start_dists, start_indices, end_dists, end_indices, len(surface_xyz))


def _get_streamline_incidence(start_dists,
                              start_indices,
                              end_dists,
                              end_indices,
                              node_count,
                              threshold=2):
    """
    returns a couple of half incidence matrices in a sparse format after
    filtering the streamlines that are far (>2mm) from their closest vertex.
    """
    # mask points that are further than the threshold from all surface coordinates
    outlier_mask = (start_dists > threshold) | (end_dists > threshold)
    _print_log('outliers located: #{} outliers ({}%, with threshold {}mm)'.format(
        sum(outlier_mask),
        (100 * sum(outlier_mask)) / len(outlier_mask),
        threshold,
    ))

    # create a sparse incidence matrix
    _print_log('creating sparse incidence matrix')
    start_dict = {}
    end_dict = {}
    indices = (i for i in range(len(outlier_mask)) if not outlier_mask[i])
    for l, i in enumerate(indices):
        start_dict[(start_indices[i], l)] = start_dict.get((start_indices[i], l), 0) + 1
        end_dict[(end_indices[i], l)] = end_dict.get((end_indices[i], l), 0) + 1

    start_inc_mat = sparse.dok_matrix(
        (
            node_count,
            (len(outlier_mask) - outlier_mask.sum())
        ),
        dtype=np.float32
    )

    for key in start_dict:
        start_inc_mat[key] = start_dict[key]

    end_inc_mat = sparse.dok_matrix(
        (
            node_count,
            (len(outlier_mask) - outlier_mask.sum())
        ),
        dtype=np.float32
    )

    for key in end_dict:
        end_inc_mat[key] = end_dict[key]

    _print_log('sparse matrix generated')

    return (start_inc_mat.tocsr(), end_inc_mat.tocsr())


def _local_geodesic_distances(max_distance, vertices, triangles):
    # distances = gdist.local_gdist_matrix(vertices.astype(np.float64), triangles.astype(np.int32), max_distance)
    distances = utils._run_in_separate_process(
        gdist.local_gdist_matrix,
        vertices.astype(np.float64),
        triangles.astype(np.int32),
        max_distance=max_distance,
    )

    # make sure maximum distance is applied
    distances[distances > max_distance] = 0
    distances = distances.minimum(distances.T)
    distances.eliminate_zeros()
    distances = distances.tolil()
    distances.setdiag(0)
    distances = distances.tocsr()
    return distances


def _local_geodesic_distances_on_surface(surface, max_distance):
    vertices = surface.darrays[0].data
    triangles = surface.darrays[1].data
    retval = _local_geodesic_distances(max_distance, vertices, triangles)
    return retval


def _trim_and_stack_local_distances(left_local_distances,
                                    right_local_distances,
                                    sample_cifti_file):
    # load a sample file to read the mapping from
    cifti = nib.load(cifti_file)

    # load the brain models from the file (first two models are the left and right cortex)
    brain_models = [x for x in cifti.header.get_index_map(1).brain_models]

    # trim left surface to cortex
    left_cortex_model = brain_models[0]
    left_cortex_indices = left_cortex_model.vertex_indices[:]
    left_cortex_local_distance = left_local_distances[left_cortex_indices, :][:, left_cortex_indices]

    # trim right surface to cortex
    right_cortex_model = brain_models[1]
    right_cortex_indices = right_cortex_model.vertex_indices[:]
    right_cortex_local_distance = right_local_distances[right_cortex_indices, :][:, right_cortex_indices]

    # concatenate local distances with diagonal stacking
    return _diagonal_stack_sparse_matrices(left_cortex_local_distance, right_cortex_local_distance)


def _get_cortical_local_distances(left_surface, right_surface, max_distance, sample_cifti_file):
    """
    This function computes the local distances on the cortical surface and returns a sparse matrix
    with dimensions equal to cortical brainordinates in the cifti file.
    """
    left_local_distances = _local_geodesic_distances_on_surface(left_surface, max_distance)
    right_local_distances = _local_geodesic_distances_on_surface(right_surface, max_distance)
    return _trim_and_stack_local_distances(left_local_distances, right_local_distances, sample_cifti_file)


def _local_distances_to_smoothing_coefficients(local_distance, sigma):
    """
    Takes a sparse local distance symmetric matrix (CSR) as input,
    Generates an assymetric coefficient sparse matrix where each
    row i, has the coefficient for smoothing a signal from node i,
    therefore, each row sum is unit (1). sigma comes from the smoothing
    variance.
    """
    # add zeros to the diagonal
    local_distance_with_diagonals = local_distance + (sparse.eye(local_distance.shape[0], dtype=local_distance.dtype).tocsr() * 0)

    # apply gaussian transform
    gaussian = -(local_distance_with_diagonals.power(2) / (2 * (sigma ** 2)))
    np.exp(gaussian.data, out=gaussian.data)

    # normalize rows of matrix
    return sklearn.preprocessing.normalize(gaussian, norm='l1')


def _get_smoothed_adjacency_from_unsmoothed_incidence(start_inc_mat, end_inc_mat, local_smoothing_coefficients):
    """
    Return a smoothed sparse adjacency matrix from the two halfs of incidence matrix.
    The smoothing is done at the network level, that is the incidence matrices are
    smoothed before creation of the adjacency.
    """
    smoothed_start_inc_mat = start_inc_mat.T.dot(local_smoothing_coefficients).T
    smoothed_end_inc_mat = end_inc_mat.T.dot(local_smoothing_coefficients).T
    A = smoothed_start_inc_mat.dot(smoothed_end_inc_mat.T)
    return A + A.T


def compute_smoothed_structural_connectivity_over_cifti_from_streamlines(
        track_file,
        left_surface_file,
        right_surface_file,
        rewarp_file,
        sample_cifti_file=_join_path(_main_dir, 'data/templates/cifti/ones.dscalar.nii'),
        sigma=_fwhm2sigma(2),
        epsilon=0.01):
    """
    This function combines previous functions to get a smoothed structural connectivity using the
    provided streamlines resulted from tractography.
    Note: As it is recommended to run tractography in the native space, a warp file is used to
          map the streamline endpoints from native space to standard space. (the space which the
          individual left and right surfaces reside in.) You may use the native surface spaces
          and not provide a warp file by rewarp_file=None instead.
    Note: Two-dimensional smoothing is done on the surface mesh provided.
    Note: A sample cifti file is used to create a mask excluding medial wall. You may use
          sample_cifti_file=None to disable surface mesh masking on medial wall. (If you chose to
          disable surface masking, make sure to mask the output high-resolution connectome to
          disregard the anatomically unfeasible connectome nodes on the medial wall.)

    Args:
        track_file (string):
            The path to a .tck file containing the streamlines generated by tractography (The
            function is tested by tck files generated by mrtrix3).
        left_surface_file (string):
            The path to a .surf.gii surface file containing the individual left white-matter surface mesh
            (tested with the 32k resolution mesh).
        right_surface_file (string):
            The path to a .surf.gii surface file containing the individual right white-matter surface mesh
            (tested with the 32k resolution mesh).
        rewarp_file (string):
            The path to a xfms file containing the warp information used to warp from standard to
            native space (The function is tested by the HCP files in
            `MNINonLinear/xfms/standard2acpc_dc.nii.gz`).
            Note: If the surface files are provided in the native space (no warp needed), just provide
                  rewarp_file=None
        sample_cifti_file (string):
            The path to a sample scalar cifti (.dscalar.nii) file to mask the medial wall from the surface
            meshes as the endpoints of the streamlines should not be on the medial wall.
            (tested with the 32k resolution mesh).
            Note: You may provide sample_cifti_file=None to disable the masking step. However, this will
                  result in a connectome with nodes on the medial wall which you must exclude later.
        sigma (float):
            The sigma value used for the smoothing (default=equivalent of 2mm FWHM).
        epsilon (float):
            The proportion of signal lost by the distance limit on the gaussian smoothing filter
            (default= 0.01 = 1%)


    Returns:
        sparse_structural_connectome:
            Nvert x Nvert sparse matrix storing the structural connectivity values of the smoothed
            connectome. (left vertices first, then right)
    """
    return _get_smoothed_adjacency_from_unsmoothed_incidence(
        *_get_streamline_incidence(
            *_get_streamline_warped_endpoints_and_dists(
                track_file,
                rewarp_file,
                left_surface_file,
                right_surface_file,
                sample_cifti_file,
            )
        ),
        _local_distances_to_smoothing_coefficients(
            _get_cortical_local_distances(
                nib.load(left_surface_file),
                nib.load(right_surface_file),
                _max_smoothing_distance(
                    sigma,
                    epsilon,
                    2  # 2 dimensional smoothing on the surface
                ),
                sample_cifti_file
            ),
            sigma
        )
    )

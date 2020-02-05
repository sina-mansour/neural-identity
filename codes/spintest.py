# Python implementation of spin-test non-parametric spatial correspondence statistic
# Author of python implementation: Sina Mansour L.
# Contact: sina.mansour.lakouraj@gmail.com

# Re-implementation of the codes from Alexander-Bloch et. al. 2018
# Original source code is available in Matlab
# https://github.com/spin-test/spin-test/

# Extra implemented features:
# Choice of two-tailed/one-tailed test (default: one tailed)

import numpy as np
import scipy.spatial as spatial
import scipy.stats as stats
import nibabel as nib


def non_parametric_spin_test(
        data_l1,
        data_r1,
        data_l2,
        data_r2,
        vertices_l,
        vertices_r,
        permutations=1000,
        random_seed=None,
        method='pearson',
        two_tailed=False):
    """
    Calculate the p-value of correlation between two surface maps based on
    the null distribution of spins of map 2

    Args:
        data_l1 (list):
            surface map 1 for the left hemisphere data (length N).
        data_r1 (list):
            surface map 1 for the right hemisphere data (length N).
        data_l2 (list):
            surface map 2 for the left hemisphere data (length N).
        data_r2 (list):
            surface map 2 for the right hemisphere data (length N).
        vertices_l (array):
            x,y,z location of each suraface vertex from left surface
            sphere mesh (shape (N,3)).
        vertices_r (array):
            x,y,z location of each suraface vertex from right surface
            sphere mesh (shape (N,3)).
        permutations (int):
            The number of permutations, default=1000.
        random_seed (int):
            Seed used for replicability in the random number generation.
            defualt=None,
        method (string):
            Choice of correspondence methods: from "pearson" and "spearman".
        two_tailed (bool):
            a flag to indicate whether the non-parametric test should be
            one- or two-tailed.
            default=False

    Returns:
        dict:
            {
                'corellation r': The Pearson's corellation of the surface maps,
                'parametric p-val': parametric p-value of the correlation in a t-test,
                'non-parametric p-val': non-parametric p-value from the spin test (main result),
                'null distribution of r values': distribution of null,
            }
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    # build kd-tree from surfaces
    ckdtree_l = spatial.cKDTree(vertices_l)
    ckdtree_r = spatial.cKDTree(vertices_r)

    # keep data1, rotate data2, and permute.
    data1 = np.concatenate([data_l1, data_r1])
    data2 = np.concatenate([data_l2, data_r2])
    mask = ~np.isnan(data1) & ~np.isnan(data2)

    # get the unrotated actual Pearson's correlation
    if method == 'spearman':
        real_rho, parametric_p = stats.spearmanr(data1[mask], data2[mask])
    else:
        real_rho, parametric_p = stats.pearsonr(data1[mask], data2[mask])

    null_rhos = []
    for p in range(permutations):
        # generate random spin
        A = np.random.normal(0, 1, (3, 3))
        TL, temp = np.linalg.qr(A)
        TL = TL.dot(np.diag(np.sign(np.diag(temp))))
        if (np.linalg.det(TL) < 0):
            TL[:, 0] = -TL[:, 0]
        I1 = np.diag([-1, 1, 1])
        TR = I1.dot(TL.dot(I1))

        # get rotated vertices
        rotated_vertices_l = vertices_l.dot(TL)
        rotated_vertices_r = vertices_r.dot(TR)

        # find nearest neighbors
        Il = ckdtree_l.query(rotated_vertices_l)[1]
        Ir = ckdtree_r.query(rotated_vertices_r)[1]

        # get rotated data
        rotated_data_l2 = [data_l2[Il[x]] for x in range(len(data_l2))]
        rotated_data_r2 = [data_r2[Ir[x]] for x in range(len(data_r2))]
        rotated_data2 = np.concatenate([rotated_data_l2, rotated_data_r2])

        # get the rotated rho for null
        mask = ~np.isnan(data1) & ~np.isnan(rotated_data2)
        if method == 'spearman':
            null_rhos.append(stats.spearmanr(data1[mask], rotated_data2[mask])[0])
        else:
            null_rhos.append(stats.pearsonr(data1[mask], rotated_data2[mask])[0])

    # compute non-parametric p-value
    if two_tailed:
        np_pval = len([x for x in null_rhos if np.abs(x) > np.abs(real_rho)]) / len(null_rhos)

    else:
        if real_rho >= 0:
            np_pval = len([x for x in null_rhos if x > real_rho]) / len(null_rhos)

        else:
            np_pval = len([x for x in null_rhos if x < real_rho]) / len(null_rhos)

    # return the statistics
    return {
        'corellation r': real_rho,
        'parametric p-val': parametric_p,
        'non-parametric p-val': np_pval,
        'null distribution of r values': null_rhos,
    }


def get_left_right_surface_data_from_cifti(cifti_file, left_surface_file, right_surface_file):
    """Get separated data from left and right hemispheres"""
    cifti = nib.load(cifti_file)
    raw_data = cifti.get_data()[0]
    brain_models = [x for x in cifti.header.get_index_map(1).brain_models]

    left_surface = nib.load(left_surface_file)
    right_surface = nib.load(right_surface_file)
    vertices_l = left_surface.darrays[0].data
    vertices_r = right_surface.darrays[0].data

    index_to_data_loc_l = {x: (i + brain_models[0].index_offset) for i, x in enumerate(brain_models[0].vertex_indices)}
    index_to_data_loc_r = {x: (i + brain_models[1].index_offset) for i, x in enumerate(brain_models[1].vertex_indices)}

    data_l_idx = [index_to_data_loc_l.get(x, np.nan) for x in range(vertices_l.shape[0])]
    data_r_idx = [index_to_data_loc_r.get(x, np.nan) for x in range(vertices_r.shape[0])]

    data_l = [x if np.isnan(x) else raw_data[x] for x in data_l_idx]
    data_r = [x if np.isnan(x) else raw_data[x] for x in data_r_idx]

    return (data_l, data_r)


def cifti_spin_test(
        cifti_file_1,
        cifti_file_2,
        left_sphere_file,
        right_sphere_file,
        permutations=1000,
        random_seed=None,
        method='pearson',
        two_tailed=False):
    """
    Read two separate cifti files (.dscalar.nii) and compute the p-value of
    the spin-test non-parametric statistic correspondence of the
    surface maps.

    Args:
        cifti_file_1 (string):
            first cifti file location.
        cifti_file_2 (string):
            second cifti file location.
        left_sphere_file (string):
            left surface sphere file location,
            default=32k mesh sphere from templates.
        right_sphere_file (string):
            right surface sphere file location,
            default=32k mesh sphere from templates.
        permutations (int):
            The number of permutations, default=1000.
        random_seed (int):
            Seed used for replicability in the random number generation,
            defualt=None.
        method (string):
            Choice of correspondence methods: from "pearson" and "spearman".
        two_tailed (bool):
            a flag to indicate whether the non-parametric test should be
            one- or two-tailed.
            default=False

    Returns:
        dict:
            {
                'corellation r': The Pearson's corellation of the surface maps,
                'parametric p-val': parametric p-value of the correlation in a t-test,
                'non-parametric p-val': non-parametric p-value from the spin test (main result),
                'null distribution of r values': distribution of null,
            }
    """
    data_l1, data_r1 = get_left_right_surface_data_from_cifti(
        cifti_file_1, left_sphere_file, right_sphere_file
    )
    data_l2, data_r2 = get_left_right_surface_data_from_cifti(
        cifti_file_2, left_sphere_file, right_sphere_file
    )

    left_sphere = nib.load(left_sphere_file)
    right_sphere = nib.load(right_sphere_file)
    vertices_l = left_sphere.darrays[0].data
    vertices_r = right_sphere.darrays[0].data

    return non_parametric_spin_test(
        data_l1,
        data_r1,
        data_l2,
        data_r2,
        vertices_l,
        vertices_r,
        permutations=permutations,
        random_seed=random_seed,
        method=method,
        two_tailed=two_tailed)

# Neural correlates of identity and behavior

---

Here, you may find the resources (code and data) used for our article on neural correlates of identity and behavior. All resources are provided as complementary to the following article:

**reference to the article will be added here ...**

---

## TLDR

This repository contains:

### Codes:

- Codes used for generation of high-resolution sparse functional connectomes:[`./codes/high_resolution.py`](https://github.com/sina-mansour/neural-identity/blob/master/codes/high_resolution.py)
- Codes used for generation of high-resolution network-smoothed structural connectomes: [`./codes/high_resolution.py`](https://github.com/sina-mansour/neural-identity/blob/master/codes/high_resolution.py)
- Implementation of **spintest** non-parametric spatial correspondence test in python: [`./codes/spintest.py`](https://github.com/sina-mansour/neural-identity/blob/master/codes/spintest.py)
- Implementation of variance component model (**VCM**) in python: [`./codes/vcm.py`](https://github.com/sina-mansour/neural-identity/blob/master/codes/vcm.py)
- Sample script we used to run mrtrix tractography: [`./codes/tractography/mrtrix_tractography.sh`](https://github.com/sina-mansour/neural-identity/blob/master/codes/tractography/mrtrix_tractography.sh)

### Data:

- Surface maps of the expressed neural uniqueness from different modalities: [`./data/uniqueness/`](https://github.com/sina-mansour/neural-identity/blob/master/data/uniqueness/)
- Individual high-resolution connectomes can be downloaded from [our public s3 bucket](https://swift.rc.nectar.org.au/v1/AUTH_ee5989f4c9184ea29012bb124cd3dff0/connectome_storage/index.html)
- Sparse mask used to sparsify functional connectomes: [`./data/sparse_mask/`](https://github.com/sina-mansour/neural-identity/blob/master/data/sparse_mask/)
- Supplementary information of our indexing and labels: [`./data/supplementary`](https://github.com/sina-mansour/neural-identity/blob/master/data/suplementary)

### Other:

- Simple [Jupyter notebooks](https://github.com/sina-mansour/neural-identity/tree/master/notebooks) to provide examples on how to use the codes (or check the [reproducible interactive binder image](https://mybinder.org/v2/gh/sina-mansour/neural-identity/master))

---

## Requirements:

---

All the code provided is written in **python 3**. The following python packages are required to use the codes:

- scipy
- numpy
- sklearn
- nibabel
- gdist

You may use a [python virtual environment](https://docs.python.org/3/library/venv.html) to install all of the required packages from [`./codes/requirements.txt`](https://github.com/sina-mansour/neural-identity/blob/master/codes/requirements.txt). The following code will take care of virtual environment setup (run in the git repository's main directory):

```bash
python3 -m venv venv/
source venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r ./codes/requirements.txt
deactivate

```

After setting up the virtual environment, you may use it by running:

```bash
source venv/bin/activate
```



---

## Codes:

---

### High-resolution functional connectivity

As mentioned in our paper, we proposed a method for sparsification of dense functional connectomes to be used in high-resolution studies. A python implementation of this method can be found in [`./codes/high_resolution.py`](https://github.com/sina-mansour/neural-identity/blob/master/codes/high_resolution.py). There is also an <span style="color:red">ipython notebook</span> with basic examples on how to use the codes.

Briefly, a sparse high-resolution mask generated from group average data is used as a sparsifier to threshold the individual high-resolution functional connectomes. This approach reduces the required memory (and hence computational load) of processing high-resolution functional connectomes.

---

### High-resolution structural connectivity

We have also proposed a method to generate smoothed high-resolution structural connectomes from tractography outputs. A python implementation of this method can be found in [`./codes/high_resolution.py`](https://github.com/sina-mansour/neural-identity/blob/master/codes/high_resolution.py). There is also an <span style="color:red">ipython notebook</span> with basic examples on how to use the codes.

The tractography input is used to locate the endpoints of the streamlines. The endpoints are warped into the standard space and mapped to closest surface vertices. The half incidence matrices (described in detail in our paper) are generated from the endpoint information. The half incidence matrices are first smoothed on the cortical surface mesh. The smoothed incidence information is used to compute the high-resolution structural connectivity on the cortical surface mesh.

---

### spintest

A python implementation of the spintest correspondence can be found in [`./codes/spintest.py`](https://github.com/sina-mansour/neural-identity/blob/master/codes/spintest.py). We have also provided an example [ipython notebook](https://github.com/sina-mansour/neural-identity/blob/master/notebooks/spintest.ipynb) on how to execute spintest. You may open this notebook using its [binder image](https://hub.gke.mybinder.org/user/sina-mansour-neural-identity-czo0ek04/notebooks/notebooks/spintest.ipynb) for a hands-on example.

This idea was proposed and implemented as a MATLAB package by the following article:

[**`Alexander-Bloch, Aaron F., et al. "On testing for spatial correspondence between maps of human brain structure and function." *Neuroimage* 178 (2018): 540-551.`**](https://doi.org/10.1016/j.neuroimage.2018.05.070)

This is a statistical test for correspondence of two spatial cortical surface maps. This test accounts for spatial correlation, thereby alleviating the overestimation of statistical significance in parametric p-values. The spin-test approach implements a spatial permutation framework by producing random rotations of the spherical representations of the cortical surface maps. This creates a null distribution to test the correspondence between the maps compared to the random overlays of maps. This approach provides the basis to account for the co-linearity caused by spatial dependence while controlling for the contralateral symmetries of the maps.

---

### VCM

A python implementation of variance component model (VCM) used in our paper to draw inference about behavior prediction performance can be found in [`./codes/vcm.py`](https://github.com/sina-mansour/neural-identity/blob/master/codes/vcm.py). We have also provided an example [ipython notebook](https://github.com/sina-mansour/neural-identity/blob/master/notebooks/VCM.ipynb) on how VCM is used. You may open this notebook using its [binder image](https://hub.gke.mybinder.org/user/sina-mansour-neural-identity-czo0ek04/notebooks/notebooks/VCM.ipynb) for a hands-on example.

This method was originally proposed and implemented as a MATLAB package in the following article:

[**`Sabuncu, Mert R., et al. "Morphometricity as a measure of the neuroanatomical signature of a trait." *Proceedings of the National Academy of Sciences* 113.39 (2016): E5749-E5756.`**](https://doi.org/10.1073/pnas.1604378113)

VCM was used to capture the association between neuroimaging measures and the behavioral characteristics.


$$
y = X \beta + a + \epsilon
$$


Where $y$ is a  $N \times 1$ vector of the selected behavioral characteristic and $N$ is the number of subjects. $X$ is the nuisance covariate matrix, $\beta$ is the fixed effect vector. The subject specific effect vector $a$ is drawn from a multivariate Gaussian distribution $a \sim N(0,\sigma_a^2 S_a)$ with a covariance matrix which is sourced from the test-test similarity matrix $S_{test-test}$. $\sigma_a^2$ is a scaling constant estimated from the data. The noise vector $\epsilon$ is drawn from a zero-mean independent normal distribution with homogeneous variance $\sigma_e^2$. The scaling constant $\sigma_a^2$ is interpreted as the degree in which a measure's similarity contributes to the observed behavioral characteristic. In other words, when $\sigma_a^2$ is large compared to $\sigma_e^2$ the individual behavioral differences are better captured by the measure similarities.

---

### Tractography scripts

We have provided a sample bash script in [`./codes/tractography/mrtrix_tractography.sh`](https://github.com/sina-mansour/neural-identity/blob/master/codes/tractography/mrtrix_tractography.sh) which was used to generate the tractography streamlines used in our study.

This script implements a Multi-Shell Multi-Tissue (MSMT) Constrained Spherical Deconvolution (CSD) Anatomically Constrained Tractography (ACT) method. First, an unsupervised method was used to estimate 3-tissue response functions for White Matter (WM), Grey Matter (GM), and Cerebro-Spinal Fluid (CSF). The fiber orientation distribution was then estimated using the MSMT CSD method. For the ACT, the 32k surface mesh was used to create a ribbon around the interface of WM and cortical GM. Probabilistic streamline tractography was then performed by 2nd order integration over fiber orientation distributions. 5 million streamlines were generated for each scan. The streamlines were uniformly seeded from the cortical GM-WM interface ribbon to create a fair chance of starting streamlines from any point on the cortical surface. A combined mask of WM and sub-cortical GM which was generated by a five-tissue-type segmentation by FSL BET was used to trim the final streamlines in order to ensure that the streamlines end at the cortical WM-GM boundary.

Note: This is just a sample script to show our exact tractography procedure. This script can be used as a starting guideline if you plan to implement the same procedure. However, the directory structure used in the script needs to modified by the user to locate the necessary initial files needed.

---

## Data

---

### Neural uniqueness maps

The localized uniqueness measures mapped for different neuroimaging modalities (structural connectivity, functional connectivity, cortical thickness, cortical curvature, sulcal depth, and myelination) and resolutions (high-resolution vertex-level maps, and atlas level maps for Glasser and Gordon atlases) are available in [`./data/uniqueness/`](https://github.com/sina-mansour/neural-identity/blob/master/data/uniqueness/) and can be opened by connectome viewer using the spec file in [`./data/uniqueness/Uniqueness.wb_spec`](https://github.com/sina-mansour/neural-identity/blob/master/data/uniqueness/Uniqueness.wb_spec).



There are 3 different sets of uniqueness maps in total:

- The high-resolution uniqueness maps are in: [`./data/uniqueness/High resolution`](https://github.com/sina-mansour/neural-identity/blob/master/data/uniqueness/High resolution)
- The atlas-level uniqueness maps generated from Glasser and Gordon atlas are in: [`./data/uniqueness/glasser`](https://github.com/sina-mansour/neural-identity/blob/master/data/uniqueness/glasser) and [`./data/uniqueness/gordon`](https://github.com/sina-mansour/neural-identity/blob/master/data/uniqueness/gordon) respectively
- The holistic approach (high-resolution maps averaged on atlases) can be found in: [`./data/uniqueness/Holistic mean`](https://github.com/sina-mansour/neural-identity/blob/master/data/uniqueness/Holistic mean)

---

### High-resolution connectomes

All the individual high-resolution connectomes produced with our pipeline could be accessed in our [publicly accessible S3 bucket](https://swift.rc.nectar.org.au/v1/AUTH_ee5989f4c9184ea29012bb124cd3dff0/connectome_storage/index.html). The directory structure is similar to the HCP directory structure. The connectomes are saved as sparse matrices that can be loaded using [scipy.sparse.load_npz](https://docs.scipy.org/doc/scipy-1.1.0/reference/generated/scipy.sparse.load_npz.html).

---

### Functional connectivity sparsifier

We utilized a sparse mask generated from the group average dense functional connectivity for sparsification of individual connectomes. The binary mask can be found in [`./data/sparse_mask/functional_sparse_mask_1%_density.npz`](https://github.com/sina-mansour/neural-identity/blob/master/data/sparse_mask/functional_sparse_mask_1%_density.npz). This is a sparse matrix that can be loaded using [scipy.sparse.load_npz](https://docs.scipy.org/doc/scipy-1.1.0/reference/generated/scipy.sparse.load_npz.html).

---

### Supplementary data

As well as the data mentioned above;

- We have provided the label indices for the two atlases used in our study.
- We have also provided the exact subject IDs comprising our **test**, and **retest** datasets

These are provided as json files located in [`./data/supplementary`](https://github.com/sina-mansour/neural-identity/blob/master/data/suplementary).
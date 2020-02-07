# Neural correlates of identity and behavior

---

Here, you may find the resources (code and data) used for our article on neural correlates of identity and behavior. All resources are provided as complementary to the following article:

**reference to the article will be added after acceptance...**

---

## TLDR

This repository contains:

### Codes:

- Codes used for generation of high-resolution sparse functional connectomes: <span style="color:red">`./codes/high_resolution.py`</span>
- Codes used for generation of high-resolution network-smoothed structural connectomes: <span style="color:red">`./codes/high_resolution.py`</span>
- Implementation of **spintest** non-parametric spatial correspondence test in python: <span style="color:red">`./codes/spintest.py`</span>
- Implementation of variance component modelling (**VCM**) in python: <span style="color:red">`./codes/vcm.py`</span>
- Sample script we used to run mrtrix tractography: <span style="color:red">`./codes/tractography/mrtrix_tractography.py`</span>

### Data:

- Surface maps of the expressed neural uniqueness from different modalities
- Samples of the high-resolution group average connectomes
- Link and information on how to access the individual high-resolution connectomes
- Sparse mask used to sparsify functional connectomes

### Other:

- Simple ipython notebooks to provide examples on how to use the provided code

---

## Requirements:

All the code provided is written in **python 3**. The following python packages are required to use the codes:

- scipy
- numpy
- sklearn
- nibabel
- gdist

You may use a [python virtual environment](https://docs.python.org/3/library/venv.html) to install all of the required packages from <span style="color:red">`./codes/requirements.txt`</span>. The following code will take care of virtual evironment setup (run in the git repository main directory):

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

---

### High-resolution structural connectivity

---

### spintest

A python implementation of the spintest correspondence can be found in `./codes/spintest.py`. There is also an <span style="color:red">ipython notebook</span> with basic examples on how to use the codes.

This idea was proposed and implemented as a MATLAB package by the following article:

**`Alexander-Bloch, Aaron F., et al. "On testing for spatial correspondence between maps of human brain structure and function." *Neuroimage* 178 (2018): 540-551.`**

This is a statistical test for correspondence of two spatial cortical surface maps. This test accounts for spatial correlation, thereby alleviating the overestimation of statistical significance in parametric p-values. The spin-test approach implements a spatial permutation framework by producing random rotations of the spherical representations of the cortical surface maps. This creates a null distribution to test the correspondence between the maps compared to the random overlays of maps. This approach provides the basis to account for the co-linearity caused by spatial dependence while controlling for the contralateral symmetries of the maps.

---

### VCM

A python implementation of variance component model (VCM) used in our paper to draw inference about behavior prediction performance can be found in `./codes/vcm.py`. This <span style="color:red">ipython notebook</span> provides a sample example on how VCM could be used.

This method was originally proposed and implemented as a MATLAB package in the following article:

**`Sabuncu, Mert R., et al. "Morphometricity as a measure of the neuroanatomical signature of a trait." *Proceedings of the National Academy of Sciences* 113.39 (2016): E5749-E5756.`**

VCM was used to capture the association between neuroimaging measures and the behavioral characteristics.


$$
y = X \beta + a + \epsilon
$$


Where $y$ is a  $N \times 1$ vector of the selected behavioral characteristic and $N$ is the number of subjects. $X$ is the nuisance covariate matrix, $\beta$ is the fixed effect vector. The subject specific effect vector $a$ is drawn from a multivariate Gaussian distribution $a \sim N(0,\sigma_a^2 S_a)$ with a covariance matrix which is sourced from the test-test similarity matrix $S_{test-test}$. $\sigma_a^2$ is a scaling constant estimated from the data. The noise vector $\epsilon$ is drawn from a zero-mean independent normal distribution with homogeneus variance $\sigma_e^2$. The scaling constant $\sigma_a^2$ is interpreted as the degree in which a measure's similarity contributes to the observed behavioral characteristic. In other words, when $\sigma_a^2$ is large compared to $\sigma_e^2$ the individual behavioral differences are better captured by the measure similarities.

---


# BrainEigenmodes
Code and data for the manuscript "[Geometric constraints on human brain function](https://www.biorxiv.org/content/10.1101/2022.10.04.510897v1)"

Uses an eigenmode-based analysis to study how the geometric structure of the brain constrains function captured by fMRI

The code also serves as a toolbox for people to calculate surface and/or volume geometric eigenmodes using their own data

## File descriptions

1. `data`: folder containing example data, parcellations, data results from the paper, template eigenmodes, and template surfaces and volumes
2. `functions_matlab`: folder containing various utility analysis and visualization MATLAB functions
3. `surface_eigenmodes.py`: Python code to calculate geometric eigenmodes of a cortical surface
4. `volume_eigenmodes.py`: Python code to calculate geometric eigenmodes of a 3D volume
5. `demo_eigenmode_calculation.sh`: Bash shell code to demonstrate how to calculate geometric eigenmodes
6. `demo_eigenmode_analysis.m`: MATLAB code to demonstrate how to calculate use eigenmodes to analyze fMRI data
7. `demo_eigenmode_visualization.m`: MATLAB code to demonstrate how to visualize eigenmodes
8. `demo_connectome_eigenmode_calculation.m`: MATLAB code to demonstrate how to calculate connectome and EDR connectome eigenmodes

## Installation

Download the repository and you're good to go.
Read the comments and documentation within each code to guide you in your usage.

## Downloading data

Due to their file sizes exceeding the limit allowed in GitHub, you will need to fill the `data/examples`, `data/examples`, `data/results`, and `data/template_eigenmodes` directories with some data that you can downloaded from this Zenodo repository (LINK TO BE UPDATED).

All empirical data are from the [Human Connectome Project](https://db.humanconnectome.org/). Please consult the link for detailed information about access and licensing.

## Dependencies

Some important aspects you need to do before running demo_eigenmode_calculation.sh, surface_eigenmodes.py, or volume_eigenmodes.py

1. Install [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) and load module (only for HPC systems).
2. Install [Connectome Workbench](https://www.humanconnectome.org/software/get-connectome-workbench) and load module (only for HPC systems).
3. Install [Gmsh](https://gmsh.info/) and load module (only for HPC systems).
4. Install the following Python libraries (e.g., via pip): [lapy](https://github.com/Deep-MI/LaPy), [nibabel](https://nipy.org/nibabel/), and [brainspace](https://brainspace.readthedocs.io/en/latest/pages/install.html)

Some important aspects you need to do before running the MATLAB-based scripts

1. Download the [gifti](https://github.com/gllmflndn/gifti) library and add to Matlab path.
2. Download the [cifti](https://github.com/Washington-University/cifti-matlab) library and add to Matlab path.
3. Download [cbrewer](https://au.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab) and add to Matlab path.
4. Download Stuart Oldham's [repository](https://github.com/StuartJO/plotSurfaceROIBoundary) for drawing ROI boundaries on a surface and add to Matlab path. 

## Compatibility

The codes have been tested on versions of Python 3.8 and MATLAB from R2019b to R2020b.

## Citation

If you use our code in your research, please cite us as follows:

[PREPRINT] J.C. Pang, K.M. Aquino, M. Oldehinkel, P.A. Robinson, B.D. Fulcher, M. Breakspear, A. Fornito, Geometric constraints on human brain function, bioRxiv (2022) (DOI: [10.1101/2022.10.04.510897](https://www.biorxiv.org/content/10.1101/2022.10.04.510897v1))

[ARTICLE] J.C. Pang, K.M. Aquino, M. Oldehinkel, P.A. Robinson, B.D. Fulcher, M. Breakspear, A. Fornito, Geometric constraints on human brain function (TO BE UPDATED)

## Further details

Please contact james.pang1@monash.edu if you need any further details.

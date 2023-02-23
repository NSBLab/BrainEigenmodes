# BrainEigenmodes
Code and data for the manuscript "[Geometric constraints on human brain function](https://www.biorxiv.org/content/10.1101/2022.10.04.510897v2)"

Uses an eigenmode-based analysis to study how the geometric structure of the brain constrains function captured by fMRI

The code also serves as a toolbox for people to calculate surface and/or volume geometric eigenmodes using their own data

## File descriptions

1. `data`: folder containing example data, parcellations, empirical data, results from the paper, figure source data, template eigenmodes, and template surfaces and volumes
2. `functions_matlab`: folder containing various utility analysis and visualization MATLAB functions
3. `surface_eigenmodes.py`: Python script to calculate geometric eigenmodes of a cortical surface
4. `volume_eigenmodes.py`: Python script to calculate geometric eigenmodes of a 3D volume
5. `demo_eigenmode_calculation.sh`: Bash shell script to demonstrate how to calculate geometric eigenmodes
6. `demo_eigenmode_analysis.m`: MATLAB script to demonstrate how to use eigenmodes to analyze fMRI data
7. `demo_eigenmode_visualization.m`: MATLAB script to demonstrate how to visualize surface and volume eigenmodes
8. `demo_connectome_eigenmode_calculation.m`: MATLAB script to demonstrate how to calculate connectome and EDR connectome eigenmodes
9. `demo_wave_model_simulation.m`: MATLAB script to demonstrate how to simulate waves on a cortical surface using eigenmodes
10. `generate_paper_figures.m`: MATLAB script to generate the main figures of the paper
11. `generate_paper_suppfigures.m`: MATLAB script to generate the supplementary figures of the paper

## Installation

Download the repository and you're good to go.
Read the comments and documentation within each code for usage guidance.

## Downloading data

Due to their file sizes exceeding the limit allowed by GitHub, you will need to fill the `data/empirical`, `data/results`, and `data/template_eigenmodes` directories with data that you can download from this Zenodo repository (LINK TO BE UPDATED). 
NOTE: Some portions of `generate_paper_figures.m` and `generate_paper_suppfigures.m` require the abovementioned data. Hence, you need them to be saved in the appropriate folders to be able to run these two scripts smoothly.

All empirical data are from the [Human Connectome Project](https://db.humanconnectome.org/). Please consult the link for detailed information about access, licensing, and terms and conditions of usage.

## Dependencies

Some important aspects you need to do before running `demo_eigenmode_calculation.sh`, `surface_eigenmodes.py`, or `volume_eigenmodes.py`

1. Install [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) and load module (only for HPC systems).
2. Install [Connectome Workbench](https://www.humanconnectome.org/software/get-connectome-workbench) and load module (only for HPC systems).
3. Install [Gmsh](https://gmsh.info/) and load module (only for HPC systems).
4. Install the following Python libraries (e.g., via pip): [lapy](https://github.com/Deep-MI/LaPy), [nibabel](https://nipy.org/nibabel/), and [brainspace](https://brainspace.readthedocs.io/en/latest/pages/install.html)
5. Make sure you also have the following common Python libraries: numpy, scipy, os, argparse, subprocess

Some of the MATLAB-based scripts depend on packages developed by others. Copies of these packages have been stored in the `functions_matlab` folder to ensure version compatibility. However, please visit their respective links below to get more details and don't forget to give them some love.

1. [gifti](https://github.com/gllmflndn/gifti)
2. [cifti](https://github.com/Washington-University/cifti-matlab)
3. [cbrewer](https://au.mathworks.com/matlabcentral/fileexchange/58350-cbrewer2?s_tid=srchtitle) (NOTE: The link is for the new version cbrewer2, but the repo uses the older version.) 
4. Stuart Oldham's [repository](https://github.com/StuartJO/plotSurfaceROIBoundary) for drawing ROI boundaries on a surface
5. Frantisek Vasa's [repository](https://github.com/frantisekvasa/rotate_parcellation) for creating rotations of a parcellated map

## Compatibility

The codes have been tested on versions of Python from 3.7 to 3.8 and versions of MATLAB from R2019b to R2020b.

## Citation

If you use our code in your research, please cite us as follows:

[PREPRINT] J.C. Pang, K.M. Aquino, M. Oldehinkel, P.A. Robinson, B.D. Fulcher, M. Breakspear, A. Fornito, Geometric constraints on human brain function, bioRxiv (2022) (DOI: [10.1101/2022.10.04.510897](https://www.biorxiv.org/content/10.1101/2022.10.04.510897v2))

[ARTICLE] J.C. Pang, K.M. Aquino, M. Oldehinkel, P.A. Robinson, B.D. Fulcher, M. Breakspear, A. Fornito, Geometric constraints on human brain function (TO BE UPDATED)

## Further details

Please contact james.pang1@monash.edu if you need any further details.

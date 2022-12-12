#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate the eigenmodes of a volume

@author: James Pang and Kevin Aquino, Monash University, 2022
"""

# Import all the libraries
from lapy import Solver, TetIO
import nibabel as nib
import numpy as np
from scipy.interpolate import griddata
import os
from argparse import ArgumentParser
import subprocess

def get_tkrvox2ras(voldim, voxres):
    """Generate transformation matrix to switch between tetrahedral and volume space.

    Parameters
    ----------
    voldim : array (1x3)
        Dimension of the volume (number of voxels in each of the 3 dimensions)
    voxres : array (!x3)
        Voxel resolution (resolution in each of the 3 dimensions)

    Returns
    ------
    T : array (4x4)
        Transformation matrix
    """

    T = np.zeros([4,4]);
    T[3,3] = 1;

    T[0,0] = -voxres[0];
    T[0,3] = voxres[0]*voldim[0]/2;

    T[1,2] = voxres[2];
    T[1,3] = -voxres[2]*voldim[2]/2;


    T[2,1] = -voxres[1];
    T[2,3] = voxres[1]*voldim[1]/2;

    return T

# Define the tetrahedral mesh using Gmsh
def make_tetra_file(nifti_input_filename):
    """Generate tetrahedral version of the ROI in the nifti file.

    Parameters
    ----------
    nifti_input_filename : str
        Filename of input volume where the relevant ROI have voxel values = 1

    Returns
    ------
    tetra_file : str
        Filename of output tetrahedral vtk file
    """

    nifti_input_file_head, nifti_input_file_tail = os.path.split(nifti_input_filename)
    nifti_input_file_main, nifti_input_file_ext = os.path.splitext(nifti_input_file_tail)

    os.system('mri_mc ' + nifti_input_filename + ' 1 ' + nifti_input_file_head + '/rh.tmp_surface.vtk')
    os.system('mv -f ' + nifti_input_file_head + '/rh.tmp_surface.vtk ' + nifti_input_filename + '.vtk')

    geo_file = nifti_input_filename + '.geo'
    tria_file = nifti_input_filename + '.vtk'
    tetra_file = nifti_input_filename + '.tetra.vtk'

    file = tria_file.rsplit('/')
    inputGeo = file[len(file)-1]
    
    with open(geo_file, 'w') as writer:
        writer.write('Mesh.Algorithm3D=4;\n')
        writer.write('Mesh.Optimize=1;\n')
        writer.write('Mesh.OptimizeNetgen=1;\n')
        writer.write('Merge "'+inputGeo+'";\n')
        writer.write('Surface Loop(1) = {1};\n')
        writer.write('Volume(1) = {1};\n')
        writer.write('Physical Volume(1) = {1};\n')

    cmd = 'gmsh -3 -o ' + tetra_file + ' ' + geo_file
    output = subprocess.check_output(cmd,shell="True")
    output = output.splitlines()

    cmd = "sed 's/double/float/g;s/UNSTRUCTURED_GRID/POLYDATA/g;s/CELLS/POLYGONS/g;/CELL_TYPES/,$d' " + tetra_file + " > " + tetra_file + "'_fixed'"
    os.system(cmd)
    os.system('mv -f ' + tetra_file + '_fixed ' + tetra_file)
    
    return tetra_file

def calc_volume(nifti_input_filename):
    """Calculate the physical volume of the ROI in the nifti file.

    Parameters
    ----------
    nifti_input_filename : str
        Filename of input volume where the relevant ROI have voxel values = 1

    Returns
    ------
    ROI_number : int
        Total number of non-zero voxels
    ROI_volume : float
        Total volume of non-zero voxels in physical dimensions   
    """

    # Load data
    ROI_data = nib.load(nifti_input_filename)
    roi_data = ROI_data.get_fdata()

    # Get voxel dimensions in mm
    voxel_dims = (ROI_data.header["pixdim"])[1:4]
    voxel_vol = np.prod(voxel_dims)

    # Compute volume
    ROI_number = np.count_nonzero(roi_data)
    ROI_volume = ROI_number * voxel_vol

    # print("Number of non-zero voxels = {}".format(ROI_number))
    # print("Volume of non-zero voxels = {} mm^3".format(ROI_volume))

    return ROI_number, ROI_volume

def normalize_vtk(tet, nifti_input_filename, normalization_type='none', normalization_factor=1):
    """Normalize tetrahedral surface.

    Parameters
    ----------
    tet : lapy compatible object
        Loaded vtk object corresponding to a surface tetrahedral mesh
    nifti_input_filename : str
        Filename of input volume where the relevant ROI have voxel values = 1
    normalization_type : str (default: 'none')
        Type of normalization
        number - normalization with respect to the total number of non-zero voxels
        volume - normalization with respect to the total volume of non-zero voxels in physical dimensions   
        constant - normalization with respect to a chosen constant
        others - no normalization
    normalization_factor : float (default: 1)
        Factor to be used in a constant normalization     

    Returns
    ------
    tet_norm : lapy compatible object
        Loaded vtk object corresponding to the normalized surface tetrahedral mesh
    """

    nifti_input_file_head, nifti_input_file_tail = os.path.split(nifti_input_filename)
    nifti_input_file_main, nifti_input_file_ext = os.path.splitext(nifti_input_file_tail)

    ROI_number, ROI_volume = calc_volume(nifti_input_filename)

    # normalization process
    tet_norm = tet
    if normalization_type == 'number':
        tet_norm.v = tet.v/(ROI_number**(1/3))
    elif normalization_type == 'volume':
        tet_norm.v = tet.v/(ROI_volume**(1/3))
    elif normalization_type == 'constant':
        tet_norm.v = tet.v/(normalization_factor**(1/3))
    else:
        pass

    # writing normalized surface to a vtk file
    if normalization_type == 'number' or normalization_type == 'volume' or normalization_type == 'constant':
        surface_output_filename = nifti_input_filename + '_norm=' + normalization_type + '.tetra.vtk'

        f = open(surface_output_filename, 'w')
        f.write('# vtk DataFile Version 2.0\n')
        f.write(nifti_input_file_tail + '\n')
        f.write('ASCII\n')
        f.write('DATASET POLYDATA\n')
        f.write('POINTS ' + str(np.shape(tet.v)[0]) + ' float\n')
        for i in range(np.shape(tet.v)[0]):
            f.write(' '.join(map(str, tet_norm.v[i, :])))
            f.write('\n')
        f.write('\n')
        f.write('POLYGONS ' + str(np.shape(tet.t)[0]) + ' ' + str(5 * np.shape(tet.t)[0]) + '\n')
        for i in range(np.shape(tet.t)[0]):
            f.write(' '.join(map(str, np.append(4, tet.t[i, :]))))
            f.write('\n')
        f.close()

    return tet_norm

def calc_eig(nifti_input_filename, output_eval_filename, output_emode_filename, num_modes, normalization_type='none', normalization_factor=1):
    """Calculate the eigenvalues and eigenmodes of the ROI volume in a nifti file.

    Parameters
    ----------
    nifti_input_filename : str
        Filename of input volume where the relevant ROI have voxel values = 1
    output_eval_filename : str  
        Filename of text file where the output eigenvalues will be stored
    output_emode_filename : str  
        Filename of text file where the output eigenmodes (in tetrahedral surface space) will be stored    
    num_modes : int
        Number of eigenmodes to be calculated
    normalization_type : str (default: 'none')
        Type of normalization
        number - normalization with respect to the total number of non-zero voxels
        volume - normalization with respect to the total volume of non-zero voxels in physical dimensions   
        constant - normalization with respect to a chosen constant
        others - no normalization
    normalization_factor : float (default: 1)
        Factor to be used in a constant normalization      

    Returns
    ------
    evals: array (num_modes x 1)
        Eigenvalues
    emodes: array (number of tetrahedral surface points x num_modes)
        Eigenmodes
    """

    # convert the ROI in the nifti file to a tetrahedral surface
    tetra_file = make_tetra_file(nifti_input_filename)

    # load tetrahedral surface (as a brainspace object)
    tetra = TetIO.import_vtk(tetra_file)

    # normalize tetrahedral surface
    tetra_norm = normalize_vtk(tetra, nifti_input_filename, normalization_type, normalization_factor)

    # calculate eigenvalues and eigenmodes
    fem = Solver(tetra_norm)
    evals, emodes = fem.eigs(k=num_modes)
    
    output_eval_file_main, output_eval_file_ext = os.path.splitext(output_eval_filename)
    output_emode_file_main, output_emode_file_ext = os.path.splitext(output_emode_filename)

    if normalization_type == 'number' or normalization_type == 'volume' or normalization_type == 'constant':
        np.savetxt(output_eval_file_main + '_norm=' + normalization_type + output_eval_file_ext, evals)
        np.savetxt(output_emode_file_main + '_norm=' + normalization_type + output_emode_file_ext, emodes)
    else:
        np.savetxt(output_eval_filename, evals)
        np.savetxt(output_emode_filename, emodes)
    
    return evals, emodes

def calc_volume_eigenmodes(nifti_input_filename, nifti_output_filename, output_eval_filename, output_emode_filename, num_modes, normalization_type='none', normalization_factor=1):
    """Main function to calculate the eigenmodes of the ROI volume in a nifti file.

    Parameters
    ----------
    nifti_input_filename : str
        Filename of input volume where the relevant ROI have voxel values = 1
    nifti_output_filename : str  
        Filename of nifti file where the output eigenmdoes (in volume space) will be stored
    output_eval_filename : str  
        Filename of text file where the output eigenvalues will be stored
    output_emode_filename : str  
        Filename of text file where the output eigenmodes (in tetrahedral surface space) will be stored    
    num_modes : int
        Number of eigenmodes to be calculated
    normalization_type : str (default: 'none')
        Type of normalization
        number - normalization with respect to the total number of non-zero voxels
        volume - normalization with respect to the total volume of non-zero voxels in physical dimensions   
        constant - normalization with respect to a chosen constant
        others - no normalization
    normalization_factor : float (default: 1)
        Factor to be used in a constant normalization
    """

    # calculate eigenvalues and eigenmodes
    evals, emodes = calc_eig(nifti_input_filename, output_eval_filename, output_emode_filename, num_modes, normalization_type, normalization_factor)


    # project eigenmodes in tetrahedral surface space into volume space

    # prepare transformation
    ROI_data = nib.load(nifti_input_filename)
    roi_data = ROI_data.get_fdata()
    inds_all = np.where(roi_data==1)
    xx = inds_all[0]
    yy = inds_all[1]
    zz = inds_all[2]

    points = np.zeros([xx.shape[0],4])
    points[:,0] = xx
    points[:,1] = yy
    points[:,2] = zz
    points[:,3] = 1

    # calculate transformation matrix
    T = get_tkrvox2ras(ROI_data.shape, ROI_data.header.get_zooms())

    # apply transformation
    points2 = np.matmul(T, np.transpose(points))

    # load tetrahedral surface
    tetra_file = nifti_input_filename + '.tetra.vtk'
    tetra = TetIO.import_vtk(tetra_file)
    points_surface = tetra.v

    # initialize nifti output array
    new_shape = np.array(roi_data.shape)
    if roi_data.ndim>3:
        new_shape[3] = num_modes
    else:
        new_shape = np.append(new_shape, num_modes)
    new_data = np.zeros(new_shape)

    # perform interpolation of eigenmodes from tetrahedral surface space to volume space
    for mode in range(0, num_modes):
        interpolated_data = griddata(points_surface, emodes[:,mode], np.transpose(points2[0:3,:]), method='linear')
        for ind in range(0, len(interpolated_data)):
            new_data[xx[ind],yy[ind],zz[ind],mode] = interpolated_data[ind]

    # save to output nifti file
    img = nib.Nifti1Image(new_data, ROI_data.affine, header=ROI_data.header)
    nib.save(img, nifti_output_filename)

    # remove all created temporary auxiliary files
    geo_file = nifti_input_filename + '.geo'
    tria_file = nifti_input_filename + '.vtk'
    if os.path.exists(geo_file):
        os.remove(geo_file)
    if os.path.exists(tria_file):
        os.remove(tria_file)

def main(raw_args=None):    
    parser = ArgumentParser(epilog="volume_eigenmodes.py -- A function to calculate the eigenmodes of an ROI volume. James Pang, Monash University, 2022 <james.pang1@monash.edu>")
    parser.add_argument("nifti_input_filename", help="An input nifti with ROI voxel values=1", metavar="volume_input.nii.gz")
    parser.add_argument("nifti_output_filename", help="An output nifti where the eigenmodes in volume space will be stored", metavar="emodes.nii.gz")
    parser.add_argument("output_eval_filename", help="An output text file where the eigenvalues will be stored", metavar="evals.txt")
    parser.add_argument("output_emode_filename", help="An output text file where the eigenmods in tetrahedral surface space will be stored", metavar="emodes.txt")
    parser.add_argument("-N", dest="num_modes", default=20, help="Number of eigenmodes to be calculated, default=20", metavar="20")
    parser.add_argument("-norm", dest="normalization_type", default='none', help="Type of normalization of tetrahedral surface", metavar="none")
    parser.add_argument("-normfactor", dest="normalization_factor", default=1, help="Value of constant normalization factor of tetrahedral surface", metavar="1")

    #--------------------    Parsing the inputs from terminal:   -------------------
    args = parser.parse_args()
    nifti_input_filename    = args.nifti_input_filename
    nifti_output_filename   = args.nifti_output_filename
    output_eval_filename    = args.output_eval_filename
    output_emode_filename   = args.output_emode_filename
    num_modes               = int(args.num_modes)
    normalization_type      = args.normalization_type
    normalization_factor    = float(args.normalization_factor)
    #-------------------------------------------------------------------------------
   
    calc_volume_eigenmodes(nifti_input_filename, nifti_output_filename, output_eval_filename, output_emode_filename, num_modes, normalization_type, normalization_factor)
   

if __name__ == '__main__':
    
    # running via commandline
    main()
    

    # running within python
    # structures = ['tha']
    # hemispheres = ['lh', 'rh']
    # num_modes = 30
    # normalization_type = 'none'
    # normalization_factor = 1

    # for structure in structures:
    #     for hemisphere in hemispheres:

    #         nifti_input_filename = 'data/template_surfaces_volumes/hcp_' + structure + '-' + hemisphere + '_thr25.nii.gz'
    #         nifti_output_filename = 'data/template_eigenmodes/hcp_' structure + '-' + hemisphere + '_emode_' + str(num_modes) + '.nii.gz'
    #         output_eval_filename = 'data/template_eigenmodes/hcp' + structure + '-' + hemisphere + '_eval_' + str(num_modes) + '.txt'
    #         output_emode_filename = 'data/template_eigenmodes/hcp' + structure + '-' + hemisphere + '_emode_' + str(num_modes) + '.txt'
            
    #         calc_volume_eigenmodes(nifti_input_filename, nifti_output_filename, output_eval_filename, output_emode_filename, num_modes, normalization_type, normalization_factor):
    
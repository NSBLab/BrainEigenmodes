#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate the eigenmodes of a cortical surface

@author: James Pang and Kevin Aquino, Monash University, 2022
"""

# Import all the libraries
from lapy import Solver, TriaIO
import numpy as np
import nibabel as nib
import brainspace.mesh as mesh
import os
from argparse import ArgumentParser

def calc_eig(tria, totalEigenNumbers):
    """Calculate the eigenvalues and eigenmodes of a surface.

    Parameters
    ----------
    tria : lapy compatible object
        Loaded vtk object corresponding to a surface triangular mesh
    totalEigenNumbers : int
        Number of eigenmodes to be calculated

    Returns
    ------
    evals : array (totalEigenNumbers x 1)
        Eigenvalues
    emodes : array (number of surface points x totalEigenNumbers)
        Eigenmodes
    """
    
    fem = Solver(tria)
    evals, emodes = fem.eigs(k=totalEigenNumbers)
    
    return evals, emodes

def create_temp_surface(surface_input, surface_output_filename):
    """Write surface to a new vtk file.

    Parameters
    ----------
    surface_input : brainspace compatible object
        Loaded vtk object corresponding to a surface triangular mesh
    surface_output_filename : str
        Filename of surface to be saved
    """

    f = open(surface_output_filename, 'w')
    f.write('# vtk DataFile Version 2.0\n')
    f.write(surface_output_filename + '\n')
    f.write('ASCII\n')
    f.write('DATASET POLYDATA\n')
    f.write('POINTS ' + str(np.shape(surface_input.Points)[0]) + ' float\n')
    for i in range(np.shape(surface_input.Points)[0]):
        f.write(' '.join(map(str, np.array(surface_input.Points[i, :]))))
        f.write('\n')
    f.write('\n')
    f.write('POLYGONS ' + str(np.shape(surface_input.polys2D)[0]) + ' ' + str(4* np.shape(surface_input.polys2D)[0]) + '\n')
    for i in range(np.shape(surface_input.polys2D)[0]):
        f.write(' '.join(map(str, np.append(3, np.array(surface_input.polys2D[i, :])))))
        f.write('\n')
    f.close()

def get_indices(surface_original, surface_new):
    """Extract indices of vertices of the two surfaces that match.

    Parameters
    ----------
    surface_original : brainspace compatible object
        Loaded vtk object corresponding to a surface triangular mesh
    surface_new : brainspace compatible object
        Loaded vtk object corresponding to a surface triangular mesh

    Returns
    ------
    indices : array
        indices of vertices
    """

    indices = np.zeros([np.shape(surface_new.Points)[0],1])
    for i in range(np.shape(surface_new.Points)[0]):
        indices[i] = np.where(np.all(np.equal(surface_new.Points[i,:],surface_original.Points), axis=1))[0][0]
    indices = indices.astype(int)
    
    return indices

def calc_surface_eigenmodes(surface_input_filename, mask_input_filename, output_eval_filename, output_emode_filename, save_cut, totalEigenNumbers):
    """Main function to calculate the eigenmodes of a cortical surface with application of a mask (e.g., to remove the medial wall).

    Parameters
    ----------
    surface_input_filename : str
        Filename of input surface
    mask_input_filename : str
        Filename of mask to be applied on the surface (e.g., cortex without medial wall, values = 1 for mask and 0 elsewhere)
    output_eval_filename : str  
        Filename of text file where the output eigenvalues will be stored
    output_emode_filename : str  
        Filename of text file where the output eigenmodes will be stored
    save_cut : boolean 
        Boolean to decide if the new surface with mask applied will be saved to a new surface file
    totalEigenNumbers : int
        Number of eigenmodes to be calculated          
    """

    # load surface (as a brainspace object)
    surface_orig = mesh.mesh_io.read_surface(surface_input_filename)
    
    # load mask
    # can be any ROI (even whole cortex)
    mask_input_file_main, mask_input_file_ext = os.path.splitext(mask_input_filename)
    if mask_input_file_ext == '.txt':
        mask = np.loadtxt(mask_input_filename)
    elif mask_input_file_ext == '.gii':
        mask = nib.load(mask_input_filename).darrays[0].data
    
    # create temporary suface based on mask
    surface_cut = mesh.mesh_operations.mask_points(surface_orig, mask)

    if save_cut == 1:
        # old method: save vtk of surface_cut and open via lapy TriaIO 
        # The writing phase of this process is very slow especially for large surfaces
        temp_cut_filename='temp_cut.vtk'
        create_temp_surface(surface_cut, temp_cut_filename)
        # load surface (as a lapy object)
        tria = TriaIO.import_vtk(temp_cut_filename)
    else:
        # new method: replace v and t of surface_orig with v and t of surface_cut
        # faster version without the need to write the vtk file
        # load surface (as a lapy object)
        tria = TriaIO.import_vtk(surface_input_filename)
        tria.v = surface_cut.Points
        tria.t = np.reshape(surface_cut.Polygons, [surface_cut.n_cells, 4])[:,1:4]

    # calculate eigenvalues and eigenmodes
    evals, emodes = calc_eig(tria, totalEigenNumbers)
    
    # get indices of vertices of surface_orig that match surface_cut
    indices = get_indices(surface_orig, surface_cut)
    
    # reshape emodes to match vertices of original surface
    emodes_reshaped = np.zeros([surface_orig.n_points,np.shape(emodes)[1]])
    for mode in range(np.shape(emodes)[1]):
        emodes_reshaped[indices,mode] = np.expand_dims(emodes[:,mode], axis=1);
        
    # save results to text files
    np.savetxt(output_eval_filename, evals)
    np.savetxt(output_emode_filename, emodes_reshaped)

    if save_cut == 0:
        if os.path.exists('temp_cut.vtk'):
            os.remove('temp_cut.vtk')

def calc_surface_eigenmodes_nomask(surface_input_filename, output_eval_filename, output_emode_filename, totalEigenNumbers):
    """Main function to calculate the eigenmodes of a cortical surface without application of a mask.

    Parameters
    ----------
    surface_input_filename : str
        Filename of input surface
    output_eval_filename : str  
        Filename of text file where the output eigenvalues will be stored
    output_emode_filename : str  
        Filename of text file where the output eigenmodes will be stored
    totalEigenNumbers : int
        Number of eigenmodes to be calculated          
    """

    # load surface (as a lapy object)
    tria = TriaIO.import_vtk(surface_input_filename)

    # calculate eigenvalues and eigenmodes
    evals, emodes = calc_eig(tria, totalEigenNumbers)
        
    # save eigenmode results
    np.savetxt(output_eval_filename, evals)
    np.savetxt(output_emode_filename, emodes)

def main(raw_args=None):    
    parser = ArgumentParser(epilog="surface_eigenmodes.py -- A function to calculate the eigenmodes of a cortical surface. James Pang, Monash University, 2022 <james.pang1@monash.edu>")
    parser.add_argument("surface_input_filename", help="An input surface in vtk format", metavar="surface_input.vtk")
    parser.add_argument("mask_input_filename", help="An input mask text or gifti file", metavar="mask.txt")
    parser.add_argument("output_eval_filename", help="An output text file where the eigenvalues will be stored", metavar="evals.txt")
    parser.add_argument("output_emode_filename", help="An output text file where the eigenmodes will be stored", metavar="emodes.txt")
    parser.add_argument("-save_cut", dest="save_cut", default=0, help="Logical value to decide whether to write the masked version of the input surface", metavar="0")
    parser.add_argument("-N", dest="totalEigenNumbers", default=20, help="Number of eigenmodes to be calculated, default=20", metavar="20")
    parser.add_argument("-is_mask", dest="is_mask", default=1, help="Logical value to decide whether to apply the mask", metavar="1")

    #--------------------    Parsing the inputs from terminal:   -------------------
    args = parser.parse_args()
    surface_input_filename   = args.surface_input_filename
    mask_input_filename      = args.mask_input_filename
    output_eval_filename     = args.output_eval_filename
    output_emode_filename    = args.output_emode_filename
    save_cut                 = int(args.save_cut)
    totalEigenNumbers        = int(args.totalEigenNumbers)
    is_mask                  = int(args.is_mask)
    #-------------------------------------------------------------------------------
   
    if is_mask == 1:
        calc_surface_eigenmodes(surface_input_filename, mask_input_filename, output_eval_filename, output_emode_filename, save_cut, totalEigenNumbers)
    else:
        calc_surface_eigenmodes_nomask(surface_input_filename, output_eval_filename, output_emode_filename, totalEigenNumbers)
    
   
if __name__ == '__main__':
    
    # running via commandline
    main()
    

    # # running within python
    # surface_interest = 'fsLR_32k'
    # structure = 'midthickness'
    # hemispheres = ['lh', 'rh']
    # totalEigenNumbers = 200
    # save_cut = 0
    
    # for hemisphere in hemispheres:
    #     print('Processing ' + hemisphere)

    #     surface_input_filename = 'data/template_surfaces_volumes/' + surface_interest + '_' + structure + '-' + hemisphere + '.vtk'
    #     mask_input_filename = 'data/template_surfaces_volumes/' + surface_interest + '_cortex-' + hemisphere + '_mask.txt'
        
    #     # with cortex mask (remove medial wall)
    #     # this is the advisable way
    #     output_eval_filename = 'data/template_eigenmodes/' + surface_interest + '_' + structure + '-' + hemisphere + '_eval_' + str(totalEigenNumbers) + '.txt'
    #     output_emode_filename = 'data/template_eigenmodes/' + surface_interest + '_' + structure + '-' + hemisphere + '_emode_' + str(totalEigenNumbers) + '.txt'

    #     calc_surface_eigenmodes(surface_input_filename, mask_input_filename, output_eval_filename, output_emode_filename, save_cut, totalEigenNumbers)
        
    #     # without cortex mask
    #     output_eval_filename = 'data/template_eigenmodes/' + 'nomask_' + surface_interest + '_' + structure + '-' + hemisphere + '_eval_' + str(totalEigenNumbers) + '.txt'
    #     output_emode_filename = 'data/template_eigenmodes/' + 'nomask_' + surface_interest + '_' + structure + '-' + hemisphere + '_emode_' + str(totalEigenNumbers) + '.txt'

    #     calc_surface_eigenmodes_nomask(surface_input_filename, output_eval_filename, output_emode_filename, totalEigenNumbers)
            
 

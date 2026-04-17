"""
Copyright (C) 2020-2026 Nico Marioni <nmarioni@seas.upenn.edu>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
""" 

# PrO-VAT: Probe-Occupiable Volume Analysis Tools
# PrO-VAT.py calculates the pore size distribution (free volume distribution, channel width distribution, etc) of the van der Waals volume of the defined system matrix from a GROMACS (xtc/trr/gro + tpr/gro) or PoreBlazer-style (xyz + dat) trajectory
# This script was specifically desgined to find the distribution of water-rich pores within a hydrated polymer system, but can be generalized to any atomic or coarse-grained system
# The output includes the Cumulative Pore Size Distribution (Cumulative PSD), Pore Size Distribution (PSD), and Free Volume Fraction (Fractional Free Volume, FFV), with optional Surface Area (SA), Tortusoity (Tau), and xyz visualizations
# Given the system matrix (e.g., polymer matrix, non-polar domain, polar domain, etc), this code can analyze the total free volumne, free volume of the largest cluster (assumed percolated), or the free volume of clusters containing (defined) solvent atoms
# This script was written based on the methods used for PoreBlazer (https://github.com/SarkisovGitHub/PoreBlazer) and is optimized for parallelized calculations over many system frames, or analysis of large (30+ nm box length) systems
#
# As written, this code reads in GROMACS trajectory or PoreBlazer xyz and dat data using MDAnalysis
# As written, this code is designed for 3D-periodic rectangular simulations
#
# When implementing this code, it is recommended to test different values of L_voxel to ensure convergence of the FFV as L_voxel decreases. Note, computation time and memory usage will grow significantly as L_voxel decreases.
# If you run into memory or extreme run times, there are debugging lines throughout the code, and several values you can change to increase or decrease memory usage.
# There are three instances where xyz files can be created to visualize 1) probe-occupiable spheres of maximum radius without overlapping the van der Waals volume of the system, 2) voxel-centers that lie within the probe-occupiable volume, and 3) voxel-center surfaces that define the Connolly or Lee-Richards surface of the probe-occupiable volume

import numpy as np
import h5py
import MDAnalysis as mda
import MDAnalysis.lib.distances as distances
from igraph import Graph
from skimage import measure
import porespy as ps

import multiprocessing as mp
import functools
import os
import time
import yaml
import argparse

# suppresses expected warning in tortuosity analysis
ps.settings.loglevel = 'ERROR'

# Set the float data type used for atom coordinates and free volume sphere radii
# Always try np.float64 first
# If memory is an issue, use np.float32 - may introduce some error due to lack of precision
float_type = np.float64

def volume_analysis(frame):
# Perform volume analysis on current frame

    # Sleep command to offset processes (limit spikes in memory usage) - no delay if N_threads = 1
    time.sleep(np.random.randint(1,6)*(frame%args.N_threads))

    with h5py.File('PrO-VAT.hdf5','r') as f:
        dset1 = f['system']; sys = dset1[frame]                                                                                                 # Position of all system atoms
        dset2 = f['sys_radii']; sys_radii = dset2[:]                                                                                            # van der Waals radii of all system atoms
        dset3 = f['solvent']; sol = dset3[frame]                                                                                                # Position of all solvent atoms
        dset3 = f['cells']; cell = dset3[frame]                                                                                                 # Size of the cell
        dset4 = f['frames']; frame_ids = dset4[:]; frame_ids = np.arange(0,len(frame_ids),1)
    N_sys = len(sys); N_sol = len(sol)
    
    # Efficiency parameters used to determine N_cube
    avg_sys_density = N_sys / cell[0] / cell[1] / cell[2]
    vol_d_inc = (4/3) * np.pi * (args.d_inc**3)

    # Track which frames are currently being processed
    print(f"Frame {frame}/{len(frame_ids)}")





    # Free volume sphere analysis
    # This part of the calculation determines the maximum size of voxel-centered free volume spheres without overlapping system atoms (the van der Waals volume), where the total volume of all spheres larger than probe_radius defines the probe-occupiable free volume of the system
    # This code will generate points at the center of voxels with side length L_voxel and grow these points into free volume spheres
    # Changing L_voxel, N_write_sph, and d_inc can reduce run time and memory usage
    vox_x = np.linspace(0, cell[0], num = np.ceil(cell[0]/args.L_voxel).astype(int), dtype=float_type); vox_x = (vox_x[:-1] + vox_x[1:])/2      # Voxel-centers in the x direction
    vox_y = np.linspace(0, cell[1], num = np.ceil(cell[1]/args.L_voxel).astype(int), dtype=float_type); vox_y = (vox_y[:-1] + vox_y[1:])/2      # Voxel-centers in the y direction
    vox_z = np.linspace(0, cell[2], num = np.ceil(cell[2]/args.L_voxel).astype(int), dtype=float_type); vox_z = (vox_z[:-1] + vox_z[1:])/2      # Voxel-centers in the z direction

    # True L_voxel in x, y, and z
    L_voxel_x = vox_x[1] - vox_x[0]
    L_voxel_y = vox_y[1] - vox_y[0]
    L_voxel_z = vox_z[1] - vox_z[0]
    # Box lengths in units of number of voxels
    l_x = len(vox_x)
    l_y = len(vox_y)
    l_z = len(vox_z)

    # Use smallest integer data types possible (without losing precision) to reduce memory usage
    indexed_type = np.min_scalar_type(np.max([l_x, l_y, l_z]) - 1)
    linear_type = np.min_scalar_type((l_x * l_y * l_z) - 1)
    signed_linear_type = np.min_scalar_type(-((l_x * l_y * l_z) - 1))
    
    if args.Voxel_dist == 'Random':
        # Add random offsets to break up the uniformity of the voxels
        vox_x = vox_x + np.random.uniform(low=-L_voxel_x/2, high=L_voxel_x/2, size=vox_x.size)
        vox_y = vox_y + np.random.uniform(low=-L_voxel_y/2, high=L_voxel_y/2, size=vox_y.size)
        vox_z = vox_z + np.random.uniform(low=-L_voxel_z/2, high=L_voxel_z/2, size=vox_z.size)

    radii_arr = np.zeros((l_x,l_y,l_z), dtype=float_type)                                                                                       # radii_arr tracks free volume sphere indices (position in array = position in voxelized system) and radius (value at that position), where we are interested in spheres of radius r >= probe_radius. All probes of r > 0 are saved for later use.

    # Divide the voxelized system into voxel cubes for efficient analysis
    N_cube = args.N_write_sph / (avg_sys_density * vol_d_inc); L_cube = N_cube**(1/3)                                                           # To improve efficiency, voxels are looped over in cubes of N_cube voxel-centers
    vox_inc = np.ceil(                                                                                                                          # vox_inc = side length of voxel cube, such that each voxel cube is about the same size
                        np.min((l_x, l_y, l_z))
                      / np.ceil(
                                  np.min((l_x, l_y, l_z))
                                / (L_cube)
                                                         )
                                                          )
    vox_inc = np.min((np.ceil(np.min((l_x, l_y, l_z)) / 2), vox_inc)).astype(int)                                                               # For a cubic cell, forces there to be a minimum of 8 cubes
    N_cube = vox_inc**3                                                                                                                         # Actual number of voxels in each voxel cube after making each cube approximately the same size
    vox_track = np.array((0,0,0), dtype=int); vox_track[0] = -vox_inc                                                                           # vox_track tracks the location of the cubes in x, y, and z compared to the position in vox_x, vox_y, and vox_z

    # Prints the number of voxels-per-cube and number of voxel cubes
    if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
        time_Spheres = time.perf_counter()
        print('##### Generating Free Volume Spheres #####')
        print(f"\nNumber of voxels-per-cube: {N_cube}")
        print(f"Number of voxel cubes: {np.ceil(l_x/vox_inc).astype(int)*np.ceil(l_y/vox_inc).astype(int)*np.ceil(l_z/vox_inc).astype(int)}")

    for x_i in np.arange(vox_inc,l_x+vox_inc,vox_inc):
        vox_track[0] += vox_inc
        if x_i > l_x: x_i = l_x

        vox_track[1] = -vox_inc
        for y_i in np.arange(vox_inc,l_y+vox_inc,vox_inc):
            vox_track[1] += vox_inc
            if y_i > l_y: y_i = l_y
    
            vox_track[2] = -vox_inc
            for z_i in np.arange(vox_inc,l_z+vox_inc,vox_inc):
                vox_track[2] += vox_inc
                if z_i > l_z: z_i = l_z

                sphere_temp = np.vstack(np.meshgrid(                                                                                            # sphere_temp contains the position of the voxel-centers within the cube of size N_cube
                                                    vox_x[vox_track[0]:x_i],
                                                    vox_y[vox_track[1]:y_i],
                                                    vox_z[vox_track[2]:z_i]
                                                                           ), dtype=float_type).reshape(3,-1).T

                # Find the approximate center of the voxel cube to find the system atoms near the voxel cube (sys_mask), where system atoms define the van der Waals volume of the system. Reduces computational cost
                center = np.array([
                                   vox_x[vox_track[0] + int((x_i - vox_track[0])/2)],
                                   vox_y[vox_track[1] + int((y_i - vox_track[1])/2)],
                                   vox_z[vox_track[2] + int((z_i - vox_track[2])/2)]
                                                                                    ], dtype=float_type)

                # To reduce the number of calculations and limit memory usage, the distance between voxel-centers and system atoms is done in steps of d_inc Angstroms
                d = 0.0                                                                                                                         # Maximum distance to calculate between every voxel-center and every system atom
                while len(sphere_temp) > 0:
                    d += args.d_inc

                    sys_mask = distances.capped_distance(center, sys, d + np.sqrt(3)*vox_inc*args.L_voxel/2 + 2*args.L_voxel, box=cell, return_distances=False)[:,1] # System atoms near the voxel cube

                    pair_arr, dist_arr = distances.capped_distance(sphere_temp, sys[sys_mask], d, box=cell)                                     # Distance between voxel-centers and system atoms

                    # Useful print command for troubleshooting memory problems
                    # Decreasing N_write_sph will reduce the number of distances generated each cycle, reducing memory usage
                    if (args.print_eff == 2) and (frame == frame_ids[-1] or args.N_threads == 1):
                        if d == args.d_inc:
                            print(f"\nVoxel block: {(vox_track/vox_inc).astype(int)}")
                        print(f"Distance, calculations, writes: {d:3.1f} {len(sphere_temp)*len(sys[sys_mask]):.1e} {len(dist_arr):.1e}")

                    if len(dist_arr) > 0:
                        dist_arr -= sys_radii[sys_mask][pair_arr[:,1]]                                                                          # Subtract radius of each system atom from the distance to get the distance from the voxel-center to the surface of the atom

                        # Fill radii_arr for all voxel-centers that contain system atoms within d distance, where the smallest distance is the radius of the free volume sphere centered on the voxel
                        index = 0; sph_save = pair_arr[0,0]
                        sphere_remove = []
                        for i,sph in enumerate(pair_arr[:,0]):
                            if (sph > sph_save) or (i+1 == len(pair_arr[:,0])):
                                r_min = np.min(dist_arr[index:i])                                                                               # Minimum distance between voxel-center and system surface

                                if r_min > 0:                                                                                                   # Sphere does not overlap the system and radius >= 0
                                    coords = np.divide(
                                                       sphere_temp[sph_save],
                                                       np.array([L_voxel_x,
                                                                 L_voxel_y,
                                                                 L_voxel_z])
                                                                            ).astype(indexed_type)
                                    radii_arr[coords[0],coords[1],coords[2]] = r_min

                                sphere_remove.append(sph_save)                                                                                  # Analysis complete, remove from future distance calculations
                                index = i; sph_save = sph
                    sphere_temp = np.delete(sphere_temp, np.array(sphere_remove), axis=0)                                                       # Remove evaluated voxel-centers
    if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1): time_Spheres = time.perf_counter() - time_Spheres
    max_radius = np.max(radii_arr); max_diameter = 2*max_radius
    del sys_mask; del sphere_temp; del pair_arr; del dist_arr; del sys; del sys_radii

    # Useful print command for troubleshooting problems: prints the number of voxel-centers within the free volume and the diameter of the largest sphere (pore)
    # Also prints the number of voxels within the system van der Waals free volume, voxels containing free volume spheres of radius r >= probe_radius, and voxels that need to be assessed whether they are in the free volume or not
    if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
        print(f"\nMaximum pore diameter: {max_diameter:.2f}")
        print(f"Number of free volume spheres (r >= probe_radius): {len(radii_arr[radii_arr >= args.probe_radius])}")
        print(f"Number of free volume voxels (r > 0): {len(radii_arr[radii_arr != 0])}")
        print(f"Time free volume spheres: {time_Spheres:.2f} s")





    # Clustering analysis
    # Only consider free volume spheres of radius r >= probe_radius that are within the desired domain.
    #   Free volume spheres outside of the desired domain are demoted from free volume spheres of radius r >= probe_radius to free volume voxels (0 < r < probe_radius)
    #       Free volume voxels are still considered in the FFV and PSD analysis.
    # NOTE: This section is the most sensitive to memory errors. Consider setting solvent_name = "" if consistently running out of memory (OOM).
    if (args.solvent_name == 'percolated') or (N_sol > 0):
        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
            time_Cluster = time.perf_counter()
            print('\n##### Performing Clustering Analysis - Percolated/Solvent-Domain #####')

        # Create an interconnected graph lattice of the voxelized system, where voxels are associated to each other through their 6 3x3x3 cube-face-center neighbors
        if args.clustering == 'Neumann':
            # Only calculate out half of the neighbors per voxel to prevent double-counting
            Neighborhood = np.array([[1, 0, 0],                                                                                                 #Neighborhood = np.array([[1, 0, 0], [-1, 0, 0],
                                     [0, 1, 0],                                                                                                 #                         [0, 1, 0], [ 0,-1, 0],
                                     [0, 0, 1]], dtype=linear_type)                                                                             #                         [0, 0, 1], [ 0, 0,-1]], dtype=linear_type)
        # Create an interconnected graph lattice of the voxelized system, where voxels are associated to each other through their 26 3x3x3 cube neighbors
        elif args.clustering == 'Moore':
            # Only calculate out half of the neighbors per voxel to prevent double-counting
            Neighborhood = np.array([[ 0,  0,  1],                                                                                              #Neighborhood = np.array([[-1, -1, -1], [-1, -1,  0], [-1, -1,  1],
                                     [ 0,  1, -1], [ 0,  1,  0], [ 0,  1,  1],                                                                  #                         [-1,  0, -1], [-1,  0,  0], [-1,  0,  1],
                                     [ 1, -1, -1], [ 1, -1,  0], [ 1, -1,  1],                                                                  #                         [-1,  1, -1], [-1,  1,  0], [-1,  1,  1],
                                     [ 1,  0, -1], [ 1,  0,  0], [ 1,  0,  1],                                                                  #                         [ 0, -1, -1], [ 0, -1,  0], [ 0, -1,  1],
                                     [ 1,  1, -1], [ 1,  1,  0], [ 1,  1,  1]], dtype=signed_linear_type)                                       #                         [ 0,  0, -1], [ 0,  0,  1], # No [ 0,  0,  0]
                                                                                                                                                #                         [ 0,  1, -1], [ 0,  1,  0], [ 0,  1,  1],
                                                                                                                                                #                         [ 1, -1, -1], [ 1, -1,  0], [ 1, -1,  1],
                                                                                                                                                #                         [ 1,  0, -1], [ 1,  0,  0], [ 1,  0,  1],
                                                                                                                                                #                         [ 1,  1, -1], [ 1,  1,  0], [ 1,  1,  1]], dtype=linear_type)
        else:
            raise ValueError(f'Clustering variable incorrectly set: {args.clustering}')
        
        # Create graph for cluster analysis, where each voxel is indexed sequentially, not in spatial (x,y,z) indices
        G = Graph(l_x * l_y * l_z, directed=False)

        # Retrieve index of all free volume spheres of radius r >= probe_radius and linearize their indices for use in the cluster graph analysis
        radii_arr = radii_arr.ravel()
        Graph_radii = radii_arr >= args.probe_radius                                                                                            # Linearized free volume sphere radii

        # For efficiency, we limit the number of edges generated per loop
        count = 0
        while count < len(radii_arr):
            count_old = count; count += min(int(args.N_edge_gen), len(radii_arr)-count_old)

            Graph_idx = np.where(Graph_radii[count_old:count])[0].astype(linear_type)                                                           # Linearized indices of free volume spheres
            idx_x, idx_y, idx_z = np.unravel_index(count_old + Graph_idx, (l_x, l_y, l_z))
            idx_x = idx_x.astype(indexed_type); idx_y = idx_y.astype(indexed_type); idx_z = idx_z.astype(indexed_type)

            # For memory efficiency, loop through each neighbor one-by-one
            for neigh in Neighborhood:
                # Retrieve linearized index of each voxel 'paired' to each free volume sphere voxel in Graph_idx
                edge_Graph_idx = (  ((idx_x + neigh[0]) % l_x) * l_y * l_z
                                  + ((idx_y + neigh[1]) % l_y) * l_z
                                  + ((idx_z + neigh[2]) % l_z)            )
                
                if edge_Graph_idx.dtype != Graph_idx.dtype:
                    edge_Graph_idx = edge_Graph_idx.astype(linear_type)

                # Add 'pairs' (edges) to graph G, where we only consider voxel edges between free volume spheres
                G.add_edges(np.stack((
                                           Graph_idx[Graph_radii[edge_Graph_idx]],
                                      edge_Graph_idx[Graph_radii[edge_Graph_idx]]
                                                                                 ), axis=1, dtype=linear_type))

        clusters = G.components()
        membership = np.array(clusters.membership, dtype=linear_type)
        cluster_ids = np.argsort(clusters.sizes())[::-1]
        G.clear(); del G; del clusters
        
        # Loop through clusters of free volume spheres in sorted order largest to smallest
        for i,id in enumerate(cluster_ids):
            # Only analyze the largest free volume sphere (radius r >= probe_radius) cluster, i.e., i = 0
            if args.solvent_name == 'percolated':
                if i == 0:
                    # Remove all free volume voxels not within the largest cluster, i.e., id(i == 0)
                    radii_arr[(membership != id) & (radii_arr >= args.probe_radius)] = args.probe_radius/2                                      # Set radii = probe_radius/2 so that these voxels are treated as free volume VOXELS and not free volume SPHERES going forward
                    break
            # Only analyze free volume sphere (radius r >= probe_radius) clusters within the solvent domain
            else:
                clust = np.where(membership == id)[0]                                                                                           # Linearized indices

                # Clusters containing 1 free volume sphere are assumed to NOT contain solvent
                #  - Significantly reduces compute time
                if len(clust) == 1:
                    radii_arr[np.isin(membership, cluster_ids[i:]) & (radii_arr >= args.probe_radius)] = args.probe_radius/2                    # Set radii = probe_radius/2 so that these voxels are treated as free volume VOXELS and not free volume SPHERES going forward
                    break

                # For efficiency, we limit the number of free volume spheres per loop to a total of N_calc_PSD distance calculations
                count = 0
                while count < len(clust):
                    count_old = count; count += min(int(args.N_calc_PSD/N_sol), len(clust)-count_old)
        
                    # Useful print command for troubleshooting memory problems
                    # Decreasing N_calc_PSD will reduce memory usage
                    if (args.print_eff == 2) and (frame == frame_ids[-1] or args.N_threads == 1) and (i == 0 or len(clust)*N_sol > args.N_calc_PSD/10):
                        if count_old == 0:
                            print(f"\nCluster size: {len(clust)}")
                        print(f"Calculations: {(count - count_old)*N_sol:.1e}")
                    
                    idx_x, idx_y, idx_z = np.unravel_index(clust[count_old:count], (l_x, l_y, l_z))                                             # Spatial indices
        
                    # Find the number of solvent atoms within probe_radius of a free volume sphere center
                    pair_arr = distances.capped_distance(
                                                         np.stack((
                                                                   vox_x[idx_x],
                                                                   vox_y[idx_y],
                                                                   vox_z[idx_z]
                                                                               ), axis=1, dtype=float_type), sol, args.probe_radius, box=cell, return_distances=False)
        
                    # If any solvent molecules are found within the free volume sphere cluster, analysis can end early
                    if len(pair_arr) != 0:
                        break
        
                # If any solvent atoms are within probe_radius of a free volume sphere, then the entire cluster is considered a part of the solvent domain
                if len(pair_arr) != 0:
                    continue
        
                # All other clusters are removed from the free volume sphere analysis.
                radii_arr[clust] = args.probe_radius/2                                                                                          # Set radii = probe_radius/2 so that these voxels are treated as free volume VOXELS and not free volume SPHERES going forward
        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1): time_Cluster = time.perf_counter() - time_Cluster
        radii_arr = radii_arr.reshape((l_x, l_y, l_z)); max_radius = np.max(radii_arr); max_diameter = 2*max_radius
        del sol; del membership; del cluster_ids

        # Useful print command for troubleshooting problems
        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
            print(f"\nMaximum pore diameter: {max_diameter:.2f}")
            if args.solvent_name == 'percolated':
                print(f"Number of free volume spheres (r >= probe_radius) within the percolated domain: {len(radii_arr[radii_arr >= args.probe_radius])}")
            else:
                print(f"Number of free volume spheres (r >= probe_radius) within the solvent domain: {len(radii_arr[radii_arr >= args.probe_radius])}")
            print(f"Time cluster: {time_Cluster:.2f} s\n")

    # Code to write coordinates and radius of each free volume sphere to a .xyz file, which can be visualized in OVITO, etc
    #     Radius = radius of the largest free volume sphere centered on the voxel-center
    if (args.print_xyz) and (frame == frame_ids[-1]):
        idx_x, idx_y, idx_z = np.where(radii_arr >= args.probe_radius)
        with open('Free_Volume_Spheres.xyz', 'w') as anaout:
            print(str(len(idx_x)), file=anaout)
            print('Properties=species:S:1:pos:R:3:Radius:R:1', file=anaout)
            for i in range(len(idx_x)):
                x = vox_x[idx_x[i]]; y = vox_y[idx_y[i]]; z = vox_z[idx_z[i]]; r = radii_arr[idx_x[i],idx_y[i],idx_z[i]]
                if args.mode == 'xyz':
                    print('X {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(x, y, z, r), file=anaout)
                else:
                    print('X {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(x - cell[0]/2, y - cell[1]/2, z - cell[2]/2, r), file=anaout)
        del idx_x; del idx_y; del idx_z
        print('Free volume sphere xyz file printed')





    # PSD/FFV analysis
    # This part of the calculation determines the free volume fraction and cumulative probe-occupiable pore size distribution, where the distribution is defined as the probability that a random point (voxel) within the free volume resides within a free volume sphere of diameter d with minimum size probe_radius
    # This code will take each voxel not within the system volume (PSD_probes) and determine 1) if it lies within the free volume (FFV), and 2) the largest free volume sphere it lies within (PSD)
    # Changing L_voxel, N_calc_PSD, and d_step can reduce run time and memory usage
    if args.PSD_FFV:
        PSD_probes = np.indices((l_x, l_y, l_z), dtype=indexed_type).reshape(3, -1).T                                                               # Indices of all

        FFV_c = 0; FFV_lr = 0; FFV_total = 0                                                                                                        # Track number of voxels within the Connolly and Lee-Richards free volume against the total number to get FFV
        d_arr = np.arange(0, args.d_max + args.d_step, args.d_step); PSD_arr = np.zeros_like(d_arr, dtype=int)                                      # d_arr is the histogram of free volume sphere sizes; PSD_arr tracks the number of instances of voxels contained within free volume spheres of size at least d

        FFV_save = np.array([[],[],[]], dtype=indexed_type).T; d_save = np.array([], dtype=indexed_type)                                            # Save voxel-centers within the free volume for surface area calculations, and the size of the largest free volume sphere containing each voxel-center for printing in Free_Volume_Voxels.xyz

        # Make sure 2*max_radius > d_max
        if max_diameter > args.d_max:
            raise ValueError(f"Largest pore diameter is greater than d_max: {max_diameter} > {args.d_max}")

        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
            time_PSD = time.perf_counter()
            print('\n##### Performing PSD/FFV Analysis #####\n')

        # Starting from the largest free volume spheres, find all free volume voxels within the desired free volume domain for FFV and PSD calulcations
        cycle = 0; err = np.inf; N_rand = int(np.ceil((l_x * l_y * l_z) * args.rand_frac)); PSD_Old = np.zeros_like(PSD_arr)
        while err > args.tol and len(PSD_probes) != 0:
            if N_rand > len(PSD_probes):
                N_rand = len(PSD_probes)

            cycle += 1
            # Track number of cycles
            if (args.print_eff == 2) and (frame == frame_ids[-1] or args.N_threads == 1):
                print(f"PSD Cycle: {cycle:5d}/{int(np.ceil(1/args.rand_frac))}")

            if args.rand_frac == 1.0:
                PSD_temp = PSD_probes[:]
                PSD_probes = np.array([])
            else:
                Rand_idx = np.random.choice(len(PSD_probes), size=N_rand, replace=False)
                PSD_temp = PSD_probes[Rand_idx]
                PSD_probes = np.delete(PSD_probes, Rand_idx, axis=0); del Rand_idx
            FFV_total += N_rand; FFV_lr += np.sum(radii_arr[PSD_temp[:,0],PSD_temp[:,1],PSD_temp[:,2]] >= args.probe_radius)
            PSD_temp = PSD_temp[radii_arr[
                                          PSD_temp[:,0],
                                          PSD_temp[:,1],
                                          PSD_temp[:,2]
                                                       ] != 0]

            for d in np.round(np.arange(args.d_max, 0, -args.d_step), decimals = 5):
            
                if d - args.d_step > max_diameter:
                    continue

                if (d < 2*args.probe_radius) or (len(PSD_temp) == 0):
                    break

                # For efficiency, we measure the distance between free volume spheres and the voxel-centers starting with the largest d_arr bin and moving down
                if (d - args.d_step)/2 < args.probe_radius:
                    idx_x, idx_y, idx_z = np.where((radii_arr <= d/2) & (radii_arr >= args.probe_radius))
                else:
                    idx_x, idx_y, idx_z = np.where((radii_arr <= d/2) & (radii_arr > (d - args.d_step)/2))
                sphere_temp = np.stack((vox_x[idx_x],vox_y[idx_y],vox_z[idx_z]), axis=1, dtype=float_type); radii_temp = radii_arr[idx_x, idx_y, idx_z] # Positions (sphere_temp) and radii (radii_temp) of free volume spheres in the current PSD bin, radius (d - d_step)/2 < r <= d/2
                if len(sphere_temp) == 0:
                    continue

                # For efficiency, we limit the number of free volume spheres per loop to a total of N_calc_PSD distance calculations
                count = 0
                while count < len(sphere_temp) and len(PSD_temp) > 0:
                    count_old = count; count += min(int(args.N_calc_PSD/len(PSD_temp)), len(sphere_temp)-count_old)

                    sph_temp = sphere_temp[count_old:count]; rad_temp = radii_temp[count_old:count]

                    pair_arr, dist_arr = distances.capped_distance(sph_temp, np.stack((                                                             # Distance between free volume spheres and voxel-centers
                                                                                       vox_x[PSD_temp[:,0]],
                                                                                       vox_y[PSD_temp[:,1]], 
                                                                                       vox_z[PSD_temp[:,2]]
                                                                                                             ), axis=1, dtype=float_type), d/2 + 0.5, box=cell)

                    # Useful print command for troubleshooting memory problems
                    # Decreasing N_calc_PSD will reduce memory usage
                    if (args.print_eff == 2) and (frame == frame_ids[-1] or args.N_threads == 1) and (len(sph_temp)*len(PSD_temp) > args.N_calc_PSD/10):
                        if count_old == 0:
                            print(f"\nDiameter: {d - args.d_step} < d <= {d}")
                        print(f"Calculations, writes: {len(sph_temp)*len(PSD_temp):.1e} {len(dist_arr):.1e}")
                    del sph_temp

                    if len(dist_arr) > 0:
                        dist_arr -= rad_temp[pair_arr[:,0]]                                                                                         # Subtract radius of each free volume sphere from the distance to get the distance from the voxel-center to the surface of the free volume sphere
                        pair_arr = np.unique(pair_arr[:,1][dist_arr < 0])                                                                           # Only consider voxel-centers that lie within the free volume sphere (adjusted distance < 0), and only count each occurence once (unique)

                        FFV_c += len(pair_arr); PSD_arr[np.where(d_arr < d)[0]] += len(pair_arr)                                                    # Voxel-centers w/n free volume sphere count towards the FFV and cumulatively towards the PSD

                        FFV_save = np.append(FFV_save, PSD_temp[pair_arr], axis=0)                                                                  # Save free volume voxel-centers for printing
                        d_save = np.append(d_save, np.zeros((len(pair_arr)), dtype=int) + int(d/args.d_step))

                        PSD_temp = np.delete(PSD_temp, pair_arr, axis=0)                                                                            # No longer consider voxel-centers that are found within a free volume sphere in future loops (prevent double-counting)

            PSD = PSD_arr / PSD_arr[0]; PSD = -(PSD[1:] - PSD[:-1]) / args.d_step
            if np.all(PSD_Old == 0):
                PSD_Old = PSD
                continue
            err = np.max(np.abs(np.divide((PSD - PSD_Old), PSD, out=np.zeros_like(PSD), where=(PSD != 0)))); PSD_Old = PSD

            if (args.print_eff == 2) and (frame == frame_ids[-1] or args.N_threads == 1):
                if args.tol > 0:
                    print(f"Maximum Error/Tolerance: {err:.1e}/{args.tol:.1e}\n")
                else:
                    print(f"Maximum Error: {err:.1e}\n")
        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1): time_PSD = time.perf_counter() - time_PSD
        del PSD_probes; del sphere_temp; del pair_arr; del dist_arr; del idx_x; del idx_y; del idx_z

        # Code to print the final FFV and PSD for the last frame analyzed
        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
            print(f"Connolly FFV: {FFV_c/FFV_total:0.3f}, {FFV_c}, {FFV_total}")
            print(f"Lee-Richards FFV: {FFV_lr/FFV_total:0.3f}, {FFV_lr}, {FFV_total}")
            print(f"\nPSD Final: {PSD_arr[0]}")
            print_string=''
            for i in PSD_arr:
                if i == 0:
                    continue
                print_string += str(np.round(i / PSD_arr[0], decimals=5)) + ' '
            print(print_string)
            print(f"Time PSD/FFV: {time_PSD:.2f} s")

        # Code to write coordinates of each voxel-center to a .xyz file, which can be visualized in OVITO
        if (args.print_xyz) and (frame == frame_ids[-1]):
            with open('Free_Volume_Voxels.xyz', 'w') as anaout:
                print(str(len(FFV_save)), file=anaout)
                print('Properties=species:S:1:pos:R:3:Radius:R:1:Alpha:R:1', file=anaout)
                for i, sph in enumerate(FFV_save):
                    x = vox_x[sph[0]]; y = vox_y[sph[1]]; z = vox_z[sph[2]]; r = args.L_voxel/2; a = d_arr[d_save[i]]
                    if args.mode == 'xyz':
                        print('X {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(x, y, z, r, a), file=anaout)
                    else:
                        print('X {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(x - cell[0]/2, y - cell[1]/2, z - cell[2]/2, r, a), file=anaout)
            print('Free volume voxel xyz file printed')
        del d_save; del vox_x; del vox_y; del vox_z
    else:
        FFV_c = -len(radii_arr.ravel()); FFV_lr = len(radii_arr[radii_arr >= args.probe_radius]); FFV_total = len(radii_arr.ravel())
        d_arr = np.arange(0, args.d_max + args.d_step, args.d_step); PSD_arr = np.zeros_like(d_arr, dtype=int) - 1
        
        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
            print(f"Lee-Richards FFV: {FFV_lr/FFV_total:0.3f}, {FFV_lr}, {FFV_total}")





    # SA analysis
    if args.Surface_area == True:
        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
            time_SA = time.perf_counter()
            print('\n##### Performing SA Analysis #####\n')

        if args.PSD_FFV:
            ######################################################
            ############### Connolly Surface Area ################
            ######################################################

            # Surface is defined by the free volume voxels
            SA_arr = np.zeros((l_x, l_y, l_z), dtype=bool); SA_arr[FFV_save[:,0], FFV_save[:,1], FFV_save[:,2]] = True                              # Create voxel lattice where free volume voxel-centers = True

            # Create a simple mesh surface around the free volume and calculate the surface area
            SA_arr = np.pad(SA_arr, pad_width = 1, mode = 'wrap')                                                                                   # Add 1 layer of wrapped coordinates around the array to properly account for periodic boundaries
            spacing = np.array([L_voxel_x, L_voxel_y, L_voxel_z])                                                                                   # Define voxel size to dimensionalize surface area calulcations

            verts_c, faces_c, _, _ = measure.marching_cubes(SA_arr, level = 0.5, spacing = spacing)                                                 # Marching cubes algorithm to create a surface mesh
            SA_c = measure.mesh_surface_area(verts_c, faces_c)                                                                                      # Calculate the surface area of the free volume
        else:
            SA_c = -1

        ######################################################
        ### Lee-Richards "Surface Accessible" Surface Area ###
        ######################################################

        # Surface defined by the *center* of the free volume *spheres* - i.e., surface-accessible free volume
        idx_x, idx_y, idx_z = np.where(radii_arr >= args.probe_radius)
        SA_arr = np.zeros((l_x, l_y, l_z), dtype=bool); SA_arr[idx_x, idx_y, idx_z] = True                                                      # Create voxel lattice where free volume sphere-centers = True

        # Create a simple mesh surface around the free volume and calculate the surface area
        SA_arr = np.pad(SA_arr, pad_width = 1, mode = 'wrap')                                                                                   # Add 1 layer of wrapped coordinates around the array to properly account for periodic boundaries
        spacing = np.array([L_voxel_x, L_voxel_y, L_voxel_z])                                                                                   # Define voxel size to dimensionalize surface area calulcations

        verts_lr, faces_lr, _, _ = measure.marching_cubes(SA_arr, level = 0.5, spacing = spacing)                                               # Marching cubes algorithm to create a surface mesh
        SA_lr = measure.mesh_surface_area(verts_lr, faces_lr)                                                                                   # Calculate the surface area of the free volume

        volume = (l_x * L_voxel_x) * (l_y * L_voxel_y) * (l_z * L_voxel_z)                                                                      # Normalize surface area to the true volume
        padded_volume = ((l_x + 2) * L_voxel_x) * ((l_y + 2) * L_voxel_y) * ((l_z + 2) * L_voxel_z)
        if args.PSD_FFV: SA_c *= (volume / padded_volume)
        SA_lr *= (volume / padded_volume)

        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1): time_SA = time.perf_counter() - time_SA

        # Code to print the surface area for the last frame analyzed
        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
            print(f"Connolly SA (A^2):  {SA_c:.2f}")
            print(f"Lee-Richards SA (A^2):  {SA_lr:.2f}")
            print(f"Time SA: {time_SA:.2f} s")

        # Code to write coordinates of each voxel-center to a .xyz file, which can be visualized in OVITO
        #     Radius = L_voxel/2
        if (args.print_xyz) and (frame == frame_ids[-1]):
            # Remove excess voxels due to padding
            verts_c_save = []
            if args.PSD_FFV:
                for i, sph in enumerate(verts_c):
                    if np.any(sph < 0 - args.L_voxel/2) or np.any(sph > cell[:3] + args.L_voxel/2):
                        continue
                    verts_c_save.append(sph)
            
            verts_lr_save = []
            for i, sph in enumerate(verts_lr):
                if np.any(sph < 0 - args.L_voxel/2) or np.any(sph > cell[:3] + args.L_voxel/2):
                    continue
                verts_lr_save.append(sph)

            with open('Free_Volume_Surface.xyz', 'w') as anaout:
                print(str(len(verts_c_save) + len(verts_lr_save)), file=anaout)
                print('Properties=species:S:1:pos:R:3:Radius:R:1', file=anaout)
                for i, sph in enumerate(verts_c_save):
                    if args.mode == 'xyz':
                        print('X {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(sph[0], sph[1], sph[2], args.L_voxel/2), file=anaout)
                    else:
                        print('X {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(sph[0] - cell[0]/2, sph[1] - cell[1]/2, sph[2] - cell[2]/2, args.L_voxel/2), file=anaout)
                for i, sph in enumerate(verts_lr_save):
                    if args.mode == 'xyz':
                        print('Y {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(sph[0], sph[1], sph[2], args.L_voxel/2), file=anaout)
                    else:
                        print('Y {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(sph[0] - cell[0]/2, sph[1] - cell[1]/2, sph[2] - cell[2]/2, args.L_voxel/2), file=anaout)
            print('Free volume surface xyz file printed')
    else:
        SA_c = 0; SA_lr = 0
    




    # Tortuosity analysis
    if args.Tortuosity == True:
        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
            time_tau = time.perf_counter()
            print('\n##### Performing Tortuosity Analysis #####\n')

        # Diffusive volume is defined by *probe-center* occupiable volume, i.e., the Lee-Richards volume
        idx_x, idx_y, idx_z = np.where(radii_arr >= args.probe_radius)
        tortuosity_arr = np.zeros((l_x, l_y, l_z), dtype=bool); tortuosity_arr[idx_x, idx_y, idx_z] = True                                      # Create voxel lattice where free volume sphere-centers = True

        try:                                                                                                                                    # Attempt tortuosity analysis across x, y, and z directions
            sim_x = ps.simulations.tortuosity_fd(tortuosity_arr, axis=0); tortuosity_x = sim_x.tortuosity                                       # Analysis fails if no percolating clusters found across that axis
            sim_y = ps.simulations.tortuosity_fd(tortuosity_arr, axis=1); tortuosity_y = sim_y.tortuosity
            sim_z = ps.simulations.tortuosity_fd(tortuosity_arr, axis=2); tortuosity_z = sim_z.tortuosity
        except Exception as e:
            if "No pores remain" in str(e):                                                                                                     # If no percolating cluster found across any axis, return -1 for failed analysis
                if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
                    print("Warning: Void space does not percolate along at least one axis. Setting tortuosity to -1.")
                tortuosity_x = -1; tortuosity_y = -1; tortuosity_z = -1
            elif "Solver failed to converge" in str(e):                                                                                         # If solver failed to converge across any axis, return -1 for failed analysis
                if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
                    print("Error: Solver failed to converge along at least one axis. Setting tortuosity to -1.")
                tortuosity_x = -1; tortuosity_y = -1; tortuosity_z = -1
            else:
                raise e
    
        tortuosity = np.mean([tortuosity_x, tortuosity_y, tortuosity_z])                                                                        # Average tortuosity across all 3 dimensions

        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1): time_tau = time.perf_counter() - time_tau

        # Code to print the surface area for the last frame analyzed
        if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
            if tortuosity == -1:
                print(f"No 1D percolated clusters found, tortuosity not measured.")
            else:
                print(f"Directional Tortuosity:  X-{tortuosity_x:.2f} Y-{tortuosity_y:.2f} Z-{tortuosity_z:.2f}")
                print(f"Average Tortuosity:  {tortuosity:.2f}")
            print(f"Time Toruosity: {time_tau:.2f} s")
    else:
        tortuosity_x = 0; tortuosity_y = 0; tortuosity_z = 0
    




    if (args.print_eff >= 1) and (frame == frame_ids[-1] or args.N_threads == 1):
        print("\n##### Summary of Calculation Times #####\n")
        print(f"Time free volume spheres: {time_Spheres:.2f} s")
        if (args.solvent_name == 'percolated') or (N_sol > 0):
            print(f"Time cluster: {time_Cluster:.2f} s")
        if args.PSD_FFV:
            print(f"Time PSD/FFV: {time_PSD:.2f} s")
        if args.Surface_area == True:
            print(f"Time SA: {time_SA:.2f} s")
        if args.Tortuosity == True:
            print(f"Time Toruosity: {time_tau:.2f} s")
    




    # Return the necessary information to complete the calculations: tortuosity_i/1000 gives the tortuosity in the ith direction, SA/100 gives the surface area, FFV_XX / FFV_total gives the probe-occupiable free volume, PSD_arr / PSD_arr[0] gives the probe-occupiable PSD
    PSD_arr = np.insert(PSD_arr, 0, FFV_total); PSD_arr = np.insert(PSD_arr, 0, FFV_lr); PSD_arr = np.insert(PSD_arr, 0, FFV_c)
    PSD_arr = np.insert(PSD_arr, 0, int(SA_lr*100)); PSD_arr = np.insert(PSD_arr, 0, int(SA_c*100))
    PSD_arr = np.insert(PSD_arr, 0, int(tortuosity_z*1000)); PSD_arr = np.insert(PSD_arr, 0, int(tortuosity_y*1000)); PSD_arr = np.insert(PSD_arr, 0, int(tortuosity_x*1000))
    return PSD_arr
    


def load_trajectory():
# loads in the trajectory and saves the necessary data to a temporary h5py .hdf5 I/O file

    if args.mode == 'xyz':
        uta = mda.Universe(args.trj_file)                        # Load in the xyz trajectory

        # Define the simulation cell from the dat file
        cell = np.zeros(6); cell[3:] = 90.0
        with open(args.top_file, 'r') as file:
            lines = file.readlines()[1]
            cell[:3] += np.array(lines.split(), dtype=float)
        print('If the following is incorrect, check your dat file format')
        print(f'XYZ cell size (A): {cell[0]} {cell[1]} {cell[2]}\n')
    else:
        uta = mda.Universe(args.top_file, args.trj_file)         # Load in the trajectory and topology

    system = uta.select_atoms(args.system_name)                  # Define the system atoms
    if args.solvent_name == 'percolated' or args.solvent_name == '':
        solvent = uta.select_atoms('not all')                    # Solvent selection is empty
    else:
        solvent = uta.select_atoms(args.solvent_name)            # Define the solvent atoms

    print("If the following is incorrect, there may be inconsistencies between your atom ID name in the topology and the Element name in the YAML file (see 'Size_arr' and 'Dummy_atoms' for more details)")
    print("\nSYSTEM ATOMS")

    # If no system atoms are detected, return error
    if len(system) == 0:
        raise ValueError("No system atoms found")

    # Remove dummy atoms from the system
    if len(Dummy_atoms) > 0 and len(system) > 0:
        for dummy in Dummy_atoms:
            if np.sum(system.names == dummy) > 0:
                print(f"Removed {np.sum(system.names == dummy)} {dummy} atoms from system analysis")
                system = system[system.names != dummy]
    
    # Create an array that tracks the radius of each system atom based on Size_array
    sys_names = system.names; sys_radii = np.zeros((len(system)), dtype=float_type); sys_count = np.zeros((len(Size_arr)), dtype=int)
    for i, name in enumerate(sys_names):
        name = str(name)
        if name in Size_arr[:,0]:
            sys_radii[i] = float(Size_arr[np.where(Size_arr[:,0] == name)[0][0],1])
            sys_count[np.where(Size_arr[:,0] == name)[0][0]] += 1
        elif name[0] in Size_arr[:,0]:
            sys_radii[i] = float(Size_arr[np.where(Size_arr[:,0] == name[0])[0][0],1])
            sys_count[np.where(Size_arr[:,0] == name[0])[0][0]] += 1
        else:
            raise ValueError(f"Missing Atom Name and Size in Size_arr: {name}")

    # Print out system atom information
    print("Element N-in-System")
    for i,j in enumerate(sys_count):
        if j > 0: print(f"{Size_arr[i,0]:>7s} {j:11d}")
    
    if len(solvent) > 0:
        print("\nSOLVENT ATOMS")

        # Remove dummy atoms from the solvent
        if len(Dummy_atoms) > 0 and len(solvent) > 0:
            for dummy in Dummy_atoms:
                if np.sum(solvent.names == dummy) > 0:
                    print(f"Removed {np.sum(solvent.names == dummy)} {dummy} atoms from solvent analysis")
                    solvent = solvent[solvent.names != dummy]
        
        # Create an array that tracks the radius of each system atom based on Size_array
        sol_names = solvent.names; sol_count = np.zeros((len(Size_arr) + 1), dtype=int)
        for name in sol_names:
            name = str(name)
            if name in Size_arr[:,0]:
                sol_count[np.where(Size_arr[:,0] == name)[0][0]] += 1
            elif name[0] in Size_arr[:,0]:
                sol_count[np.where(Size_arr[:,0] == name[0])[0][0]] += 1
            else:
                sol_count[-1] += 1

        # Print out solvent atom information
        print("Element N-in-System")
        for i,j in enumerate(sol_count):
            if j > 0:
                if i == len(sol_count) - 1:
                    print(f"{"Other":>7s} {j:11d}")
                    print("\nElement 'Other' means the atom is not defined as part of the van der Waals volume list.")
                    print("This is not an error and does not impact the code output.")
                else:
                    print(f"{Size_arr[i,0]:>7s} {j:11d}")

    # Define the system times/frames to be calculated over
    print()
    if args.mode == 'xyz' or '.gro' in args.trj_file:
        frame_ids = np.array([0], dtype=int)                     # .xyz and .gro trajectories are assumed to contained only a single frame
    else:
        if args.t_min == -1:    args.t_min    = uta.trajectory[0].time
        if args.t_max == -1:    args.t_max    = uta.trajectory[-1].time
        if args.N_frames == -1: args.N_frames = args.N_threads

        dt = np.round((uta.trajectory[1].time - uta.trajectory[0].time),3)

        if args.N_frames == 1:                                   # If only analyzing one frame, analyze the final frame
            frame_ids = np.array([int((args.t_max - uta.trajectory[0].time)/dt)], dtype=int)
        else:
            frame_ids = np.linspace(int((args.t_min - uta.trajectory[0].time)/dt), int((args.t_max - uta.trajectory[0].time)/dt), args.N_frames, dtype=int)
            print(f"Timestep: ~{dt*(frame_ids[1] - frame_ids[0])} ps")
    print(f"Number of frames: {len(frame_ids)}")

    # Load in the necessary data: "system" atom positions, "solvent" atom positions, cell dimensions
    r_system = np.zeros((len(frame_ids), len(system), 3), dtype=float_type)
    r_solvent = np.zeros((len(frame_ids), len(solvent), 3), dtype=float_type)
    cells = np.zeros((len(frame_ids), 6))
    for i, frame in enumerate(frame_ids):
        ts = uta.trajectory[frame]

        r_system[i] = system.positions
        r_solvent[i] = solvent.positions

        if args.mode == 'xyz':
            cells[i] = cell
        else:
            cell = ts.dimensions
            cells[i] = cell

    # Save necessary information to a temporary .hdf5 file for later use in the calculation
    with h5py.File('PrO-VAT.hdf5','w') as f:
        dset1 = f.create_dataset("system", data=r_system, dtype=float_type)
        dset2 = f.create_dataset("sys_radii", data = sys_radii, dtype=float_type)
        dset3 = f.create_dataset("solvent", data=r_solvent, dtype=float_type)
        dset4 = f.create_dataset("cells", data = cells)
        dset5 = f.create_dataset("frames", data = frame_ids)



def main():
    # Load in the trajectory file and exit the code to purge the memory before multiprocessing
    if not os.path.exists('PrO-VAT.hdf5'):
        print('Loading trajectory data\n')
        try:
            load_trajectory()
        except ValueError as e:
            print(f"ERROR - {e}")
        print('\nTrajectory loaded, terminating process. Run again to perform analysis')
        exit()

    with h5py.File('PrO-VAT.hdf5','r') as f:
        dset1 = f['frames']; frame_ids = dset1[:]
    frame_ids = np.arange(0,len(frame_ids),1)

    print("PSD Analysis\n")
    # Perform the analysis using multiprocessing
    try:
        pool = mp.Pool(processes=args.N_threads)
        func = functools.partial(volume_analysis)
        out_arr = pool.map(func, list(frame_ids))
        pool.close()
        pool.join()
        out_arr = np.array(out_arr)
    except ValueError as e:
        print(f"ERROR - {e}")
    
    # Return the average and standard deviation (over the frames processed) of the tortuosity
    tortuosity_x = out_arr[:,0]/1000; tortuosity_y = out_arr[:,1]/1000; tortuosity_z = out_arr[:,2]/1000
    if args.Tortuosity:
        if np.any(tortuosity_x == -1):
            tortuosity = np.array([-1, -1, -1, -1, -1, -1])
        else:
            tortuosity = np.array([np.mean(tortuosity_x), np.std(tortuosity_x), np.mean(tortuosity_y), np.std(tortuosity_y), np.mean(tortuosity_z), np.std(tortuosity_z)])

        with open('Tau.dat', 'w') as anaout:
            print("# Tortuosity Std - 0.0, 1.0, 2.0 = X, Y, and Z direction - value of -1 denotes a failed tortuosity analysis on at least 1 frame", file=anaout)
            print('0.0 {:10.5f} {:10.5f}'.format(tortuosity[0], tortuosity[1]), file=anaout)
            print('1.0 {:10.5f} {:10.5f}'.format(tortuosity[2], tortuosity[3]), file=anaout)
            print('2.0 {:10.5f} {:10.5f}'.format(tortuosity[4], tortuosity[5]), file=anaout)

    # Return the average and standard deviation (over the frames processed) of the surface area
    SA_c = out_arr[:,3]/100; SA_lr = out_arr[:,4]/100; SA = np.array([np.mean(SA_c), np.std(SA_c), np.mean(SA_lr), np.std(SA_lr)])
    if args.Surface_area:
        with open('SA.dat', 'w') as anaout:
            print("# SA (A^2) Std - 0.0 = Connolly, 1.0 = Lee-Richards", file=anaout)
            print('0.0 {:15.5f} {:10.5f}'.format(SA[0], SA[1]), file=anaout)
            print('1.0 {:15.5f} {:10.5f}'.format(SA[2], SA[3]), file=anaout)

    # Return the average and standard deviation (over the frames processed) of the probe-occupiable fractional free volume
    FFV = out_arr[:,5:8]; FFV_c = FFV[:,0] / FFV[:,2]; FFV_lr = FFV[:,1] / FFV[:,2]; FFV = np.array([np.mean(FFV_c), np.std(FFV_c),np.mean(FFV_lr), np.std(FFV_lr)])
    with open('FFV.dat', 'w') as anaout:
        print("# FFV Std - 0.0 = Connolly, 1.0 = Lee-Richards", file=anaout)
        print('0.0 {:10.5f} {:10.5f}'.format(FFV[0], FFV[1]), file=anaout)
        print('1.0 {:10.5f} {:10.5f}'.format(FFV[2], FFV[3]), file=anaout)

    if args.PSD_FFV:
        # Return the average and standard deviation (over the frames processed) of the probe-occupiable pore size ditribution
        d_arr = np.arange(0, args.d_max + args.d_step, args.d_step)
        PSD_all = out_arr[:,8:]; PSD_all = np.divide(PSD_all.T, PSD_all[:,0], dtype=float).T
        PSD_Cumulative = np.array([np.mean(PSD_all, axis=0), np.std(PSD_all, axis = 0)])
        # PSD is the negative derivative of the cumulative sum
        PSD = np.array([np.mean(-(PSD_all[:,1:] - PSD_all[:,:len(d_arr)-1])/(d_arr[1:] - d_arr[:len(d_arr)-1]), axis=0), np.std(-(PSD_all[:,1:] - PSD_all[:,:len(d_arr)-1])/(d_arr[1:] - d_arr[:len(d_arr)-1]), axis=0)])

        with open('Cumulative_PSD.dat', 'w') as anaout:
            print("# d (A) Cumulative_PSD Std", file=anaout)
            for i in range(len(PSD_Cumulative[0,:])):
                print(' {:10.5f} {:10.5f} {:10.5f}'.format(np.round(d_arr[i], decimals=3), PSD_Cumulative[0,i], PSD_Cumulative[1,i]), file=anaout)

        with open('PSD.dat', 'w') as anaout:
            print("# d (A) PSD Std", file=anaout) 
            for i in range(len(PSD[0,:])):
                    if i == 0:
                        print(' {:10.5f} {:10.5f} {:10.5f}'.format(np.round(d_arr[i], decimals=3), 0.0, 0.0), file=anaout)
                    else:
                        print(' {:10.5f} {:10.5f} {:10.5f}'.format(np.round(d_arr[i], decimals=3), PSD[0,i-1], PSD[1,i-1]), file=anaout)

    # Deletes the temporary .hdf5 file
    os.remove('PrO-VAT.hdf5')



def readable_file(path):
    # Check if a path exists and is a file.
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"The file '{path}' does not exist.")
    elif not os.access(path, os.R_OK):
        raise argparse.ArgumentTypeError(f"The file '{path}' is not readable.")
    return path

def string2bool(input):
    # Convert string input to boolean
    if input == 'True':
        return True
    elif input == 'False':
        return False



def loadArgs():
    # Define parser for YAML config file
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument('yaml_file', type = readable_file)

    # Add helpful error message if YAML file is not provided
    try:
        args, remaining_argv = config_parser.parse_known_args()
    except:
        print('\nERROR - Please provide a YAML config file, e.g., python3 /{Path}/PrO-VAT.py /{Path}/config.yaml')
        exit()

    # Load the YAML data
    with open(args.yaml_file, 'r') as f:
        config = yaml.safe_load(f)

    # Define parser for inputs + help menu
    parser = argparse.ArgumentParser(description="PrO-VAT: Probe-Occupiable Volume Analysis Tools")

    # Define subparsers
    subparsers = parser.add_subparsers(dest="mode", help = "Input file mode", required=True)

    #################################
    #################################
    ## SUBPARSER 1: For .xyz files ##
    #################################
    #################################
    xyz_parser = subparsers.add_parser('xyz', help = "Process PoreBlazer-style xyz + dat trajectory files")
       ###################################
       ## GROUP 1: Required input files ##
       ###################################
    xyz_files = xyz_parser.add_argument_group('Required input files')
    xyz_files.add_argument('trj_file', type = readable_file,
                           help = "Path to xyz file")
    xyz_files.add_argument('top_file', type = readable_file,
                           help = "Path to dat file")
       ##########################################
       ## GROUP 2: Frame selection and threads ##
       ##########################################
    xyz_frames = xyz_parser.add_argument_group('Frame selection and threads')
    xyz_frames.add_argument('-n', '--N_frames', type = int, default = 1, choices = [1],
                            help = "Number of frames to analyze [Locked to 1 frame for xyz analysis]")
    xyz_frames.add_argument('-t', '--N_threads', type = int, default = 1, choices = [1],
                            help = "Number of threads for parallelization [Locked to 1 thread for xyz analysis]")
       ###########################################
       ## GROUP 3: MDAnalysis selection strings ##
       ###########################################
    xyz_selection = xyz_parser.add_argument_group('MDAnalysis selection strings')
    xyz_selection.add_argument('-m', '--system_name', type = str, default = config['system_name'],
                               help = "MDAnalysis selection string defining the system matrix, e.g., 'all', 'moltype MOL', 'resname PEO', 'resname SOL LI CL' [default = YAML; Typically 'all' for PoreBlazer-style xyz + dat input.]")
    xyz_selection.add_argument('-s', '--solvent_name', type = str, default = config['solvent_name'],
                               help = "MDAnalysis selection string defining the solvent matrix, e.g., '', 'percolated', 'resname SOL LI CL' [default = YAML; Typically '' or 'percolated' for PoreBlazer-style xyz + dat input.]")
       ##################################
       ## GROUP 4: Important variables ##
       ##################################
    vars = xyz_parser.add_argument_group('Important variables')
    vars.add_argument('-L', '--L_voxel', type = float, default = config['L_voxel'],
                      help = "Voxel side length (A) [default = YAML]")
    vars.add_argument('-r', '--probe_radius', type = float, default = config['probe_radius'],
                      help = "Probe radius (A) [default = YAML]")
    vars.add_argument('--d_max', type = float, default = config['d_max'],
                      help = "Max PSD diameter (A) [default = YAML]")
    vars.add_argument('--d_step', type = float, default = config['d_step'],
                      help = "PSD bin size (A) [default = YAML]")
    vars.add_argument('--Voxel_dist', type = str, choices = ['Uniform', 'Random'], default = config['Voxel_dist'],
                      help = "Voxel distribution setting [default = YAML; Locked to 'Uniform' or 'Random']")
    vars.add_argument('--PSD_FFV', type = string2bool, choices = [True, False], default = config['PSD_FFV'],
                      help = "Pore size distribution and free volume fraction calculation setting [default = YAML; Locked to True or False]")
    vars.add_argument('--Surface_area', type = string2bool, choices = [True, False], default = config['Surface_area'],
                      help = "Surface area calculation setting; Requires --Voxel_dist 'Uniform' and --tol -1 [default = YAML; Locked to True or False]")
    vars.add_argument('--Tortuosity', type = string2bool, choices = [True, False], default = config['Tortuosity'],
                      help = "Tortuosity calculation setting; Requires --Voxel_dist 'Uniform' and --tol -1 [default = YAML; Locked to True or False]")
       ###################################################
       ## GROUP 5: Terminal printing and xyz generation ##
       ###################################################
    printing = xyz_parser.add_argument_group('Terminal printing and xyz generation')
    printing.add_argument('--print_eff', type = int, choices = [0, 1, 2], default = config['print_eff'],
                          help = "Level of printing [default = YAML; Locked to 0, 1, or 2]")
    printing.add_argument('--print_xyz', type = string2bool, choices = [True, False], default = config['print_xyz'],
                          help = "xyz visulatization flag [default = YAML; Locked to True or False]")
       ####################################
       ## GROUP 6: Efficiency parameters ##
       ####################################
    efficiency = xyz_parser.add_argument_group('Efficiency parameters - see YAML description for more details [default = YAML]')
    efficiency.add_argument('--clustering', type = str, choices = ['Neumann', 'Moore'], default = config['clustering'],)
    efficiency.add_argument('--N_write_sph', type = float, default = config['N_write_sph'],)
    efficiency.add_argument('--d_inc', type = float, default = config['d_inc'],)
    efficiency.add_argument('--N_calc_PSD', type = float, default = config['N_calc_PSD'],)
    efficiency.add_argument('--N_edge_gen', type = float, default = config['N_edge_gen'])
    efficiency.add_argument('--tol', type = float, default = config['tol'],)
    efficiency.add_argument('--rand_frac', type = float, default = config['rand_frac'],)

    ########################################################
    ########################################################
    ## SUBPARSER 2: For trajectory files (xtc, trr, etc.) ##
    ########################################################
    ########################################################
    traj_parser = subparsers.add_parser('gmx', help = "Process GROMACS trajectory files")
       ###################################
       ## GROUP 1: Required input files ##
       ###################################
    trj_files = traj_parser.add_argument_group('Required input files')
    trj_files.add_argument('trj_file', type = readable_file,
                            help = "Path to xtc/trr/gro file")
    trj_files.add_argument('top_file', type = readable_file,
                            help = "Path to tpr/gro file")
       ##########################################
       ## GROUP 2: Frame selection and threads ##
       ##########################################
    traj_frames = traj_parser.add_argument_group('Frame selection and threads')
    traj_frames.add_argument('-b', '--t_min', type = float, default = config['t_min'],
                             help = "Start time (ps) [default = YAML]")
    traj_frames.add_argument('-e', '--t_max', type = float, default = config['t_max'],
                             help = "End time (ps) [default = YAML]")
    traj_frames.add_argument('-n', '--N_frames', type = int, default = config['N_frames'],
                             help = "Number of frames to analyze [default = YAML]")
    traj_frames.add_argument('-t', '--N_threads', type = int, default = config['N_threads'],
                             help = "Number of threads for parallelization [default = YAML]")
       ###########################################
       ## GROUP 3: MDAnalysis selection strings ##
       ###########################################
    traj_selection = traj_parser.add_argument_group('MDAnalysis selection strings')
    traj_selection.add_argument('-m', '--system_name', type = str, default = config['system_name'],
                                help = "MDAnalysis selection string defining the system matrix, e.g., 'moltype MOL', 'resname PEO', 'resname SOL LI CL' [default = YAML]")
    traj_selection.add_argument('-s', '--solvent_name', type = str, default = config['solvent_name'],
                                help = "MDAnalysis selection string defining the solvent matrix, e.g., '', 'percolated', 'resname SOL LI CL' [default = YAML]") 
       ##################################
       ## GROUP 4: Important variables ##
       ##################################
    vars = traj_parser.add_argument_group('Important variables')
    vars.add_argument('-L', '--L_voxel', type = float, default = config['L_voxel'],
                      help = "Voxel side length (A) [default = YAML]")
    vars.add_argument('-r', '--probe_radius', type = float, default = config['probe_radius'],
                      help = "Probe radius (A) [default = YAML]")
    vars.add_argument('--d_max', type = float, default = config['d_max'],
                      help = "Max PSD diameter (A) [default = YAML]")
    vars.add_argument('--d_step', type = float, default = config['d_step'],
                      help = "PSD bin size (A) [default = YAML]")
    vars.add_argument('--Voxel_dist', type = str, choices = ['Uniform', 'Random'], default = config['Voxel_dist'],
                      help = "Voxel distribution setting [default = YAML; Locked to 'Uniform' or 'Random']")
    vars.add_argument('--PSD_FFV', type = string2bool, choices = [True, False], default = config['PSD_FFV'],
                      help = "Pore size distribution and free volume fraction calculation setting [default = YAML; Locked to True or False]")
    vars.add_argument('--Surface_area', type = string2bool, choices = [True, False], default = config['Surface_area'],
                      help = "Surface area calculation setting; Requires --Voxel_dist 'Uniform' and --tol -1 [default = YAML; Locked to True or False]")
    vars.add_argument('--Tortuosity', type = string2bool, choices = [True, False], default = config['Tortuosity'],
                      help = "Tortuosity calculation setting; Requires --Voxel_dist 'Uniform' and --tol -1 [default = YAML; Locked to True or False]")
       ###################################################
       ## GROUP 5: Terminal printing and xyz generation ##
       ###################################################
    printing = traj_parser.add_argument_group('Terminal printing and xyz generation')
    printing.add_argument('--print_eff', type = int, choices = [0, 1, 2], default = config['print_eff'],
                          help = "Level of printing [default = YAML; Locked to 0, 1, or 2]")
    printing.add_argument('--print_xyz', type = string2bool, choices = [True, False], default = config['print_xyz'],
                          help = "xyz visulatization flag [default = YAML; Locked to True or False]")
       ####################################
       ## GROUP 6: Efficiency parameters ##
       ####################################
    efficiency = traj_parser.add_argument_group('Efficiency parameters - see YAML description for more details [default = YAML]')
    efficiency.add_argument('--clustering', type = str, choices = ['Neumann', 'Moore'], default = config['clustering'],)
    efficiency.add_argument('--N_write_sph', type = float, default = config['N_write_sph'],)
    efficiency.add_argument('--d_inc', type = float, default = config['d_inc'],)
    efficiency.add_argument('--N_calc_PSD', type = float, default = config['N_calc_PSD'],)
    efficiency.add_argument('--N_edge_gen', type = float, default = config['N_edge_gen'])
    efficiency.add_argument('--tol', type = float, default = config['tol'],)
    efficiency.add_argument('--rand_frac', type = float, default = config['rand_frac'],)

    # Define args
    args = parser.parse_args(remaining_argv)

    # If a .gro is provided instead of .xtc or .trr, N_threads = N_frames = 1
    if args.mode == 'gmx' and '.gro' in args.trj_file:
        if args.N_threads != 1:                        parser.error("gro file inputs require N_threads = 1")
        if args.N_frames != 1 and args.N_frames != -1: parser.error("gro file inputs require N_frames = 1")
    
    # --Voxel_dist 'Uniform' and --tol -1 are required for SA calculations
    if args.Surface_area == True:
        if args.Voxel_dist != 'Uniform':               parser.error("SA calculation requires --Voxel_dist 'Uniform'")
        if args.tol != -1:                             parser.error("SA calculation requires --tol -1")
    
    # --Voxel_dist 'Uniform' and --tol -1 are required for Tau calculations
    if args.Tortuosity == True:
        if args.Voxel_dist != 'Uniform':               parser.error("Tortuosity calculation requires --Voxel_dist 'Uniform'")
        if args.tol != -1:                             parser.error("Tortuosity calculation requires --tol -1")

    # Define data arrays from YAML
    Size_arr = np.array(config['Size_arr'], dtype=object)
    Dummy_atoms = np.array(config['Dummy_atoms'])

    return args, Size_arr, Dummy_atoms


if __name__ == "__main__":
    args, Size_arr, Dummy_atoms = loadArgs()

    print('########################################')
    print('########### Input Parameters ###########')
    print('########################################\n')
    for key, value in vars(args).items():
        if 'mode' in key:          print(  '    ############# Mode #############')
        elif 'trj_file' in key:    print('\n    ### Files and Run Parameters ###')
        elif 'system_name' in key: print('\n    ######## System/Solvent ########')
        elif 'L_voxel' in key:     print('\n    ########## Variables ###########')
        elif 'clustering' in key:  print('\n    #### Efficiency Parameters #####')

        if '_calc' in key or '_write' in key or 'target_' in key or '_gen' in key: print(f"    {key:18}: {value:.0e}")
        else:                                                                      print(f"    {key:18}: {value}")
    print('\n########################################')
    print(  '########################################')
    print(  '########################################\n')


    main()


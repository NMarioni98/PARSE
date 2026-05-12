#!/bin/bash

######################################################################################
########## Running PrO-VAT using an xyz + dat file (PoreBlazer-style input) ##########
######################################################################################

# Load in the trajectory data creating a PrO-VAT.hdf5 data file
python3 ../../../PrO-VAT.py ../../config.yaml xyz pore.xyz input.dat -L 0.25 -r 0.1 --d_max 5 --d_step 0.2 --Voxel_dist 'Uniform'

# Run the PSD analysis
python3 ../../../PrO-VAT.py ../../config.yaml xyz pore.xyz input.dat -L 0.25 -r 0.1 --d_max 5 --d_step 0.2 --Voxel_dist 'Uniform'

# For a xyz file input, PrO-VAT.py relies on less input variables - several inputs are "Locked"
# PrO-VAT.py takes the following inputs
#    python3 {PATH}/PrO-VAT.py {Path}/{YAML config file} xyz {PATH}/{input .XYZ file} {PATH}/{input .DAT file} {Optional arguments}
#     - All files, inclduing PrO-VAT.py, can be in the current folder, e.g., 'polymer_matrix.xyz', or a different folder, e.g., '{Path-to-file}/polymer_matrix.xyz'
#     - All output files are generated in the current working directory
#     - All input parameters are defined within the config.yaml file
#       - See below for important parameters
#     - There are 2 ways to change the parameters
#       1) Edit the config.yaml file that is called into PrO-VAT.py
#       2) Overwrite the default yaml parameter, e.g., adding "--Voxel_dist 'Uniform'" to the PrO-VAT.py command overwrites the "Voxel_dist" parameter in the yaml file, etc

# Important considerations when running PrO-VAT.py:
#   --system_name:          MDAnalysis selection string defining the system matrix, e.g., ''all', 'moltype MOL', 'resname PEO', 'resname SOL LI CL'. Typically 'all' for "xyz" mode.
#   --solvent_name:         solvent_name is either '' to probe the entire van der Waals free volume of the provided atoms, 'percolated' to probe the largest (assumed percolated) free volume cluster, or an MDAnalysis selection string to only probe free volume clusters containg solvent atoms.
#   --L_voxel:              defines the approximate size of the voxels the system is broken down in to. Smaller voxels take exponentially longer to analyze. Typically 0.5-1.0 angstroms.
#   --probe_radius:         defines the size of the probe, where the minimum pore size analyzed has a diameter of 2*probe_radius. Smaller values take exponentially longer to analyze. Typically 1.4-1.575 for a "water molecule" probe.
#   --Voxel_dist:           create 'Uniform' or pseudo-"Random" voxel positions. Typically "Random", but visualizing the system using print_xyz is clearer with 'Uniform.
#   --PSD_FFV:              calculate the PSD and Connolly FFV of the solvent matrix.
#   --Surface_area:         calculate the Connolly (--PSD_FFV True) and Lee-Richards surface area of the solvent matrix. Requires --Voxel_dist 'Uniform' and --tol -1.
#   --Tortuosity:           calculate the 1D diffusional tortuosity of the solvent matrix. Requires --Voxel_dist 'Uniform' and --tol -1. Memory intensive on large systems or small --L_voxel.
#   --d_max and d_step:     defines the binning for the PSD. It may be useful to change d_step to achieve smoother profiles. Typically 50.0 and 0.25-0.50, respectively.
#   --print_eff:            defines how much information is printing while PrO-VAT.py is running. Typically 1, but 2 is useful for trouble shooting memory errors or significant slowdowns in compute time.
#   --print_xyz:            defines whether PrO-VAT.py generates xyz files to visualize the probe-occupiable volume that is analyzed. Typically False to conserve harddrive space (these xyz can get large for small L_voxel or large simulation boxes).
#   --N_repeats:            Number of times to analyze each frame. Typically 1.
#   --N_threads:            number of threads.
#   - NOTE: If more frames are available, always prioritize increasing --N_frames over --N_repeats. --N_repeats > 1 can smooth out single-frame analyses and reduce artifacts due to the cubic nature of voxels.
#   - NOTE: system_name, t_min, t_max, N_frames are "Locked" for xyz inputs

# - Large systems or small --L_voxel can run very slow. Setting tol > 0 can significantly speed up the calculation without losing significant accuracy.
#   --tol:            set the tolerance for measuring the PSD. --tol -1 means that all voxels are analyzed. tol > 0 means that the PSD calculation will end early when the largest error in the PSD is less than that value. Typically -1.
#   --rand_frac:      Fraction of voxels to analyze each cycle during PSD and FFV calculation. rand_frac defaults to a value of 1 if --tol -1 or --rand_frac >= 0.5. Typically 0.01.
#                     If tol > 0, then 2*rand_frac*N_voxels are analyzed at a minimum for PSD and FFV, e.g., --rand_frac 0.01 means that a minimum of 2% of all voxels (2 cycles) are analyzed

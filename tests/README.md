# PARSE Examples

All tests performed on a Windows laptop with an 8-core 8-thread Intel(R) Core(TM) Ultra 7 258V CPU and 32 GB of LPDDR5X-8533MT/s RAM

### Tests and Examples
 - **/tests/xyz/Example_\*/:** example analyses on PoreBlazer-style xyz/dat trajectory input for PARSE
   - Example_{Cylindrical_Pore, Rectangular_Pore, 2D_Channel} contain ideal pore and channel geometries
     - PARSE and PoreBlazer complete calculations in approx. 1 min.
   - Example_AEM contains an anion exchange membrane (*p*5CNMe3 - *λ* = 10) from: https://doi.org/10.1021/acs.macromol.5c01789
     - PARSE completes calculations in approx. 40 s, PoreBlazer completes calculations in approx. 5 min.
   - Example_GRF contains a Gaussian random field reconstruction of a cation exchange membrane (*p*5PhSH - *Y* = 69, *λ* = 9) from: Wang, L.; Kronenberger, S.; Marioni, N.; Frischknecht, A.L.; Jayaraman, A.; Winey, K.I. *In Preparation* **2026**.
     - PARSE completes calculations in approx. 25 min., PoreBlazer fails to complete calculations within 24 hr.
 - **/tests/trj/Example_\*/:** example analyses on molecular dynamics trajectory inputs for PARSE
   - Example_CEM contains a cation exchange membrane (*p*5PhSH - *Y* = 70, *λ* = 9) from: https://doi.org/10.1021/jacsau.5c00218
     - Example using GROMACS gro/tpr/xtc files
     - PARSE completes calculations over 1 frame in approx. 40 s and over 24 frames in approx. 3.5 min., PoreBlazer completes calculations over 1 frame in approx. 5 min.
   - Adding Soon: Example_AEM contains an anion exchange membrane...
     - Example using LAMMPS data/dcd file...
 - **NOTE:** It is recommended to average results over many different frames and several independent repeats for the best results. These systems just serve as simple, fast to analyze examples of using PARSE.

 ## Example files and folders
 - **/tests/{xyz/trj}/Example_\*/run.sh:** bash script describing how to perform the analysis for the example systems
   - See the file for more details on running PARSE
 - xyz input files:
   - **/tests/xyz/Example_\*/\*.xyz:** xyz file that defines the atoms making up the system matrix, i.e., the solvent domain has been deleted to probe the solvent-phase PSD, etc
   - **/tests/xyz/\*/input.dat:** input file that defines the box size for PARSE
     - See the file for more details on creating and formatting the file
 - trj input files:
   - **/tests/trj/Example_CEM/md.gro:** GROMACS gro file defining the atom positions and limited topology information
   - **/tests/trj/Example_CEM/md.xtc:** GROMACS xtc file defining the atom positions over 24 frames with a timestep of 2500 ps (t_start = 242500 ps, t_end = 300000 ps)
   - **/tests/trj/Example_CEM/md.tpr:** GROMACS tpr file defining the system topology
   - **/tests/trj/Example_AEM/......:** LAMMPS data file....
 - **/tests/{xyz/trj}/\*/output.txt:** an example of PARSE's output when it is run as shown in run.sh
   - NOTE: PARSE must be run twice. First to generate a data file, then to calculate the PSD, etc
 - **/tests/{xyz/trj}/Example_\*/Example_output_files/:** contains example files generated when PARSE is run as shown in run.sh
   - PSD.dat: pore size distribution (PSD, or free volume distribution, channel width distribution, etc)
   - Cumulative_PSD.dat: cumulative PSD, where the PSD is the negative derivative of this profile
   - PSD_Plot.xlsx: Excel plot of the cumulative PSD and PSD plus comparisons to PoreBlazer or multi-frame analyses where applicable
     - **NOTE**: PSDs are recalculated from the cumulative distribution with the diameters in units of nm. Results from PoreBlazer are shifted to be right-justified and the derivative is recalculated to match PARSE's output.
     - **NOTE**: Some deviations between PARSE and PoreBlazer are expected due to differences in the voxelization strategy and methodology. However, we observe that the shape and location of the PSD curves are in good agreement
   - FFV.dat: Connolly and Lee-Richards volume fraction (FFV, free volume fraction, water volume fraction, porosity, etc)
   - SA.dat: a simple marching-cubes mesh surface area calculation ([scikit-image](https://scikit-image.org/)) of the Connolly and Lee-Richards pore surface
     - **NOTE:** The SA calculation requires --Voxel_dist 'Uniform' and --tol -1 (see "Surface_area" in config.yaml)
   - Tau.dat: 1D diffusional tortuosity of the percolated domain in the X, Y, and Z direction using a simple Fickian diffusion algorithm ([PoreSpy](https://porespy.org/))
     - Tortuosity does not account for PBCs and applies to the Lee-Richards volume
     - **NOTE:** Tortuosity calculation is memory intensive. Larger --L_voxel may be needed for large systems, i.e., /tests/xyz/Example_GRF/
     - **NOTE:** The tortuosity calculation requires --Voxel_dist 'Uniform' and --tol -1 (see "Tortuosity" in config.yaml)
   - \*.xyz: xyz files to visualize the free volume probed by PARSE using OVITO
     - Free_Volume_Spheres visualizes the free volume *spheres* of maximum radius --probe_radius (defined in config.yaml) which make up the probe-occupiable free volume
     - Free_Volume_Voxels visualizes the free volume *voxels* of side length --L_voxel (defined in config.yaml) which make up the probe-occupiable free volume
       - "Alpha" values represent the bin of the largest free volume *sphere* that contains each *voxel*, enabling visualization of voxels based on their associated bin in PSD
     - Free_Volume_Surface visualizes the free volume *voxels* of side length --L_voxel (defined in config.yaml) which defines the surface of the probe-occupiable free volume
       - "X" particles define the Connolly surface while "Y" particles define the Lee-Richards surface
       - **NOTE:** voxels are centered on the voxel face at the surface
   - Example.ovito: OVITO visualization state showing the above. Created using [OVITO(-basic) version 3.7.12](https://www.ovito.org/download_history/)
     - In the top right, you can select the different xyz files (pipelines) and turn on and off the Particles (under "Visual elements")
     - If GROMACS gro and PARSE output xyz boxes do not overlap in OVITO, uncheck "Center simulation box on coordinate origin" under the "Gromacs reader" for the gro file pipeline
 - **/tests/{xyz/trj}/Example_\*/PoreBlazer/:** contains example files generated when [PoreBlazer v4.0](https://github.com/SarkisovGitHub/PoreBlazer) is run with inputs comparable to PARSE's inputs in run.sh

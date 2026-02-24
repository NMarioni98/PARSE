# PrO-VAT (Probe-Occupiable Volume Analysis Tools)

Written by: Nico Marioni, nmarioni@seas.upenn.edu
 - Developed using Python 3.12.X
   - Packages: numpy 2.3.3+, h5py 3.14.0+, MDAnalysis 2.9.0+, igraph 1.0.0, scikit-image 0.25.0+, PyYAML 6.0.3
 - Based on **/PyAnalysis/analysis_PSD_Voxel.py** (https://github.com/Ganesan-Group-Codes-and-Analysis/PolyAnalysis) by the same author

PrO-VAT.py calculates the pore size distribution (free volume distribution, channel width distribution, etc) of the van der Waals volume of the defined system matrix from a GROMACS (gro/xtc/trr + tpr) or PoreBlazer-style (xyz + dat) trajectory. This script was specifically desgined to find the distribution of water-rich pores within a hydrated polymer system, but can be generalized to any atomic or coarse-grained system. The output includes the Cumulative Pore Size Distribution (Cumulative PSD), Pore Size Distribution (PSD), Free Volume Fraction (Fractional Free Volume, FFV), and Surface Area (SA) with optional xyz visualizations. This script was written based on the methods used for PoreBlazer v4.0 (https://github.com/SarkisovGitHub/PoreBlazer, https://doi.org/10.1021/acs.chemmater.0c03575) and is optimized for parallelized calculations over many system frames, or analysis of large (30+ nm box length) systems.

## Getting started

PrO-VAT.py requires the following inputs:
 - ```python3 PrO-VAT.py {YAML} {Mode} {Trajectory} {Topology} {Optional arguments}```
   - **YAML:** "config.yaml", configuration file containing default PrO-VAT.py inputs
     - ```python3 PrO-VAT.py -h``` for more information
   - **Mode:** "xyz" or "gmx" for PoreBlazer-style or GROMACS trajectory input, respectively
     - ```python3 PrO-VAT.py {YAML} -h``` for more information
   - **Trajectory:** xyz or gro/xtc/trr file input for "xyz" or "gmx" mode, respectively
   - **Topology:** dat or tpr file input for "xyz" or "gmx" mode, respectively
   - **Optional arguments:** additional (optional) arguments can be added to overwrite the default settings defined in {**YAML**}
     - e.g., "-r 1.4" or "--probe_radius 1.4"
     - ```python3 PrO-VAT.py {YAML} {Mode} -h``` for more information
 - **NOTE:** PrO-VAT.py must be run twice: First to generate a PrO-VAT.hdf5 run file, Second to perform the analysis.
   - If PrO-VAT.py successfully runs the second time, it will delete the PrO-VAT.hdf5 file. However, make sure to delete this file and rebuild it if you change the inputs for PrO-VAT.py (yaml file or arguments) in between analysis attempts.

## Repository Contents

### Files
 - **PrO-VAT.py:** python analysis script
 - **config.yaml:** yaml file containing default inputs for PrO-VAT.py
   - See the file for more details on each input
 - **run.sh:** bash script describing how to perform the analysis for two example systems
   - See the file for more details on running PrO-VAT.py

### Examples
 - **/xyz/:** example analysis on a PoreBlazer-style xyz/dat trajectory input for PrO-VAT.py
   - Example system contains an anion exchange membrane from: https://doi.org/10.1021/acs.macromol.5c01789
     - *p*5CNMe3 - *λ* = 10
 - **/gmx/:** example analysis on a GROMACS gro/tpr (similarly xtc, trr) input for PrO-VAT.py
   - Example system contains a cation exchange membrane from: https://doi.org/10.1021/jacsau.5c00218
     - *p*5PhSH - *Y* = 70, *λ* = 6

### Example files and folders
 - xyz input files:
   - **/xyz/polymer_matrix.xyz:** XYZ file that defines the atoms making up the polymer matrix, i.e., the domain of interest has been deleted
   - **/xyz/input.dat:** input file that defines the box size for PrO-VAT
     - See the file for more details on creating and formatting the file
 - gmx input files:
   - **/gmx/md.gro:** GROMACS gro file defining the atom positions
   - **/gmx/md.tpr:** GROMACS tpr file defining the atom topology
   - **/gmx/topol.top:** GROMACS top file defining the atom topology - not used by PrO-VAT.py, but human readable
 - **/{xyz/gmx}/output.txt:** an example of PrO-VAT's output when it is run as shown in run.sh
   - NOTE: PrO-VAT must be run twice. First to generate a data file, then to calculate the PSD, etc
 - **/{xyz/gmx}/Example_output_files/:** contains example files generated when PrO-VAT.py is run as shown in run.sh
   - PSD.dat: pore size distribution (PSD, or free volume distribution, channel width distribution, etc)
   - Cumulative_PSD.dat: cumulative PSD, where the PSD is the derivative of this profile
   - FFV.dat: fractional free volume
   - SA.dat: a simple marching-cubes mesh surface area calculation of the pore surface
   - {}.xyz: xyz files to visualize the free volume probed by PrO-VAT using OVITO
     - Free_Volume_Spheres visualizes the free volume *spheres* of maximum radius R which make up the probed free volume
     - Free_Volume_Voxels visualizes the free volume *voxels* of side length L_voxel (defined in "config.yaml") which make up the probed free volume
     - Free_Volume_Surface visualizes the free volume *voxels* of side length L_voxel (defined in "config.yaml") which defines the surface of the probed free volume
      - **NOTE:** voxels are shifted by L_voxel/2 so it is centered on the voxel face at the surface
      - **NOTE:** The SA calculation assumes the Free_Volume_Voxels are "Uniform" not "Random" (see "Voxel_dist" in "congif.yaml")
   - Example.ovito: Ovito visualization state showing the above. I recommend installing Ovito(-basic) version 3.7.12, https://www.ovito.org/download_history/
     - In the top right, you can select the different xyz files (pipelines) and turn on and off the Particles (under "Visual elements")

# Local Python Install using WSL (Windows Subsystem for Linux)

## WSL

### Install
 - Open Windows PowerShell and run ```wsl --install ubuntu```
 - Restart computer and start up wsl (may take some time)
   - If wsl instantly closes upon opening it, run ```wsl --install ubuntu``` again

### Using WSL
 - WSL mounts a Linux subsystem in Windows
 - If you open WSL, you will be in the '/home/{username}' directory
   - Your windows files are located at '/mnt/c/Users/{username}'
   - If you want to access Box, you need to mount is separately. Run ```sudo mkdir /mnt/box/``` once, then add the following to ~/.bashrc:
     ```
     if mount | grep /mnt/box > /dev/null; then
       sudo mount -t drvfs 'C:\Users\{username}\Box\' /mnt/box/
     ```
    - This does require you to input your password the first time you start a wsl instance.
    - You can also add a cd command, aliases, etc to ~/.bashrc to automate certain tasks every time a terminal is started.
 - Recommended settings
   - In the windows WSL terminal, open settings
     - Interaction > Automatically copy selection to clipboard = On
     - Defaults > Advanced > Bell notification style = None
     - Ubuntu > Advanced > Bell notification style = None


### Installing python
 - open wsl
```
cd ~
sudo apt update && sudo apt upgrade
sudo apt install python3-full
mkdir -p ~/.venvs
python3 -m venv ~/.venvs/myPython
source ~/.venvs/myPython/bin/activate
python3 -m pip install PyYAML h5py MDAnalysis igraph scikit-image
```
 - Must ```source ~/.venvs/myPython/bin/activate``` to start python in virtual environment
   - Add this command to ~/.bashrc if you want this to happen automatically when you open wsl

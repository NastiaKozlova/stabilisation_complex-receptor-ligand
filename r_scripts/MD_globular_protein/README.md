# MD_globular_protein

# Introduction

This group of scripts supports preparing and analyzing MD simulation of globular proteins.
MD simulation are run using only NAMD. 
You should provide your own NAMD copy and CHARMM force fields.
Scripts were tested using CHARMM 36 force fields  and NAMD 2.12 and NAMD 2.14.
Scripts were tested using Debian 10.9.
You should copy structure of your protein to the _*start/structure*_ directory.
You should copy your toppar folder with force fields to the _*start/*_ directory.

Main script from with program will run is _*r_scripts/main_script.R*_
"r_scripts/prepare_to_stabilisation_MD.R" -- prepare all scpipts and structures 
to run MD simulations. MD simulations are contained 0.1ns minimisation, 0.3ns heating,
2 ns eqvilibration and 125 ns productive simulation(5\*25ns). 
You can alter time by changing _*num_din*_

After this step you can find your namd script _*MD\_globular\_protein/r\_scripts/namd\_script.txt*_
This script have to be altered to accommodate different configurations of computers. 

# List of dependenses

R pacages:
1. bio3d
2. dssp
3. cowlpot
4. ggplot2
5. dplyr

Programs:
1. VMD
2. NAMD
3. R

Extra modules for vmd: orient and la101psx
https://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/

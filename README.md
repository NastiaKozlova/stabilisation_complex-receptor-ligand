# stabilisation_complex-receptor-ligand

# MD_globular_protein

# Introduction

This group of scripts supports preparing and analyzing MD simulations of the complex receptor(globular proteins)-ligand.
Scripts are separate to the four parts:
1. Protein structure prediction by homoloigy
2. Protein structure stabilization by MD simulation using NAMD
3. Receptor-ligand complex constraction
4. Receptor-ligand complex structure stabilization by MD simulation using NAMD

#Protein structure prediction

To find rigth homology  
MD simulations are run using only NAMD. 
You should provide your own NAMD copy and CHARMM force fields.
Scripts were tested using CHARMM 36 force fields and NAMD 2.12 and NAMD 2.14.
Scripts were tested using Debian 10.9.
You should copy the structure of your protein to the _*start/structure*_ directory.
You should copy your toppar folder with force fields to the _*start/*_ directory.

The main script from with program will run is _*r_scripts/main_script.R*_
"r_scripts/prepare_to_stabilisation_MD.R" -- prepare all scripts and structures 
to run MD simulations. MD simulations are contained 0.1ns minimisation, 0.3ns heating,
2 ns eqvilibration and 100 ns productive simulation(100\*1ns). 
You can alter time by changing _*num_din*_

After this step you can find your namd script _*MD\_globular\_protein/r\_scripts/namd\_script.txt*_
This script has to be altered to accommodate different configurations of computers. 

# List of dependencies

R pacages:
1. bio3d install.packages("bio3d")
2. cowlpot install.packages("cowlpot")
3. ggplot2 install.packages("ggplot2")
4. dplyr install.packages("dplyr")
5. httr install.packages("httr")

Programs:
1. VMD
2. NAMD
3. R
4. dssp 
5. openbabel

Extra modules for VMD: orient and la101psx
https://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/

How to install essential packages for Linux and Debian

sudo apt-get install libgdal-dev libssl-dev
sudo apt-get install muscle
sudo apt-get install openbabel

Where is 2 points where is essential to exit from the R environment
1. to run namd to stabilise structure of the receptor
2. to run namd to stabilise structure of receptor-ligand complex
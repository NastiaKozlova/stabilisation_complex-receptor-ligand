# Preparation of stable receptor-ligand structure 

# MD_globular_protein

# Introduction

This group of scripts supports preparing and analyzing MD simulations of the complex receptor(globular proteins)-ligand.
Scripts are separate to the four parts:
1. Protein structure prediction by homoloigy
2. Protein structure stabilization by MD simulation using NAMD
3. Receptor-ligand complex constraction
4. Receptor-ligand complex structure stabilization by MD simulation using NAMD

The main script from with program will run is _*r_scripts/master_script.R*_

MD simulations are run using only NAMD. 
You should provide your own NAMD copy and CHARMM force fields.
Scripts were tested using CHARMM 36 force fields and NAMD 2.12 and NAMD 2.14.
Scripts were tested using Debian 10.9.
You should copy your toppar folder with force fields to the stabilisation_complex-receptor-ligand/start/*_ directory and add there all additional fource fields.
This project support running multiple protein structures at the same time. In this README ProteinName - is one of the proteins name. you shold not use spaces, "_", "-". in the protein and ligand, and center snames. This symblos are confliting with scripts. Scripts will work wrong. 

# Protein structure prediction, optional for uncomplite proteins

Protein sequences should be in and _*stabilisation_complex-receptor-ligand/start/structure_prediction/sequence/ProteinName.fasta*_ 

To find rigth homology model you should have sthuctures to chouse from put into _*/home/nastia/projects/current/stabilisation_complex-receptor-ligand/start/structure_prediction/pdb/ProteinName*_
You should check most appropriate strutures chousen by scripts from stabilisation_complex-receptor-ligand/predicted/ProteinName folder, and predict full structure based on chousen structure by homology prediction method of your chouse. We recomend Robetta. 

# Protein structure stablilisation using MD simulations (NAMD), needed for relaxing protein structure before docking

You should put predicted or your good known strcutre into _*stabilisation_complex-receptor-ligand/start/MD_stabilisation/ProteinName.pdb*_
"stabilisation_complex-receptor-ligand/r_scripts/MD_globular_protein/r_scripts/prepare_to_stabilisation_MD.R" -- prepare all scripts and structures to run MD simulations. 
MD simulations are contained 0.1ns minimisation, 0.3ns heating, 2 ns eqvilibration and 100 ns productive simulation(100\*1ns).
You can adjast length of simulation by editing _*num_din*_ parameter. _*num_din*_ =time in ns

After this step you can find your namd script _*MD\_globular\_protein/r\_scripts/namd\_script.txt*_
This script has to be altered to accommodate different configurations of computers. 

# Prediction of structure receptor-ligand complex this docking AutoDock, needed for putting ligand roughly near real coordinates of ligands in ligand-receptor complex.

You should put all of ligand structures in pdb form into _*stabilisation_complex-receptor-ligand/start/docking/docking_first/ligand_start/*_

# Receptor-ligand structure stablilisation using MD simulations (NAMD)

This step needed to find correct plase and interactions between receptor and ligand.

The main script from with program will run is _*r_scripts/master_script.R*_
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

Extra modules for VMD: orient and la101psx shold be downloaded and untar to _*stabilisation_complex-receptor-ligand/programs/*_
https://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/

autodock vina into _*stabilisation_complex-receptor-ligand/programs/autodock_vina_1_1_2_linux_x86*_ from http://vina.scripps.edu/download.html
MGLTools into _*stabilisation_complex-receptor-ligand/programs/MGLTools-1.5.7*_ from https://ccsb.scripps.edu/mgltools/downloads/, and you shold chouse this wersion of the ptogram _*mgltools_Linux-x86_64_1.5.7_Install (Linux 64 GUI installer 109Mb)*_ or rename folder atherthords to _*stabilisation_complex-receptor-ligand/programs/MGLTools-1.5.7*_

How to install essential packages for Linux and Debian

sudo apt-get install libgdal-dev libssl-dev
sudo apt-get install muscle
sudo apt-get install openbabel

Where is 3 points where is essential to exit from the R environment
1. to predict structure of the protein using Robetta or other programm
2. to run namd to stabilise structure of the receptor
3. to run namd to stabilise structure of receptor-ligand complex

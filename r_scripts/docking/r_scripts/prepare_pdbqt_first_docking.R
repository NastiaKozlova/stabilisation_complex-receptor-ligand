for (i in 1:nrow(df_ligand)) {
  system(command = paste0("obabel ",part_name,"ligand_start/",df_ligand$ligand[i], ".pdb -O ",part_name,"ligand/",df_ligand$ligand[i], ".pdbqt"),ignore.stdout=T,wait = T)
}
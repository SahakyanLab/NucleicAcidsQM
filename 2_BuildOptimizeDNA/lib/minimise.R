minimise <- function(prmtop, coor, maxcyc, ncyc, target,
                     igb, output.path) {
  
  # Minimisation wrapper for sander (imin = 1)
  # ncyc    <int>   Number of steepest descent minimization
  # maxcyc  <int>   Total number of minimisation.
  #                 steepest descent + conjugate gradient
  # igb     <int>   Minimise with implicit solvation model.
  #                 0: in vacuo
  #                 1: standard generalized Born method
  # ntb     <int>   1: Use constant volume periodic boundaries.
  # restrain <bool> ntr=1: Turn on Cartesian restraints
  # restraint_wt <numeric> force constant for restraint
  # restraintmask <character> Restraint atoms/residues
  # drms  <numeric> Convergence criterion. Default is 10−4 kcal/(mol.Å).
  
  condition = paste0("min",
                     "_igb_", igb,
                     "_step_", maxcyc,
                     "_optimize_", target)
  
  if(target=="all_atoms"){restrain_atoms<-"!@*"}
  if(target=="no_atoms"){restrain_atoms<-"@*"}
  if(target=="only_H"){restrain_atoms<-"!@H="}
  if(target=="other_than_H"){restrain_atoms<-"@H="} 
  
  
  writeLines(paste0(
    'minimisation setup\n',
    ' &cntrl\n',
    '  imin = 1,\n',
    if(maxcyc=="FULL") {paste0('  ntmin = 0,\n','  maxcyc = 100000,\n')}
    else{paste0('  maxcyc =', maxcyc,',\n','  ncyc =', ncyc,',\n')},     
    '  ntb = 0,', '\n',
    '  igb = ', igb, ',\n',
    '  cut = 100', ',\n',
    '  ntr = 1, ',
    '  restraint_wt = 10, ',
    '  restraintmask = "', restrain_atoms, '",\n',
    ' /'
    
  ), paste0(output.path, "/", condition,".in"))
  
  cmd <- paste0("sander -O -i ", output.path, condition,".in ",
                         " -o ", output.path, condition, ".out ",
                         " -c ", coor,
                         " -p ", prmtop,
                         " -r ", output.path, condition, ".ncrst",
                         " -inf ", output.path, condition, "_mdinfo",
                         " -ref ", coor)
  
  system(cmd)

  # output pdb file after minimization   
  min_pdb <- paste0(output.path, condition,".pdb")
  min_coor <-paste0(output.path, condition, ".ncrst")
  system(paste0("ambpdb -p ", prmtop, " -c ", min_coor, " > ", min_pdb))
  
}

genTopCoor <- function(dna.pdb, neutralise=FALSE, solvate=FALSE,
                       output.path="./", output.name=NULL) {
  # dna.pdb   <character>   A path to DNA pdb file
  # solvate   <bool>        Default is FALSE.
  
  if (is.null(output.name)) {
    output.name <- gsub("\\.pdb$", "", basename(dna.pdb))
  }
  output.name <- paste0(output.path, "/", output.name)
  
  writeLines(paste0(
    'source leaprc.DNA.OL15\n',
    if(solvate | neutralise) 'source leaprc.water.tip3p\n',
    'dna = loadpdb ', dna.pdb, '\n',
    if(neutralise) 'addions dna Na+ 0\n',
    if(solvate) 'solvatebox dna TIP3PBOX 8.0\n',
    'saveamberparm dna ', output.name, '.prmtop ', output.name, '.rst7\n',
    'quit'
    
  ), paste0(output.path, "tleap.in"))
  
  cmd <- paste0("tleap -f ", output.path, "/tleap.in")
  system(cmd)
  file.rename("leap.log", paste0(output.path, "/", "leap.log"))

}
    

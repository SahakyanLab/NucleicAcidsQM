################################################################################
#
# AIM: Build and optimize K-mer DNA
#
# TASKS
#   1) To generate all sequences of K-mer DNA
#   2) To build DNA structure 
#   3) To optimise the DNA structure
#   4) To calculate RMSD
#

## Configuration ###############################################################
options(warn = 2)

## Task  #######################################################################
toDo = NULL

toDo$generateKmers = F
Kmer.num = 8

toDo$initDNA      = F

toDo$optDNA       = F
opt.igb = 1
opt.maxcyc = 5000 
opt.ncyc = 2500
opt.atoms = "all_atoms"           

toDo$calculateRmsd = F
rmsd.atom = "onlyP"
rmsd.align = TRUE
rmsd.targetpdb = "min_igb_1_step_5000_optimize_all_atoms.pdb"


## input conformation ##########################################################
conformation = "DNA_B_Handbook1999"
strand       = "duplex"
neutralise   = FALSE
solvate      = FALSE


# location of output folder
output.path = paste0("data/", conformation,"_",strand)
if (neutralise) {output.path = paste0(output.path, "_neutrized")}
if (solvate) {output.path = paste0(output.path, "_solvated")}


# AmberTools path
ambertools.path <- "/opt/anaconda3/envs/AmberTools21"
Sys.setenv(PATH = paste0(ambertools.path, "/bin:", Sys.getenv("PATH")))
Sys.setenv(AMBERHOME = ambertools.path)


## Dependant functions #########################################################

source("lib/GenerateKmers.R")
source("lib/buildDNA.R")
source("lib/genTopCoor.R")
source("lib/minimise.R")
source("lib/calculateRMSD.R")

## Directory setup #############################################################
dir.create(output.path, recursive = TRUE, showWarnings = FALSE)

## Main Code ###################################################################

#-------------------------------------------------------------------------------
if (toDo$generateKmers) {
  sequences_arrays=GenerateKmers(Kmer.num)
  file.create("sequences.txt")
  write(sequences_arrays,"sequences.txt",sep="\n",append=T)
}

#-------------------------------------------------------------------------------
if (toDo$initDNA) {
  
  input.seqs = "sequences.txt"
  dna.seqs <- readLines(input.seqs)
  
  for (dna.seq in dna.seqs) {
    output <- paste0(output.path, "/", dna.seq, "/")
    dir.create(output, showWarnings = FALSE)
    buildDNA(conformation, strand, dna.seq, output)
    genTopCoor(paste0(output, dna.seq, ".pdb"), neutralise, solvate, 
               output.path = output, output.name = "dna")
  }
  
}


#-------------------------------------------------------------------------------
if (toDo$optDNA) {

  input.seqs = "sequences.txt"
  dna.seqs <- readLines(input.seqs)
  
  for (dna.seq in dna.seqs) {
    output <- paste0(output.path, "/", dna.seq, "/")
    minimise(paste0(output, "dna.prmtop"),
             paste0(output, "dna.rst7"),
             maxcyc = opt.maxcyc,
             ncyc   = opt.ncyc,
             target = opt.atoms,
             igb = opt.igb,
             output)
  }
}


#-------------------------------------------------------------------------------
if (toDo$calculateRmsd) {
  library(bio3d)
  
  input.seqs = "sequences.txt"
  dna.seqs <- readLines(input.seqs)
  
  reference_pdb=paste0(output.path,"/",dna.seqs,"/",dna.seqs,".pdb")
  target_pdb=paste0(output.path,"/",dna.seqs,"/",rmsd.targetpdb)
  rmsd=mapply(rmsd_bio3d, target_pdb, reference_pdb, rmsd.atom, rmsd.align)
  
  file_names=paste0(output.path,"/",
                    "rmsd_",rmsd.atom,
                    "_align_",rmsd.align,
                    "_for_",gsub(".pdb", "",rmsd.targetpdb),".txt")
  file.create(file_names)
  writeLines(paste0("seq, rmsd"), file_names)
  write(paste0(dna.seqs, ", ", rmsd), file_names,sep="\n",append=T)  
}

################################################################################
#
# AIM: To calculate theoretical UV-Vis spectra
#
# TASKS
#   1) To build DNA structure of UV octamers using nab program
#   2) To optimise the DNA structure
#   3) To calculate UV-Vis spectra
#

## Configuration ###############################################################
options(warn = 2)

## Task  #######################################################################
toDo = NULL

toDo$initDNA                     = T

toDo$optDNA                      = T
opt.igb = 1
opt.maxcyc = 500 
opt.ncyc = 250 
opt.atoms = "all_atoms"           

# input conformation
conformation = "DNA_B_Handbook1999"
strand       = "duplex"
neutralise   = FALSE
solvate      = FALSE

# input sequences
input.seqs = "/Users/kairi/Desktop/2021_12_17_Generate_Kmer/sequences.txt"

# location of output folder
output.path = paste0("data/", conformation,"_",strand)

if (neutralise) {output.path = paste0(output.path, "_neutrized")}
if (solvate) {output.path = paste0(output.path, "_solvated")}


# AmberTools path
ambertools.path <- "/opt/anaconda3/envs/AmberTools21"
Sys.setenv(PATH = paste0(ambertools.path, "/bin:", Sys.getenv("PATH")))
Sys.setenv(AMBERHOME = ambertools.path)


## Dependant functions #########################################################

source("lib/buildDNA.R")
source("lib/genTopCoor.R")
source("lib/minimise.R")


## Directory setup #############################################################
dir.create(output.path, recursive = TRUE, showWarnings = FALSE)

## Main Code ###################################################################

#-------------------------------------------------------------------------------
if (toDo$initDNA) {
  
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

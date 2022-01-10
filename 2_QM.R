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

toDo$QMcalculation= F
QM.hamiltonian="PM6-DH+"
QM.calculation="1SCF"
QM.COSMO=TRUE
QM.charge="twice_of_sequence"

toDo$QMsummary= F
summary.target="QM_PM6-DH+_1SCF_COSMO_from_min_igb_1_step_10_optimize_all_atoms.arc"

toDo$QMRmsd = T
QM_rmsd.atom = "onlyP"
QM_rmsd.align = TRUE
QM_rmsd.targetpdb = "QM_PM6-DH+_1SCF_COSMO_from_min_igb_1_step_10_optimize_all_atoms.pdb"

## input #######################################################################
inputfolder = "data/DNA_B_Handbook1999_duplex_3mer"
inputseqs = "sequences.txt"
inputpdb = "min_igb_1_step_10_optimize_all_atoms.pdb"

#output location is same as inputfolder

## MOPAC path ###################################################################
mopac.path <- "/Users/kairi/Desktop/MOPAC2016_for_Macintosh/"
Sys.setenv(PATH = paste0(Sys.getenv("PATH"),":",mopac.path))


## Dependant functions #########################################################

source("lib/QMminimise.R")
source("lib/QMsummary.R")
source("lib/calculateRMSD.R")

## Main Code ###################################################################

#-------------------------------------------------------------------------------
if (toDo$QMcalculation) {
  dna.seqs <- readLines(inputseqs)
  
  for (dna.seq in dna.seqs) {
    QMminimise(inputfolder  = inputfolder,
               sequence     = dna.seq,
               inputfile    = inputpdb,
               hamiltonian  = QM.hamiltonian,
               calculation  = QM.calculation,
               COSMO        = QM.COSMO,
               charge       = QM.charge)
  }
}

#-------------------------------------------------------------------------------
if (toDo$QMsummary) {
  file_names=paste0("summary_of_",gsub(".arc","",summary.target),".txt")
  file.create(paste0(inputfolder,"/",file_names))
  writeLines(paste0("SEQUENCE,", 
                    "HEAT_OF_FORMATION(KCAL/MOL),",
                    "DIPOLE(DEBYE),",
                    "IONIZATION_POTENTIAL(EV),",
                    "HOMO(EV),",
                    "LUMO(EV)"),
             paste0(inputfolder,"/",file_names))

  dna.seqs <- readLines(inputseqs)
  input_arcs=paste0(inputfolder, "/", dna.seqs, "/", summary.target)
  HFvalues=lapply(input_arcs,extract_HF)
  DIvalues=lapply(input_arcs,extract_DI)
  IPvalues=lapply(input_arcs,extract_IP)
  HOvalues=lapply(input_arcs,extract_HOMO)
  LUvalues=lapply(input_arcs,extract_LUMO)
  
  write(paste0(dna.seqs, ", ",
               HFvalues, ", ",
               DIvalues, ", ",
               IPvalues, ", ",
               HOvalues, ", ",
               LUvalues),
        paste0(inputfolder,"/",file_names),sep="\n",append=T)
  
}


#-------------------------------------------------------------------------------
if (toDo$QMRmsd) {
  library(bio3d)
  
  input.seqs = "sequences.txt"
  dna.seqs <- readLines(input.seqs)
  
  QM_reference_pdb=paste0(inputfolder,"/",dna.seqs,"/",dna.seqs,".pdb")
  QM_target_pdb=paste0(inputfolder,"/",dna.seqs,"/",QM_rmsd.targetpdb)
  QM_rmsd=mapply(rmsd_bio3d, QM_target_pdb, QM_reference_pdb, QM_rmsd.atom, QM_rmsd.align)
  
  file_names=paste0(inputfolder,"/",
                    "rmsd_",QM_rmsd.atom,
                    "_align_",QM_rmsd.align,
                    "_for_",gsub(".pdb", "",QM_rmsd.targetpdb),".txt")
  file.create(file_names)
  writeLines(paste0("seq, rmsd"), file_names)
  write(paste0(dna.seqs, ", ", QM_rmsd), file_names,sep="\n",append=T)  
}
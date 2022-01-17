################################################################################
#
# AIM: QM calculation of K-mer DNA
#
# TASKS
#   1) To calculate K-mer DNA by MOPAC
#   2) To summarize the results of QM calculation 
#   3) To calculate RMSD
#

## Configuration ###############################################################
options(warn = 2)

## Task  #######################################################################
toDo = NULL

toDo$QM = F
QM.hami="PM6-DH+"
QM.cal="1SCF"
QM.COSMO="RSOLV=1.3 EPS=78.4 NSPA=42"
QM.charge="twice_of_sequence"
QM.targetpdb = "min_igb_1_step_5000_optimize_all_atoms.pdb"


toDo$SubQM = T
SubQM.conf="del_mid_2nd_strand_opt"
SubQM.hami=QM.hami
SubQM.COSMO=QM.COSMO
SubQM.charge="depend_on_conf"
SubQM.targetpdb="QM_PM6-DH+_1SCF_COSMO_from_min_igb_1_step_5000_optimize_all_atoms.pdb"


toDo$QMsummary= F
QMsummary.target="QM_PM6-DH+_1SCF_COSMO_from_min_igb_1_step_5000_optimize_all_atoms.arc"


toDo$QMrmsd = F
QMrmsd.atom = "onlyP"
QMrmsd.align = TRUE
QMrmsd.targetpdb = "QM_PM6-DH+_1SCF_COSMO_from_min_igb_1_step_5000_optimize_all_atoms.pdb"

## input #######################################################################
inputfolder = "data/DNA_B_Handbook1999_duplex_7mer"
inputseqs = "sequences.txt"

#output location is same as inputfolder

## MOPAC path ###################################################################
mopac.path <- "/Users/kairimasuda/Desktop/MOPAC2016_for_Macintosh/"
Sys.setenv(PATH = paste0(Sys.getenv("PATH"),":",mopac.path))


## Dependant functions #########################################################

source("lib/QMminimise.R")
source("lib/SubQMcal.R")
source("lib/QMsummary.R")
source("lib/calculateRMSD.R")

## Main Code ###################################################################

#-------------------------------------------------------------------------------
if (toDo$QM) {
  dna.seqs <- readLines(inputseqs)
  
  step=0
  for (dna.seq in dna.seqs) {
    step=step+1
    print(paste0(step,' ',dna.seq))
    QMminimise(inputfolder  = inputfolder,
               sequence     = dna.seq,
               inputfile    = QM.targetpdb,
               hamiltonian  = QM.hami,
               calculation  = QM.cal,
               COSMO        = QM.COSMO,
               charge       = QM.charge)
  }
}

#-------------------------------------------------------------------------------
if (toDo$SubQM) {
  dna.seqs <- readLines(inputseqs)
  
  step=0
  for (dna.seq in dna.seqs) {
    step=step+1
    print(paste0(step,' ',dna.seq))
    SubQMcal(
      inputfolder  = inputfolder,
      sequence     = dna.seq,
      inputfile    = SubQM.targetpdb,
      inputconf    = SubQM.conf,
      hamiltonian  = SubQM.hami,
      COSMO        = SubQM.COSMO,
      charge       = SubQM.charge
    )
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
if (toDo$QMrmsd) {
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
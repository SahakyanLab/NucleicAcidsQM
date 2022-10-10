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

## input conformation ##########################################################
conformation = "DNA_Z_Handbook1999"
neutralise   = FALSE
solvate      = FALSE
sequences = c("readfile","sequence.txt") # c("generateKmer",7) 

## Task  #######################################################################
toDo = NULL

#step1 Build model
toDo$initDNA = F
Bld.outpdb="init.pdb"

#step2 MM and/or QM optimization
toDo$MMopt = F
MMopt.igb  = 1
MMopt.maxcyc = 5000
MMopt.ncyc = 2500
MMopt.atoms = "all_atoms"
MMopt.inpdb=Bld.outpdb
MMopt.outpdb=paste0("MMopt",
                    "_igb_", MMopt.igb,
                    "_step_", MMopt.maxcyc,
                    "_target_", MMopt.atoms,".pdb")

toDo$QMopt = F
QMopt.HF="RHF"
QMopt.hami="PM6-DH+"
QMopt.norm=10.0
QMopt.COSMO="RSOLV=1.3 EPS=78.4 NSPA=42"
QMopt.charge="number_of_phosphate"
QMopt.other="MMOK THREADS=4"
QMopt.inpdb = MMopt.outpdb
QMopt.outpdb=paste0('QMopt_',QMopt.hami,"_COSMO_",
                    "from_", gsub(".pdb","",QMopt.inpdb),".pdb")


#step3 Subsequent single point calculation or analyze
toDo$SubQM = F
SubQM.conf="del_mid_2nd_strand_sp" #"del_mid_2nd_strand_sp"
SubQM.HF=QMopt.HF
SubQM.hami=QMopt.hami
SubQM.COSMO=QMopt.COSMO
SubQM.charge="number_of_phosphate"
SubQM.other=paste0(QMopt.other," ","DISP MULLIK PRTCHAR ALLVEC GRAPHF LARGE BONDS SUPER")  # LOCAL BOND SUPER GRAPHF 
SubQM.inpdb=QMopt.outpdb
SubQM.outpdb=paste0('SubQM_',SubQM.conf,
                   "_from_",gsub(".pdb","",SubQM.inpdb),".pdb")

toDo$QMsummary= F
QMsum.inarc=paste0(gsub(".pdb",".arc",SubQM.outpdb))
QMsum.outfile=paste0("summary_of_",gsub(".arc","",QMsum.inarc),".txt")

toDo$QMcharge= T
QMcharge.inout=paste0(gsub(".pdb",".out",SubQM.outpdb))
QMcharge.outfile=paste0("Mulliken_of_",gsub(".arc","",QMsum.inarc),".txt")

toDo$ESPcharge= F
ESPcharge.inarc=paste0(gsub(".pdb",".arc",SubQM.outpdb))
ESPcharge.outfile=paste0("ESPcharge_of_",gsub(".arc","",QMsum.inarc),".txt")

toDo$mecha= F
mecha.inpdb=QMopt.outpdb
mecha.outfile=paste0("Various_mechanical_parameter_of_",gsub(".pdb","",mecha.inpdb))

toDo$mecha_comp= F
mecha_comp.inpdb=QMopt.outpdb
mecha_comp.outfile=paste0("Various_mechanical_parameter_of_comple_",gsub(".pdb","",mecha.inpdb))

toDo$Rmsd = F
Rmsd.atom = "nonH"
Rmsd.align = TRUE
Rmsd.refpdb1 = Bld.outpdb
Rmsd.refpdb2 = QMopt.outpdb
Rmsd.outfile = paste0("rmsd_",Rmsd.atom,
                     "_align_",Rmsd.align,
                     "_for_",gsub(".pdb", "",Rmsd.refpdb2),".txt")

toDo$DNAshape=F
DNAshape.outfile="DNAshape.txt"

toDo$DNAshape_comp=F
DNAshape_comp.outfile="DNAshape_comp.txt"

# location of output folder
output.path = paste0("data/", conformation)
if (neutralise) {output.path = paste0(output.path, "_neutrized")}
if (solvate) {output.path = paste0(output.path, "_solvated")}


## Dependant functions #########################################################

source("lib/generateKmers.R")
source("lib/buildDNA.R")
source("lib/genTopCoor.R")
source("lib/MMopt.R")
source("lib/QMopt.R")
source("lib/SubQM.R")
source("lib/charge.R")
source("lib/mechanical.R")
source("lib/QMsummary.R")
source("lib/calculateRMSD.R")
#source("lib/getDNAshape.R")
library(bio3d)
#library(DNAshapeR)

## AmberTools path
#ambertools.path <- "/opt/anaconda3/envs/AmberTools21"
#Sys.setenv(PATH = paste0(ambertools.path, "/bin:", Sys.getenv("PATH")))
#Sys.setenv(AMBERHOME = ambertools.path)

## MOPAC path ###################################################################
#mopac.path <- "/Users/kairimasuda/Desktop/MOPAC2016_for_Macintosh/"
#mopac.path <- "/home/k/kamasuda/Desktop/MOPAC2016_for_Linux_64_bit"
mopac.path <-"/home/imm/clab0629/Desktop/MOPAC2016_for_Linux_64_bit"
Sys.setenv(PATH = paste0(Sys.getenv("PATH"),":",mopac.path))
Sys.setenv(LD_LIBRARY_PATH = paste0(Sys.getenv("LD_LIBRARY_PATH"),":",mopac.path))

## curve+ path ###################################################################
#curve.path <- "/Users/kairimasuda/Desktop/curves+/"
curve.path <- "/home/imm/clab0629/Desktop/curves+/"
#mopac.path <- "/home/k/kamasuda/Desktop/MOPAC2016_for_Linux_64_bit"
Sys.setenv(PATH = paste0(Sys.getenv("PATH"),":",curve.path))

## Directory setup #############################################################
dir.create(output.path, recursive = TRUE, showWarnings = FALSE)

## Main Code ###################################################################

#step0 read sequence------------------------------------------------------------
if (sequences[1]=="readfile") {
  file.copy(sequences[2],"temp.txt")
  write("\n","temp.txt", append=TRUE)
  dna.seqs <- readLines("temp.txt")
  dna.seqs <-dna.seqs[dna.seqs !=""]
  invisible(file.remove("temp.txt"))
}else if (sequences[1]=="generateKmer"){
  sequences_arrays=GenerateKmers(as.numeric(sequences[2]))
  filename=paste0(sequences[2],"mer_all.txt")
  file.create(filename)
  write(sequences_arrays,filename,sep="\n",append=T)
  dna.seqs <- readLines(filename)
}

for (dna.seq in dna.seqs) {
output <- paste0(output.path, "/", dna.seq, "/")
dir.create(output, showWarnings = FALSE)
}

#step1 build model--------------------------------------------------------------
if (toDo$initDNA) {
  for (dna.seq in dna.seqs) {
    output <- paste0(output.path, "/", dna.seq, "/")
    buildDNA(conformation, dna.seq, output, Bld.outpdb)
    genTopCoor(paste0(output, Bld.outpdb), neutralise, solvate,
               output.path = output, output.name = "dna")
  }
}

#step2-1 MM optimization--------------------------------------------------------
if (toDo$MMopt) {
  step=0
  for (dna.seq in dna.seqs) {
    step=step+1
    print(paste0(step," ",dna.seq))
    output <- paste0(output.path, "/", dna.seq, "/")
    MMopt(paste0(output, "dna.prmtop"),
             paste0(output, "dna.rst7"),
             maxcyc = MMopt.maxcyc,
             ncyc   = MMopt.ncyc,
             target = MMopt.atoms,
             igb = MMopt.igb,
             output,
             MMopt.outpdb)
  }
}

#step2-2 QM optimization--------------------------------------------------------
if (toDo$QMopt) {
  step=0
  for (dna.seq in dna.seqs) {
    step=step+1
    print(paste0(step,' ',dna.seq))
    
    QMopt(outputfolder = output.path,
          sequence     = dna.seq,
          infile       = QMopt.inpdb,
          outfile      = QMopt.outpdb,
          HF           = QMopt.HF,
          hamiltonian  = QMopt.hami,
          norm         = QMopt.norm,
          COSMO        = QMopt.COSMO,
          charge       = QMopt.charge,
          other        = QMopt.other)
  }
}


#step 3 Subsequent calculation------------------------------------------------
if (toDo$SubQM) {
  step=0
  for (dna.seq in dna.seqs) {
    step=step+1
    print(paste0(step,' ',dna.seq))
    
    SubQM(outputfolder = output.path,
          sequence     = dna.seq,
          infile       = SubQM.inpdb,
          outfile      = SubQM.outpdb,
          inputconf    = SubQM.conf,
          HF           = SubQM.HF,
          hamiltonian  = SubQM.hami,
          COSMO        = SubQM.COSMO,
          charge       = SubQM.charge,
          other        = SubQM.other
    )
  }
}

#step4-1 Make a summary of QM calculation---------------------------------------
if (toDo$QMsummary) {
  file.create(paste0(output.path,"/",QMsum.outfile))
  writeLines(paste0("seq,", 
                    "Ehof,",
                    "Dipole,",
                    "IP,",
                    "Ehomo,",
                    "Elumo"),
             paste0(output.path,"/",QMsum.outfile))
  
  input_arcs=paste0(output.path, "/", dna.seqs, "/", QMsum.inarc)
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
        paste0(output.path,"/",QMsum.outfile),sep="\n",append=T)
}


#step4-2 Make a summary of charge--------------------------------------------
if (toDo$QMcharge) {
  Mulliken_values<-c()
  Mulliken_values=lapply(dna.seqs,extract_Mulliken,
                         output.path=output.path,QMcharge.inout=QMcharge.inout)

  file.create(paste0(output.path,"/",QMcharge.outfile))
  writeLines("seq,res_num,res_name,atom_name,element,charge,population",
             paste0(output.path,"/",QMcharge.outfile), sep = "\n")
  write(unlist(Mulliken_values),paste0(output.path,"/",QMcharge.outfile), append=T)
  Mulliken_values<-c()
}

if (toDo$ESPcharge) {
  ESP_values<-c()
  ESP_values[dna.seqs]=lapply(dna.seqs,extract_ESP,
                              output.path=output.path,ESPcharge.inout=ESPcharge.inarc)
  file.create(paste0(output.path,"/",ESPcharge.outfile))
  writeLines("seq,res_num,res_name,atom_name,element,charge,population",
             paste0(output.path,"/",ESPcharge.outfile), sep = "\n")
  write(unlist(ESP_values),paste0(output.path,"/",ESPcharge.outfile),sep="\n",append=T) 
  ESP_values<-c()
}

#step4-3 Make a summary of mechanical parameters--------------------------------
if (toDo$mecha) {
  lapply(dna.seqs, curve_plus, 
         output.path=output.path, 
         mecha.inpdb=mecha.inpdb,
         mecha.outfile=mecha.outfile,
         curve.path=curve.path)
  file.create(paste0(output.path,"/",mecha.outfile,".txt"))
  
  # make axis pareameter summary  
  axis_param<-c()
  axis_param=lapply(dna.seqs, BPAxis_read, output.path=output.path, mecha.outfile=mecha.outfile)
  
  write(paste0("seq, res_num1, res_name1, res_num2, res_name2, Xdisp, Ydisp, Inclin, Tip, Ax-bend"),
             paste0(output.path,"/",mecha.outfile,".txt"), sep = "\n",append=T)
  write(unlist(axis_param), paste0(output.path,"/",mecha.outfile,".txt"),sep="\n",append=T)    

  # make Intra BP pareameter summary  
  intra_param<-c()
  intra_param=lapply(dna.seqs, intraBP_read, output.path=output.path, mecha.outfile=mecha.outfile)

  write(paste0(""),paste0(output.path,"/",mecha.outfile,".txt"), sep = "\n",append=T)  
  write(paste0("seq, res_num1, res_name1, res_num2, res_name2, Shear, Stretch, Stagger, Buckle, Propel, Opening"),
             paste0(output.path,"/",mecha.outfile,".txt"), sep = "\n",append=T)
  write(unlist(intra_param), paste0(output.path,"/",mecha.outfile,".txt"),sep="\n",append=T)   

  # make Inter BP pareameter summary  
  inter_param<-c()
  inter_param=lapply(dna.seqs,interBP_read, output.path=output.path, mecha.outfile=mecha.outfile)
  
  write(paste0(""),paste0(output.path,"/",mecha.outfile,".txt"), sep = "\n",append=T)  
  write(paste0("seq, res_num1, res_name1, res_num2, res_name2, Shift, Slide, Rise, Tilt, Roll, Twist, H-Ris, H-Twi"),
             paste0(output.path,"/",mecha.outfile,".txt"), sep = "\n",append=T)
  write(unlist(inter_param), paste0(output.path,"/",mecha.outfile,".txt"),sep="\n",append=T)   

  # make backbone pareameter summary
  backbone_param<-c()
  backbone_param=lapply(dna.seqs,backbone_read, output.path=output.path, mecha.outfile=mecha.outfile)
  
  write(paste0(""),paste0(output.path,"/",mecha.outfile,".txt"), sep = "\n",append=T)  
  write(paste0("seq, res_num1, res_name, Alpha, Beta, Gamma, Delta, Epsil, Zeta, Chi, Phase, Ampli, Puckr"),
             paste0(output.path,"/",mecha.outfile,".txt"), sep = "\n",append=T)
  write(unlist(backbone_param), paste0(output.path,"/",mecha.outfile,".txt"),sep="\n",append=T)
  }

#step4-3 Make a summary of mechanical parameters for complementary sequece -----
if (toDo$mecha_comp) {
  lapply(dna.seqs, curve_plus_comp, 
         output.path=output.path,
         mecha.inpdb=mecha_comp.inpdb,
         mecha.outfile=mecha_comp.outfile,
         curve.path=curve.path)
  file.create(paste0(output.path,"/",mecha_comp.outfile,".txt"))
  
  # make axis pareameter summary  
  axis_param<-c()
  axis_param=lapply(dna.seqs, BPAxis_read, output.path=output.path, mecha.outfile=mecha_comp.outfile)
  
  write(paste0("seq, res_num1, res_name1, res_num2, res_name2, Xdisp, Ydisp, Inclin, Tip, Ax-bend"),
        paste0(output.path,"/",mecha_comp.outfile,".txt"), sep = "\n",append=T)
  write(unlist(axis_param), paste0(output.path,"/",mecha_comp.outfile,".txt"),sep="\n",append=T)    

  # make Intra BP pareameter summary  
  intra_param<-c()
  intra_param=lapply(dna.seqs, intraBP_read, output.path=output.path, mecha.outfile=mecha_comp.outfile)
  
  write(paste0(""),paste0(output.path,"/",mecha_comp.outfile,".txt"), sep = "\n",append=T)  
  write(paste0("seq, res_num1, res_name1, res_num2, res_name2, Shear, Stretch, Stagger, Buckle, Propel, Opening"),
        paste0(output.path,"/",mecha_comp.outfile,".txt"), sep = "\n",append=T)
  write(unlist(intra_param), paste0(output.path,"/",mecha_comp.outfile,".txt"),sep="\n",append=T)   

  # make Inter BP pareameter summary  
  inter_param<-c()
  inter_param=lapply(dna.seqs,interBP_read, output.path=output.path, mecha.outfile=mecha_comp.outfile)
  
  write(paste0(""),paste0(output.path,"/",mecha_comp.outfile,".txt"), sep = "\n",append=T)  
  write(paste0("seq, res_num1, res_name1, res_num2, res_name2, Shift, Slide, Rise, Tilt, Roll, Twist, H-Ris, H-Twi"),
        paste0(output.path,"/",mecha_comp.outfile,".txt"), sep = "\n",append=T)
  write(unlist(inter_param), paste0(output.path,"/",mecha_comp.outfile,".txt"),sep="\n",append=T)   

  # make backbone pareameter summary
  backbone_param<-c()
  backbone_param=lapply(dna.seqs, backbone_read, output.path=output.path, mecha.outfile=mecha_comp.outfile)
  
  write(paste0(""),paste0(output.path,"/",mecha_comp.outfile,".txt"), sep = "\n",append=T)  
  write(paste0("seq, res_num1, res_name, Alpha, Beta, Gamma, Delta, Epsil, Zeta, Chi, Phase, Ampli, Puckr"),
        paste0(output.path,"/",mecha_comp.outfile,".txt"), sep = "\n",append=T)
  write(unlist(backbone_param), paste0(output.path,"/",mecha_comp.outfile,".txt"),sep="\n",append=T)
}


#step4-4 Calculation of rmsd----------------------------------------------------
if (toDo$Rmsd) {
  reference_pdb=paste0(output.path,"/",dna.seqs,"/",Rmsd.refpdb1)
  target_pdb=paste0(output.path,"/",dna.seqs,"/",Rmsd.refpdb2)
  rmsd=mapply(rmsd_bio3d, target_pdb, reference_pdb, Rmsd.atom, Rmsd.align)

  file_names=paste0(output.path,"/",Rmsd.outfile)
  file.create(file_names)
  writeLines(paste0("seq, rmsd"), file_names)
  write(paste0(dna.seqs, ", ", rmsd), file_names,sep="\n",append=T)
}

#step4-5 Calculation of DNAshapeR-----------------------------------------------
if (toDo$DNAshape) {
  mecha=vector()
  mecha[dna.seqs]<-lapply(dna.seqs, getDNAshape)
  MGW    <-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["MGW"]),collapse=","))})
  HelT   <-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["HelT"]),collapse=","))})
  ProT   <-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["ProT"]),collapse=","))})
  Roll   <-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Roll"]),collapse=","))})
  EP     <-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["EP"]),collapse=","))})
  Stretch<-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Stretch"]),collapse=","))})
  Tilt   <-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Tilt"]),collapse=","))})
  Buckle <-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Buckle"]),collapse=","))})
  Shear  <-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Shear"]),collapse=","))})
  Opening<-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Opening"]),collapse=","))})
  Rise   <-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Rise"]),collapse=","))})
  Shift  <-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Shift"]),collapse=","))})
  Stagger<-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Stagger"]),collapse=","))})
  Slide  <-lapply(dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Slide"]),collapse=","))})
  
  file_name=paste0(output.path,"/",DNAshape.outfile)
  max_res_num=max(nchar(dna.seqs))
  writeLines(paste0("seq,type,",paste(1:max_res_num,collapse = ",")), file_name, sep = "\n")
  write(paste0(dna.seqs,",MGW,",unlist(MGW)),file_name,sep = "\n",append=T)
  write(paste0(dna.seqs,",ProT,",unlist(ProT)),file_name,sep = "\n",append=T)
  write(paste0(dna.seqs,",Stretch,",unlist(Stretch)),file_name,sep = "\n",append=T)
  write(paste0(dna.seqs,",Buckle,",unlist(Buckle)),file_name,sep = "\n",append=T)
  write(paste0(dna.seqs,",Shear,",unlist(Shear)),file_name,sep = "\n",append=T)
  write(paste0(dna.seqs,",Opening,",unlist(Opening)),file_name,sep = "\n",append=T)
  write(paste0(dna.seqs,",Stagger,",unlist(Stagger)),file_name,sep = "\n",append=T)

  write(paste0(""),file_name, sep = "\n",append=T)  
  write(paste0("seq,type,",paste(1:(max_res_num-1),collapse = ",")), file_name, sep = "\n",append=T)
  write(paste0(dna.seqs,",HelT,",unlist(HelT)),file_name,sep = "\n",append=T)
  write(paste0(dna.seqs,",Roll,",unlist(Roll)),file_name,sep = "\n",append=T)
  write(paste0(dna.seqs,",Tilt,",unlist(Tilt)),file_name,sep = "\n",append=T)
  write(paste0(dna.seqs,",Rise,",unlist(Rise)),file_name,sep = "\n",append=T)
  write(paste0(dna.seqs,",Shift,",unlist(Shift)),file_name,sep = "\n",append=T)  
  write(paste0(dna.seqs,",Slide,",unlist(Slide)),file_name,sep = "\n",append=T)
  
  delfiles <- dir(path="./",pattern="temp.fa")
  file.remove(delfiles)
}


if (toDo$DNAshape_comp) {
  generate_complement<-function(seq){
    print(seq)
    len = nchar(seq)
    seq=unlist(strsplit(seq,""))
    wc_seq<-rep(NA,len)
    pair<-function(x){
      if(x=="A"){return("T")}
      if(x=="T"){return("A")}
      if(x=="G"){return("C")}
      if(x=="C"){return("G")}
    }
    wc_seq=unlist(lapply(seq,pair))
    wc_seq=rev(wc_seq)
    wc_seq=paste(wc_seq,collapse="")
    return(wc_seq)
  }
  comp_dna.seqs=unlist(lapply(dna.seqs,generate_complement))
  mecha=vector()
  mecha[comp_dna.seqs]<-lapply(comp_dna.seqs, getDNAshape)
  MGW    <-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["MGW"]),collapse=","))})
  HelT   <-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["HelT"]),collapse=","))})
  ProT   <-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["ProT"]),collapse=","))})
  Roll   <-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Roll"]),collapse=","))})
  EP     <-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["EP"]),collapse=","))})
  Stretch<-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Stretch"]),collapse=","))})
  Tilt   <-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Tilt"]),collapse=","))})
  Buckle <-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Buckle"]),collapse=","))})
  Shear  <-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Shear"]),collapse=","))})
  Opening<-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Opening"]),collapse=","))})
  Rise   <-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Rise"]),collapse=","))})
  Shift  <-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Shift"]),collapse=","))})
  Stagger<-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Stagger"]),collapse=","))})
  Slide  <-lapply(comp_dna.seqs, function(x,aaa=mecha){return(paste(unlist(aaa[[x]]["Slide"]),collapse=","))})
  
  file_name=paste0(output.path,"/",DNAshape_comp.outfile)
  max_res_num=max(nchar(comp_dna.seqs))
  writeLines(paste0("seq,type,",paste(1:max_res_num,collapse = ",")), file_name, sep = "\n")
  write(paste0(comp_dna.seqs,",MGW,",unlist(MGW)),file_name,sep = "\n",append=T)
  write(paste0(comp_dna.seqs,",ProT,",unlist(ProT)),file_name,sep = "\n",append=T)
  write(paste0(comp_dna.seqs,",Stretch,",unlist(Stretch)),file_name,sep = "\n",append=T)
  write(paste0(comp_dna.seqs,",Buckle,",unlist(Buckle)),file_name,sep = "\n",append=T)
  write(paste0(comp_dna.seqs,",Shear,",unlist(Shear)),file_name,sep = "\n",append=T)
  write(paste0(comp_dna.seqs,",Opening,",unlist(Opening)),file_name,sep = "\n",append=T)
  write(paste0(comp_dna.seqs,",Stagger,",unlist(Stagger)),file_name,sep = "\n",append=T)
  
  write(paste0(""),file_name, sep = "\n",append=T)  
  write(paste0("seq,type,",paste(1:(max_res_num-1),collapse = ",")), file_name, sep = "\n",append=T)
  write(paste0(comp_dna.seqs,",HelT,",unlist(HelT)),file_name,sep = "\n",append=T)
  write(paste0(comp_dna.seqs,",Roll,",unlist(Roll)),file_name,sep = "\n",append=T)
  write(paste0(comp_dna.seqs,",Tilt,",unlist(Tilt)),file_name,sep = "\n",append=T)
  write(paste0(comp_dna.seqs,",Rise,",unlist(Rise)),file_name,sep = "\n",append=T)
  write(paste0(comp_dna.seqs,",Shift,",unlist(Shift)),file_name,sep = "\n",append=T)  
  write(paste0(comp_dna.seqs,",Slide,",unlist(Slide)),file_name,sep = "\n",append=T)
  
  delfiles <- dir(path="./",pattern="temp.fa")
  file.remove(delfiles)
}


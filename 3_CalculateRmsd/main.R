############################################################
### Input: reference_pdb=[ref1.pdb,ref2.pdb,...],        ###
###        target_pdb=[target1.pdb, target2.pdb,...]     ###
############################################################

data_folder="/Users/kairi/Desktop/2021_12_17_DNAPhysBasisUV/data"
input.seqs = "/Users/kairi/Desktop/2021_12_17_Generate_Kmer/sequences.txt"

conformation="DNA_B_Handbook1999_duplex"
dna.seqs <- readLines(input.seqs)
condtion="min_igb_1_step_500_optimize_all_atoms"
rmsd_target="other_than_H"

reference_pdb=paste0(data_folder,"/",conformation,"/",dna.seqs,"/",dna.seqs,".pdb")
target_pdb=paste0(data_folder,"/",conformation,"/",dna.seqs,"/",condtion,".pdb")

############################################################
### Calculate rmsd of (ref1,target1), (ref2,target2),... ###
############################################################
source("lib/calculate_rmsd.R")
rmsd=mapply(rmsd_bio3d, target_pdb, reference_pdb, rmsd_target)

file.create("rmsd.txt")
writeLines(paste0("file, rmsd_",rmsd_target) ,"./rmsd.txt")
write(paste0(names(rmsd), ",", rmsd),"rmsd.txt",sep="\n",append=T)

############################################################
### Visualization                                        ###
############################################################
barplot(main=paste0(conformation,"\n",condtion),
        xlab="sequence",
        ylab=paste0("rmsd (",rmsd_target,")"),
        rmsd,
        names.arg=dna.seqs,
        ylim=c(0,1.0),
        las=2
)






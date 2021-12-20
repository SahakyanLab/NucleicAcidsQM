# DNAPhysBasisUV

Research and development code to unravel the underlying physical chemistry of variable UV-damage sensitivity of DNA.


#  Automation in R console
file.create("temp.R")
test <- readLines("main.R")
writeLines(test[25:91],"temp.R")

for (i in list(500, 1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,"FULL")) {
  for (j in list(0,1)) {
    for (k in list("only_H", "other_than_H", "all_atoms")){
    toDo = NULL
    toDo$initDNA = T
    toDo$optDNA  = T
    if (i=="FULL"){
    opt.maxcyc = i
    opt.ncyc = NULL
    } else{
    opt.maxcyc = i
    opt.ncyc = i/4
    }
    opt.igb = j
    opt.atoms = k    
    source("temp.R")
    }
  }
}

file.remove("temp.R")

# Calculate rmsd

rmsd_conf<-"DNA_B_Handbook1999_duplex"
rmsd_seq  <- "atgcagta"
rmsd_igbs  <- c(0,1)
rmsd_steps <- c("500", "1000","2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000","FULL")
rmsd_opts<- c("only_H", "other_than_H", "all_atoms" )

file.create(paste0("./data/rmsd_",rmsd_conf,".txt"))
writeLines("seq, igb, step, optimize, rmsd",paste0("./data/rmsd_",rmsd_conf,".txt"))
 
for (rmsd_igb in rmsd_igbs) {
  for (rmsd_step in rmsd_steps) {
    for (rmsd_opt in rmsd_opts){
      condition<-paste0("min","_igb_",rmsd_igb,"_step_",rmsd_step,"_optimize_",rmsd_opt)
      nab.in <- paste0(
               'molecule mol1, mol2;
                float r;
                file out;
                mol1=getpdb("./',rmsd_conf,'/',rmsd_seq,'/',rmsd_seq,'.pdb");
                mol2=getpdb("./',rmsd_conf,'/',rmsd_seq,'/', condition,'.pdb");
                rmsd(mol1, "::",mol2,"::",r);
                out=fopen("./rmsd_',rmsd_conf,'.txt","a");
                fprintf(out,"',rmsd_seq,",",rmsd_igb,",",rmsd_step,",",rmsd_opt,",",' %f\\n",r);'
                )
    
      writeLines(nab.in, "./data/rmsd.nab")
      system(paste0("cd ./data; nab rmsd.nab && ./rmsd"))
    }
  }
}


# visualization
rmsd_conf<-"DNA_B_Handbook1999_duplex"
my_data<-read.delim(paste0("./data/rmsd_",rmsd_conf,".txt"),
                    sep=",",header = T)

par(mfrow = c(1,2),oma=c(0,0,2,0))
vac=subset(my_data, igb==0)
counts_vac<-xtabs(vac$rmsd~vac$optimize + vac$step,)
counts_vac<-counts_vac[,c(500,1000,10000,"FULL")]
counts_vac<-counts_vac[c("only_H","other_than_H","all_atoms"),]

barplot(main="vacuum",
        xlab="step",
        ylab="rmsd",
        counts_vac,
        beside=T,
        ylim=c(0,10.0),
        col=c("green","red","blue")
        )

legend("topleft",
       legend=c("opt only H","opt other than H","opt all"),
        col=c("green","red","blue"),
        lty=1,cex=0.8,lwd=5,bty="n")

sol=subset(my_data, igb==1)
counts_sol<-xtabs(sol$rmsd~sol$optimize + sol$step,)
counts_sol<-counts_sol[,c(500,1000,10000,"FULL")]
counts_sol<-counts_sol[c("only_H","other_than_H","all_atoms"),]

barplot(main="implicit solvent",
        counts_sol,
        beside=T,
        ylim=c(0,10.0),
        col=c("green","red","blue")
        )

mtext(rmsd_conf,outer=TRUE,cex=1)

dev.copy(png,filename=paste0("./data/rmsd_",rmsd_conf,".png"),width=600, height=400)
dev.off()


### visualize in detail
rmsd_conf<-"DNA_B_Handbook1999_duplex"
par(mfrow = c(1,2),oma=c(0,0,2,0))
my_data<-read.delim(paste0("./data/rmsd_",rmsd_conf,".txt"),
                    sep=",",header = T)

vac=subset(my_data, igb==0&step!="FULL")
counts_vac<-xtabs(vac$rmsd~vac$optimize + vac$step,)
counts_vac<-counts_vac[,c("500","1000","2000","3000","4000","5000","6000","7000","8000","9000","10000")]
counts_vac<-counts_vac[c("only_H","other_than_H","all_atoms"),]

barplot(main="vacuum",
        xlab="step",
        ylab="rmsd",
        counts_vac,
        beside=T,
        ylim=c(0,5.0),
        col=c("green","red","blue")
        )

legend("topleft",
       legend=c("opt only H","opt other than H","opt all"),
        col=c("green","red","blue"),
        lty=1,cex=0.8,lwd=5,bty="n")

sol=subset(my_data, igb==1&step!="FULL")
counts_sol<-xtabs(sol$rmsd~sol$optimize + sol$step,)
counts_sol<-counts_sol[,c("500","1000","2000","3000","4000","5000","6000","7000","8000","9000","10000")]
counts_sol<-counts_sol[c("only_H","other_than_H","all_atoms"),]

barplot(main="implicit solvent",
        counts_sol,
        beside=T,
        ylim=c(0,5.0),
        col=c("green","red","blue")
        )
        
mtext(rmsd_conf,outer=TRUE,cex=1)

dev.copy(png,filename=paste0("./data/rmsd_detailed_",rmsd_conf,".png"),width=600, height=400)
dev.off()


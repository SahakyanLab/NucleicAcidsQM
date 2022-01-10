extract_HF<-function(input_data){
  summarys<-readLines(input_data)
  HFline=grep("HEAT OF FORMATION       =", summarys)
  HFvalue=gsub(" ","",substr(summarys[HFline],36,53))   
  return(HFvalue)
}

extract_DI<-function(input_data){
  summarys<-readLines(input_data)
  DIline=grep("DIPOLE                  =", summarys)
  DIvalue=gsub(" ","",substr(summarys[DIline],36,53))
  return(DIvalue)
  
}

extract_IP<-function(input_data){
  summarys<-readLines(input_data)
  IPline=grep("IONIZATION POTENTIAL    =", summarys)
  IPvalue=gsub(" ","",substr(summarys[IPline],36,54)) 
  return(IPvalue)
  
}

extract_HOMO<-function(input_data){
  summarys<-readLines(input_data)
  HLline=grep("HOMO LUMO ENERGIES \\(EV\\) =", summarys)
  HOvalue=gsub(" ","",substr(summarys[HLline],36,51))
  return(HOvalue)
}

extract_LUMO<-function(input_data){
  summarys<-readLines(input_data)
  HLline=grep("HOMO LUMO ENERGIES \\(EV\\) =", summarys)
  LUvalue=gsub(" ","",substring(summarys[HLline],52))
  return(LUvalue)
}












# if (toDo$QMsummary) {
#   file_names=paste0("summary_of_",summary.target,".txt")
#   file.create(paste0(outputfolder,"/",file_names))
#   writeLines(paste0("SEQUENCE, 
#                     HEAT OF FORMATION, 
#                     DIPOLE, 
#                     IONIZATION POTENTIAL, 
#                     HOMO, 
#                     LUMO"), 
#              file_names)  
#   
#   dna.seqs <- readLines(inputseqs)
#   
#   for (dna.seq in dna.seqs) {
#     summarys<-readLines(paste0(inputfolder, "/", dna.seq, "/QM.arc"))
#       HFline=grep("HEAT OF FORMATION       =", summarys)
#       HFvalue=gsub(" ","",substr(summarys[HFline],36,53))
#       
#       DIline=grep("DIPOLE                  =", summarys)
#       DIvalue=gsub(" ","",substr(summarys[DIline],36,53))
#       
#       IPline=grep("IONIZATION POTENTIAL    =", summarys)
#       IPvalue=gsub(" ","",substr(summarys[IPline],36,54))
#       
#       HLline=grep("HOMO LUMO ENERGIES \\(EV\\) =", summarys)
#       HOvalue=gsub(" ","",substr(summarys[HLline],36,51))
#       LUvalue=gsub(" ","",substring(summarys[HLline],52))
#       
#       write(paste0(dna.seq, ", ", 
#                    HFvalue, ", ",
#                    DIvalue, ", ", 
#                    IPvalue, ", ", 
#                    HOvalue, ", ", 
#                    LUvalue), 
#             file_names,sep="\n",append=T)  
#       
#   }
# }


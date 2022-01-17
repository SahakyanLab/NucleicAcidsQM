QMminimise <- function(inputfolder,sequence,inputfile,
                       hamiltonian,calculation,COSMO,charge) {
  
  
  output_name=paste0('QM_',hamiltonian,"_",
                           calculation,"_",
             if(COSMO != F) {paste0("COSMO_")},
                          "from_",
                           gsub(".pdb","",inputfile)
  )

  if(charge=="twice_of_sequence"){char_num=-2*nchar(sequence)+2}
  
  writeLines(paste0(
    'RHF MMOK PDBOUT GRAPHF THREADS=4', ' ',
    hamiltonian, ' ',
   if(calculation=="OPT") {paste0('GNORM=8.0', ' ')},
   if(calculation=="1SCF") {paste0('1SCF', ' ')},
   'CHARGE=',char_num, ' ',
   if(COSMO != F) {paste0(COSMO, ' ')},
   'GEO_DAT="',inputfile,'"' 
     ), paste0(inputfolder,'/',sequence,'/',output_name,'.mop'))
      
  cmd <- paste0("cd ",inputfolder,'/',sequence,";",
                "MOPAC2016.exe ", output_name,'.mop')
  system(cmd)

}




# # write geometry data
# temp <- readLines(paste0(inputfolder,'/',sequence,'/',inputfile))
# temp<-temp[- grep("TER", temp)]
# temp<-temp[- grep("END", temp)]  
# for(i in grep("ATOM      1",temp):(length(temp))){
#   temp_vec=unlist(strsplit(temp[i], "[ ]+"))
#   line= paste0(temp_vec[11]," ",temp_vec[6]," 1 ",temp_vec[7]," 1 ",temp_vec[8]," 1 ")
#   write(line,paste0(inputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)
# }



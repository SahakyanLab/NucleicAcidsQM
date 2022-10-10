QMopt <- function(outputfolder, sequence, infile, outfile,
                  HF,hamiltonian,norm,COSMO,charge,other) {
  
  
  output_name=gsub(".pdb","",outfile)

  if(charge=="number_of_phosphate"){char_num=-2*nchar(sequence)+2}
  if(!is.na(suppressWarnings(as.numeric(charge)))){char_num=as.numeric(charge)}
  
  writeLines(paste0(
    HF, ' ',
    hamiltonian, ' ',
    'GNORM=',norm, ' ',
    'CHARGE=',char_num, ' ',
    if(COSMO != F) {paste0(COSMO, ' ')},
    other, ' ',
    'PDB PDBOUT','\n','\n'
     ), paste0(outputfolder,'/',sequence,'/',output_name,'.mop'))
      
  temp <- readLines(paste0(outputfolder,'/',sequence,'/',infile))
  if(length(grep("TER", temp))>0){temp<-temp[- grep("TER", temp)]}
  if(length(grep("END", temp))>0){temp<-temp[- grep("END", temp)]}

  for(i in grep("ATOM      1",temp):(length(temp))){
    line= temp[i]
    write(line,paste0(outputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)
  }
  
  cmd <- paste0("cd ",outputfolder,'/',sequence,";",
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



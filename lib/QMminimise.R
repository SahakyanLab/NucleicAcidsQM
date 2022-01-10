QMminimise <- function(inputfolder,sequence,inputfile,
                       hamiltonian,calculation,COSMO,charge) {
  
  if(charge=="twice_of_sequence"){char_num=-2*nchar(sequence)}
  
  
  output_name=paste0('QM_',hamiltonian,"_",
                           calculation,"_",
             if(COSMO==T) {paste0("COSMO_")},
                          "from_",
                           gsub(".pdb","",inputfile)
  )
  
  writeLines(paste0(
    'RHF', ' ',
    'MMOK', ' ',
    'PDBOUT', ' ',
    'GRAPHF', ' ',
    'THREADS=4', ' ',
    hamiltonian, ' ',
    calculation, ' ',
   'CHARGE=',char_num, ' ',
   if(COSMO==T) {paste0('RSOLV=1.3 EPS=78.4 NSPA=42 ')},
    'GEO_DAT="',inputfile,'"' 
     ), paste0(inputfolder,'/',sequence,'/',output_name,'.mop'))
  
  cmd <- paste0("cd ",inputfolder,'/',sequence,";",
                "MOPAC2016.exe ", output_name,'.mop')
  system("cd")
  system(cmd)

}



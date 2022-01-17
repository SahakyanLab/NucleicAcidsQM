SubQMcal <- function(inputfolder,sequence,inputfile,inputconf,
                       hamiltonian,COSMO,charge){
  
  output_name=paste0('SubQM_',inputconf,
                     "_from_",gsub(".pdb","",inputfile))
  
  #-----------------------------------------------------------------------------
  if(inputconf=="duplex_1scf"){
    
    if(charge=="depend_on_conf"){char_num=-2*nchar(sequence)+2}
    
  writeLines(paste0(
    'RHF MMOK PDB PDBOUT GRAPHF THREADS=4 1SCF ',
    hamiltonian, ' ',
   'CHARGE=',char_num, ' ',
   if(COSMO != F) {paste0(COSMO, ' ')},
    'GEO_DAT="',inputfile,'"' 
     ), paste0(inputfolder,'/',sequence,'/',output_name,'.mop'))
  }
  
  #-----------------------------------------------------------------------------
  if(inputconf=="1st_strand_1scf"){
    
    if(charge=="depend_on_conf"){char_num=-nchar(sequence)+1}
    
    writeLines(paste0(
      'RHF MMOK XYZ PDBOUT GRAPHF THREADS=4 1SCF', ' ',
      hamiltonian, ' ',
      'CHARGE=',char_num, ' ',
      if(COSMO != F) {paste0(COSMO, ' ')},'\n','\n'
    ), paste0(inputfolder,'/',sequence,'/',output_name,'.mop'))
    
    temp <- readLines(paste0(inputfolder,'/',sequence,'/',inputfile))
    temp<-temp[- grep("TER", temp)]
    temp<-temp[- grep("END", temp)]
    
    for(i in grep("ATOM      1",temp):(length(temp))){
      temp_vec=unlist(strsplit(temp[i], "[ ]+"))
      if(as.numeric(temp_vec[6])<=nchar(sequence)){
        line= paste0(temp_vec[13]," ",temp_vec[7]," 0 ",temp_vec[8]," 0 ",temp_vec[9]," 0 ")
        write(line,paste0(inputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)
      }
    }
  }
  
  #-----------------------------------------------------------------------------
  if(inputconf=="2nd_strand_1scf"){
    
    if(charge=="depend_on_conf"){char_num=-nchar(sequence)+1}
    
    writeLines(paste0(
      'RHF MMOK XYZ PDBOUT GRAPHF THREADS=4 1SCF', ' ',
      hamiltonian, ' ',
      'CHARGE=',char_num, ' ',
      if(COSMO != F) {paste0(COSMO, ' ')},'\n','\n'
    ), paste0(inputfolder,'/',sequence,'/',output_name,'.mop'))
    
    temp <- readLines(paste0(inputfolder,'/',sequence,'/',inputfile))
    temp<-temp[- grep("TER", temp)]
    temp<-temp[- grep("END", temp)]  
    
    for(i in grep("ATOM      1",temp):(length(temp))){
      temp_vec=unlist(strsplit(temp[i], "[ ]+"))
      if(as.numeric(temp_vec[6])>nchar(sequence)){
        line= paste0(temp_vec[13]," ",temp_vec[7]," 0 ",temp_vec[8]," 0 ",temp_vec[9]," 0 ")
        write(line,paste0(inputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)
      }
    }
  }
  
  #-----------------------------------------------------------------------------
  if(inputconf=="del_mid_1st_strand_opt"){
    
    if(charge=="depend_on_conf"){char_num=-2*nchar(sequence)+2}
    
    writeLines(paste0(
      'RHF MMOK XYZ PDBOUT GRAPHF THREADS=4 GNORM=8.0', ' ',
      hamiltonian, ' ',
      'CHARGE=',char_num, ' ',
      if(COSMO != F) {paste0(COSMO, ' ')},'\n','\n'
    ), paste0(inputfolder,'/',sequence,'/',output_name,'.mop'))
    
    # extract a base located on the mid of the first strand. 
    # After extraction, C1' is capped by hydrogen.
  
    temp <- readLines(paste0(inputfolder,'/',sequence,'/',inputfile))
    temp<-temp[- grep("TER", temp)]
    temp<-temp[- grep("END", temp)] 
    
    cen1<-c(0,0,0) # A position of C1'
    cen2<-c(0,0,0) # A position of N1 or N9 of base. 
    cen3<-c(0,0,0) # A vector N1-C1' or N9-C1'. Hydrogen will be put at C1'+(N1-C1')/2 or C1'+(N9-C1')/2
    
    #For example, temp[1]=("ATOM      1 HO5' DA5 A   1      -2.734   2.425   8.612  1.00  0.35      PROT H")
    #Thus, temp_vec=("ATOM","1","HO5'","DA5","A","1","-2.734","2.425","8.612","1.00","0.35","PROT","H")    
    #Especially, temp_vec[3]=atom type, temp_vec[4]=residue type, temp_vec[6]=residue number,
    #            temp_vec[7-9]=atom coordinates, temp_vec[13]=atom class.
    
    for(i in grep("ATOM      1",temp):(length(temp))){
      temp_vec=unlist(strsplit(temp[i], "[ ]+"))
      if(as.numeric(temp_vec[6])==round(nchar(sequence)/2)){ #When in the mid residue
        if( grepl("\'", temp_vec[3])){ # backbone atoms other than P,OP1,OP2 include "'" such as C1', HO5'
          line= paste0(temp_vec[13]," ",temp_vec[7]," 0 ",temp_vec[8]," 0 ",temp_vec[9]," 0 ")
          write(line,paste0(inputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)     
          if(grepl("C1\'", temp_vec[3])){      
            cen1<-c(as.numeric(temp_vec[7]),as.numeric(temp_vec[8]),as.numeric(temp_vec[9]))
          } 
        } else { # base atoms do not include "'" such as N1 and N9
          if((temp_vec[3]=="N9" && grepl("DA", temp_vec[4]))
            ||(temp_vec[3]=="N9" && grepl("DG", temp_vec[4]))
            ||(temp_vec[3]=="N1" && grepl("DC", temp_vec[4]))
            ||(temp_vec[3]=="N1" && grepl("DT", temp_vec[4]))){   
            cen2<-c(as.numeric(temp_vec[7]),as.numeric(temp_vec[8]),as.numeric(temp_vec[9]))
            cen3<-cen1+(cen2-cen1)/2
            line= paste0("H"," ",cen3[1]," 1 ",cen3[2]," 1 ",cen3[3]," 1 ")
            write(line,paste0(inputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)   
          } else if( temp_vec[3]=="P" || temp_vec[3]=="OP1" || temp_vec[3]=="OP2"){
            line= paste0(temp_vec[13]," ",temp_vec[7]," 0 ",temp_vec[8]," 0 ",temp_vec[9]," 0 ")
            write(line,paste0(inputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)    
          }
        }
      } else { #When not in the mid residue
      line= paste0(temp_vec[13]," ",temp_vec[7]," 0 ",temp_vec[8]," 0 ",temp_vec[9]," 0 ")
      write(line,paste0(inputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)
      }
    }
  }

  #-----------------------------------------------------------------------------
  if(inputconf=="del_mid_2nd_strand_opt"){
    
    if(charge=="depend_on_conf"){char_num=-2*nchar(sequence)+2}
    
    writeLines(paste0(
      'RHF MMOK XYZ PDBOUT GRAPHF THREADS=4 GNORM=8.0', ' ',
      hamiltonian, ' ',
      'CHARGE=',char_num, ' ',
      if(COSMO != F) {paste0(COSMO, ' ')},'\n','\n'
    ), paste0(inputfolder,'/',sequence,'/',output_name,'.mop'))
    
    temp <- readLines(paste0(inputfolder,'/',sequence,'/',inputfile))
    temp<-temp[- grep("TER", temp)]
    temp<-temp[- grep("END", temp)]  

    cen1<-c(0,0,0)
    cen2<-c(0,0,0)
    cen3<-c(0,0,0)
    
    for(i in grep("ATOM      1",temp):(length(temp))){
      temp_vec=unlist(strsplit(temp[i], "[ ]+"))
      if(as.numeric(temp_vec[6])==(nchar(sequence)+round(nchar(sequence)/2))){ # the resiude number of mid of the second strand = 1st strand + 2nd strand/2  
        if( grepl("\'", temp_vec[3])){
          line= paste0(temp_vec[13]," ",temp_vec[7]," 0 ",temp_vec[8]," 0 ",temp_vec[9]," 0 ")
          write(line,paste0(inputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)   
          if(grepl("C1\'", temp_vec[3])){      
            cen1<-c(as.numeric(temp_vec[7]),as.numeric(temp_vec[8]),as.numeric(temp_vec[9]))
          }
        } else {
          if((temp_vec[3]=="N9" && grepl("DA", temp_vec[4]))
             ||(temp_vec[3]=="N9" && grepl("DG", temp_vec[4]))
             ||(temp_vec[3]=="N1" && grepl("DC", temp_vec[4]))
             ||(temp_vec[3]=="N1" && grepl("DT", temp_vec[4]))){   
            cen2<-c(as.numeric(temp_vec[7]),as.numeric(temp_vec[8]),as.numeric(temp_vec[9]))
            cen3<-cen1+(cen2-cen1)/2
            line= paste0("H"," ",cen3[1]," 1 ",cen3[2]," 1 ",cen3[3]," 1 ")
            write(line,paste0(inputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)   
          }else if( temp_vec[3]=="P" || temp_vec[3]=="OP1" || temp_vec[3]=="OP2"){
            line= paste0(temp_vec[13]," ",temp_vec[7]," 0 ",temp_vec[8]," 0 ",temp_vec[9]," 0 ")
            write(line,paste0(inputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)    
          }
        }
      } else {
        line= paste0(temp_vec[13]," ",temp_vec[7]," 0 ",temp_vec[8]," 0 ",temp_vec[9]," 0 ")
        write(line,paste0(inputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)
      }
    }
  }  
  
  #-------------------------------------------------------------------------------
  
  cmd <- paste0("cd ",inputfolder,'/',sequence,";",
                "MOPAC2016.exe ", output_name,'.mop')
  system(cmd)

}



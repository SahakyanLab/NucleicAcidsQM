SubQM <- function(outputfolder,sequence,infile, outfile, inputconf,
                  HF,hamiltonian,COSMO,charge,other){
  
  output_name=gsub(".pdb","",outfile)
  
  #-----------------------------------------------------------------------------
  if(inputconf=="duplex_sp"){
    
    if(charge=="number_of_phosphate"){char_num=-2*nchar(sequence)+2}
    if(!is.na(suppressWarnings(as.numeric(charge)))){char_num=as.numeric(charge)}
    
    writeLines(paste0(
      HF, ' ',
      hamiltonian, ' ',
      '1SCF', ' ',
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
    
  }
  
  #-----------------------------------------------------------------------------
  if(inputconf=="1st_strand_sp"){
    
    if(charge=="number_of_phosphate"){char_num=-nchar(sequence)+1}
    if(!is.na(suppressWarnings(as.numeric(charge)))){char_num=as.numeric(charge)}
    
    writeLines(paste0(
      HF, ' ',
      hamiltonian, ' ',
      '1SCF', ' ',
      'CHARGE=',char_num, ' ',
      if(COSMO != F) {paste0(COSMO, ' ')},
      other, ' ',
      'PDB PDBOUT','\n','\n'
    ), paste0(outputfolder,'/',sequence,'/',output_name,'.mop'))
    
    temp <- readLines(paste0(outputfolder,'/',sequence,'/',infile))
    if(length(grep("TER", temp))>0){temp<-temp[- grep("TER", temp)]}
    if(length(grep("END", temp))>0){temp<-temp[- grep("END", temp)]}
    
    for(i in grep("ATOM      1",temp):(length(temp))){
      ResNum<-substr(temp[i],23,26)
      if(as.numeric(ResNum)<=nchar(sequence)){
        line= temp[i]
        write(line,paste0(outputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)
      }
    }
  }
  
  # temp<-"ATOM      1 HO5' DA5 A   1      -2.734   2.425   8.612  1.00  0.35      PROT H"
  # Header<-substr(temp,1,4)
  # AtomID<-substr(temp,7,11)
  # AtomName<-substr(temp,13,16)
  # AlterLoc<-substr(temp,17,17)
  # ResName<-substr(temp,18,20)
  # Chain<-substr(temp,22,22)
  # ResNum<-substr(temp,23,26)
  # InserRes<-substr(temp,27,27)
  # PosX<-substr(temp,31,38)
  # PosY<-substr(temp,39,46)
  # PosZ<-substr(temp,47,54)
  # Occu<-substr(temp,55,60)
  # TempFac<-substr(temp,61,66)
  # SegIden<-substr(temp,73,76)
  # Element<-substr(temp,77,78)
  # paste0(Header,strrep(" ",2),AtomID,strrep(" ",1),AtomName,AlterLoc,ResName,strrep(" ",1),
  #        Chain,ResNum,InserRes,strrep(" ",3),PosX,PosY,PosZ,Occu,TempFac,strrep(" ",6),SegIden,Element)
  
  
  
  #-----------------------------------------------------------------------------
  if(inputconf=="2nd_strand_sp"){
    
    if(charge=="number_of_phosphate"){char_num=-nchar(sequence)+1}
    if(!is.na(suppressWarnings(as.numeric(charge)))){char_num=as.numeric(charge)}
    
    writeLines(paste0(
      HF, ' ',
      hamiltonian, ' ',
      '1SCF', ' ',
      'CHARGE=',char_num, ' ',
      if(COSMO != F) {paste0(COSMO, ' ')},
      other, ' ',
      'PDB PDBOUT','\n','\n'
    ), paste0(outputfolder,'/',sequence,'/',output_name,'.mop'))
    
    temp <- readLines(paste0(outputfolder,'/',sequence,'/',infile))
    if(length(grep("TER", temp))>0){temp<-temp[- grep("TER", temp)]}
    if(length(grep("END", temp))>0){temp<-temp[- grep("END", temp)]}
    
    for(i in grep("ATOM      1",temp):(length(temp))){
      ResNum<-substr(temp[i],23,26)
      if(as.numeric(ResNum)>nchar(sequence)){
        line= temp[i]
        write(line,paste0(outputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)
      }
    }
  }
  
  #-----------------------------------------------------------------------------
  if(inputconf=="del_mid_1st_strand_sp"){
    
    if(charge=="number_of_phosphate"){char_num=-2*nchar(sequence)+2}
    if(!is.na(suppressWarnings(as.numeric(charge)))){char_num=as.numeric(charge)}
    
    writeLines(paste0(
      HF, ' ',
      hamiltonian, ' ',
      'GNORM=1.0', ' ', # the deleted base is replaced by H. Only this H is optimized.
      'CHARGE=',char_num, ' ',
      if(COSMO != F) {paste0(COSMO, ' ')},
      other, ' ',
      'PDBOUT','\n','\n'
    ), paste0(outputfolder,'/',sequence,'/',output_name,'.mop'))
    
    # extract a base located on the mid of the first strand. 
    # After extraction, C1' is capped by hydrogen.
    
    temp <- readLines(paste0(outputfolder,'/',sequence,'/',infile))
    if(length(grep("TER", temp))>0){temp<-temp[- grep("TER", temp)]}
    if(length(grep("END", temp))>0){temp<-temp[- grep("END", temp)]}
    
    cen1<-c(0,0,0) # A position of C1'
    cen2<-c(0,0,0) # A position of N1 or N9 of base. 
    cen3<-c(0,0,0) # A vector N1-C1' or N9-C1'. Hydrogen will be put at C1'+(N1-C1')/3 or C1'+(N9-C1')/3
    
    #For example, temp[1]=("ATOM      1 HO5' DA5 A   1      -2.734   2.425   8.612  1.00  0.35      PROT H")
    
    for(i in grep("ATOM      1",temp):(length(temp))){#For format, please refer to http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
      Header<-substr(temp[i],1,4)
      AtomID<-substr(temp[i],7,11)
      AtomName<-substr(temp[i],13,16)
      AlterLoc<-substr(temp[i],17,17)
      ResName<-substr(temp[i],18,20)
      Chain<-substr(temp[i],22,22)
      ResNum<-substr(temp[i],23,26)
      InserRes<-substr(temp[i],27,27)
      PosX<-substr(temp[i],31,38)
      PosY<-substr(temp[i],39,46)
      PosZ<-substr(temp[i],47,54)
      Occu<-substr(temp[i],55,60)
      TempFac<-substr(temp[i],61,66)
      SegIden<-substr(temp[i],73,76)
      Element<-substr(temp[i],77,78)
      if(as.numeric(ResNum)==round(nchar(sequence)/2)){ #When in the mid residue
        if( grepl("\'", AtomName)){ # backbone atoms other than P,OP1,OP2 include "'" such as C1', HO5'
          line= paste0(Element,"(ATOM  ",AtomID," ",AtomName,AlterLoc,"THF"," ",Chain,ResNum,")"," ", 
                       PosX," 0 ",PosY," 0 ",PosZ," 0 ")  #Refer to http://openmopac.net/manual/atom_numbers.html
          write(line,paste0(outputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)     
          if(trimws(AtomName)=="C1\'"){      
            cen1<-c(as.numeric(PosX),as.numeric(PosY),as.numeric(PosZ))
          } 
        } else { # base atoms do not include "'" such as N1 and N9
          if((trimws(AtomName)=="N9" && grepl("DA", ResName))
             ||(trimws(AtomName)=="N9" && grepl("DG", ResName))
             ||(trimws(AtomName)=="N1" && grepl("DC", ResName))
             ||(trimws(AtomName)=="N1" && grepl("DT", ResName))){   
            cen2<-c(as.numeric(PosX),as.numeric(PosY),as.numeric(PosZ))
            cen3<-cen1+(cen2-cen1)/2
            line= paste0("H","(ATOM  ",AtomID," ","H1\'\'",AlterLoc,"THF"," ",Chain,ResNum,")"," ", #For the atom name, we refered the supplemental information of E. Bignon et al. Sci. Rep., vol. 10, 17314 (2020).
                         cen3[1]," 1 ",cen3[2]," 1 ",cen3[3]," 1 ")
            write(line,paste0(outputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)   
          } else if( trimws(AtomName)=="P" || trimws(AtomName)=="OP1" || trimws(AtomName)=="OP2"){
            line= paste0(Element,"(ATOM  ",AtomID," ",AtomName,AlterLoc,"THF"," ",Chain,ResNum,")"," ",
                         PosX," 0 ",PosY," 0 ",PosZ," 0 ")
            write(line,paste0(outputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)    
          }
        }
      } else { #When not in the mid residue
        line= paste0(Element,"(ATOM  ",AtomID," ",AtomName,AlterLoc,ResName," ",Chain,ResNum,")"," ",
                     PosX," 0 ",PosY," 0 ",PosZ," 0 ")
        write(line,paste0(outputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)
      }
    }
  }   
  
  #-----------------------------------------------------------------------------
  if(inputconf=="del_mid_2nd_strand_sp"){
    
    if(charge=="number_of_phosphate"){char_num=-2*nchar(sequence)+2}
    if(!is.na(suppressWarnings(as.numeric(charge)))){char_num=as.numeric(charge)}
    
    writeLines(paste0(
      HF, ' ',
      hamiltonian, ' ',
      'GNORM=1.0', ' ', # the deleted base is replaced by H. Only this H is optimized.
      'CHARGE=',char_num, ' ',
      if(COSMO != F) {paste0(COSMO, ' ')},
      other, ' ',
      'PDBOUT','\n','\n'
    ), paste0(outputfolder,'/',sequence,'/',output_name,'.mop'))
    
    temp <- readLines(paste0(outputfolder,'/',sequence,'/',infile))
    if(length(grep("TER", temp))>0){temp<-temp[- grep("TER", temp)]}
    if(length(grep("END", temp))>0){temp<-temp[- grep("END", temp)]} 
    
    cen1<-c(0,0,0)
    cen2<-c(0,0,0)
    cen3<-c(0,0,0)
    
    for(i in grep("ATOM      1",temp):(length(temp))){#For format, please refer to http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
      Header<-substr(temp[i],1,4)
      AtomID<-substr(temp[i],7,11)
      AtomName<-substr(temp[i],13,16)
      AlterLoc<-substr(temp[i],17,17)
      ResName<-substr(temp[i],18,20)
      Chain<-substr(temp[i],22,22)
      ResNum<-substr(temp[i],23,26)
      InserRes<-substr(temp[i],27,27)
      PosX<-substr(temp[i],31,38)
      PosY<-substr(temp[i],39,46)
      PosZ<-substr(temp[i],47,54)
      Occu<-substr(temp[i],55,60)
      TempFac<-substr(temp[i],61,66)
      SegIden<-substr(temp[i],73,76)
      Element<-substr(temp[i],77,78)
      if(as.numeric(ResNum)==(nchar(sequence)+round(nchar(sequence)/2))){ # the resiude number of mid of the second strand = 1st strand + 2nd strand/2  
        if( grepl("\'", AtomName)){
          line= paste0(Element,"(ATOM  ",AtomID," ",AtomName,AlterLoc,"THF"," ",Chain,ResNum,")"," ", #and http://openmopac.net/manual/atom_numbers.html
                       PosX," 0 ",PosY," 0 ",PosZ," 0 ")
          write(line,paste0(outputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)   
          if(trimws(AtomName)=="C1\'"){      
            cen1<-c(as.numeric(PosX),as.numeric(PosY),as.numeric(PosZ))
          }
        } else {
          if((trimws(AtomName)=="N9" && grepl("DA", ResName))
             ||(trimws(AtomName)=="N9" && grepl("DG", ResName))
             ||(trimws(AtomName)=="N1" && grepl("DC", ResName))
             ||(trimws(AtomName)=="N1" && grepl("DT", ResName))){   
            cen2<-c(as.numeric(PosX),as.numeric(PosY),as.numeric(PosZ))
            cen3<-cen1+(cen2-cen1)/2
            line= paste0("H","(ATOM  ",AtomID," ","H1\'\'",AlterLoc,"THF"," ",Chain,ResNum,")"," ",  #For the atom name, we refered the supplemental information of E. Bignon et al. Sci. Rep., vol. 10, 17314 (2020).
                         cen3[1]," 1 ",cen3[2]," 1 ",cen3[3]," 1 ")
            write(line,paste0(outputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)   
          }else if( trimws(AtomName)=="P" || trimws(AtomName)=="OP1" || trimws(AtomName)=="OP2"){
            line= paste0(Element,"(ATOM  ",AtomID," ",AtomName,AlterLoc,"THF"," ",Chain,ResNum,")"," ",
                         PosX," 0 ",PosY," 0 ",PosZ," 0 ")
            write(line,paste0(outputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)    
          }
        }
      } else {
        line= paste0(Element,"(ATOM  ",AtomID," ",AtomName,AlterLoc,ResName," ",Chain,ResNum,")"," ",
                     PosX," 0 ",PosY," 0 ",PosZ," 0 ")
        write(line,paste0(outputfolder,'/',sequence,'/',output_name,'.mop'),append=TRUE)
      }
    }
  }  
  
  #-------------------------------------------------------------------------------
  
  cmd <- paste0("cd ",outputfolder,'/',sequence,";",
                "MOPAC2016.exe ", output_name,'.mop')
  system(cmd)
  
}



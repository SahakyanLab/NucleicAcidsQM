################################################################################
extract_Mulliken<-function(output.path,dna.seqs,QMcharge.inout){
  ### extract atom information from arc file ###
  input_arc=paste0(output.path, "/", dna.seqs, "/", gsub(".out",".arc",QMcharge.inout))
  arc_summary<-readLines(input_arc)
  num_of_atom=as.numeric(substr(arc_summary[8], start = nchar(arc_summary[8])-11, stop = nchar(arc_summary[8])-5))
  start_col_arc=length(arc_summary)-num_of_atom
  
  num<-start_col_arc[1]-1
  atom_num_array=c(1:num_of_atom)
  
  Element  <- c()
  AtomID   <- c()
  AtomName <- c()
  ResName  <- c()
  ResNum   <- c()
  
  Element[atom_num_array]  <- trimws(substr(arc_summary[num+atom_num_array],1,3))
  #AtomID[atom_num_array]   <- trimws(substr(arc_summary[num+atom_num_array],11,15))
  AtomName[atom_num_array] <- trimws(substr(arc_summary[num+atom_num_array],17,20))
  ResName[atom_num_array]  <- trimws(substr(arc_summary[num+atom_num_array],22,24))
  ResNum[atom_num_array]   <- trimws(substr(arc_summary[num+atom_num_array],28,30))
  
  ### extract charge information from out file ###
  input_out=paste0(output.path, "/", dna.seqs, "/", QMcharge.inout)
  out_summary <- system(paste0("tail -n ",num_of_atom+20 ," ",input_out), intern = TRUE)
  start_col_out=grep("NO.  ATOM   POPULATION      CHARGE", out_summary)
  
  Mullik_char<-c()
  Mullik_pop<-c()
  
  Mullik_char[atom_num_array]=substr(out_summary[start_col_out+atom_num_array],30,41)
  Mullik_pop[atom_num_array]=substr(out_summary[start_col_out+atom_num_array],16,27)
  data=paste(rep(dna.seqs, num_of_atom),ResNum,ResName,AtomName,Element,Mullik_char,Mullik_pop, sep=",")
  
  print(paste0("Mullik"," ",dna.seqs))
  return(data)    
}



################################################################################
extract_ESP<-function(output.path, dna.seqs, ESPcharge.inout){
  input_out=paste0(output.path, "/", dna.seqs, "/", ESPcharge.inout)
  arc_summary<-readLines(input_out)
  num_of_atom=as.numeric(substr(arc_summary[8], start = nchar(arc_summary[8])-11, stop = nchar(arc_summary[8])-5))
  start_col_arc=length(arc_summary)-num_of_atom
  
  num<-start_col_arc[1]-1
  atom_num_array=c(1:num_of_atom)
  
  Element=c()  
  AtomID=c()
  AtomName=c()
  ResName=c()
  ResNum=c()
  Charge=c()
  Population=c()
  
  Element[atom_num_array]  <- trimws(substr(arc_summary[num+atom_num_array],1,3))
  AtomID[atom_num_array]   <- trimws(substr(arc_summary[num+atom_num_array],11,15))
  AtomName[atom_num_array] <- trimws(substr(arc_summary[num+atom_num_array],17,20))
  ResName[atom_num_array]  <- trimws(substr(arc_summary[num+atom_num_array],22,24))
  ResNum[atom_num_array]   <- trimws(substr(arc_summary[num+atom_num_array],28,30))
  Charge[atom_num_array]   <- trimws(substr(arc_summary[num+atom_num_array],100,106))
  
  cal_pop<-function(x, Charge, Element){
    if(trimws(Element[x])=="H") return(1-as.numeric(Charge[x]))          
    if(trimws(Element[x])=="C") return(4-as.numeric(Charge[x]))          
    if(trimws(Element[x])=="N") return(5-as.numeric(Charge[x]))          
    if(trimws(Element[x])=="O") return(6-as.numeric(Charge[x]))          
    if(trimws(Element[x])=="P") return(7-as.numeric(Charge[x])) 
  }
  
  Population[atom_num_array]   <- lapply(atom_num_array, cal_pop, Charge=Charge, Element=Element)
  
  data=paste(rep(dna.seqs, num_of_atom),ResNum,ResName, AtomName,Element,Charge,Population, sep=",")
  print(paste0("ESP"," ",dna.seqs))
  return(data)  
}



################################################################################
read_charge<-function(table, mid_res_num){
  table=filter(table, res_num==mid_res_num)
  
  phos_atoms=c("P","OP1","OP2")
  sugar_atoms=c("O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'","H1'","C3'","H3'","C2'","H2'","H2''","O3'")
  phos=filter(table,atom_name %in% phos_atoms) %>% group_by(seq) %>%summarise(phosphate = mean(charge))
  sugar=filter(table,atom_name %in% sugar_atoms) %>% group_by(seq) %>%summarise(sugar = mean(charge))
  base=filter(table,!atom_name %in% append(phos_atoms,sugar_atoms)) %>% group_by(seq) %>%summarise(base = mean(charge))
  
  table$atom_name_new <- paste0(table$res_name,"_", table$atom_name)
  table=select(table,seq, atom_name_new, charge)
  table=pivot_wider(table, names_from = "atom_name_new", values_from = "charge")
  
  if(ncol(table)>1){
    table=left_join(table, phos, by = "seq")
    table=left_join(table, sugar, by = "seq")
    table=left_join(table, base, by = "seq")
  }
  return(table)    
}

read_population<-function(table, mid_res_num){
  table=filter(table, res_num==mid_res_num)
  
  phos_atoms=c("P","OP1","OP2")
  sugar_atoms=c("O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'","H1'","C3'","H3'","C2'","H2'","H2''","O3'")
  phos=filter(table,atom_name %in% phos_atoms) %>% group_by(seq) %>%summarise(phosphate = mean(population))
  sugar=filter(table,atom_name %in% sugar_atoms) %>% group_by(seq) %>%summarise(sugar = mean(population))
  base=filter(table,!atom_name %in% append(phos_atoms,sugar_atoms)) %>% group_by(seq) %>%summarise(base = mean(population))
  
  table$atom_name_new <- paste0(table$res_name,"_", table$atom_name)
  table=select(table,seq, atom_name_new, population)
  table=pivot_wider(table, names_from = "atom_name_new", values_from = "population")
  
  if(ncol(table)>1){
    table=left_join(table, phos, by = "seq")
    table=left_join(table, sugar, by = "seq")
    table=left_join(table, base, by = "seq")
  }
  
  return(table)    
}

read_curves<-function(curves_data, mid_res_num){
  space_num=which(readLines(curves_data)=="")
  
  BPaxis_table = read_csv(curves_data,n_max=(space_num[1]-2))
  BPaxis_table=filter(BPaxis_table, res_num1==mid_res_num)
  BPaxis_table=select(BPaxis_table, seq,Xdisp,Ydisp,Inclin,Tip,`Ax-bend`)
  
  intraBP_table = read_csv(curves_data,skip=space_num[1],n_max=(space_num[2]-space_num[1]-2))
  intraBP_table=filter(intraBP_table, res_num1==mid_res_num)
  intraBP_table=select(intraBP_table, seq,Shear,Stretch,Stagger,Buckle,Propel,Opening)
  
  interBP_table = read_csv(curves_data,skip=space_num[2],n_max=(space_num[3]-space_num[2]-2))
  interBP_table_1=filter(interBP_table, res_num1==mid_res_num)
  interBP_table_1=select(interBP_table_1, seq, Shift, Slide, Rise, Tilt, Roll, Twist, `H-Ris`, `H-Twi`)
  interBP_table_1=rename_with(interBP_table_1, function(x){paste0(x,"_to_3")}, .cols =c(Shift, Slide, Rise, Tilt, Roll, Twist, `H-Ris`, `H-Twi`) )
  interBP_table_2=filter(interBP_table, res_num2==mid_res_num)
  interBP_table_2=select(interBP_table_2, seq, Shift, Slide, Rise, Tilt, Roll, Twist, `H-Ris`, `H-Twi`)
  interBP_table_2=rename_with(interBP_table_2, function(x){paste0(x,"_from_5")}, .cols =c(Shift, Slide, Rise, Tilt, Roll, Twist, `H-Ris`, `H-Twi`) )
  
  Backbone_table = read_csv(curves_data,skip=space_num[3])
  Backbone_table=filter(Backbone_table, res_num==mid_res_num)
  Backbone_table=select(Backbone_table, seq, Alpha, Beta, Gamma, Delta, Epsil, Zeta, Chi, Phase, Ampli)
  
  mechanical_table=right_join(BPaxis_table,intraBP_table, by = "seq")
  mechanical_table=right_join(mechanical_table, interBP_table_1, by = "seq")
  mechanical_table=right_join(mechanical_table, interBP_table_2, by = "seq")
  mechanical_table=right_join(mechanical_table, Backbone_table, by = "seq")
  mechanical_table[2:ncol(mechanical_table)]=mechanical_table[2:ncol(mechanical_table)] %>% mutate(across(where(is.character), as.double))
  
  return(mechanical_table)    
}

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

generate_res_num<-function(num, k){
  if(num<=k){
    return(num+k)
  }
  if(num>k){
    return(num-k)
  }
}

################################################################################
add_comp=function(ori_data1,ori_data2,seq_convert, anti_res_num){
  data1 = read_csv(ori_data1)
  data2 = read_csv(ori_data2)
  comp_seq=seq_convert[data2$seq]
  comp_num=anti_res_num[data2$res_num]
  data2 =mutate(data2,seq=unlist(comp_seq))
  data2 =mutate(data2,res_num=unlist(comp_num))
  table=bind_rows(data1,data2)
  table=distinct(table ,seq, res_num, res_name, atom_name, element, charge, population)
  return(table)
}

wide=function(table1){
  table1$atom_name_new <- paste0(table1$res_name,"_", table1$atom_name)
  table1=select(table1,seq, atom_name_new, charge)
  table1=pivot_wider(table1, names_from = "atom_name_new", values_from = "charge") 
  return(table1)
}

make_raw_table<-function(table,seq_table_raw,j,kmer_num,flag,type,type_name){
  table_sense=filter(table, res_num==j)
  if(nrow(table_sense)>0){
    table_sense$atom_name_new <- paste0(table_sense$res_name,"_", table_sense$atom_name)
    table_sense=select(table_sense, seq, atom_name_new, starts_with(type))
    table_sense=pivot_wider(table_sense, names_from = "atom_name_new", values_from = type) 
    table_sense=rename_with(table_sense, function(x){paste0(flag,"_strandPlus_",j,"_",x,".",type_name)}, .cols =c(-seq) )
    seq_table_raw=left_join(seq_table_raw, table_sense, by = "seq")
  }
  
  table_anti=filter(table, res_num==kmer_num*2-j+1)
  if(nrow(table_anti)>0){
    table_anti$atom_name_new <- paste0(table_anti$res_name,"_", table_anti$atom_name)
    table_anti=select(table_anti, seq, atom_name_new, starts_with(type))
    table_anti=pivot_wider(table_anti, names_from = "atom_name_new", values_from = type)     
    table_anti=rename_with(table_anti, function(x){paste0(flag,"_strandMinus_",j,"_",x,".",type_name)}, .cols =c(-seq) )
    seq_table_raw=left_join(seq_table_raw, table_anti, by = "seq")
  }
  return(seq_table_raw) 
}


make_cal_table<-function(table,seq_table_cal,j,kmer_num,flag, type, type_name){
  cal_values=function(table1,type){
    phos_atoms=c("P","OP1","OP2")
    sugar_atoms=c("HO5'","O5'","C5'","H5'","H5''","C4'","H4'","O4'","C1'","H1'","C3'","H3'","C2'","H2'","H2''","O3'","HO3'")
    
    table=distinct(select(table1,seq))

    phos_mean=filter(table1,atom_name %in% phos_atoms) %>% group_by(seq) %>%summarise_at(vars(type), funs(phos_mean = mean))
    phos_max =filter(table1,atom_name %in% phos_atoms) %>% group_by(seq) %>%summarise_at(vars(type), funs(phos_max = max))
    phos_min =filter(table1,atom_name %in% phos_atoms) %>% group_by(seq) %>%summarise_at(vars(type), funs(phos_min = min))
    table=left_join(table, phos_mean, by='seq') %>%left_join(., phos_max, by='seq')%>%left_join(., phos_min, by='seq')

    sugar_mean=filter(table1,atom_name %in% sugar_atoms) %>% group_by(seq) %>%summarise_at(vars(type), funs(sugar_mean = mean))
    sugar_max =filter(table1,atom_name %in% sugar_atoms) %>% group_by(seq) %>%summarise_at(vars(type), funs(sugar_max = max))
    sugar_min =filter(table1,atom_name %in% sugar_atoms) %>% group_by(seq) %>%summarise_at(vars(type), funs(sugar_min = min))
    table=left_join(table, sugar_mean, by='seq') %>%left_join(., sugar_max, by='seq')%>%left_join(., sugar_min, by='seq') 
    
    base_mean=filter(table1,!atom_name %in% append(phos_atoms,sugar_atoms)) %>% group_by(seq) %>%summarise_at(vars(type), funs(base_mean = mean))
    base_max =filter(table1,!atom_name %in% append(phos_atoms,sugar_atoms)) %>% group_by(seq) %>%summarise_at(vars(type), funs(base_max = max))
    base_min =filter(table1,!atom_name %in% append(phos_atoms,sugar_atoms)) %>% group_by(seq) %>%summarise_at(vars(type), funs(base_min = min))
    table=left_join(table, base_mean, by='seq')%>%left_join(.,base_max, by='seq') %>%left_join(.,base_min, by='seq') 
    
    return(table)
  }
  
  table_sense=filter(table, res_num==j)
  if(nrow(table_sense)>0){
    table_sense=cal_values(table_sense,type)
    table_sense=rename_with(table_sense, function(x){paste0(flag,"_strandPlus_",j,"_",x,".",type_name)}, .cols =c(-seq) )
    seq_table_cal=left_join(seq_table_cal, table_sense, by = "seq")
  }
  
  table_anti=filter(table, res_num==kmer_num*2-j+1)
  if(nrow(table_anti)>0){
    table_anti=cal_values(table_anti,type)
    table_anti=rename_with(table_anti, function(x){paste0(flag,"_strandMinus_",j,"_",x,".",type_name)}, .cols =c(-seq) )
    seq_table_cal=left_join(seq_table_cal, table_anti, by = "seq")
  }
  return(seq_table_cal) 
}

################################################################################
add_comp_summary=function(ori_data1,ori_data2,seq_convert){
  data1 = read_csv(ori_data1)
  data2 = read_csv(ori_data2)
  comp_seq=seq_convert[data1$seq]
  data2 =mutate(data2,seq=unlist(comp_seq))
  table=bind_rows(data1,data2)
  table=distinct(table,seq, .keep_all = TRUE)
  return(table)
}

################################################################################
add_comp_curves1=function(data1, data2, seq_convert, anti_res_num){
  comp_seq=seq_convert[data2$seq]
  comp_num=anti_res_num[data2$res_num1]
  data2=mutate(data2,seq=unlist(comp_seq))
  data2=mutate(data2,res_num1=unlist(comp_num))
  table=bind_rows(data1,data2)
  table=distinct(table ,across(everything()))
  return(table)
}

add_comp_curves2=function(data1, data2, seq_convert, anti_res_num){
  comp_seq=seq_convert[data2$seq]
  comp_num1=anti_res_num[data2$res_num1]
  comp_num2=anti_res_num[data2$res_num2]
  data2 =mutate(data2,seq=unlist(comp_seq))
  data2 =mutate(data2,res_num1=unlist(comp_num1))
  data2 =mutate(data2,res_num2=unlist(comp_num2))
  table=bind_rows(data1, data2)  
  table=distinct(table ,across(everything()))
  return(table)
}

add_comp_curves3=function(data1, data2,seq_convert, anti_res_num){
  comp_seq=seq_convert[data2$seq]
  comp_num=anti_res_num[data2$res_num1]
  data2=mutate(data2,seq=unlist(comp_seq))
  data2=mutate(data2,res_num1=unlist(comp_num))
  table=bind_rows(data1,data2)
  table=distinct(table ,across(everything()))
  return(table)
}

################################################################################
GenerateKmers_with_comlement<-function(x){
  
  ### generate all k-mers ###
  bases<-c("A","C","G","T")
  kmers<-c("A","C","G","T")
  count = 0
  
  while(count < x-1 ){
    kmers_temp=c()
    for (kmer in kmers){
      for(base in bases){
        kmers_temp[length(kmers_temp)+1]=paste0(kmer,base)
      }
    }
    kmers=kmers_temp
    count = count+1
  }
  return(kmers)
}

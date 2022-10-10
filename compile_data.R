library(readr)
library(dplyr)
library(tidyr)
source("lib/read_datas.R")

seq_data="sequence.txt"
kmer_num=7

MullikCharge_add = F
MullikCharge_res_nums=c(1,2,3,4,5,6,7)
MullikCharge_data=c("B_ds"     = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/Mulliken_of_SubQM_duplex_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    #"B_ssPlus" = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/Mulliken_of_SubQM_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    #"B_ssMinus"= "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/Mulliken_of_SubQM_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    #"B_del1"   = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/Mulliken_of_SubQM_del_mid_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    #"B_del2"   = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/Mulliken_of_SubQM_del_mid_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    "A_ds"     = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/Mulliken_of_SubQM_duplex_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    #"A_ssPlus" = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/Mulliken_of_SubQM_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    #"A_ssMinus"= "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/Mulliken_of_SubQM_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    #"A_del1"   = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/Mulliken_of_SubQM_del_mid_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    #"A_del2"   = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/Mulliken_of_SubQM_del_mid_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    "Z_ds"     = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/Mulliken_of_SubQM_duplex_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt")#,
                    #"Z_ssPlus" = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/Mulliken_of_SubQM_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    #"Z_ssMinus"= "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/Mulliken_of_SubQM_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    #"Z_del1"   = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/Mulliken_of_SubQM_del_mid_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                    #"Z_del2"   = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/Mulliken_of_SubQM_del_mid_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt")
MullikCharge_raw_output="./compiled_data/Mullik_Charge_raw.txt"
MullikCharge_cal_output="./compiled_data/Mullik_Charge.txt"


MullikPop_add = F
MullikPop_res_nums=MullikCharge_res_nums
MullikPop_data=MullikCharge_data
MullikPop_raw_output="./compiled_data/Mullik_Density_raw.txt"
MullikPop_cal_output="./compiled_data/Mullik_Density.txt"


ESPCharge_add = F
ESPCharge_res_nums=c(1,2,3,4,5,6,7)
ESPCharge_data=c("B_ds"      = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/ESPcharge_of_SubQM_duplex_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 #"B_ssPlus"  = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/ESPcharge_of_SubQM_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 #"B_ssMinus" = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/ESPcharge_of_SubQM_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 #"B_del1"    = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/ESPcharge_of_SubQM_del_mid_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 #"B_del2"    = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/ESPcharge_of_SubQM_del_mid_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 "A_ds"      = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/ESPcharge_of_SubQM_duplex_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 #"A_ssPlus"  = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/ESPcharge_of_SubQM_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 #"A_ssMinus" = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/ESPcharge_of_SubQM_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 #"A_del1"    = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/ESPcharge_of_SubQM_del_mid_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 #"A_del2"    = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/ESPcharge_of_SubQM_del_mid_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 "Z_ds"      = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/ESPcharge_of_SubQM_duplex_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt")#,
                 #"Z_ssPlus"  = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/ESPcharge_of_SubQM_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 #"Z_ssMinus" = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/ESPcharge_of_SubQM_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 #"Z_del1"    = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/ESPcharge_of_SubQM_del_mid_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
                 #"Z_del2"    = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/ESPcharge_of_SubQM_del_mid_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt")

ESPCharge_raw_output="./compiled_data/ESP_Charge_raw.txt"
ESPCharge_cal_output="./compiled_data/ESP_Charge.txt"


ESPPop_add = F
ESPPop_res_nums=ESPCharge_res_nums
ESPPop_data=ESPCharge_data
ESPPop_raw_output="./compiled_data/ESP_Density_raw.txt"
ESPPop_cal_output="./compiled_data/ESP_Density.txt"

Curves_add = F
Curves_res_nums=c(1,2,3,4,5,6,7)
Curves_data=c("B_ds"     = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/Various_mechanical_parameter_of_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
              "B_ds_comp"= "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/Various_mechanical_parameter_of_comple_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
              "A_ds"     = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/Various_mechanical_parameter_of_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
              "A_ds_comp"= "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/Various_mechanical_parameter_of_comple_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
              "Z_ds"     = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/Various_mechanical_parameter_of_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
              "Z_ds_comp"= "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/Various_mechanical_parameter_of_comple_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt")
Curves_output="./compiled_data/Curves.txt"


DNAshape_add = F
DNAshape_res_nums=c(1,2,3,4,5,6,7)
DNAshape_data=c("B_ds" = "./data/DNA_B_Handbook1999/DNAshape.txt",
                "B_ds_comp" = "./data/DNA_B_Handbook1999/DNAshape_comp.txt")
DNAshape_output="./1_compiled_data/DNAshape_raw.txt"


Summary_add = T
Summary_data=c("B_ds"     = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/summary_of_SubQM_duplex_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "B_ssPlus" = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/summary_of_SubQM_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "B_ssMinus"= "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/summary_of_SubQM_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "B_del1"   = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/summary_of_SubQM_del_mid_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "B_del2"   = "../2022_3_25_BDNA_QMopt/data/DNA_B_Handbook1999/summary_of_SubQM_del_mid_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "A_ds"     = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/summary_of_SubQM_duplex_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "A_ssPlus" = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/summary_of_SubQM_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "A_ssMinus"= "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/summary_of_SubQM_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "A_del1"   = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/summary_of_SubQM_del_mid_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "A_del2"   = "../2022_3_25_ADNA_QMopt/data/DNA_A_Handbook1999/summary_of_SubQM_del_mid_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "Z_ds"     = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/summary_of_SubQM_duplex_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "Z_ssPlus" = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/summary_of_SubQM_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "Z_ssMinus"= "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/summary_of_SubQM_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "Z_del1"   = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/summary_of_SubQM_del_mid_1st_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt",
               "Z_del2"   = "../2022_3_25_ZDNA_QMopt/data/DNA_Z_Handbook1999/summary_of_SubQM_del_mid_2nd_strand_sp_from_QMopt_PM6-DH+_COSMO_from_MMopt_igb_1_step_5000_target_all_atoms.txt")
summary_raw_output="./compiled_data/energy.txt"
summary_cal_output="./compiled_data/denergy.txt"
#summary2_cal_output="./compiled_data/summary2_cal.txt" #Summary between A, B, and Z conformations



## Add complementary sequence data##############################################
ALLkmer=rename(as_tibble(GenerateKmers_with_comlement(7)), "seq"=value)
seq_convert=c()
seq_convert[ALLkmer$seq]=lapply(ALLkmer$seq,generate_complement)

seq_table1=read_csv(seq_data, col_names = "seq")
comp_seq=seq_convert[seq_table1$seq]
seq_table2=mutate(seq_table1,seq=unlist(comp_seq))
seq_table=bind_rows(seq_table1,seq_table2)
seq_table=distinct(seq_table,seq,.keep_all = TRUE)
seq_table=seq_table[order(seq_table$seq),]

anti_res_num=c()
anti_res_num[c(1:(kmer_num*2))]=lapply(c(1:(kmer_num*2)),generate_res_num,k=kmer_num)

## Read Mulliken charges and add table #########################################
if (MullikCharge_add) {
  data=MullikCharge_data
  raw_output=MullikCharge_raw_output
  cal_output=MullikCharge_cal_output
  res_nums=MullikCharge_res_nums
  type="charge"
  type_name="MullikCharge"
  
  seq_table_raw=seq_table
  seq_table_cal=seq_table
  for(i in 1:length(data)){
    ### Add complementary sequence data ###
    flag=names(data[i])
    if(flag=="B_ds")     {table=add_comp(data["B_ds"],     data["B_ds"],     seq_convert, anti_res_num)}
    if(flag=="B_ssPlus") {table=add_comp(data["B_ssPlus"], data["B_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="B_ssMinus"){table=add_comp(data["B_ssMinus"],data["B_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="B_del1")   {table=add_comp(data["B_del1"],   data["B_del2"],   seq_convert, anti_res_num)}
    if(flag=="B_del2")   {table=add_comp(data["B_del2"],   data["B_del1"],   seq_convert, anti_res_num)}
    if(flag=="A_ds")     {table=add_comp(data["A_ds"],     data["A_ds"],     seq_convert, anti_res_num)}
    if(flag=="A_ssPlus") {table=add_comp(data["A_ssPlus"], data["A_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="A_ssMinus"){table=add_comp(data["A_ssMinus"],data["A_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="A_del1")   {table=add_comp(data["A_del1"],   data["A_del2"],   seq_convert, anti_res_num)}
    if(flag=="A_del2")   {table=add_comp(data["A_del2"],   data["A_del1"],   seq_convert, anti_res_num)}
    if(flag=="Z_ds")     {table=add_comp(data["Z_ds"],     data["Z_ds"],     seq_convert, anti_res_num)}
    if(flag=="Z_ssPlus") {table=add_comp(data["Z_ssPlus"], data["Z_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="Z_ssMinus"){table=add_comp(data["Z_ssMinus"],data["Z_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="Z_del1")   {table=add_comp(data["Z_del1"],   data["Z_del2"],   seq_convert, anti_res_num)}
    if(flag=="Z_del2")   {table=add_comp(data["Z_del2"],   data["Z_del1"],   seq_convert, anti_res_num)}
    
    for (j in res_nums){
      seq_table_raw=make_raw_table(table,seq_table_raw,j,kmer_num, flag, type, type_name)
      seq_table_cal=make_cal_table(table,seq_table_cal,j,kmer_num, flag, type, type_name)
    }    
  }
  write.table(seq_table_raw, raw_output, sep = ",", row.names = FALSE)
  write.table(seq_table_cal, cal_output, sep = ",", row.names = FALSE)
}

## Read Mulliken population and add table #########################################
if (MullikPop_add) {
  data=MullikPop_data
  raw_output=MullikPop_raw_output
  cal_output=MullikPop_cal_output
  res_nums=MullikPop_res_nums
  type="population"
  type_name="MullikDensity"
  
  seq_table_raw=seq_table
  seq_table_cal=seq_table
  for(i in 1:length(data)){
    ### Add complementary sequence data ###
    flag=names(data[i])
    if(flag=="B_ds")     {table=add_comp(data["B_ds"],     data["B_ds"],     seq_convert, anti_res_num)}
    if(flag=="B_ssPlus") {table=add_comp(data["B_ssPlus"], data["B_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="B_ssMinus"){table=add_comp(data["B_ssMinus"],data["B_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="B_del1")   {table=add_comp(data["B_del1"],   data["B_del2"],   seq_convert, anti_res_num)}
    if(flag=="B_del2")   {table=add_comp(data["B_del2"],   data["B_del1"],   seq_convert, anti_res_num)}
    if(flag=="A_ds")     {table=add_comp(data["A_ds"],     data["A_ds"],     seq_convert, anti_res_num)}
    if(flag=="A_ssPlus") {table=add_comp(data["A_ssPlus"], data["A_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="A_ssMinus"){table=add_comp(data["A_ssMinus"],data["A_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="A_del1")   {table=add_comp(data["A_del1"],   data["A_del2"],   seq_convert, anti_res_num)}
    if(flag=="A_del2")   {table=add_comp(data["A_del2"],   data["A_del1"],   seq_convert, anti_res_num)}
    if(flag=="Z_ds")     {table=add_comp(data["Z_ds"],     data["Z_ds"],     seq_convert, anti_res_num)}
    if(flag=="Z_ssPlus") {table=add_comp(data["Z_ssPlus"], data["Z_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="Z_ssMinus"){table=add_comp(data["Z_ssMinus"],data["Z_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="Z_del1")   {table=add_comp(data["Z_del1"],   data["Z_del2"],   seq_convert, anti_res_num)}
    if(flag=="Z_del2")   {table=add_comp(data["Z_del2"],   data["Z_del1"],   seq_convert, anti_res_num)}
    
    for (j in res_nums){
      seq_table_raw=make_raw_table(table,seq_table_raw,j,kmer_num, flag, type, type_name)
      seq_table_cal=make_cal_table(table,seq_table_cal,j,kmer_num, flag, type, type_name)
    }    
  }
  write.table(seq_table_raw, raw_output, sep = ",", row.names = FALSE)
  write.table(seq_table_cal, cal_output, sep = ",", row.names = FALSE)
}

## Read ESP charges and add table #########################################
if (ESPCharge_add) {
  data=ESPCharge_data
  raw_output=ESPCharge_raw_output
  cal_output=ESPCharge_cal_output
  res_nums=ESPCharge_res_nums
  type="charge"
  type_name="ESPCharge"
  
  seq_table_raw=seq_table
  seq_table_cal=seq_table
  for(i in 1:length(data)){
    ### Add complementary sequence data ###
    flag=names(data[i])
    if(flag=="B_ds")     {table=add_comp(data["B_ds"],     data["B_ds"],     seq_convert, anti_res_num)}
    if(flag=="B_ssPlus") {table=add_comp(data["B_ssPlus"], data["B_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="B_ssMinus"){table=add_comp(data["B_ssMinus"],data["B_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="B_del1")   {table=add_comp(data["B_del1"],   data["B_del2"],   seq_convert, anti_res_num)}
    if(flag=="B_del2")   {table=add_comp(data["B_del2"],   data["B_del1"],   seq_convert, anti_res_num)}
    if(flag=="A_ds")     {table=add_comp(data["A_ds"],     data["A_ds"],     seq_convert, anti_res_num)}
    if(flag=="A_ssPlus") {table=add_comp(data["A_ssPlus"], data["A_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="A_ssMinus"){table=add_comp(data["A_ssMinus"],data["A_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="A_del1")   {table=add_comp(data["A_del1"],   data["A_del2"],   seq_convert, anti_res_num)}
    if(flag=="A_del2")   {table=add_comp(data["A_del2"],   data["A_del1"],   seq_convert, anti_res_num)}
    if(flag=="Z_ds")     {table=add_comp(data["Z_ds"],     data["Z_ds"],     seq_convert, anti_res_num)}
    if(flag=="Z_ssPlus") {table=add_comp(data["Z_ssPlus"], data["Z_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="Z_ssMinus"){table=add_comp(data["Z_ssMinus"],data["Z_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="Z_del1")   {table=add_comp(data["Z_del1"],   data["Z_del2"],   seq_convert, anti_res_num)}
    if(flag=="Z_del2")   {table=add_comp(data["Z_del2"],   data["Z_del1"],   seq_convert, anti_res_num)}
    
    for (j in res_nums){
      seq_table_raw=make_raw_table(table,seq_table_raw,j,kmer_num, flag, type, type_name)
      seq_table_cal=make_cal_table(table,seq_table_cal,j,kmer_num, flag, type, type_name)
    }    
  }
  write.table(seq_table_raw, raw_output, sep = ",", row.names = FALSE)
  write.table(seq_table_cal, cal_output, sep = ",", row.names = FALSE)
}

## Read ESP population and add table #########################################
if (ESPPop_add) {
  data=ESPPop_data
  raw_output=ESPPop_raw_output
  cal_output=ESPPop_cal_output
  res_nums=ESPPop_res_nums
  type="population"
  type_name="ESPDensity"
  
  seq_table_raw=seq_table
  seq_table_cal=seq_table
  for(i in 1:length(data)){
    ### Add complementary sequence data ###
    flag=names(data[i])
    if(flag=="B_ds")     {table=add_comp(data["B_ds"],     data["B_ds"],     seq_convert, anti_res_num)}
    if(flag=="B_ssPlus") {table=add_comp(data["B_ssPlus"], data["B_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="B_ssMinus"){table=add_comp(data["B_ssMinus"],data["B_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="B_del1")   {table=add_comp(data["B_del1"],   data["B_del2"],   seq_convert, anti_res_num)}
    if(flag=="B_del2")   {table=add_comp(data["B_del2"],   data["B_del1"],   seq_convert, anti_res_num)}
    if(flag=="A_ds")     {table=add_comp(data["A_ds"],     data["A_ds"],     seq_convert, anti_res_num)}
    if(flag=="A_ssPlus") {table=add_comp(data["A_ssPlus"], data["A_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="A_ssMinus"){table=add_comp(data["A_ssMinus"],data["A_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="A_del1")   {table=add_comp(data["A_del1"],   data["A_del2"],   seq_convert, anti_res_num)}
    if(flag=="A_del2")   {table=add_comp(data["A_del2"],   data["A_del1"],   seq_convert, anti_res_num)}
    if(flag=="Z_ds")     {table=add_comp(data["Z_ds"],     data["Z_ds"],     seq_convert, anti_res_num)}
    if(flag=="Z_ssPlus") {table=add_comp(data["Z_ssPlus"], data["Z_ssMinus"],seq_convert, anti_res_num)}
    if(flag=="Z_ssMinus"){table=add_comp(data["Z_ssMinus"],data["Z_ssPlus"], seq_convert, anti_res_num)}
    if(flag=="Z_del1")   {table=add_comp(data["Z_del1"],   data["Z_del2"],   seq_convert, anti_res_num)}
    if(flag=="Z_del2")   {table=add_comp(data["Z_del2"],   data["Z_del1"],   seq_convert, anti_res_num)}
    
    for (j in res_nums){
      seq_table_raw=make_raw_table(table,seq_table_raw,j,kmer_num, flag, type, type_name)
      seq_table_cal=make_cal_table(table,seq_table_cal,j,kmer_num, flag, type, type_name)
    }    
  }
  write.table(seq_table_raw, raw_output, sep = ",", row.names = FALSE)
  write.table(seq_table_cal, cal_output, sep = ",", row.names = FALSE)
}


### Mechanical table ###########################################################
if (Curves_add){
  test<-names(Curves_data)
　flags=test[-which(test %in% c("B_ds_comp","A_ds_comp","Z_ds_comp"))] 
  for (flag_i in flags){
    if(flag_i=="B_ds"){flag_j="B_ds_comp"}
    if(flag_i=="A_ds"){flag_j="A_ds_comp"}
    if(flag_i=="Z_ds"){flag_j="Z_ds_comp"}
    
    ### Read data and add complementary data ###
    space_num1=which(readLines(Curves_data[flag_i])=="")
    space_num2=which(readLines(Curves_data[flag_j])=="")
  
    BPaxis_table1 = read_csv(Curves_data[flag_i],n_max=(space_num1[1]-2), na = c("---", "NA"))
    BPaxis_table2 = read_csv(Curves_data[flag_j],n_max=(space_num2[1]-2), na = c("---", "NA"))
    BPaxis_table  = add_comp_curves1(BPaxis_table1, BPaxis_table2, seq_convert, anti_res_num)
  
    intraBP_table1 = read_csv(Curves_data[flag_i],skip=space_num1[1],n_max=(space_num1[2]-space_num1[1]-2))
    intraBP_table2 = read_csv(Curves_data[flag_j],skip=space_num2[1],n_max=(space_num2[2]-space_num2[1]-2))
    intraBP_table  = add_comp_curves1(intraBP_table1,intraBP_table2, seq_convert, anti_res_num) 
  
    interBP_table1 = read_csv(Curves_data[flag_i],skip=space_num1[2],n_max=(space_num1[3]-space_num1[2]-2))
    interBP_table2 = read_csv(Curves_data[flag_j],skip=space_num2[2],n_max=(space_num2[3]-space_num2[2]-2))
    interBP_table  = add_comp_curves2(interBP_table1, interBP_table2, seq_convert, anti_res_num)
  
    Backbone_table1 = read_csv(Curves_data[flag_i],skip=space_num1[3], na = c("----", "NA"),)
    Backbone_table2 = read_csv(Curves_data[flag_j],skip=space_num2[3], na = c("----", "NA"),)
    Backbone_table = add_comp_curves3(Backbone_table1, Backbone_table2, seq_convert, anti_res_num)
  
  
    ### make tables ###
    for (j in Curves_res_nums){
      BPaxis_table_1=filter(BPaxis_table, res_num1==j)
      BPaxis_table_1=select(BPaxis_table_1, seq,Xdisp,Ydisp,Inclin,Tip,`Ax-bend`)
      BPaxis_table_1=rename_with(BPaxis_table_1, function(x){paste0(names(Curves_data[flag_i]),"_strandPlus_",j,".Curves_",x)}, .cols =c(-seq) )
    
      intraBP_table_1=filter(intraBP_table, res_num1==j)
      intraBP_table_1=select(intraBP_table_1, seq,Shear,Stretch,Stagger,Buckle,Propel,Opening) 
      intraBP_table_1=rename_with(intraBP_table_1, function(x){paste0(names(Curves_data[flag_i]),"_strandPlus_",j,".Curves_",x)}, .cols =c(-seq) )
    
      interBP_table_1=filter(interBP_table, res_num1==j)
      interBP_table_1=select(interBP_table_1, seq, Shift, Slide, Rise, Tilt, Roll, Twist, `H-Ris`, `H-Twi`)
      interBP_table_1=rename_with(interBP_table_1, function(x){paste0(names(Curves_data[flag_i]),"_strandPlus_",j,".Curves_",x,"_to_3")}, .cols =c(Shift, Slide, Rise, Tilt, Roll, Twist, `H-Ris`, `H-Twi`) )
      interBP_table_2=filter(interBP_table, res_num2==j)
      interBP_table_2=select(interBP_table_2, seq, Shift, Slide, Rise, Tilt, Roll, Twist, `H-Ris`, `H-Twi`)
      interBP_table_2=rename_with(interBP_table_2, function(x){paste0(names(Curves_data[flag_i]),"_strandPlus_",j,".Curves_",x,"_from_5")}, .cols =c(Shift, Slide, Rise, Tilt, Roll, Twist, `H-Ris`, `H-Twi`) )
    
      Backbone_table1=filter(Backbone_table, res_num1==j)
      Backbone_table1=select(Backbone_table1, seq, Alpha, Beta, Gamma, Delta, Epsil, Zeta, Chi, Phase, Ampli)
      Backbone_table1=rename_with(Backbone_table1, function(x){paste0(names(Curves_data[flag_i]),"_strandPlus_",j,".Curves_",x)}, .cols =c(-seq) )
      Backbone_table2=filter(Backbone_table, res_num1==kmer_num*2-j+1)
      Backbone_table2=select(Backbone_table2, seq, Alpha, Beta, Gamma, Delta, Epsil, Zeta, Chi, Phase, Ampli)
      Backbone_table2=rename_with(Backbone_table2, function(x){paste0(names(Curves_data[flag_i]),"_strandMinus_",j,".Curves_",x)}, .cols =c(-seq) )
      ### combine and write table ###
      seq_table=left_join(seq_table, BPaxis_table_1,  by = "seq")
      seq_table=left_join(seq_table, intraBP_table_1, by = "seq")
      seq_table=left_join(seq_table, interBP_table_1, by = "seq")
      seq_table=left_join(seq_table, interBP_table_2, by = "seq")
      seq_table=left_join(seq_table, Backbone_table1, by = "seq")
      seq_table=left_join(seq_table, Backbone_table2, by = "seq")
     }
   }
  write.table(seq_table, Curves_output, sep = ",",row.names = FALSE)
}


if (DNAshape_add) {
    space_num1=which(readLines(DNAshape_data["B_ds"])=="")
    space_num2=which(readLines(DNAshape_data["B_ds_comp"])=="")
    
    intra1 = read_csv(DNAshape_data["B_ds"],n_max=(space_num1[1]-2))
    intra2 = read_csv(DNAshape_data["B_ds_comp"],n_max=(space_num2[1]-2))
    intra=bind_rows(intra1, intra2)  
    intra=distinct(intra,across(everything()))
      
    inter1 = read_csv(DNAshape_data["B_ds"],skip=space_num1[1])
    inter2 = read_csv(DNAshape_data["B_ds_comp"],skip=space_num2[1])
    inter=bind_rows(inter1, inter2)  
    inter=distinct(inter,across(everything()))
    
    for (j in DNAshape_res_nums){
      intra_1=select(intra, seq, type, as.character(j))
      intra_1=pivot_wider(intra_1, names_from = "type", values_from = as.character(j))
      intra_1=rename_with(intra_1, function(x){paste0(names(DNAshape_data["B_ds"]),"_",j,".DNAshapeR_",x)}, .cols =c(-seq) )
      seq_table=left_join(seq_table, intra_1, by = "seq")
      
      if(j>1){
        inter_1=select(inter, seq, type, (as.character(j-1)))
        inter_1=pivot_wider(inter_1, names_from = "type", values_from = (as.character(j-1)))
        inter_1=rename_with(inter_1, function(x){paste0(names(DNAshape_data["B_ds"]),"_",j,".DNAshapeR_",x,"_from_5")}, .cols =c(-seq) )
        seq_table=left_join(seq_table, inter_1, by = "seq")
      }
      
      if(j<kmer_num){
        inter_2=select(inter, seq, type, as.character(j))
        inter_2=pivot_wider(inter_2, names_from = "type", values_from = as.character(j))
        inter_2=rename_with(inter_2, function(x){paste0(names(DNAshape_data["B_ds"]),"_",j,".DNAshapeR_",x,"_to_3")}, .cols =c(-seq) )
        seq_table=left_join(seq_table, inter_2, by = "seq")
      }
  }
  write.table(seq_table, DNAshape_output, sep = ",",row.names = FALSE)
}



## Read summary table ###################################################
if (Summary_add){
  for(i in 1:length(Summary_data)){
    ### Add complementary sequence data ###
    flag_i=names(Summary_data[i])
    if(flag_i=="B_ds")     {flag_j="B_ds"}
    if(flag_i=="B_ssPlus") {flag_j="B_ssMinus"}
    if(flag_i=="B_ssMinus"){flag_j="B_ssPlus"}
    if(flag_i=="B_del1")   {flag_j="B_del2"}
    if(flag_i=="B_del2")   {flag_j="B_del1"}
    if(flag_i=="A_ds")     {flag_j="A_ds"}
    if(flag_i=="A_ssPlus") {flag_j="A_ssMinus"}
    if(flag_i=="A_ssMinus"){flag_j="A_ssPlus"}
    if(flag_i=="A_del1")   {flag_j="A_del2"}
    if(flag_i=="A_del2")   {flag_j="A_del1"}
    if(flag_i=="Z_ds")     {flag_j="Z_ds"}
    if(flag_i=="Z_ssPlus") {flag_j="Z_ssMinus"}
    if(flag_i=="Z_ssMinus"){flag_j="Z_ssPlus"}
    if(flag_i=="Z_del1")   {flag_j="Z_del2"}
    if(flag_i=="Z_del2")   {flag_j="Z_del1"}
    table=add_comp_summary(Summary_data[flag_i],Summary_data[flag_j],seq_convert)
    table=rename_with(table, function(x){paste0(names(Summary_data[flag_i]),".",x)}, .cols =c(-seq) )
    seq_table=left_join(seq_table, table, by = "seq")
  }
  
  data_ds=select(seq_table, seq, 
                 B_ds.Ehof, B_ds.Dipole,   B_ds.IP,   B_ds.Ehomo,   B_ds.Elumo,
                 A_ds.Ehof, A_ds.Dipole,   A_ds.IP,   A_ds.Ehomo,   A_ds.Elumo,
                 Z_ds.Ehof, Z_ds.Dipole,   Z_ds.IP,   Z_ds.Ehomo,   Z_ds.Elumo)
  write.table(data_ds, summary_raw_output, sep = ",",row.names = FALSE)
  

  ### Calculation and write table ###
  data=seq_table
  ###
  data=mutate(data, B_ds.dEhof_ds_ss    =B_ds.Ehof  -(B_ssPlus.Ehof  +B_ssMinus.Ehof))
  data=mutate(data, B_ds.dDipole_ds_ss  =B_ds.Dipole-(B_ssPlus.Dipole+B_ssMinus.Dipole))
  data=mutate(data, B_ds.dIP_ds_ss      =B_ds.IP    -(B_ssPlus.IP    +B_ssMinus.IP))
  data=mutate(data, B_ds.dEhomo_ds_ss   =B_ds.Ehomo -(B_ssPlus.Ehomo +B_ssMinus.Ehomo))
  data=mutate(data, B_ds.dElumo_ds_ss   =B_ds.Elumo -(B_ssPlus.Elumo +B_ssMinus.Elumo))
  data=mutate(data, B_ds.dEhof_ds_del1  =B_ds.Ehof  - B_del1.Ehof)
  data=mutate(data, B_ds.dDipole_ds_del1=B_ds.Dipole- B_del1.Dipole)
  data=mutate(data, B_ds.dIP_ds_del1    =B_ds.IP    - B_del1.IP)
  data=mutate(data, B_ds.dEhomo_ds_del1 =B_ds.Ehomo - B_del1.Ehomo)
  data=mutate(data, B_ds.dElumo_ds_del1 =B_ds.Elumo - B_del1.Elumo)
  data=mutate(data, B_ds.dEhof_ds_del2  =B_ds.Ehof  - B_del2.Ehof)
  data=mutate(data, B_ds.dDipole_ds_del2=B_ds.Dipole- B_del2.Dipole)
  data=mutate(data, B_ds.dIP_ds_del2    =B_ds.IP    - B_del2.IP)
  data=mutate(data, B_ds.dEhomo_ds_del2 =B_ds.Ehomo - B_del2.Ehomo)
  data=mutate(data, B_ds.dElumo_ds_del2 =B_ds.Elumo - B_del2.Elumo)

  ###
  data=mutate(data, A_ds.dEhof_ds_ss    =A_ds.Ehof  -(A_ssPlus.Ehof  + A_ssMinus.Ehof))
  data=mutate(data, A_ds.dDipole_ds_ss  =A_ds.Dipole-(A_ssPlus.Dipole+ A_ssMinus.Dipole))
  data=mutate(data, A_ds.dIP_ds_ss      =A_ds.IP    -(A_ssPlus.IP    + A_ssMinus.IP))
  data=mutate(data, A_ds.dEhomo_ds_ss   =A_ds.Ehomo -(A_ssPlus.Ehomo + A_ssMinus.Ehomo))
  data=mutate(data, A_ds.dElumo_ds_ss   =A_ds.Elumo -(A_ssPlus.Elumo + A_ssMinus.Elumo))
  data=mutate(data, A_ds.dEhof_ds_del1  =A_ds.Ehof  - A_del1.Ehof)
  data=mutate(data, A_ds.dDipole_ds_del1=A_ds.Dipole- A_del1.Dipole)
  data=mutate(data, A_ds.dIP_ds_del1    =A_ds.IP    - A_del1.IP)
  data=mutate(data, A_ds.dEhomo_ds_del1 =A_ds.Ehomo - A_del1.Ehomo)
  data=mutate(data, A_ds.dElumo_ds_del1 =A_ds.Elumo - A_del1.Elumo)
  data=mutate(data, A_ds.dEhof_ds_del2  =A_ds.Ehof  - A_del2.Ehof)
  data=mutate(data, A_ds.dDipole_ds_del2=A_ds.Dipole- A_del2.Dipole)
  data=mutate(data, A_ds.dIP_ds_del2    =A_ds.IP    - A_del2.IP)
  data=mutate(data, A_ds.dEhomo_ds_del2 =A_ds.Ehomo - A_del2.Ehomo)
  data=mutate(data, A_ds.dElumo_ds_del2 =A_ds.Elumo - A_del2.Elumo)

  ###
  data=mutate(data, Z_ds.dEhof_ds_ss    =Z_ds.Ehof  -(Z_ssPlus.Ehof  + Z_ssMinus.Ehof))
  data=mutate(data, Z_ds.dDipole_ds_ss  =Z_ds.Dipole-(Z_ssPlus.Dipole+ Z_ssMinus.Dipole))
  data=mutate(data, Z_ds.dIP_ds_ss      =Z_ds.IP    -(Z_ssPlus.IP    + Z_ssMinus.IP))
  data=mutate(data, Z_ds.dEhomo_ds_ss   =Z_ds.Ehomo -(Z_ssPlus.Ehomo + Z_ssMinus.Ehomo))
  data=mutate(data, Z_ds.dElumo_ds_ss   =Z_ds.Elumo -(Z_ssPlus.Elumo + Z_ssMinus.Elumo))
  data=mutate(data, Z_ds.dEhof_ds_del1  =Z_ds.Ehof  - Z_del1.Ehof)
  data=mutate(data, Z_ds.dDipole_ds_del1=Z_ds.Dipole- Z_del1.Dipole)
  data=mutate(data, Z_ds.dIP_ds_del1    =Z_ds.IP    - Z_del1.IP)
  data=mutate(data, Z_ds.dEhomo_ds_del1 =Z_ds.Ehomo - Z_del1.Ehomo)
  data=mutate(data, Z_ds.dElumo_ds_del1 =Z_ds.Elumo - Z_del1.Elumo)
  data=mutate(data, Z_ds.dEhof_ds_del2  =Z_ds.Ehof  - Z_del2.Ehof)
  data=mutate(data, Z_ds.dDipole_ds_del2=Z_ds.Dipole- Z_del2.Dipole)
  data=mutate(data, Z_ds.dIP_ds_del2    =Z_ds.IP    - Z_del2.IP)
  data=mutate(data, Z_ds.dEhomo_ds_del2 =Z_ds.Ehomo - Z_del2.Ehomo)
  data=mutate(data, Z_ds.dElumo_ds_del2 =Z_ds.Elumo - Z_del2.Elumo)
  
  # ### Between A,B, and Z conformations
  # data=mutate(data, B_ds.dEhof_A_ds  = B_ds.Ehof  - A_ds.Ehof)
  # data=mutate(data, B_ds.dDipole_A_ds= B_ds.Dipole- A_ds.Dipole)
  # data=mutate(data, B_ds.dIP_A_ds    = B_ds.IP    - A_ds.IP)
  # data=mutate(data, B_ds.dEhomo_A_ds = B_ds.Ehomo - A_ds.Ehomo)
  # data=mutate(data, B_ds.dElumo_A_ds = B_ds.Elumo - A_ds.Elumo)
  # data=mutate(data, B_ds.dEhof_Z_ds  = B_ds.Ehof  - Z_ds.Ehof)
  # data=mutate(data, B_ds.dDipole_Z_ds= B_ds.Dipole- Z_ds.Dipole)
  # data=mutate(data, B_ds.dIP_Z_ds    = B_ds.IP    - Z_ds.IP)
  # data=mutate(data, B_ds.dEhomo_Z_ds = B_ds.Ehomo - Z_ds.Ehomo)
  # data=mutate(data, B_ds.dElumo_Z_ds = B_ds.Elumo - Z_ds.Elumo)
  # data=mutate(data, A_ds.dEhof_Z_ds  = A_ds.Ehof  - Z_ds.Ehof)
  # data=mutate(data, A_ds.dDipole_Z_ds= A_ds.Dipole- Z_ds.Dipole)
  # data=mutate(data, A_ds.dIP_Z_ds    = A_ds.IP    - Z_ds.IP)
  # data=mutate(data, A_ds.dEhomo_Z_ds = A_ds.Ehomo - Z_ds.Ehomo)
  # data=mutate(data, A_ds.dElumo_Z_ds = A_ds.Elumo - Z_ds.Elumo)
  
  ###
  data1=select(data, seq, 
              B_ds.dEhof_ds_ss,    B_ds.dIP_ds_ss,   B_ds.dEhomo_ds_ss,   B_ds.dElumo_ds_ss,   #B_ds.dDipole_ds_ss,
              B_ds.dEhof_ds_del1,  B_ds.dIP_ds_del1, B_ds.dEhomo_ds_del1, B_ds.dElumo_ds_del1, #B_ds.dDipole_ds_del1,
              B_ds.dEhof_ds_del2,  B_ds.dIP_ds_del2, B_ds.dEhomo_ds_del2, B_ds.dElumo_ds_del2, #B_ds.dDipole_ds_del2,
              A_ds.dEhof_ds_ss,    A_ds.dIP_ds_ss,   A_ds.dEhomo_ds_ss,   A_ds.dElumo_ds_ss,   #A_ds.dDipole_ds_ss, 
              A_ds.dEhof_ds_del1,  A_ds.dIP_ds_del1, A_ds.dEhomo_ds_del1, A_ds.dElumo_ds_del1, #A_ds.dDipole_ds_del1,
              A_ds.dEhof_ds_del2,  A_ds.dIP_ds_del2, A_ds.dEhomo_ds_del2, A_ds.dElumo_ds_del2, #A_ds.dDipole_ds_del2,
              Z_ds.dEhof_ds_ss,    Z_ds.dIP_ds_ss,   Z_ds.dEhomo_ds_ss,   Z_ds.dElumo_ds_ss,   #Z_ds.dDipole_ds_ss, 
              Z_ds.dEhof_ds_del1,  Z_ds.dIP_ds_del1, Z_ds.dEhomo_ds_del1, Z_ds.dElumo_ds_del1, #Z_ds.dDipole_ds_del1,
              Z_ds.dEhof_ds_del2,  Z_ds.dIP_ds_del2, Z_ds.dEhomo_ds_del2, Z_ds.dElumo_ds_del2)#, #Z_ds.dDipole_ds_del2)
  write.table(data1, summary_cal_output, sep = ",",row.names = FALSE)
  
  # ###  
  # data2=select(data, seq, 
  #             B_ds.dEhof_A_ds,   B_ds.dDipole_A_ds,   B_ds.dIP_A_ds,   B_ds.dEhomo_A_ds,   B_ds.dElumo_A_ds,
  #             B_ds.dEhof_Z_ds,   B_ds.dDipole_Z_ds,   B_ds.dIP_Z_ds,   B_ds.dEhomo_Z_ds,   B_ds.dElumo_Z_ds,
  #             A_ds.dEhof_Z_ds,   A_ds.dDipole_Z_ds,   A_ds.dIP_Z_ds,   A_ds.dEhomo_Z_ds,   A_ds.dElumo_Z_ds)  
  # write.table(data2, summary2_cal_output, sep = ",",row.names = FALSE)  
  
}


################################################################################
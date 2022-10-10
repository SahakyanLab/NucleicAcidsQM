###### function for generating a complement sequence ######
generate_complement<-function(seq){
  len = nchar(seq)
  seq=unlist(strsplit(seq,""))
  wc_seq<-rep(NA,len)
  pair<-function(x){
    if(x=="A"){return("T")}
    if(x=="C"){return("G")}
    if(x=="G"){return("C")}
    if(x=="T"){return("A")}

  }
  wc_seq=unlist(lapply(seq,pair))
  wc_seq=rev(wc_seq)
  wc_seq=paste(wc_seq,collapse="")
  return(wc_seq)  
}

###### function for searching a sequence number ######
search_seq_num<-function(seq){
  len = nchar(seq)
  seq_num=1
  for (i in 1:len){
    base=substr(seq,i,i)
    if(base=="A"){seq_num=seq_num+4^(len-(i-1))*0/4}
    if(base=="C"){seq_num=seq_num+4^(len-(i-1))*1/4}
    if(base=="G"){seq_num=seq_num+4^(len-(i-1))*2/4}
    if(base=="T"){seq_num=seq_num+4^(len-(i-1))*3/4}
  }
  return(seq_num)  
}

###### function for judging redundant sequence numbers ######
redundant_judge<-function(X,Y){
  if(X>Y) {return (X)}
}

###### main function  ######
GenerateKmers<-function(x){
  
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
  
  seq_num=c(1:length(kmers))
  
  ### generate all complement k-mers ###
  comp_kmers=c(NA,length(kmers))
  comp_kmers=unlist(lapply(kmers,generate_complement))
  comp_seq_num=unlist(lapply(comp_kmers,search_seq_num))
  
  ### subtract redundant sequences ###
  redundant_seq_num=unlist(mapply(redundant_judge,seq_num,comp_seq_num))
  kmers_subtracted=kmers[-c(redundant_seq_num)]
  
  return(kmers_subtracted)
}
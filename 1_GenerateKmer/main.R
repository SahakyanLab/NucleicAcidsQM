###input###
k=3
###make sequences of Kmers###
source("lib/GenerateKmers.R")
sequences_arrays=GenerateKmers(k)

### write sequences ###
file.create("sequences.txt")
write(sequences_arrays,"sequences.txt",sep="\n",append=T)








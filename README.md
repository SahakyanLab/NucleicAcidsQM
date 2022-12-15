# NucleicAcidsQM
Research and development code for MM and QM modelling on nucleic acids.

(1) Run main.R  
a) DNA model generation according to sequence.txt.  
b) Molecular mechanics (MM) and semi-quantum mechanical (QM) optimization.  
c) Single point calculation and Curves+ calculation for electronic and geometric features.  
d) Summary of features are put at ./data directory.  

(2) Run compile_data.R  
Extract and complie features obtained from the above calculation.  
Obtained datas are put at ./compiled_data directory.

(3) Note for BDNA
The following sequences, the initial position of hydrogen cap for the center deletaion state is C1'+(N1-C1')/2 or C1'+(N9-C1')/2 instead of C1'+(N1-C1')/3 or C1'+(N9-C1')/3. (in ./lib/SubQM cen3<-cen1+(cen2-cen1)/3 to cen3<-cen1+(cen2-cen1)/2)

del1
GTCTACA
GATCCAA
GAGTCAA
CTTCTGC
CGCTCCG
AGTCGGG
ACAATCG
AATCTAT
AATCCGA

del2
TAACTGA
CAGGGTA     
ACGCCTG
ACAGTTC
AACAGCG

The following sequences, the initial position of hydrogen cap for the center deletaion state is C1'+(N1-C1')/2.5 or C1'+(N9-C1')/2.5 instead of C1'+(N1-C1')/3 or C1'+(N9-C1')/3. (in ./lib/SubQM cen3<-cen1+(cen2-cen1)/3 to cen3<-cen1+(cen2-cen1)/2.5)

del2
ACTACCT


(4) Note for ADNA
The following sequences, the initial position of hydrogen cap for the center deletaion state is C1'+(N1-C1')/2 or C1'+(N9-C1')/2 instead of C1'+(N1-C1')/3 or C1'+(N9-C1')/3.

del1
CGATGTC
ATGTCCC
ACGTCGA

del2 
CAAGTCC
AGAGGTG
ACGACGC

The following sequences, the initial position of hydrogen cap for the center deletaion state is C1'+(N1-C1')/2.5 or C1'+(N9-C1')/2.5 instead of C1'+(N1-C1')/3 or C1'+(N9-C1')/3.

del2 
GTTGCAA



(5) Note for ZDNA
The following sequences, the initial position of hydrogen cap for the center deletaion state is C1'+(N1-C1')/2 or C1'+(N9-C1')/2 instead of C1'+(N1-C1')/3 or C1'+(N9-C1')/3. (in ./lib/SubQM)

del1
GGGCTGA
CTACTGA
CGACCAG
ATCTCGC
ACCCTTG

del2
GATATAC
AGTTCTA
CGTATTA

The following sequences, the initial position of hydrogen cap for the center deletaion state is C1'+(N1-C1')/2.5 or C1'+(N9-C1')/2.5 instead of C1'+(N1-C1')/3 or C1'+(N9-C1')/3. (in ./lib/SubQM)
del2
CAAATGA
ATGAGTA
ATGAGCC
GTTATGA 



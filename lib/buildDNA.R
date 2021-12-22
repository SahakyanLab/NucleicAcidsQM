buildDNA <- function(input.conformation, strand, seq, output.path="./") {

if (input.conformation =="DNA_A") {
  nab.in <- paste0(
      
 'molecule m;
  m = fd_helix( "adna", "', seq, '", "dna" );
  putpdb( "', seq, '.pdb", m, "-wwpdb");'
      
  )
    
  writeLines(nab.in, paste0(output.path, "/fd_helix.nab"))
    
  system(paste0("cd ", output.path,"; nab fd_helix.nab && ./fd_helix"))
  
  if (strand =="first_strand") {
    temp <- readLines(paste0(output.path,dna.seq,".pdb"))
    writeLines(temp[1:grep("TER",temp)[1]],paste0(output.path,dna.seq,".pdb"))
  } else if (strand =="second_strand") {
    temp <- readLines(paste0(output.path,dna.seq,".pdb"))
    writeLines(temp[(grep("TER",temp)[1]+1):grep("TER",temp)[2]],paste0(output.path,dna.seq,".pdb"))
  } 
  
}##end of conformation =="DNA_A"  

if (input.conformation =="DNA_A_Handbook1999") {
##################################################################################
### A-DNA by Oxford handbook of nucleic acid structure (1999) by S. Neidle.    ###
### Template pdb was made by 3DNA (25 repeating of sequence: ccccccccccccccccccccccccc )###
### Hydrogens were attached by addH of Chimera.                                ###
### Furthermore, we modified H2',H2'',H5',H5'' to H2'1,H2'2,H5'1,H5'2          ###
##################################################################################
    
template.in <- paste0(
      "ATOM   7472  P   DC  A 250       2.808   0.717  -8.564  1.00  1.00           P
ATOM   7473  O1P DC  A 250       4.023   0.637  -9.405  1.00  1.00           O
ATOM   7474  O2P DC  A 250       2.791  -0.154  -7.366  1.00  1.00           O
ATOM   7475  O5' DC  A 250       2.554   2.234  -8.128  1.00  1.00           O
ATOM   7476  C5' DC  A 250       2.260   3.212  -9.143  1.00  1.00           C
ATOM   7477  C4' DC  A 250       1.931   4.543  -8.502  1.00  1.00           C
ATOM   7478  O4' DC  A 250       0.592   4.469  -7.931  1.00  1.00           O
ATOM   7479  C3' DC  A 250       2.813   4.950  -7.319  1.00  1.00           C
ATOM   7480  O3' DC  A 250       4.063   5.472  -7.750  1.00  1.00           O
ATOM   7481  C2' DC  A 250       1.881   5.966  -6.664  1.00  1.00           C
ATOM   7482  C1' DC  A 250       0.550   5.225  -6.731  1.00  1.00           C
ATOM   7483  N1  DC  A 250       0.320   4.302  -5.583  1.00  1.00           N
ATOM   7484  C2  DC  A 250      -0.227   4.837  -4.421  1.00  1.00           C
ATOM   7485  O2  DC  A 250      -0.496   6.044  -4.388  1.00  1.00           O
ATOM   7486  N3  DC  A 250      -0.444   4.013  -3.363  1.00  1.00           N
ATOM   7487  C4  DC  A 250      -0.139   2.712  -3.439  1.00  1.00           C
ATOM   7488  N4  DC  A 250      -0.371   1.950  -2.381  1.00  1.00           N
ATOM   7489  C5  DC  A 250       0.424   2.141  -4.627  1.00  1.00           C
ATOM   7490  C6  DC  A 250       0.634   2.982  -5.671  1.00  1.00           C
ATOM   7491 H5'1 DC  A 250       1.407   2.875  -9.732  1.00  0.00           H
ATOM   7492 H5'2 DC  A 250       3.126   3.329  -9.794  1.00  0.00           H
ATOM   7493  H4' DC  A 250       1.958   5.326  -9.260  1.00  0.00           H
ATOM   7494  H3' DC  A 250       2.969   4.099  -6.655  1.00  0.00           H
ATOM   7495 H2'1 DC  A 250       1.848   6.899  -7.226  1.00  0.00           H
ATOM   7496 H2'2 DC  A 250       2.172   6.147  -5.629  1.00  0.00           H
ATOM   7497  H1' DC  A 250      -0.265   5.947  -6.785  1.00  0.00           H
ATOM   7498  H5  DC  A 250       0.669   1.091  -4.688  1.00  0.00           H
ATOM   7499  H6  DC  A 250       1.057   2.597  -6.587  1.00  0.00           H
ATOM   7500  H41 DC  A 250      -0.768   2.354  -1.545  1.00  0.00           H
ATOM   7501  H42 DC  A 250      -0.151   0.965  -2.410  1.00  0.00           H
TER
ATOM  23255  P   DG  B 751      -5.462   5.230   6.840  1.00  1.00           P
ATOM  23256  O1P DG  B 751      -6.677   5.617   7.592  1.00  1.00           O
ATOM  23257  O2P DG  B 751      -5.445   3.851   6.302  1.00  1.00           O
ATOM  23258  O5' DG  B 751      -5.208   6.273   5.655  1.00  1.00           O
ATOM  23259  C5' DG  B 751      -4.914   7.644   5.981  1.00  1.00           C
ATOM  23260  C4' DG  B 751      -4.585   8.418   4.723  1.00  1.00           C
ATOM  23261  O4' DG  B 751      -3.246   8.047   4.282  1.00  1.00           O
ATOM  23262  C3' DG  B 751      -5.467   8.122   3.508  1.00  1.00           C
ATOM  23263  O3' DG  B 751      -6.717   8.795   3.589  1.00  1.00           O
ATOM  23264  C2' DG  B 751      -4.535   8.624   2.407  1.00  1.00           C
ATOM  23265  C1' DG  B 751      -3.204   8.036   2.865  1.00  1.00           C
ATOM  23266  N9  DG  B 751      -2.974   6.641   2.396  1.00  1.00           N
ATOM  23267  C8  DG  B 751      -3.226   5.451   3.046  1.00  1.00           C
ATOM  23268  N7  DG  B 751      -2.906   4.386   2.351  1.00  1.00           N
ATOM  23269  C5  DG  B 751      -2.406   4.907   1.159  1.00  1.00           C
ATOM  23270  C6  DG  B 751      -1.901   4.240   0.012  1.00  1.00           C
ATOM  23271  O6  DG  B 751      -1.789   3.034  -0.191  1.00  1.00           O
ATOM  23272  N1  DG  B 751      -1.499   5.152  -0.974  1.00  1.00           N
ATOM  23273  C2  DG  B 751      -1.575   6.527  -0.866  1.00  1.00           C
ATOM  23274  N2  DG  B 751      -1.139   7.223  -1.919  1.00  1.00           N
ATOM  23275  N3  DG  B 751      -2.049   7.151   0.209  1.00  1.00           N
ATOM  23276  C4  DG  B 751      -2.444   6.281   1.176  1.00  1.00           C
ATOM  23277 H5'1 DG  B 751      -4.061   7.678   6.659  1.00  0.00           H
ATOM  23278 H5'2 DG  B 751      -5.780   8.094   6.466  1.00  0.00           H
ATOM  23279  H4' DG  B 751      -4.612   9.486   4.938  1.00  0.00           H
ATOM  23280  H3' DG  B 751      -5.623   7.048   3.409  1.00  0.00           H
ATOM  23281 H2'1 DG  B 751      -4.502   9.713   2.376  1.00  0.00           H
ATOM  23282 H2'2 DG  B 751      -4.825   8.217   1.438  1.00  0.00           H
ATOM  23283  H1' DG  B 751      -2.389   8.673   2.522  1.00  0.00           H
ATOM  23284  H8  DG  B 751      -3.649   5.401   4.038  1.00  0.00           H
ATOM  23285  H1  DG  B 751      -1.124   4.773  -1.832  1.00  0.00           H
ATOM  23286  H21 DG  B 751      -0.781   6.740  -2.730  1.00  0.00           H
ATOM  23287  H22 DG  B 751      -1.167   8.232  -1.904  1.00  0.00           H
TER                                                                             
")

writeLines(template.in, paste0(output.path, "/template.pdb"))
    
nab.in <- paste0(
'molecule mol_temp, mol1, mol2, mol_test;
residue sense_res, anti_res;
matrix mat_offset, mat_rise;
string res_lib;
string seq, anti_seq, sense_res_name, anti_res_name;
string loup[ hashed ],phos_suger_atoms_array[18];
int i, j, seq_length;
float rise, twist, total_rise, total_twist;
point pts1[dynamic], pts2[dynamic];
point temp_sense_5cap_pos, temp_sense_3cap_pos, temp_anti_5cap_pos,  temp_anti_3cap_pos;

//###################
//###### Input ######
//###################
seq="', seq, '"; 
res_lib="all_nucleic94.lib";
anti_seq = wc_complement( seq, res_lib, "dna" );

seq_length = length( seq );
loup["g"] = "G"; loup["a"] = "A"; loup["t"] = "T"; loup["c"] = "C";
loup["G"] = "G"; loup["A"] = "A"; loup["T"] = "T"; loup["C"] = "C";

mol_temp = getpdb("./template.pdb");

phos_suger_atoms_array[1]="P";
phos_suger_atoms_array[2]="O1P";
phos_suger_atoms_array[3]="O2P";
phos_suger_atoms_array[4]="O5\'";
phos_suger_atoms_array[5]="C5\'";
phos_suger_atoms_array[6]="H5\'1";
phos_suger_atoms_array[7]="H5\'2";
phos_suger_atoms_array[8]="C4\'";
phos_suger_atoms_array[9]="H4\'";
phos_suger_atoms_array[10]="O4\'";
phos_suger_atoms_array[11]="C1\'";
phos_suger_atoms_array[12]="H1\'";
phos_suger_atoms_array[13]="C3\'";
phos_suger_atoms_array[14]="H3\'";
phos_suger_atoms_array[15]="C2\'";
phos_suger_atoms_array[16]="H2\'1";
phos_suger_atoms_array[17]="H2\'2";
phos_suger_atoms_array[18]="O3\'";


temp_sense_5cap_pos.x = 3.307;  temp_sense_5cap_pos.y = 1.705;  temp_sense_5cap_pos.z =-8.401; //Cap positions obtained by addH of UCSF Chimera to templates.
temp_sense_3cap_pos.x = 4.589;  temp_sense_3cap_pos.y = 5.716;  temp_sense_3cap_pos.z =-6.985;
temp_anti_5cap_pos.x=-5.975;    temp_anti_5cap_pos.y = 5.982;    temp_anti_5cap_pos.z =  6.154;   
temp_anti_3cap_pos.x=-7.243;    temp_anti_3cap_pos.y = 8.588;    temp_anti_3cap_pos.z = 2.813;





//######################################
//###### handle the first residue ######
//######################################
//read residue
sense_res_name = "D"+loup[substr(seq, 1 ,1 )]+"5";
anti_res_name = "D"+loup[substr(anti_seq, 1 ,1 )]+"3";
sense_res = getresidue( sense_res_name, res_lib );
anti_res = getresidue( anti_res_name, res_lib );
mol1 = wc_basepair( sense_res, anti_res );


//adjust base pair position to template base pair position
if (loup[substr(seq, 1 ,1 )]=="A" || loup[substr(seq, 1 ,1 )]=="G"){
    setframe(1, mol1, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9", "sense:1:C8");
} else if (loup[substr(seq, 1 ,1 )]=="T" || loup[substr(seq, 1 ,1 )]=="C"){ 
    setframe(1, mol1, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1", "sense:1:C6");
}

setframe(1, mol_temp, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
alignframe(mol1, mol_temp);


//adjust backbone position to template backbone position
allocate pts1[3*mol1.natoms];
allocate pts2[3*mol1.natoms];

for(j = 1; j <= 18; j = j + 1){
   setxyz_from_mol(mol_temp,"1:1:"+phos_suger_atoms_array[j], pts1);
   setmol_from_xyz(mol1,"sense:1:"+phos_suger_atoms_array[j], pts1);
   setxyz_from_mol(mol_temp,"2:1:"+phos_suger_atoms_array[j], pts2);
   setmol_from_xyz(mol1,"anti:1:"+phos_suger_atoms_array[j], pts2);
   }

setmol_from_xyz(mol1, "sense:1:H5T", temp_sense_5cap_pos);
setmol_from_xyz(mol1, "anti:1:H3T",  temp_anti_3cap_pos);


//adjust rise and twist   
rise= 0.0;
twist=-0.0;

total_rise=rise;
total_twist=twist;


//######################################
//###### handle the second~ residue ######
//######################################
for(i = 2; i <= seq_length-1; i = i + 1){
  //read residue
  sense_res_name = "D"+loup[substr(seq, i ,1 )];
  anti_res_name = "D"+loup[substr(anti_seq, i ,1 )];
  sense_res = getresidue( sense_res_name, res_lib );
  anti_res = getresidue( anti_res_name, res_lib );
  mol2 = wc_basepair( sense_res, anti_res );

  if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
      setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
  } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){ 
      setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
  }

  setframe(1, mol_temp, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
  alignframe(mol2, mol_temp);

  allocate pts1[3*mol2.natoms];
  allocate pts2[3*mol2.natoms];

  for(j = 1; j <= 18; j = j + 1){
      setxyz_from_mol(mol_temp,"1:1:"+phos_suger_atoms_array[j], pts1);
      setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
      setxyz_from_mol(mol_temp,"2:1:"+phos_suger_atoms_array[j], pts2);
      setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
   }
         
   rise  = 2.548; 
   twist = -32.7;
   total_rise=total_rise+rise;
   total_twist=total_twist+twist;
 

  mat_rise=newtransform(total_rise, 0, 0, total_twist, 0, 0);
  transformmol(mat_rise, mol2, NULL);
  mergestr(mol1,"sense","last",mol2,"sense","first");
  connectres(mol1, "sense", i-1, "O3\'", i, "P");
  mergestr(mol1,"anti","first",mol2,"anti","last");
  connectres(mol1, "anti", 1, "O3\'", 2, "P");
  freemolecule(mol2);
}


//######################################
//###### handle the final residue ######
//######################################

i = seq_length;

if( i > 1 ){
  //read residue
  sense_res_name = "D"+loup[substr(seq, i ,1 )]+"3";
  anti_res_name = "D"+loup[substr(anti_seq, i ,1 )]+"5";
  sense_res = getresidue( sense_res_name, res_lib );
  anti_res = getresidue( anti_res_name, res_lib );
  mol2 = wc_basepair( sense_res, anti_res );

  
   if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
       setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
   } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){ 
       setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
   }

   setframe(1, mol_temp, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
   alignframe(mol2, mol_temp);

   allocate pts1[3*mol2.natoms];
   allocate pts2[3*mol2.natoms];

   for(j = 1; j <= 18; j = j + 1){
      setxyz_from_mol(mol_temp,"1:1:"+phos_suger_atoms_array[j], pts1);
      setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
      setxyz_from_mol(mol_temp,"2:1:"+phos_suger_atoms_array[j], pts2);
      setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
   }
         
   setmol_from_xyz(mol2, "sense:1:H3T", temp_sense_3cap_pos);
   setmol_from_xyz(mol2, "anti:1:H5T",  temp_anti_5cap_pos);

   rise  = 2.548; 
   twist = -32.7;
   total_rise=total_rise+rise;
   total_twist=total_twist+twist;


  mat_rise=newtransform(total_rise, 0, 0, total_twist, 0, 0);
  transformmol(mat_rise, mol2, NULL);
  mergestr(mol1,"sense","last",mol2,"sense","first");
  connectres(mol1, "sense", i-1, "O3\'", i, "P");
  mergestr(mol1,"anti","first",mol2,"anti","last");
  connectres(mol1, "anti", 1, "O3\'", 2, "P");
  freemolecule(mol2);
}


putpdb("', seq, '.pdb",mol1);'

)
    
writeLines(nab.in, paste0(output.path, "/model_build.nab"))
    
system(paste0("cd ", output.path,"; nab model_build.nab && ./model_build"))
    
file.remove(paste0(output.path, "/template.pdb"))

if (strand =="first_strand") {
  temp <- readLines(paste0(output.path,dna.seq,".pdb"))
  writeLines(temp[1:grep("TER",temp)[1]],paste0(output.path,dna.seq,".pdb"))
} else if (strand =="second_strand") {
  temp <- readLines(paste0(output.path,dna.seq,".pdb"))
  writeLines(temp[(grep("TER",temp)[1]+1):grep("TER",temp)[2]],paste0(output.path,dna.seq,".pdb"))
} 
    
}##end of conformation =="DNA_A_Handbook1999"   
  
if (input.conformation =="DNA_B") {
  nab.in <- paste0(

  'molecule m;
  m = fd_helix( "abdna", "', seq, '", "dna" );
  putpdb( "', seq, '.pdb", m, "-wwpdb");'

  )

  writeLines(nab.in, paste0(output.path, "/fd_helix.nab"))

  system(paste0("cd ", output.path,"; nab fd_helix.nab && ./fd_helix"))
  
  if (strand =="first_strand") {
    temp <- readLines(paste0(output.path,dna.seq,".pdb"))
    writeLines(temp[1:grep("TER",temp)[1]],paste0(output.path,dna.seq,".pdb"))
  } else if (strand =="second_strand") {
    temp <- readLines(paste0(output.path,dna.seq,".pdb"))
    writeLines(temp[(grep("TER",temp)[1]+1):grep("TER",temp)[2]],paste0(output.path,dna.seq,".pdb"))
  } 

}##end of conformation =="B"

if (input.conformation =="DNA_B_Handbook1999") {
##################################################################################
### B-DNA by Oxford handbook of nucleic acid structure (1999) by S. Neidle.    ###
### Template pdb was made by 3DNA (1 repeating of sequence: cccccccccccccccccccccccccc )###
### Hydrogens were attached by addH of Chimera.                                ###
### Furthermore, we modified H2',H2'',H5',H5'' to H2'1,H2'2,H5'1,H5'2          ###
##################################################################################
  
  template.in <- paste0(
    "ATOM    362  P   DC  A  13      -3.573   2.823   8.772  1.00  1.00           P
ATOM    363  O1P DC  A  13      -4.329   3.122  10.009  1.00  1.00           O
ATOM    364  O2P DC  A  13      -2.465   3.748   8.449  1.00  1.00           O
ATOM    365  O5' DC  A  13      -3.019   1.322   8.830  1.00  1.00           O
ATOM    366  C5' DC  A  13      -3.870   0.247   8.391  1.00  1.00           C
ATOM    367  C4' DC  A  13      -3.079  -0.723   7.536  1.00  1.00           C
ATOM    368  O4' DC  A  13      -3.076  -0.349   6.128  1.00  1.00           O
ATOM    369  C3' DC  A  13      -1.597  -0.860   7.883  1.00  1.00           C
ATOM    370  O3' DC  A  13      -1.205  -2.216   7.705  1.00  1.00           O
ATOM    371  C2' DC  A  13      -0.874  -0.005   6.842  1.00  1.00           C
ATOM    372  C1' DC  A  13      -1.740  -0.406   5.653  1.00  1.00           C
ATOM    373  N1  DC  A  13      -1.624   0.516   4.489  1.00  1.00           N
ATOM    374  C2  DC  A  13      -1.715  -0.034   3.212  1.00  1.00           C
ATOM    375  O2  DC  A  13      -1.887  -1.252   3.099  1.00  1.00           O
ATOM    376  N3  DC  A  13      -1.610   0.791   2.138  1.00  1.00           N
ATOM    377  C4  DC  A  13      -1.424   2.106   2.304  1.00  1.00           C
ATOM    378  N4  DC  A  13      -1.328   2.869   1.225  1.00  1.00           N
ATOM    379  C5  DC  A  13      -1.327   2.692   3.608  1.00  1.00           C
ATOM    380  C6  DC  A  13      -1.433   1.852   4.668  1.00  1.00           C
ATOM    381 H5'1 DC  A  13      -4.266  -0.278   9.260  1.00  0.00           H
ATOM    382 H5'2 DC  A  13      -4.695   0.654   7.807  1.00  0.00           H
ATOM    383  H4' DC  A  13      -3.540  -1.707   7.624  1.00  0.00           H
ATOM    384  H3' DC  A  13      -1.391  -0.515   8.896  1.00  0.00           H
ATOM    385 H2'1 DC  A  13       0.170  -0.293   6.716  1.00  0.00           H
ATOM    386 H2'2 DC  A  13      -0.970   1.058   7.064  1.00  0.00           H
ATOM    387  H1' DC  A  13      -1.499  -1.424   5.345  1.00  0.00           H
ATOM    388  H5  DC  A  13      -1.176   3.753   3.741  1.00  0.00           H
ATOM    389  H6  DC  A  13      -1.365   2.249   5.670  1.00  0.00           H
ATOM    390  H41 DC  A  13      -1.396   2.455   0.306  1.00  0.00           H
ATOM    391  H42 DC  A  13      -1.187   3.864   1.321  1.00  0.00           H
TER
ATOM   1214  P   DG  B  40       0.121  -2.145  -8.983  1.00  1.00           P
ATOM   1215  O1P DG  B  40       0.878  -2.530 -10.195  1.00  1.00           O
ATOM   1216  O2P DG  B  40      -0.987  -1.187  -9.186  1.00  1.00           O
ATOM   1217  O5' DG  B  40      -0.431  -3.458  -8.253  1.00  1.00           O
ATOM   1218  C5' DG  B  40       0.420  -4.148  -7.319  1.00  1.00           C
ATOM   1219  C4' DG  B  40      -0.370  -4.535  -6.085  1.00  1.00           C
ATOM   1220  O4' DG  B  40      -0.373  -3.484  -5.076  1.00  1.00           O
ATOM   1221  C3' DG  B  40      -1.852  -4.833  -6.311  1.00  1.00           C
ATOM   1222  O3' DG  B  40      -2.243  -5.900  -5.454  1.00  1.00           O
ATOM   1223  C2' DG  B  40      -2.576  -3.563  -5.865  1.00  1.00           C
ATOM   1224  C1' DG  B  40      -1.709  -3.287  -4.640  1.00  1.00           C
ATOM   1225  N9  DG  B  40      -1.827  -1.894  -4.124  1.00  1.00           N
ATOM   1226  C8  DG  B  40      -2.020  -0.721  -4.820  1.00  1.00           C
ATOM   1227  N7  DG  B  40      -2.080   0.346  -4.061  1.00  1.00           N
ATOM   1228  C5  DG  B  40      -1.916  -0.155  -2.771  1.00  1.00           C
ATOM   1229  C6  DG  B  40      -1.891   0.522  -1.524  1.00  1.00           C
ATOM   1230  O6  DG  B  40      -2.010   1.722  -1.297  1.00  1.00           O
ATOM   1231  N1  DG  B  40      -1.700  -0.372  -0.460  1.00  1.00           N
ATOM   1232  C2  DG  B  40      -1.553  -1.740  -0.585  1.00  1.00           C
ATOM   1233  N2  DG  B  40      -1.381  -2.418   0.552  1.00  1.00           N
ATOM   1234  N3  DG  B  40      -1.576  -2.374  -1.756  1.00  1.00           N
ATOM   1235  C4  DG  B  40      -1.761  -1.521  -2.798  1.00  1.00           C
ATOM   1236 H5'1 DG  B  40       1.243  -3.494  -7.030  1.00  0.00           H
ATOM   1237 H5'2 DG  B  40       0.820  -5.047  -7.789  1.00  0.00           H
ATOM   1238  H4' DG  B  40       0.092  -5.421  -5.650  1.00  0.00           H
ATOM   1239  H3' DG  B  40      -2.058  -5.064  -7.356  1.00  0.00           H
ATOM   1240 H2'1 DG  B  40      -2.518  -2.766  -6.607  1.00  0.00           H
ATOM   1241 H2'2 DG  B  40      -3.608  -3.771  -5.584  1.00  0.00           H
ATOM   1242  H1' DG  B  40      -1.949  -3.997  -3.849  1.00  0.00           H
ATOM   1243  H8  DG  B  40      -2.114  -0.684  -5.895  1.00  0.00           H
ATOM   1244  H1  DG  B  40      -1.667   0.015   0.472  1.00  0.00           H
ATOM   1245  H21 DG  B  40      -1.267  -3.421   0.529  1.00  0.00           H
ATOM   1246  H22 DG  B  40      -1.365  -1.928   1.435  1.00  0.00           H
TER                                                                             
")
  
writeLines(template.in, paste0(output.path, "/template.pdb"))
  
nab.in <- paste0(
    'molecule mol_temp, mol1, mol2, mol_test;
residue sense_res, anti_res;
matrix mat_offset, mat_rise;
string res_lib;
string seq, anti_seq, sense_res_name, anti_res_name;
string loup[ hashed ],phos_suger_atoms_array[18];
int i, j, seq_length;
float rise, twist, total_rise, total_twist;
point pts1[dynamic], pts2[dynamic];
point temp_sense_5cap_pos, temp_sense_3cap_pos, temp_anti_5cap_pos,  temp_anti_3cap_pos;

//###################
//###### Input ######
//###################
seq="', seq, '"; 
res_lib="all_nucleic94.lib";
anti_seq = wc_complement( seq, res_lib, "dna" );

seq_length = length( seq );
loup["g"] = "G"; loup["a"] = "A"; loup["t"] = "T"; loup["c"] = "C";
loup["G"] = "G"; loup["A"] = "A"; loup["T"] = "T"; loup["C"] = "C";

mol_temp = getpdb("./template.pdb");

phos_suger_atoms_array[1]="P";
phos_suger_atoms_array[2]="O1P";
phos_suger_atoms_array[3]="O2P";
phos_suger_atoms_array[4]="O5\'";
phos_suger_atoms_array[5]="C5\'";
phos_suger_atoms_array[6]="H5\'1";
phos_suger_atoms_array[7]="H5\'2";
phos_suger_atoms_array[8]="C4\'";
phos_suger_atoms_array[9]="H4\'";
phos_suger_atoms_array[10]="O4\'";
phos_suger_atoms_array[11]="C1\'";
phos_suger_atoms_array[12]="H1\'";
phos_suger_atoms_array[13]="C3\'";
phos_suger_atoms_array[14]="H3\'";
phos_suger_atoms_array[15]="C2\'";
phos_suger_atoms_array[16]="H2\'1";
phos_suger_atoms_array[17]="H2\'2";
phos_suger_atoms_array[18]="O3\'";


temp_sense_5cap_pos.x =-3.559;  temp_sense_5cap_pos.y = 2.049;  temp_sense_5cap_pos.z = 9.147; //Cap positions obtained by addH of UCSF Chimera to templates.
temp_sense_3cap_pos.x =-1.651;  temp_sense_3cap_pos.y =-2.766;  temp_sense_3cap_pos.z = 8.353;
temp_anti_5cap_pos.x= 0.109;    temp_anti_5cap_pos.y =-3.000;    temp_anti_5cap_pos.z = -8.901;   
temp_anti_3cap_pos.x=-1.796;    temp_anti_3cap_pos.y =-6.706;    temp_anti_3cap_pos.z = -5.722;





//######################################
//###### handle the first residue ######
//######################################
//read residue
sense_res_name = "D"+loup[substr(seq, 1 ,1 )]+"5";
anti_res_name = "D"+loup[substr(anti_seq, 1 ,1 )]+"3";
sense_res = getresidue( sense_res_name, res_lib );
anti_res = getresidue( anti_res_name, res_lib );
mol1 = wc_basepair( sense_res, anti_res );


//adjust base pair position to template base pair position
if (loup[substr(seq, 1 ,1 )]=="A" || loup[substr(seq, 1 ,1 )]=="G"){
    setframe(1, mol1, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9", "sense:1:C8");
} else if (loup[substr(seq, 1 ,1 )]=="T" || loup[substr(seq, 1 ,1 )]=="C"){ 
    setframe(1, mol1, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1", "sense:1:C6");
}

setframe(1, mol_temp, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
alignframe(mol1, mol_temp);


//adjust backbone position to template backbone position
allocate pts1[3*mol1.natoms];
allocate pts2[3*mol1.natoms];

for(j = 1; j <= 18; j = j + 1){
   setxyz_from_mol(mol_temp,"1:1:"+phos_suger_atoms_array[j], pts1);
   setmol_from_xyz(mol1,"sense:1:"+phos_suger_atoms_array[j], pts1);
   setxyz_from_mol(mol_temp,"2:1:"+phos_suger_atoms_array[j], pts2);
   setmol_from_xyz(mol1,"anti:1:"+phos_suger_atoms_array[j], pts2);
   }

setmol_from_xyz(mol1, "sense:1:H5T", temp_sense_5cap_pos);
setmol_from_xyz(mol1, "anti:1:H3T",  temp_anti_3cap_pos);


//adjust rise and twist   
rise= 0.0;
twist=-0.0;

total_rise=rise;
total_twist=twist;


//######################################
//###### handle the second~ residue ######
//######################################
for(i = 2; i <= seq_length-1; i = i + 1){
  //read residue
  sense_res_name = "D"+loup[substr(seq, i ,1 )];
  anti_res_name = "D"+loup[substr(anti_seq, i ,1 )];
  sense_res = getresidue( sense_res_name, res_lib );
  anti_res = getresidue( anti_res_name, res_lib );
  mol2 = wc_basepair( sense_res, anti_res );

  if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
      setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
  } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){ 
      setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
  }

  setframe(1, mol_temp, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
  alignframe(mol2, mol_temp);

  allocate pts1[3*mol2.natoms];
  allocate pts2[3*mol2.natoms];

  for(j = 1; j <= 18; j = j + 1){
      setxyz_from_mol(mol_temp,"1:1:"+phos_suger_atoms_array[j], pts1);
      setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
      setxyz_from_mol(mol_temp,"2:1:"+phos_suger_atoms_array[j], pts2);
      setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
   }
         
   rise  = 3.375; 
   twist = -36.0;
   total_rise=total_rise+rise;
   total_twist=total_twist+twist;
 

  mat_rise=newtransform(total_rise, 0, 0, total_twist, 0, 0);
  transformmol(mat_rise, mol2, NULL);
  mergestr(mol1,"sense","last",mol2,"sense","first");
  connectres(mol1, "sense", i-1, "O3\'", i, "P");
  mergestr(mol1,"anti","first",mol2,"anti","last");
  connectres(mol1, "anti", 1, "O3\'", 2, "P");
  freemolecule(mol2);
}


//######################################
//###### handle the final residue ######
//######################################

i = seq_length;

if( i > 1 ){
  //read residue
  sense_res_name = "D"+loup[substr(seq, i ,1 )]+"3";
  anti_res_name = "D"+loup[substr(anti_seq, i ,1 )]+"5";
  sense_res = getresidue( sense_res_name, res_lib );
  anti_res = getresidue( anti_res_name, res_lib );
  mol2 = wc_basepair( sense_res, anti_res );

  
   if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
       setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
   } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){ 
       setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
   }

   setframe(1, mol_temp, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
   alignframe(mol2, mol_temp);

   allocate pts1[3*mol2.natoms];
   allocate pts2[3*mol2.natoms];

   for(j = 1; j <= 18; j = j + 1){
      setxyz_from_mol(mol_temp,"1:1:"+phos_suger_atoms_array[j], pts1);
      setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
      setxyz_from_mol(mol_temp,"2:1:"+phos_suger_atoms_array[j], pts2);
      setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
   }
         
   setmol_from_xyz(mol2, "sense:1:H3T", temp_sense_3cap_pos);
   setmol_from_xyz(mol2, "anti:1:H5T",  temp_anti_5cap_pos);

   rise  = 3.375; 
   twist = -36.0;
   total_rise=total_rise+rise;
   total_twist=total_twist+twist;


  mat_rise=newtransform(total_rise, 0, 0, total_twist, 0, 0);
  transformmol(mat_rise, mol2, NULL);
  mergestr(mol1,"sense","last",mol2,"sense","first");
  connectres(mol1, "sense", i-1, "O3\'", i, "P");
  mergestr(mol1,"anti","first",mol2,"anti","last");
  connectres(mol1, "anti", 1, "O3\'", 2, "P");
  freemolecule(mol2);
}


putpdb("', seq, '.pdb",mol1);'

  )
  
writeLines(nab.in, paste0(output.path, "/model_build.nab"))
  
system(paste0("cd ", output.path,"; nab model_build.nab && ./model_build"))
  
file.remove(paste0(output.path, "/template.pdb"))
  
if (strand =="first_strand") {
  temp <- readLines(paste0(output.path,dna.seq,".pdb"))
  writeLines(temp[1:grep("TER",temp)[1]],paste0(output.path,dna.seq,".pdb"))
} else if (strand =="second_strand") {
  temp <- readLines(paste0(output.path,dna.seq,".pdb"))
  writeLines(temp[(grep("TER",temp)[1]+1):grep("TER",temp)[2]],paste0(output.path,dna.seq,".pdb"))
} 
  
}##end of conformation =="DNA_B_Handbook1999"     

if (input.conformation =="DNA_Z1_wang_1981") {
##################################################################################
### Z1-dna by Wang et al. Science, vol. 211, pp.171-176 (1981).                ###
### Template pdb can be found at https://modelarchive.org/doi/10.5452/ma-cz2oa ###
### Note that a figure at the left in this site is not related to Z-DNA.       ###
##################################################################################
  
template_CG.in <- paste0(
"ATOM      1  P     C A   1       6.780   2.740  -3.700  1.00  0.00           P  
ATOM      2  O1P   C A   1       6.390   2.010  -4.890  1.00  0.00           O  
ATOM      3  O2P   C A   1       8.170   2.680  -3.220  1.00  0.00           O  
ATOM      4  O5'   C A   1       5.790   2.330  -2.510  1.00  0.00           O 
ATOM      5  C5'   C A   1       6.240   1.500  -1.450  1.00  0.00           C  
ATOM      6  C4'   C A   1       5.150   0.510  -1.120  1.00  0.00           C  
ATOM      7  O4'   C A   1       3.990   1.230  -0.650  1.00  0.00           O  
ATOM      8  C3'   C A   1       4.660  -0.350  -2.320  1.00  0.00           C  
ATOM      9  O3'   C A   1       4.480  -1.700  -1.890  1.00  0.00           O  
ATOM     10  C2'   C A   1       3.330   0.250  -2.730  1.00  0.00           C  
ATOM     11  C1'   C A   1       2.830   0.790  -1.360  1.00  0.00           C  
ATOM     12  N1    C A   1       1.920   1.960  -1.460  1.00  0.00           N  
ATOM     13  C2    C A   1       0.570   1.680  -1.700  1.00  0.00           C  
ATOM     14  O2    C A   1       0.170   0.520  -1.820  1.00  0.00           O  
ATOM     15  N3    C A   1      -0.300   2.720  -1.800  1.00  0.00           N  
ATOM     16  C4    C A   1       0.120   3.990  -1.670  1.00  0.00           C  
ATOM     17  N4    C A   1      -0.770   4.980  -1.780  1.00  0.00           N  
ATOM     18  C5    C A   1       1.510   4.280  -1.420  1.00  0.00           C  
ATOM     19  C6    C A   1       2.350   3.250  -1.330  1.00  0.00           C  
ATOM     20 H5'1   C A   1       6.450   2.090  -0.660  1.00  0.00           H  
ATOM     21 H5'2   C A   1       7.070   1.030  -1.750  1.00  0.00           H  
ATOM     22  H4'   C A   1       5.530  -0.070  -0.400  1.00  0.00           H  
ATOM     23  H3'   C A   1       5.300  -0.420  -3.080  1.00  0.00           H  
ATOM     24 H2'1   C A   1       3.460   0.980  -3.400  1.00  0.00           H  
ATOM     25 H2'2   C A   1       2.720  -0.450  -3.110  1.00  0.00           H  
ATOM     26  H1'   C A   1       2.330   0.080  -0.850  1.00  0.00           H  
ATOM     27  H41   C A   1      -0.700   5.900  -1.720  1.00  0.00           H  
ATOM     28  H42   C A   1      -1.500   4.810  -1.910  1.00  0.00           H  
ATOM     29  H5    C A   1       1.840   5.220  -1.320  1.00  0.00           H  
ATOM     30  H6    C A   1       3.320   3.420  -1.160  1.00  0.00           H  
ATOM    725  P     G B  12      -5.610  -2.810   2.030  1.00  0.00           P  
ATOM    726  O1P   G B  12      -6.360  -2.510   3.270  1.00  0.00           O  
ATOM    727  O2P   G B  12      -5.020  -4.120   1.840  1.00  0.00           O  
ATOM    728  O5'   G B  12      -6.570  -2.520   0.780  1.00  0.00           O  
ATOM    729  C5'   G B  12      -5.990  -2.330  -0.510  1.00  0.00           C  
ATOM    730  C4'   G B  12      -7.120  -2.110  -1.510  1.00  0.00           C  
ATOM    731  O4'   G B  12      -7.710  -0.820  -1.320  1.00  0.00           O  
ATOM    732  C3'   G B  12      -6.650  -2.140  -3.000  1.00  0.00           C  
ATOM    733  O3'   G B  12      -6.950  -3.440  -3.530  1.00  0.00           O  
ATOM    734  C2'   G B  12      -7.420  -1.070  -3.710  1.00  0.00           C  
ATOM    735  C1'   G B  12      -7.810  -0.090  -2.560  1.00  0.00           C  
ATOM    736  N9    G B  12      -6.960   1.050  -2.460  1.00  0.00           N  
ATOM    737  C8    G B  12      -7.410   2.320  -2.410  1.00  0.00           C  
ATOM    738  N7    G B  12      -6.440   3.210  -2.300  1.00  0.00           N  
ATOM    739  C5    G B  12      -5.290   2.450  -2.280  1.00  0.00           C  
ATOM    740  C6    G B  12      -3.930   2.850  -2.180  1.00  0.00           C  
ATOM    741  O6    G B  12      -3.500   4.020  -2.090  1.00  0.00           O  
ATOM    742  N1    G B  12      -3.050   1.790  -2.190  1.00  0.00           N  
ATOM    743  C2    G B  12      -3.430   0.490  -2.290  1.00  0.00           C  
ATOM    744  N2    G B  12      -2.510  -0.410  -2.290  1.00  0.00           N  
ATOM    745  N3    G B  12      -4.710   0.100  -2.380  1.00  0.00           N  
ATOM    746  C4    G B  12      -5.570   1.120  -2.370  1.00  0.00           C  
ATOM    747 H5'1   G B  12      -5.460  -3.140  -0.750  1.00  0.00           H  
ATOM    748 H5'2   G B  12      -5.390  -1.530  -0.470  1.00  0.00           H  
ATOM    749  H4'   G B  12      -7.800  -2.810  -1.320  1.00  0.00           H  
ATOM    750  H3'   G B  12      -5.670  -2.020  -3.150  1.00  0.00           H  
ATOM    751 H2'1   G B  12      -6.850  -0.630  -4.400  1.00  0.00           H  
ATOM    752 H2'2   G B  12      -8.230  -1.470  -4.150  1.00  0.00           H  
ATOM    753  H1'   G B  12      -8.720   0.300  -2.700  1.00  0.00           H  
ATOM    754  H3'   G B  12      -8.380   2.560  -2.440  1.00  0.00           H  
ATOM    755  H1    G B  12      -2.070   1.980  -2.130  1.00  0.00           H  
ATOM    756  H21   G B  12      -2.730  -1.290  -2.350  1.00  0.00           H  
ATOM    757  H22   G B  12      -1.810  -0.010  -2.230  1.00  0.00           H  
TER     758        G B  12                                                      
MASTER       81    0    0    0    0    0    0   24  756    2    0    2          
END                                                                             
")
writeLines(template_CG.in, paste0(output.path, "/template_CG.pdb"))
    
template_GC.in <- paste0(
  "ATOM     31  P     G A   2       5.610  -2.810  -2.030  1.00  0.00           P  
ATOM     32  O1P   G A   2       6.360  -2.510  -3.270  1.00  0.00           O  
ATOM     33  O2P   G A   2       5.020  -4.120  -1.840  1.00  0.00           O  
ATOM     34  O5'   G A   2       6.570  -2.520  -0.780  1.00  0.00           O  
ATOM     35  C5'   G A   2       5.990  -2.330   0.510  1.00  0.00           C  
ATOM     36  C4'   G A   2       7.120  -2.110   1.510  1.00  0.00           C  
ATOM     37  O4'   G A   2       7.710  -0.820   1.320  1.00  0.00           O  
ATOM     38  C3'   G A   2       6.650  -2.140   3.000  1.00  0.00           C  
ATOM     39  O3'   G A   2       6.950  -3.440   3.530  1.00  0.00           O  
ATOM     40  C2'   G A   2       7.420  -1.070   3.710  1.00  0.00           C  
ATOM     41  C1'   G A   2       7.810  -0.090   2.560  1.00  0.00           C  
ATOM     42  N9    G A   2       6.960   1.050   2.460  1.00  0.00           N  
ATOM     43  C8    G A   2       7.410   2.320   2.410  1.00  0.00           C  
ATOM     44  N7    G A   2       6.440   3.210   2.300  1.00  0.00           N  
ATOM     45  C5    G A   2       5.290   2.450   2.280  1.00  0.00           C  
ATOM     46  C6    G A   2       3.930   2.850   2.180  1.00  0.00           C  
ATOM     47  O6    G A   2       3.500   4.020   2.090  1.00  0.00           O  
ATOM     48  N1    G A   2       3.050   1.790   2.190  1.00  0.00           N  
ATOM     49  C2    G A   2       3.430   0.490   2.290  1.00  0.00           C  
ATOM     50  N2    G A   2       2.510  -0.410   2.290  1.00  0.00           N  
ATOM     51  N3    G A   2       4.710   0.100   2.380  1.00  0.00           N  
ATOM     52  C4    G A   2       5.570   1.120   2.370  1.00  0.00           C  
ATOM     53 H5'1   G A   2       5.460  -3.140   0.750  1.00  0.00           H  
ATOM     54 H5'2   G A   2       5.390  -1.530   0.470  1.00  0.00           H  
ATOM     55  H4'   G A   2       7.800  -2.810   1.320  1.00  0.00           H  
ATOM     56  H3'   G A   2       5.670  -2.020   3.150  1.00  0.00           H  
ATOM     57 H2'1   G A   2       6.850  -0.630   4.400  1.00  0.00           H  
ATOM     58 H2'2   G A   2       8.230  -1.470   4.150  1.00  0.00           H  
ATOM     59  H1'   G A   2       8.720   0.300   2.700  1.00  0.00           H  
ATOM     60  H3'   G A   2       8.380   2.560   2.440  1.00  0.00           H  
ATOM     61  H1    G A   2       2.070   1.980   2.130  1.00  0.00           H  
ATOM     62  H21   G A   2       2.730  -1.290   2.350  1.00  0.00           H  
ATOM     63  H22   G A   2       1.810  -0.010   2.230  1.00  0.00           H  
ATOM    695  P     C B  11      -6.780   2.740   3.700  1.00  0.00           P  
ATOM    696  O1P   C B  11      -8.170   2.680   3.220  1.00  0.00           O  
ATOM    697  O2P   C B  11      -6.390   2.010   4.890  1.00  0.00           O  
ATOM    698  O5'   C B  11      -5.790   2.330   2.510  1.00  0.00           O  
ATOM    699  C5'   C B  11      -6.240   1.500   1.450  1.00  0.00           C  
ATOM    700  C4'   C B  11      -5.150   0.510   1.120  1.00  0.00           C  
ATOM    701  O4'   C B  11      -3.990   1.230   0.650  1.00  0.00           O  
ATOM    702  C3'   C B  11      -4.660  -0.350   2.320  1.00  0.00           C  
ATOM    703  O3'   C B  11      -4.480  -1.700   1.890  1.00  0.00           O  
ATOM    704  C2'   C B  11      -3.330   0.250   2.730  1.00  0.00           C  
ATOM    705  C1'   C B  11      -2.830   0.790   1.360  1.00  0.00           C  
ATOM    706  N1    C B  11      -1.920   1.960   1.460  1.00  0.00           N  
ATOM    707  C2    C B  11      -0.570   1.680   1.700  1.00  0.00           C  
ATOM    708  O2    C B  11      -0.170   0.520   1.820  1.00  0.00           O  
ATOM    709  N3    C B  11       0.300   2.720   1.800  1.00  0.00           N  
ATOM    710  C4    C B  11      -0.120   3.990   1.670  1.00  0.00           C  
ATOM    711  N4    C B  11       0.770   4.980   1.780  1.00  0.00           N  
ATOM    712  C5    C B  11      -1.510   4.280   1.420  1.00  0.00           C  
ATOM    713  C6    C B  11      -2.350   3.250   1.330  1.00  0.00           C  
ATOM    714 H5'1   C B  11      -6.450   2.090   0.660  1.00  0.00           H  
ATOM    715 H5'2   C B  11      -7.070   1.030   1.750  1.00  0.00           H  
ATOM    716  H4'   C B  11      -5.530  -0.070   0.400  1.00  0.00           H  
ATOM    717  H3'   C B  11      -5.300  -0.420   3.080  1.00  0.00           H  
ATOM    718 H2'1   C B  11      -3.460   0.980   3.400  1.00  0.00           H  
ATOM    719 H2'2   C B  11      -2.720  -0.450   3.110  1.00  0.00           H  
ATOM    720  H1'   C B  11      -2.330   0.080   0.850  1.00  0.00           H  
ATOM    721  H41   C B  11       0.700   5.900   1.720  1.00  0.00           H  
ATOM    722  H42   C B  11       1.500   4.810   1.910  1.00  0.00           H  
ATOM    723  H5    C B  11      -1.840   5.220   1.320  1.00  0.00           H  
ATOM    724  H6    C B  11      -3.320   3.420   1.160  1.00  0.00           H                                                   
MASTER       81    0    0    0    0    0    0   24  756    2    0    2          
END                                                                                                                                                         
")
writeLines(template_GC.in, paste0(output.path, "/template_GC.pdb"))

nab.in <- paste0(
'molecule mol_zig, mol_zag, mol1, mol2, mol_test;
residue sense_res, anti_res;
matrix mat_offset, mat_rise;
string res_lib;
string seq, anti_seq, sense_res_name, anti_res_name;
string loup[ hashed ],phos_suger_atoms_array[18];
int i, j, seq_length;
float rise, twist, total_rise, total_twist;
point pts1[dynamic], pts2[dynamic];
point zig_sense_5cap_pos, zig_sense_3cap_pos, zig_anti_5cap_pos,  zig_anti_3cap_pos;
point zag_sense_5cap_pos, zag_sense_3cap_pos, zag_anti_5cap_pos,  zag_anti_3cap_pos;

//###################
//###### Input ######
//###################
seq="',seq,'"; 
res_lib="all_nucleic94.lib";
anti_seq = wc_complement( seq, res_lib, "dna" );

seq_length = length( seq );
loup["g"] = "G"; loup["a"] = "A"; loup["t"] = "T"; loup["c"] = "C";
loup["G"] = "G"; loup["A"] = "A"; loup["T"] = "T"; loup["C"] = "C";

mol_zig = getpdb("./template_CG.pdb");
mol_zag = getpdb("./template_GC.pdb");

phos_suger_atoms_array[1]="P";
phos_suger_atoms_array[2]="O1P";
phos_suger_atoms_array[3]="O2P";
phos_suger_atoms_array[4]="O5\'";
phos_suger_atoms_array[5]="C5\'";
phos_suger_atoms_array[6]="H5\'1";
phos_suger_atoms_array[7]="H5\'2";
phos_suger_atoms_array[8]="C4\'";
phos_suger_atoms_array[9]="H4\'";
phos_suger_atoms_array[10]="O4\'";
phos_suger_atoms_array[11]="C1\'";
phos_suger_atoms_array[12]="H1\'";
phos_suger_atoms_array[13]="C3\'";
phos_suger_atoms_array[14]="H3\'";
phos_suger_atoms_array[15]="C2\'";
phos_suger_atoms_array[16]="H2\'1";
phos_suger_atoms_array[17]="H2\'2";
phos_suger_atoms_array[18]="O3\'";


zig_sense_5cap_pos.x= 6.474;  zig_sense_5cap_pos.y= 2.966;  zig_sense_5cap_pos.z=-2.731;   //Cap positions obtained by addH of UCSF Chimera to templates.
zig_sense_3cap_pos.x= 4.178;  zig_sense_3cap_pos.y=-2.234;  zig_sense_3cap_pos.z=-2.628;
zig_anti_5cap_pos.x =-5.874;  zig_anti_5cap_pos.y =-2.661;  zig_anti_5cap_pos.z = 1.427;
zig_anti_3cap_pos.x =-6.670;  zig_anti_3cap_pos.y =-3.485;  zig_anti_3cap_pos.z =-4.447;

zag_sense_5cap_pos.x= 5.874;  zag_sense_5cap_pos.y=-2.661;  zag_sense_5cap_pos.z=-1.427;   
zag_sense_3cap_pos.x= 6.670;  zag_sense_3cap_pos.y=-3.485;  zag_sense_3cap_pos.z= 4.447;
zag_anti_5cap_pos.x =-6.474;  zag_anti_5cap_pos.y = 2.966;  zag_anti_5cap_pos.z = 2.731;
zag_anti_3cap_pos.x =-4.178;  zag_anti_3cap_pos.y =-2.234;  zag_anti_3cap_pos.z = 2.628;


//######################################
//###### handle the first residue ######
//######################################
//read residue
sense_res_name = "D"+loup[substr(seq, 1 ,1 )]+"5";
anti_res_name = "D"+loup[substr(anti_seq, 1 ,1 )]+"3";
sense_res = getresidue( sense_res_name, res_lib );
anti_res = getresidue( anti_res_name, res_lib );
mol1 = wc_basepair( sense_res, anti_res );


//adjust base pair position to template base pair position
if (loup[substr(seq, 1 ,1 )]=="A" || loup[substr(seq, 1 ,1 )]=="G"){
    setframe(1, mol1, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9", "sense:1:C8");
} else if (loup[substr(seq, 1 ,1 )]=="T" || loup[substr(seq, 1 ,1 )]=="C"){ 
    setframe(1, mol1, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1", "sense:1:C6");
}

setframe(1, mol_zig, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
alignframe(mol1, mol_zig);


//adjust backbone position to template backbone position
allocate pts1[3*mol1.natoms];
allocate pts2[3*mol1.natoms];

for(j = 1; j <= 18; j = j + 1){
   setxyz_from_mol(mol_zig,"1:1:"+phos_suger_atoms_array[j], pts1);
   setmol_from_xyz(mol1,"sense:1:"+phos_suger_atoms_array[j], pts1);
   setxyz_from_mol(mol_zig,"2:1:"+phos_suger_atoms_array[j], pts2);
   setmol_from_xyz(mol1,"anti:1:"+phos_suger_atoms_array[j], pts2);
   }

setmol_from_xyz(mol1, "sense:1:H5T", zig_sense_5cap_pos);
setmol_from_xyz(mol1, "anti:1:H3T",  zig_anti_3cap_pos);


//adjust rise and twist   
rise= 0.0;
twist=-0.0;

total_rise=rise;
total_twist=twist;


//######################################
//###### handle the second~ residue ######
//######################################
for(i = 2; i <= seq_length-1; i = i + 1){
  //read residue
  sense_res_name = "D"+loup[substr(seq, i ,1 )];
  anti_res_name = "D"+loup[substr(anti_seq, i ,1 )];
  sense_res = getresidue( sense_res_name, res_lib );
  anti_res = getresidue( anti_res_name, res_lib );
  mol2 = wc_basepair( sense_res, anti_res );

  if (i%2==0){
      if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
          setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
      } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){
          setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
      }

      setframe(1, mol_zag, "1:1:N9", "1:1:N9", "2:1:N1", "1:1:N9","1:1:C8");
      alignframe(mol2, mol_zag);

      allocate pts1[3*mol2.natoms];
      allocate pts2[3*mol2.natoms];

      for(j = 1; j <= 18; j = j + 1){
        setxyz_from_mol(mol_zag,"1:1:"+phos_suger_atoms_array[j], pts1);
        setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
        setxyz_from_mol(mol_zag,"2:1:"+phos_suger_atoms_array[j], pts2);
        setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
      }

      rise= 0.0;
      twist=0.0;
      total_rise=total_rise+rise;
      total_twist=total_twist+twist;


  } else if (i%2!=0){
        if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
            setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
        } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){ 
            setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
        }

        setframe(1, mol_zig, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
        alignframe(mol2, mol_zig);

        allocate pts1[3*mol2.natoms];
        allocate pts2[3*mol2.natoms];

        for(j = 1; j <= 18; j = j + 1){
           setxyz_from_mol(mol_zig,"1:1:"+phos_suger_atoms_array[j], pts1);
           setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
           setxyz_from_mol(mol_zig,"2:1:"+phos_suger_atoms_array[j], pts2);
           setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
        }
         
        rise= 7.43; //=44.58/6
        twist=-60.0;
        total_rise=total_rise+rise;
        total_twist=total_twist+twist;
  }

  mat_rise=newtransform(0, 0, total_rise, 0, 0, total_twist);
  transformmol(mat_rise, mol2, NULL);
  mergestr(mol1,"sense","last",mol2,"sense","first");
  connectres(mol1, "sense", i-1, "O3\'", i, "P");
  mergestr(mol1,"anti","first",mol2,"anti","last");
  connectres(mol1, "anti", 1, "O3\'", 2, "P");
  freemolecule(mol2);
}


//######################################
//###### handle the final residue ######
//######################################

i = seq_length;

if( i > 1 ){
  //read residue
  sense_res_name = "D"+loup[substr(seq, i ,1 )]+"3";
  anti_res_name = "D"+loup[substr(anti_seq, i ,1 )]+"5";
  sense_res = getresidue( sense_res_name, res_lib );
  anti_res = getresidue( anti_res_name, res_lib );
  mol2 = wc_basepair( sense_res, anti_res );

  if (i%2==0){
      if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
          setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
      } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){
          setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
      }

      setframe(1, mol_zag, "1:1:N9", "1:1:N9", "2:1:N1", "1:1:N9","1:1:C8");
      alignframe(mol2, mol_zag);

      allocate pts1[3*mol2.natoms];
      allocate pts2[3*mol2.natoms];

      for(j = 1; j <= 18; j = j + 1){
        setxyz_from_mol(mol_zag,"1:1:"+phos_suger_atoms_array[j], pts1);
        setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
        setxyz_from_mol(mol_zag,"2:1:"+phos_suger_atoms_array[j], pts2);
        setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
      }

      setmol_from_xyz(mol2, "sense:1:H3T", zag_sense_3cap_pos);
      setmol_from_xyz(mol2, "anti:1:H5T",  zag_anti_5cap_pos);

      rise= 0.0;
      twist=0.0;
      total_rise=total_rise+rise;
      total_twist=total_twist+twist;


  } else if (i%2!=0){
        if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
            setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
        } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){ 
            setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
        }

        setframe(1, mol_zig, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
        alignframe(mol2, mol_zig);

        allocate pts1[3*mol2.natoms];
        allocate pts2[3*mol2.natoms];

        for(j = 1; j <= 18; j = j + 1){
           setxyz_from_mol(mol_zig,"1:1:"+phos_suger_atoms_array[j], pts1);
           setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
           setxyz_from_mol(mol_zig,"2:1:"+phos_suger_atoms_array[j], pts2);
           setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
        }
         
        setmol_from_xyz(mol2, "sense:1:H3T", zig_sense_3cap_pos);
        setmol_from_xyz(mol2, "anti:1:H5T",  zig_anti_5cap_pos);

        rise= 7.43; //=44.58/6
        twist=-60.0;
        total_rise=total_rise+rise;
        total_twist=total_twist+twist;
  }

  mat_rise=newtransform(0, 0, total_rise, 0, 0, total_twist);
  transformmol(mat_rise, mol2, NULL);
  mergestr(mol1,"sense","last",mol2,"sense","first");
  connectres(mol1, "sense", i-1, "O3\'", i, "P");
  mergestr(mol1,"anti","first",mol2,"anti","last");
  connectres(mol1, "anti", 1, "O3\'", 2, "P");
  freemolecule(mol2);
}

putpdb("',seq,'.pdb",mol1);'
  
)


writeLines(nab.in, paste0(output.path, "/fd_helix.nab"))

system(paste0("cd ", output.path,"; nab fd_helix.nab && ./fd_helix"))

file.remove(paste0(output.path, "/template_CG.pdb"))
file.remove(paste0(output.path, "/template_GC.pdb"))

if (strand =="first_strand") {
  temp <- readLines(paste0(output.path,dna.seq,".pdb"))
  writeLines(temp[1:grep("TER",temp)[1]],paste0(output.path,dna.seq,".pdb"))
} else if (strand =="second_strand") {
  temp <- readLines(paste0(output.path,dna.seq,".pdb"))
  writeLines(temp[(grep("TER",temp)[1]+1):grep("TER",temp)[2]],paste0(output.path,dna.seq,".pdb"))
} 

}##end of conformation =="Z1_wang_1981"  

if (input.conformation =="DNA_Z_Handbook1999") {
##################################################################################
### Z-DNA by Oxford handbook of nucleic acid structure (1999) by S. Neidle.    ###
### Template pdb was made by 3DNA (100 repeating )                             ###
### Hydrogens were attached by addH of Chimera.                                ###
### Furthermore, we modified H2',H2'',H5',H5'' to H2'1,H2'2,H5'1,H5'2          ###
##################################################################################
  
  template_CG.in <- paste0(
    "ATOM   3122  P   DC  A 100      -4.484  -6.673  -3.014  1.00  1.00           P
ATOM   3123  O1P DC  A 100      -5.678  -6.317  -2.216  1.00  1.00           O
ATOM   3124  O2P DC  A 100      -3.991  -8.061  -2.864  1.00  1.00           O
ATOM   3125  O5' DC  A 100      -3.299  -5.646  -2.695  1.00  1.00           O
ATOM   3126  C5' DC  A 100      -2.207  -6.077  -1.861  1.00  1.00           C
ATOM   3127  C4' DC  A 100      -1.853  -4.988  -0.867  1.00  1.00           C
ATOM   3128  O4' DC  A 100      -1.255  -3.856  -1.561  1.00  1.00           O
ATOM   3129  C3' DC  A 100      -3.022  -4.417  -0.061  1.00  1.00           C
ATOM   3130  O3' DC  A 100      -2.614  -4.260   1.292  1.00  1.00           O
ATOM   3131  C2' DC  A 100      -3.258  -3.038  -0.674  1.00  1.00           C
ATOM   3132  C1' DC  A 100      -1.850  -2.649  -1.114  1.00  1.00           C
ATOM   3133  N1  DC  A 100      -1.823  -1.690  -2.221  1.00  1.00           N
ATOM   3134  C2  DC  A 100      -1.831  -0.335  -1.906  1.00  1.00           C
ATOM   3135  O2  DC  A 100      -1.860   0.001  -0.716  1.00  1.00           O
ATOM   3136  N3  DC  A 100      -1.806   0.570  -2.918  1.00  1.00           N
ATOM   3137  C4  DC  A 100      -1.776   0.166  -4.194  1.00  1.00           C
ATOM   3138  N4  DC  A 100      -1.753   1.087  -5.145  1.00  1.00           N
ATOM   3139  C5  DC  A 100      -1.768  -1.225  -4.540  1.00  1.00           C
ATOM   3140  C6  DC  A 100      -1.792  -2.113  -3.513  1.00  1.00           C
ATOM   3141 H5'1 DC  A 100      -1.339  -6.292  -2.485  1.00  0.00           H
ATOM   3142 H5'2 DC  A 100      -2.497  -6.979  -1.322  1.00  0.00           H
ATOM   3143  H4' DC  A 100      -1.118  -5.388  -0.168  1.00  0.00           H
ATOM   3144  H3' DC  A 100      -3.908  -5.047  -0.141  1.00  0.00           H
ATOM   3145 H2'1 DC  A 100      -3.651  -2.338   0.063  1.00  0.00           H
ATOM   3146 H2'2 DC  A 100      -3.925  -3.107  -1.533  1.00  0.00           H
ATOM   3147  H1' DC  A 100      -1.292  -2.252  -0.266  1.00  0.00           H
ATOM   3148  H5  DC  A 100      -1.744  -1.551  -5.569  1.00  0.00           H
ATOM   3149  H6  DC  A 100      -1.786  -3.172  -3.726  1.00  0.00           H
ATOM   3150  H41 DC  A 100      -1.758   2.067  -4.899  1.00  0.00           H
ATOM   3151  H42 DC  A 100      -1.730   0.810  -6.116  1.00  0.00           H
TER
ATOM   9455  P   DG  B 301       2.970   5.408   2.348  1.00  1.00           P
ATOM   9456  O1P DG  B 301       2.848   4.873   3.722  1.00  1.00           O
ATOM   9457  O2P DG  B 301       4.275   6.014   2.000  1.00  1.00           O
ATOM   9458  O5' DG  B 301       1.798   6.463   2.077  1.00  1.00           O
ATOM   9459  C5' DG  B 301       0.432   6.007   2.090  1.00  1.00           C
ATOM   9460  C4' DG  B 301      -0.504   7.172   1.837  1.00  1.00           C
ATOM   9461  O4' DG  B 301      -0.334   7.610   0.456  1.00  1.00           O
ATOM   9462  C3' DG  B 301      -1.997   6.867   1.969  1.00  1.00           C
ATOM   9463  O3' DG  B 301      -2.490   7.185   3.264  1.00  1.00           O
ATOM   9464  C2' DG  B 301      -2.679   7.691   0.878  1.00  1.00           C
ATOM   9465  C1' DG  B 301      -1.571   7.981  -0.129  1.00  1.00           C
ATOM   9466  N9  DG  B 301      -1.716   7.231  -1.410  1.00  1.00           N
ATOM   9467  C8  DG  B 301      -1.877   7.716  -2.689  1.00  1.00           C
ATOM   9468  N7  DG  B 301      -1.976   6.782  -3.605  1.00  1.00           N
ATOM   9469  C5  DG  B 301      -1.872   5.595  -2.880  1.00  1.00           C
ATOM   9470  C6  DG  B 301      -1.907   4.249  -3.325  1.00  1.00           C
ATOM   9471  O6  DG  B 301      -2.040   3.818  -4.467  1.00  1.00           O
ATOM   9472  N1  DG  B 301      -1.766   3.353  -2.255  1.00  1.00           N
ATOM   9473  C2  DG  B 301      -1.611   3.712  -0.931  1.00  1.00           C
ATOM   9474  N2  DG  B 301      -1.492   2.708  -0.060  1.00  1.00           N
ATOM   9475  N3  DG  B 301      -1.578   4.977  -0.516  1.00  1.00           N
ATOM   9476  C4  DG  B 301      -1.714   5.859  -1.540  1.00  1.00           C
ATOM   9477 H5'1 DG  B 301       0.294   5.257   1.311  1.00  0.00           H
ATOM   9478 H5'2 DG  B 301       0.207   5.566   3.061  1.00  0.00           H
ATOM   9479  H4' DG  B 301      -0.245   7.992   2.507  1.00  0.00           H
ATOM   9480  H3' DG  B 301      -2.164   5.808   1.773  1.00  0.00           H
ATOM   9481 H2'1 DG  B 301      -3.484   7.122   0.412  1.00  0.00           H
ATOM   9482 H2'2 DG  B 301      -3.065   8.622   1.293  1.00  0.00           H
ATOM   9483  H1' DG  B 301      -1.556   9.050  -0.343  1.00  0.00           H
ATOM   9484  H8  DG  B 301      -1.918   8.771  -2.917  1.00  0.00           H
ATOM   9485  H1  DG  B 301      -1.779   2.366  -2.470  1.00  0.00           H
ATOM   9486  H21 DG  B 301      -1.518   1.752  -0.385  1.00  0.00           H
ATOM   9487  H22 DG  B 301      -1.376   2.904   0.924  1.00  0.00           H
TER                                                                             
")

writeLines(template_CG.in, paste0(output.path, "/template_CG.pdb"))
  
template_GC.in <- paste0(
    "ATOM   3152  P   DG  A 101      -2.970  -5.408   2.348  1.00  1.00           P
ATOM   3153  O1P DG  A 101      -2.848  -4.873   3.722  1.00  1.00           O
ATOM   3154  O2P DG  A 101      -4.275  -6.014   2.000  1.00  1.00           O
ATOM   3155  O5' DG  A 101      -1.798  -6.463   2.077  1.00  1.00           O
ATOM   3156  C5' DG  A 101      -0.432  -6.007   2.090  1.00  1.00           C
ATOM   3157  C4' DG  A 101       0.504  -7.172   1.837  1.00  1.00           C
ATOM   3158  O4' DG  A 101       0.334  -7.610   0.456  1.00  1.00           O
ATOM   3159  C3' DG  A 101       1.997  -6.867   1.969  1.00  1.00           C
ATOM   3160  O3' DG  A 101       2.490  -7.185   3.264  1.00  1.00           O
ATOM   3161  C2' DG  A 101       2.679  -7.691   0.878  1.00  1.00           C
ATOM   3162  C1' DG  A 101       1.571  -7.981  -0.129  1.00  1.00           C
ATOM   3163  N9  DG  A 101       1.716  -7.231  -1.410  1.00  1.00           N
ATOM   3164  C8  DG  A 101       1.877  -7.716  -2.689  1.00  1.00           C
ATOM   3165  N7  DG  A 101       1.976  -6.782  -3.605  1.00  1.00           N
ATOM   3166  C5  DG  A 101       1.872  -5.595  -2.880  1.00  1.00           C
ATOM   3167  C6  DG  A 101       1.907  -4.249  -3.325  1.00  1.00           C
ATOM   3168  O6  DG  A 101       2.040  -3.818  -4.467  1.00  1.00           O
ATOM   3169  N1  DG  A 101       1.766  -3.353  -2.255  1.00  1.00           N
ATOM   3170  C2  DG  A 101       1.611  -3.712  -0.931  1.00  1.00           C
ATOM   3171  N2  DG  A 101       1.492  -2.708  -0.060  1.00  1.00           N
ATOM   3172  N3  DG  A 101       1.578  -4.977  -0.516  1.00  1.00           N
ATOM   3173  C4  DG  A 101       1.714  -5.859  -1.540  1.00  1.00           C
ATOM   3174 H5'1 DG  A 101      -0.294  -5.257   1.311  1.00  0.00           H
ATOM   3175 H5'2 DG  A 101      -0.207  -5.566   3.061  1.00  0.00           H
ATOM   3176  H4' DG  A 101       0.245  -7.992   2.507  1.00  0.00           H
ATOM   3177  H3' DG  A 101       2.164  -5.808   1.773  1.00  0.00           H
ATOM   3178 H2'1 DG  A 101       3.484  -7.122   0.412  1.00  0.00           H
ATOM   3179 H2'2 DG  A 101       3.065  -8.622   1.293  1.00  0.00           H
ATOM   3180  H1' DG  A 101       1.556  -9.050  -0.343  1.00  0.00           H
ATOM   3181  H8  DG  A 101       1.918  -8.771  -2.917  1.00  0.00           H
ATOM   3182  H1  DG  A 101       1.779  -2.366  -2.470  1.00  0.00           H
ATOM   3183  H21 DG  A 101       1.518  -1.752  -0.385  1.00  0.00           H
ATOM   3184  H22 DG  A 101       1.376  -2.904   0.924  1.00  0.00           H
TER
ATOM   9425  P   DC  B 300       4.484   6.673  -3.014  1.00  1.00           P
ATOM   9426  O1P DC  B 300       5.678   6.317  -2.216  1.00  1.00           O
ATOM   9427  O2P DC  B 300       3.991   8.061  -2.864  1.00  1.00           O
ATOM   9428  O5' DC  B 300       3.299   5.646  -2.695  1.00  1.00           O
ATOM   9429  C5' DC  B 300       2.207   6.077  -1.861  1.00  1.00           C
ATOM   9430  C4' DC  B 300       1.853   4.988  -0.867  1.00  1.00           C
ATOM   9431  O4' DC  B 300       1.255   3.856  -1.561  1.00  1.00           O
ATOM   9432  C3' DC  B 300       3.022   4.417  -0.061  1.00  1.00           C
ATOM   9433  O3' DC  B 300       2.614   4.260   1.292  1.00  1.00           O
ATOM   9434  C2' DC  B 300       3.258   3.038  -0.674  1.00  1.00           C
ATOM   9435  C1' DC  B 300       1.850   2.649  -1.114  1.00  1.00           C
ATOM   9436  N1  DC  B 300       1.823   1.690  -2.221  1.00  1.00           N
ATOM   9437  C2  DC  B 300       1.831   0.335  -1.906  1.00  1.00           C
ATOM   9438  O2  DC  B 300       1.860  -0.001  -0.716  1.00  1.00           O
ATOM   9439  N3  DC  B 300       1.806  -0.570  -2.918  1.00  1.00           N
ATOM   9440  C4  DC  B 300       1.776  -0.166  -4.194  1.00  1.00           C
ATOM   9441  N4  DC  B 300       1.753  -1.087  -5.145  1.00  1.00           N
ATOM   9442  C5  DC  B 300       1.768   1.225  -4.540  1.00  1.00           C
ATOM   9443  C6  DC  B 300       1.792   2.113  -3.513  1.00  1.00           C
ATOM   9444 H5'1 DC  B 300       2.499   6.978  -1.321  1.00  0.00           H
ATOM   9445 H5'2 DC  B 300       1.339   6.294  -2.484  1.00  0.00           H
ATOM   9446  H4' DC  B 300       1.118   5.388  -0.168  1.00  0.00           H
ATOM   9447  H3' DC  B 300       3.908   5.047  -0.141  1.00  0.00           H
ATOM   9448 H2'1 DC  B 300       3.651   2.338   0.063  1.00  0.00           H
ATOM   9449 H2'2 DC  B 300       3.925   3.107  -1.533  1.00  0.00           H
ATOM   9450  H1' DC  B 300       1.292   2.252  -0.266  1.00  0.00           H
ATOM   9451  H5  DC  B 300       1.744   1.551  -5.569  1.00  0.00           H
ATOM   9452  H6  DC  B 300       1.786   3.172  -3.726  1.00  0.00           H
ATOM   9453  H41 DC  B 300       1.758  -2.067  -4.899  1.00  0.00           H
ATOM   9454  H42 DC  B 300       1.730  -0.810  -6.116  1.00  0.00           H
TER                                                                                                                                                         
")

writeLines(template_GC.in, paste0(output.path, "/template_GC.pdb"))
  
nab.in <- paste0(
    'molecule mol_zig, mol_zag, mol1, mol2, mol_test;
residue sense_res, anti_res;
matrix mat_offset, mat_rise;
string res_lib;
string seq, anti_seq, sense_res_name, anti_res_name;
string loup[ hashed ],phos_suger_atoms_array[18];
int i, j, seq_length;
float rise, twist, total_rise, total_twist;
point pts1[dynamic], pts2[dynamic];
point zig_sense_5cap_pos, zig_sense_3cap_pos, zig_anti_5cap_pos,  zig_anti_3cap_pos;
point zag_sense_5cap_pos, zag_sense_3cap_pos, zag_anti_5cap_pos,  zag_anti_3cap_pos;

//###################
//###### Input ######
//###################
seq="',seq,'"; 
res_lib="all_nucleic94.lib";
anti_seq = wc_complement( seq, res_lib, "dna" );

seq_length = length( seq );
loup["g"] = "G"; loup["a"] = "A"; loup["t"] = "T"; loup["c"] = "C";
loup["G"] = "G"; loup["A"] = "A"; loup["T"] = "T"; loup["C"] = "C";

mol_zig = getpdb("./template_CG.pdb");
mol_zag = getpdb("./template_GC.pdb");

phos_suger_atoms_array[1]="P";
phos_suger_atoms_array[2]="O1P";
phos_suger_atoms_array[3]="O2P";
phos_suger_atoms_array[4]="O5\'";
phos_suger_atoms_array[5]="C5\'";
phos_suger_atoms_array[6]="H5\'1";
phos_suger_atoms_array[7]="H5\'2";
phos_suger_atoms_array[8]="C4\'";
phos_suger_atoms_array[9]="H4\'";
phos_suger_atoms_array[10]="O4\'";
phos_suger_atoms_array[11]="C1\'";
phos_suger_atoms_array[12]="H1\'";
phos_suger_atoms_array[13]="C3\'";
phos_suger_atoms_array[14]="H3\'";
phos_suger_atoms_array[15]="C2\'";
phos_suger_atoms_array[16]="H2\'1";
phos_suger_atoms_array[17]="H2\'2";
phos_suger_atoms_array[18]="O3\'";


zig_sense_5cap_pos.x =-3.517;  zig_sense_5cap_pos.y =-6.340;  zig_sense_5cap_pos.z =-3.321; //Cap positions obtained by addH of UCSF Chimera to templates.
zig_sense_3cap_pos.x =-3.341;  zig_sense_3cap_pos.y =-3.902;  zig_sense_3cap_pos.z = 1.806;
zig_anti_5cap_pos.x= 2.385;    zig_anti_5cap_pos.y= 5.721;    zig_anti_5cap_pos.z= 2.238;   
zig_anti_3cap_pos.x=-3.427;    zig_anti_3cap_pos.y= 6.981;    zig_anti_3cap_pos.z= 3.309;


zag_sense_5cap_pos.x=-2.385;  zag_sense_5cap_pos.y=-5.721;  zag_sense_5cap_pos.z= 2.238;   
zag_sense_3cap_pos.x= 3.427;  zag_sense_3cap_pos.y=-6.981;  zag_sense_3cap_pos.z= 3.309;
zag_anti_5cap_pos.x = 3.517;  zag_anti_5cap_pos.y = 6.340;  zag_anti_5cap_pos.z =-3.321;
zag_anti_3cap_pos.x = 3.341;  zag_anti_3cap_pos.y = 3.902;  zag_anti_3cap_pos.z = 1.806;


//######################################
//###### handle the first residue ######
//######################################
//read residue
sense_res_name = "D"+loup[substr(seq, 1 ,1 )]+"5";
anti_res_name = "D"+loup[substr(anti_seq, 1 ,1 )]+"3";
sense_res = getresidue( sense_res_name, res_lib );
anti_res = getresidue( anti_res_name, res_lib );
mol1 = wc_basepair( sense_res, anti_res );


//adjust base pair position to template base pair position
if (loup[substr(seq, 1 ,1 )]=="A" || loup[substr(seq, 1 ,1 )]=="G"){
    setframe(1, mol1, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9", "sense:1:C8");
} else if (loup[substr(seq, 1 ,1 )]=="T" || loup[substr(seq, 1 ,1 )]=="C"){ 
    setframe(1, mol1, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1", "sense:1:C6");
}

setframe(1, mol_zig, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
alignframe(mol1, mol_zig);


//adjust backbone position to template backbone position
allocate pts1[3*mol1.natoms];
allocate pts2[3*mol1.natoms];

for(j = 1; j <= 18; j = j + 1){
   setxyz_from_mol(mol_zig,"1:1:"+phos_suger_atoms_array[j], pts1);
   setmol_from_xyz(mol1,"sense:1:"+phos_suger_atoms_array[j], pts1);
   setxyz_from_mol(mol_zig,"2:1:"+phos_suger_atoms_array[j], pts2);
   setmol_from_xyz(mol1,"anti:1:"+phos_suger_atoms_array[j], pts2);
   }

setmol_from_xyz(mol1, "sense:1:H5T", zig_sense_5cap_pos);
setmol_from_xyz(mol1, "anti:1:H3T",  zig_anti_3cap_pos);


//adjust rise and twist   
rise= 0.0;
twist=-0.0;

total_rise=rise;
total_twist=twist;


//######################################
//###### handle the second~ residue ######
//######################################
for(i = 2; i <= seq_length-1; i = i + 1){
  //read residue
  sense_res_name = "D"+loup[substr(seq, i ,1 )];
  anti_res_name = "D"+loup[substr(anti_seq, i ,1 )];
  sense_res = getresidue( sense_res_name, res_lib );
  anti_res = getresidue( anti_res_name, res_lib );
  mol2 = wc_basepair( sense_res, anti_res );

  if (i%2==0){
      if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
          setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
      } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){
          setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
      }

      setframe(1, mol_zag, "1:1:N9", "1:1:N9", "2:1:N1", "1:1:N9","1:1:C8");
      alignframe(mol2, mol_zag);

      allocate pts1[3*mol2.natoms];
      allocate pts2[3*mol2.natoms];

      for(j = 1; j <= 18; j = j + 1){
        setxyz_from_mol(mol_zag,"1:1:"+phos_suger_atoms_array[j], pts1);
        setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
        setxyz_from_mol(mol_zag,"2:1:"+phos_suger_atoms_array[j], pts2);
        setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
      }

      rise= 0.0;
      twist=0.0;
      total_rise=total_rise+rise;
      total_twist=total_twist+twist;


  } else if (i%2!=0){
        if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
            setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
        } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){ 
            setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
        }

        setframe(1, mol_zig, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
        alignframe(mol2, mol_zig);

        allocate pts1[3*mol2.natoms];
        allocate pts2[3*mol2.natoms];

        for(j = 1; j <= 18; j = j + 1){
           setxyz_from_mol(mol_zig,"1:1:"+phos_suger_atoms_array[j], pts1);
           setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
           setxyz_from_mol(mol_zig,"2:1:"+phos_suger_atoms_array[j], pts2);
           setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
        }
         
        rise= 7.25; 
        twist=60.0;
        total_rise=total_rise+rise;
        total_twist=total_twist+twist;
  }

  mat_rise=newtransform(total_rise, 0, 0, total_twist, 0, 0);
  transformmol(mat_rise, mol2, NULL);
  mergestr(mol1,"sense","last",mol2,"sense","first");
  connectres(mol1, "sense", i-1, "O3\'", i, "P");
  mergestr(mol1,"anti","first",mol2,"anti","last");
  connectres(mol1, "anti", 1, "O3\'", 2, "P");
  freemolecule(mol2);
}


//######################################
//###### handle the final residue ######
//######################################

i = seq_length;

if( i > 1 ){
  //read residue
  sense_res_name = "D"+loup[substr(seq, i ,1 )]+"3";
  anti_res_name = "D"+loup[substr(anti_seq, i ,1 )]+"5";
  sense_res = getresidue( sense_res_name, res_lib );
  anti_res = getresidue( anti_res_name, res_lib );
  mol2 = wc_basepair( sense_res, anti_res );

  if (i%2==0){
      if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
          setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
      } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){
          setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
      }

      setframe(1, mol_zag, "1:1:N9", "1:1:N9", "2:1:N1", "1:1:N9","1:1:C8");
      alignframe(mol2, mol_zag);

      allocate pts1[3*mol2.natoms];
      allocate pts2[3*mol2.natoms];

      for(j = 1; j <= 18; j = j + 1){
        setxyz_from_mol(mol_zag,"1:1:"+phos_suger_atoms_array[j], pts1);
        setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
        setxyz_from_mol(mol_zag,"2:1:"+phos_suger_atoms_array[j], pts2);
        setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
      }

      setmol_from_xyz(mol2, "sense:1:H3T", zag_sense_3cap_pos);
      setmol_from_xyz(mol2, "anti:1:H5T",  zag_anti_5cap_pos);

      rise= 0.0;
      twist=0.0;
      total_rise=total_rise+rise;
      total_twist=total_twist+twist;


  } else if (i%2!=0){
        if (loup[substr(seq, i ,1 )]=="A" || loup[substr(seq, i ,1 )]=="G"){
            setframe(1, mol2, "sense:1:N9", "sense:1:N9", "anti:1:N1", "sense:1:N9","sense:1:C8");
        } else if (loup[substr(seq, i ,1 )]=="T" || loup[substr(seq, i ,1 )]=="C"){ 
            setframe(1, mol2, "sense:1:N1", "sense:1:N1", "anti:1:N9", "sense:1:N1","sense:1:C6");
        }

        setframe(1, mol_zig, "1:1:N1", "1:1:N1", "2:1:N9", "1:1:N1","1:1:C6");
        alignframe(mol2, mol_zig);

        allocate pts1[3*mol2.natoms];
        allocate pts2[3*mol2.natoms];

        for(j = 1; j <= 18; j = j + 1){
           setxyz_from_mol(mol_zig,"1:1:"+phos_suger_atoms_array[j], pts1);
           setmol_from_xyz(mol2,"sense:1:"+phos_suger_atoms_array[j], pts1);
           setxyz_from_mol(mol_zig,"2:1:"+phos_suger_atoms_array[j], pts2);
           setmol_from_xyz(mol2,"anti:1:"+phos_suger_atoms_array[j], pts2);
        }
         
        setmol_from_xyz(mol2, "sense:1:H3T", zig_sense_3cap_pos);
        setmol_from_xyz(mol2, "anti:1:H5T",  zig_anti_5cap_pos);

        rise= 7.25; 
        twist=60.0;
        total_rise=total_rise+rise;
        total_twist=total_twist+twist;
  }

  mat_rise=newtransform(total_rise, 0, 0, total_twist, 0, 0);
  transformmol(mat_rise, mol2, NULL);
  mergestr(mol1,"sense","last",mol2,"sense","first");
  connectres(mol1, "sense", i-1, "O3\'", i, "P");
  mergestr(mol1,"anti","first",mol2,"anti","last");
  connectres(mol1, "anti", 1, "O3\'", 2, "P");
  freemolecule(mol2);
}

putpdb("',seq,'.pdb",mol1);'

)
  

writeLines(nab.in, paste0(output.path, "/build_structure.nab"))
system(paste0("cd ", output.path,"; nab build_structure.nab && ./build_structure"))
  
file.remove(paste0(output.path, "/template_CG.pdb"))
file.remove(paste0(output.path, "/template_GC.pdb"))
  
if (strand =="first_strand") {
 temp <- readLines(paste0(output.path,dna.seq,".pdb"))
 writeLines(temp[1:grep("TER",temp)[1]],paste0(output.path,dna.seq,".pdb"))
} else if (strand =="second_strand") {
 temp <- readLines(paste0(output.path,dna.seq,".pdb"))
 writeLines(temp[(grep("TER",temp)[1]+1):grep("TER",temp)[2]],paste0(output.path,dna.seq,".pdb"))
} 
  
}##end of conformation =="DNA_Z_Handbook1999"

}##end of function buildDNA  

###### By bio3d ###
rmsd_bio3d<-function(x,y,z){
  library(bio3d)
  target<-read.pdb(x)
  ref<-read.pdb(y)
  if(z=="other_than_H"){rmsd_atoms<-atom.select(ref,"noh")}
  if(z=="all_atoms"){rmsd_atoms<-NULL}
  return(rmsd(a=target,b=ref,a.inds=rmsd_atoms,b.inds=rmsd_atoms))
}

###### By AmberTools path (much slower than bio3d) ###
 rmsd_nab<-function(x,y,z){
    ambertools.path <- "/opt/anaconda3/envs/AmberTools21"
    Sys.setenv(PATH = paste0(ambertools.path, "/bin:", Sys.getenv("PATH")))
    Sys.setenv(AMBERHOME = ambertools.path)
    
    if(z=="all_atoms") {target="::"}
    if(z=="other_than_H") {target="::[^H]*"}
    
    nab.in <- paste0(
   'molecule mol1, mol2;
    float r;
    mol1=getpdb("',x,'");
    mol2=getpdb("',y,'");
    rmsd(mol1, "',target,'",mol2,"',target,'",r);
    printf("%f",r);')
    writeLines(nab.in, "./rmsd.nab")
    system(paste0("nab rmsd.nab && ./rmsd"),intern=T)
 }
################################################################################
curve_plus <- function(output.path, sequence, mecha.inpdb, mecha.outfile, curve.path){
  
  system(paste0("cd ",output.path,'/',sequence,";",'rm ',mecha.outfile,'*.*',"\n"))
  writeLines(paste0(
    curve.path,'Cur+<<!',"\n",
    ' &inp file=',mecha.inpdb,"\n",
    ' lis=',mecha.outfile,"\n",
    ' lib=',curve.path,"standard,","\n",
    ' test = .t.,',"\n",
    ' &end',"\n",
    '2 1 -1 0 0',"\n",
    '1:',nchar(sequence),"\n",
    2*nchar(sequence),':',nchar(sequence)+1,"\n",
    '!'
  ), paste0(output.path,'/',sequence,'/',mecha.outfile,'.sh'))
  
  cmd <- paste0("cd ",output.path,'/',sequence,";",
                "chmod +x ", mecha.outfile,'.sh',";",
                "./", mecha.outfile,'.sh')
  
  print(paste0("Curve"," ",sequence))
  system(cmd)
}

################################################################################
curve_plus_comp <- function(output.path, sequence, mecha.inpdb, mecha.outfile, curve.path){
  
  system(paste0("cd ",output.path,'/',sequence,";",'rm ',mecha.outfile,'*.*',"\n"))
  writeLines(paste0(
    curve.path,'Cur+<<!',"\n",
    ' &inp file=',mecha.inpdb,"\n",
    ' lis=',mecha.outfile,"\n",
    ' lib=',curve.path,"standard,","\n",
    ' test = .t.,',"\n",
    ' &end',"\n",
    '2 1 -1 0 0',"\n",
    nchar(sequence)+1,':',2*nchar(sequence),"\n",
    nchar(sequence),':',"1","\n",
    '!'
  ), paste0(output.path,'/',sequence,'/',mecha.outfile,'.sh'))
  
  cmd <- paste0("cd ",output.path,'/',sequence,";",
                "chmod +x ", mecha.outfile,'.sh',";",
                "./", mecha.outfile,'.sh')
  
  print(paste0("Curve"," ",sequence))
  system(cmd)
}

################################################################################
BPAxis_read <- function(output.path, sequence, mecha.outfile){
  
  temp <- readLines(paste0(output.path,'/',sequence,'/',mecha.outfile,".lis"))
  # read (A) BP axis Parameters
  low= match("  (A) BP-Axis        Xdisp   Ydisp   Inclin    Tip  Ax-bend", temp)
  axis_param<-c()
  low=low+2
  i=1
  while(trimws(temp[low])!=""){
    num<-substr(temp[low],1,5)
    resnam1<-substr(temp[low],7,9)
    resnum1<-substr(temp[low],10,12)
    resnam2<-substr(temp[low],14,15)
    resnum2<-substr(temp[low],16,18)
    Xdisp <-substr(temp[low],20,26)
    Ydisp  <-substr(temp[low],28,34)
    Inclin <-substr(temp[low],36,42)
    Tip <-substr(temp[low],44,50)
    Axbend <-substr(temp[low],52,58)

    axis_param[as.numeric(i)]=
      paste(sequence, resnum1, resnam1, resnum2, resnam2, Xdisp, Ydisp, Inclin, Tip, Axbend, sep=",")
    low=low+1
    i=i+1
  }
  print(paste0("BPaxis"," ",sequence))
  return(axis_param)
}


################################################################################
intraBP_read <- function(output.path, sequence, mecha.outfile){
  temp <- readLines(paste0(output.path,'/',sequence,'/',mecha.outfile,".lis"))
  # read (B) Intra BP Parameters
  low= match("  Strands 1-2       Shear  Stretch Stagger  Buckle  Propel Opening", temp)
  intra_param<-c()
  low=low+2
  i=1
  while(trimws(temp[low])!=""){
    num<-substr(temp[low],1,5)
    resnam1<-substr(temp[low],7,9)
    resnum1<-substr(temp[low],10,12)
    resnam2<-substr(temp[low],14,15)
    resnum2<-substr(temp[low],16,18)
    shear <-substr(temp[low],20,26)
    stretch  <-substr(temp[low],28,34)
    stagger <-substr(temp[low],36,42)
    buckle <-substr(temp[low],44,50)
    propel <-substr(temp[low],52,58)
    opening <-substr(temp[low],60,66)
    
    intra_param[as.numeric(i)]=
      paste(sequence, resnum1, resnam1, resnum2, resnam2, shear, stretch, stagger, buckle, propel, opening, sep=",")
    low=low+1
    i=i+1
  }
  print(paste0("Intra"," ",sequence))
  return(intra_param)
}

################################################################################
interBP_read <- function(output.path, sequence, mecha.outfile){
  temp <- readLines(paste0(output.path,'/',sequence,'/',mecha.outfile,".lis"))
  # read (c) Inter BP Parameters
  low= match("  (C) Inter-BP       Shift   Slide    Rise    Tilt    Roll   Twist   H-Ris   H-Twi", temp)
  inter_param<-c()
  low=low+2
  i=1
  while(trimws(temp[low])!=""){
    num<-substr(temp[low],1,5)
    resnam1<-substr(temp[low],7,9)
    resnum1<-substr(temp[low],10,12)
    resnam2<-substr(temp[low],14,15)
    resnum2<-substr(temp[low],16,18)
    shift <-substr(temp[low],20,26)
    slide  <-substr(temp[low],28,34)
    rise <-substr(temp[low],36,42)
    tilt <-substr(temp[low],44,50)
    roll <-substr(temp[low],52,58)
    twist <-substr(temp[low],60,66)
    hris <-substr(temp[low],68,74)
    htwi <-substr(temp[low],76,82)
    
    inter_param[as.numeric(i)]=
      paste(sequence, resnum1, resnam1, resnum2, resnam2, shift, slide, rise, tilt, roll, twist, hris, htwi, sep=",")
    low=low+1
    i=i+1
  }
  
  print(paste0("Inter"," ",sequence))
  return(inter_param)
}


################################################################################
backbone_read <- function(output.path, sequence, mecha.outfile){
  
  temp <- readLines(paste0(output.path,'/',sequence,'/',mecha.outfile,".lis"))
  # read (D) Backbone Parameters
  low= grep("Strand 1     Alpha  Beta   Gamma  Delta  Epsil  Zeta   Chi    Phase  Ampli  Puckr ", temp)
  backbone_param<-c()
  low=low+2
  while(trimws(temp[low])!=""){
    resnam<-substr(temp[low],7,9)
    resnum<-substr(temp[low],10,12)
    alpha <-substr(temp[low],16,21)
    beta  <-substr(temp[low],23,28)
    gamma <-substr(temp[low],30,35)
    delta <-substr(temp[low],37,42)
    epsil <-substr(temp[low],44,49)
    zeta  <-substr(temp[low],51,56)
    chi   <-substr(temp[low],58,63)
    phase <-substr(temp[low],65,70)
    ampli <-substr(temp[low],72,77)
    puckr <-substr(temp[low],79,84)
    backbone_param[as.numeric(resnum)]=
      paste(sequence, resnum, resnam, alpha, beta, gamma, delta, epsil,
            zeta, chi, phase, ampli, puckr, sep=",")
    low=low+1
  }
  
  low= grep("Strand 2     Alpha  Beta   Gamma  Delta  Epsil  Zeta   Chi    Phase  Ampli  Puckr ", temp)
  #backbone_param<-c()
  low=low+2
  while(trimws(temp[low])!=""){
    resnam<-substr(temp[low],7,9)
    resnum<-substr(temp[low],10,12)
    alpha <-substr(temp[low],16,21)
    beta  <-substr(temp[low],23,28)
    gamma <-substr(temp[low],30,35)
    delta <-substr(temp[low],37,42)
    epsil <-substr(temp[low],44,49)
    zeta  <-substr(temp[low],51,56)
    chi   <-substr(temp[low],58,63)
    phase <-substr(temp[low],65,70)
    ampli <-substr(temp[low],72,77)
    puckr <-substr(temp[low],79,84)
    backbone_param[as.numeric(resnum)]=
      paste(sequence, resnum, resnam, alpha, beta, gamma, delta, epsil,
            zeta,chi,phase,ampli,puckr, sep=",")
    low=low+1
  }  
  print(paste0("Backbone"," ",sequence))
  return(backbone_param)
}


  
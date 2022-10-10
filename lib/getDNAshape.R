library(DNAshapeR)

getDNAshape<-function(sequence){
  file_name=paste0("temp.fa")
  file.create(file_name)
  write(sequence, file_name)
  target<-c('MGW','HelT','ProT','Roll','EP','Stretch','Tilt','Buckle','Shear','Opening','Rise','Shift','Stagger','Slide')
  pred <- getShape(file_name,shapeType = target)
  print(sequence)
  return(pred)
}







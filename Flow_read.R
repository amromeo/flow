olwd<-getwd()
setwd("H:\\Desktop\\ResidentsAutopsy\\Hematology Lab\\Database")
allfiles<-list.files(all.files=TRUE, pattern=".fso")
if (file.exists("filelist.txt")){
  oldfiles<-readLines("filelist.txt")
  newfiles<-setdiff(allfiles,oldfiles)
} else {
  newfiles<-allfiles
}

for(file in newfiles){
  tempdata<-scan(file,character(0),skipNul =TRUE,sep="/")
  olddata<-rbind(olddata,tempdata)
  rm(tempdata)
}
--------------------------
  
  junk<-scan("251619.FSO",character(0),skipNul =TRUE,sep="/")  
write.csv(x=junk,file="junk3.csv" )

write(allfiles,file = "filelist.txt",sep = "\t")
olddata<-mdstest[11:20,]

for(i in 1:10){
  write.csv(x=mdstest[i,],file=paste("file",i,".csv",sep=""))
}


?readFCSheader

-------------------------
library(purrr)
files<-dir(path='D:')  
fso<-readBin(readfile,what="raw",2000)
#replace all 'null' values with a space
fso[fso==00]<-as.raw(32)
#convert to characters
fsonew<-rawToChar(fso)
match<-regmatches(fsonew, regexpr('Sysmex_RET%/.+Sysmex_RET', fsonew,perl=T))
ret<-substr(match,13,17)
grepl("     ",ret)
match<-regmatches(fsonew, regexpr('SMNO/.+/*SYS', fsonew,perl=T))
barcode<-substr(match,6,16)


#Generating FCS files from fso files



fso<-readBin(readfile,what="raw",100000000)
#replace all 'null' values with a space
fso[fso==00]<-as.raw(32)
#convert to characters
fsonew<-rawToChar(fso)
#find the infiidvidual FCS files in the file
pos = gregexpr('FCS', fsonew)
pos<-pos[[1]]

#get the name of the test so that the fcs files can get good names
testid = gregexpr('fso', readfile)
testid<-testid[[1]]


#write the individual flow files
if(length(pos)==1){
  first<-pos[1]
  writeBin(substr(fsonew, first, nchar(fsonew)),paste(substr(readfile,testid-7,testid-2),"_1.fcs",sep=""))
}
if(length(pos)==2){
  first<-pos[1]
  second<-pos[2]
  writeBin(substr(fsonew, first, second-1),paste(substr(readfile,testid-7,testid-2),"_1.fcs",sep=""))
  writeBin(substr(fsonew, second, nchar(fsonew)),paste(substr(readfile,testid-7,testid-2),"_2.fcs",sep=""))
} 
if(length(pos)==3){
  first<-pos[1]
  second<-pos[2]
  third<-pos[3]
  writeBin(substr(fsonew, first, second-1),paste(substr(readfile,testid-7,testid-2),"_1.fcs",sep=""))
  writeBin(substr(fsonew, second, third-1),paste(substr(readfile,testid-7,testid-2),"_2.fcs",sep=""))
  writeBin(substr(fsonew, third, nchar(fsonew)),paste(substr(readfile,testid-7,testid-2),"_3.fcs",sep=""))
} 
if(length(pos)==4){
  first<-pos[1]
  second<-pos[2]
  third<-pos[3]
  fourth<-pos[4]
  writeBin(substr(fsonew, first, second-1),paste(substr(readfile,testid-7,testid-2),"_1.fcs",sep=""))
  writeBin(substr(fsonew, second, third-1),paste(substr(readfile,testid-7,testid-2),"_2.fcs",sep=""))
  writeBin(substr(fsonew, third, fourth),paste(substr(readfile,testid-7,testid-2),"_3.fcs",sep=""))
  writeBin(substr(fsonew, fourth, nchar(fsonew)),paste(substr(readfile,testid-7,testid-2),"_4.fcs",sep=""))
} 


-----------------------
  
  #start to flow
  library(flowCore)
library(flowViz)
library(ggplot2)
file.name <- choose.files()
x <- read.FCS(file.name, transformation=FALSE,emptyValue=FALSE )
plot(x)
plot(x,c('SFL','FSC'))

#In order to generate a scatterplot that looks like a regular flow plot you need to add a value "d" with the color value
y <- as.data.frame(exprs(x))
d = densCols(y$SFL, y$SSC, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
y<-(data.frame(y,d=d))
p <- ggplot(y) +
  geom_point(aes(SSC, SFL, col = d), size = 2,alpha=.5) +
  scale_color_identity() +
  theme_bw()+
  xlim(0,250)

#Scatterplot with right colors without changing the datafram
p <- ggplot(y) + geom_hex(aes(SFL, FSC), bins = 300) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6)))

#create flowset
y<-read.flowSet(files=file.name, names.keyword="SAMPLE ID", emptyValue = FALSE)

#create gate around presumed retics
reticgate<-polygonGate(filterId = "Retics",cbind("FSC"=c(25,40,40,250,250),"SFL"=c(50,150,500,500,75)))
rbcgate<-polygonGate(filterId = "RBC",cbind("FSC"=c(25,20,250,250),"SFL"=c(50,0,10,75)))
filters1 <- filters(list(reticgate, rbcgate))
xyplot(data=x,FSC~SFL, filter=filters1, smooth=FALSE, outline=TRUE,stat=TRUE,margin=TRUE)
summary(filter(x,rbcgate))
summary(filter(x,reticgate))


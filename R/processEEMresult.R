# frame_files <- lapply(sys.frames(), function(x) x$ofile)
# frame_files <- Filter(Negate(is.null), frame_files)
# setwd(dirname(frame_files[[length(frame_files)]]))
# 
# eemFile<-"../coadExp_hallmark.eem"
# expFile<-"../coadExp.tab"
# #eemFile<-"../testExp_testGs.eem"
# #expFile<-"../testExp.tab"
# Pcut <- 3
# treeCut <- 0.75
# minClusterSize<-0
maxClusterCount<-10

outdir<- paste(gsub(".eem","",eemFile),"posteem",sep=".")

tmp<-read.table(expFile, sep="\t")
s<-as.vector(t(tmp[1,2:ncol(tmp)]))
g<-as.vector(tmp[2:nrow(tmp),1])
tmp<-tmp[-1,]
tmp<-tmp[,-1]
E<-matrix(as.numeric(as.matrix(tmp)), ncol=length(s), nrow=length(g))
colnames(E)<-s
rownames(E)<-g

tmp<-strsplit(scan(eemFile, what="character", sep="\n"), "\t")
EEM<-list()
ids<-NULL
P<-NULL
for( i  in  1:length(tmp)){
 ids<-c(ids, tmp[[i]][1]) 
 p<-as.numeric(tmp[[i]][2])
 P<-c(P,p)
 stmp<-strsplit(tmp[[i]][4],"[A-Za-z]+=", perl=T)[[1]][6]
 mtmp<-strsplit(tmp[[i]][4],"[A-Za-z]+=", perl=T)[[1]][7]
 seedGenes<-strsplit(gsub('\\[|\\]| ',"",stmp,perl=T), ",")[[1]]
 moduleGenes<-strsplit(gsub('\\[|\\]| ',"",mtmp,perl=T), ",")[[1]]
  EEM<-c(EEM, list(c(list(p), list(seedGenes),  list(moduleGenes))))
}
names(EEM)<-ids
names(P)<-ids
getId<-function(x){names(EEM)[x]}
getP<-function(x){ EEM[[x]][[1]]}
getSeedGenes<-function(x){ EEM[[x]][[2]]}
getModuleGenes<-function(x){ EEM[[x]][[3]]}

library("gplots")
plotGeneHeatmap<-function(x){
  moduleG<-getModuleGenes(x)
  seedG<-getSeedGenes(x)
  
  seedE<-t(scale(t(E[seedG,])))
  moduleE<-seedE[moduleG,]
  moduleAct<-apply(moduleE,2,mean)
  sampleIds<-names(rev(sort(moduleAct)))
  moduleE<-moduleE[,sampleIds]
  seedE<-seedE[,sampleIds]
  moduleAct<-moduleAct[sampleIds]
  
  coln<-100
  col<-my.col(coln)
  #M<-max(moduleAct)
  #m<-min(moduleAct)
  M<-max(seedE)
  m<-min(seedE)
  moduleAct2<-as.integer((coln-1)*(moduleAct-m)/(M-m))+1
  sampCol<-col[moduleAct2]
  
  geneCol<-rep("gray90",nrow(seedE))
  names(geneCol)<-rownames(seedE)
  geneCol[moduleG]<-"red"
  
  heatmap(seedE, col=col, scale="none",
          hclustfun=function(x){hclust(x, method = "ward.D2")},
          #labRow=NA,labCol=NA,
          Colv=NA, 
          RowSideColors=geneCol,
          ColSideColors=sampCol)
}
getModuleAct<-function(x){
  moduleG<-getModuleGenes(x)
  seedG<-getSeedGenes(x)
  seedE<-t(scale(t(E[seedG,])))
  moduleE<-seedE[moduleG,]
  moduleAct<-apply(moduleE,2,mean)
  return(moduleAct)
}

my.var<-function(x) var(x)*(length(x)-1)/length(x)

my.scale <- function(x){
  if(is.vector(x)){
    (x-mean(x))/sqrt(my.var(x))
  }else{
    v<-sqrt(apply(x, 2, my.var))
    scale(x, center=TRUE, scale=v)
  }
}

my.dist <- function(x){
  y<-t(as.matrix(my.scale((t(x)))))
  d<-dist(y)
  (d^2)/(2*ncol(y))
}
#red1 <-  "red"
#blue1 <-  "blue"
red1 <- (colorpanel(20, low = "blue", high = "red", mid = "white"))[19]
blue1 <- (colorpanel(20, low = "blue", high = "red", mid = "white"))[2]
white <- "white"

my.col <- function(coln=100){
  colorpanel(coln, low = blue1, high = red1, mid = white)
}

red <-  "firebrick1"
blue <-  "dodgerblue"
orange <- "darkorange"
yellow <- "gold"
green <-  "darkolivegreen4"

my.col2 <- function(n){
  if(n==2){
    return( c(blue, red))
  }
  if(n==3){     
    return(c(blue,yellow,red))
  }
  if(n==4){     
    return(c(blue, green, orange, red))
  }
  if(n==5){     
    return(c(blue, green, yellow, orange, red))
  }
  if(n==6){
    return(colorpanel(7, low = blue, high = red, mid = yellow)[c(1:3,5:7)])
  }
  for(i in n:(n+1000)){
    tmp<-unique(colorpanel(i, low = blue, high = red, mid = yellow))
    if(length(tmp) == n){
      return(tmp)
    }
    if(length(tmp) == n+1){
      return(tmp[-((n/2)+1)])
    }
  }
  return(NULL);
}
mapColor <- function(B){
  lev<-rev(sort(as.numeric(levels(factor(B)))))
  C2<-NULL
  C <-  my.col2(length(lev))
  for(i in 1:length(B)){
    if(is.na(B[i])){
      C2<-c(C2,"White")
      next
    }
    for(j in 1:length(lev)){
      if(lev[j] == B[i]){
        C2<-c(C2, C[j])
      }
    }
  }
  C2
}
my.hclust<-function(D){
  #hclust(D,method="ward")
  hclust(D, method="complete")
}
collapse.matrix<-function(E,cluster){
  E2<-NULL	
  for(i in 1:length(cluster)){
    if(length(cluster[[i]])==1){
      E2<-rbind(E2, E[cluster[[i]],])
    }else{
      E2<-rbind(E2, apply(E[cluster[[i]],],2,mean))
    }
  }
  rownames(E2)<-names(cluster)
  return(E2)
}
printMatrix<-function(x,file){
  write(t(as.matrix(c("", colnames(x)))), file,  append=F, sep="\t", ncolumns=ncol(x)+1)
  write.table(x, file, quote=F, col.names=F, append=T, sep="\t")
}
printGmt<-function(x,file){
  for( i in 1:length(x)){	
    write(c(names(x)[i], length(x[[i]]), x[[i]] ), file, ncolumns=(length(x[[i]])+2), append=(if(i==1){F}else{T}), sep="\t")
  }
}

ids2<-ids[P>Pcut]
if(length(ids2)==0){
  quit()
}
dir.create(outdir)

for(i in 1:length(ids2)){
  #print(i)
  outfile <- paste(ids2[i], "pdf", sep=".")
  pdf(paste(outdir,outfile, sep="/"))
  plotGeneHeatmap(ids2[i])
  dev.off()
}

if(length(ids2)==1){
  quit()
}

M<-NULL
for(i in 1:length(ids2)){
  #print(i)
  M<-rbind(M,getModuleAct(ids2[i])) 
}
rownames(M)<-ids2
M<-t(scale(t(M)))
outfile <- paste("module", "pdf", sep=".")
pdf(paste(outdir,outfile, sep="/"))
heatmap(M, col=my.col(), scale="none",
         hclustfun=function(x){hclust(x, method = "ward.D2")})
dev.off()


E<-M
E<-t(my.scale(t(E)))
outfile<-paste(outdir,"module", sep="/")
printMatrix(E, paste(outfile, ".tab", sep=""))

D<-my.dist(E)
h<-my.hclust(D)
maxHeight <- max(h$height)

c<-cutree(h,h=1-treeCut)

if(length(levels(factor(c)))>maxClusterCount){
  c<-cutree(h,maxClusterCount)	
}

names(c)<-rownames(E)

# get cluster members
clusterId<-NULL	 
cluster<-NULL
clusterSize<-NULL
tmp<-as.numeric(levels(factor(c)))
for(i in 1:length(tmp)){	
  member<-rownames(E)[c==tmp[i]]
  cluster<-c(cluster, list(member))
  clusterId<-c(clusterId, paste("cluster", i, sep=""))
  clusterSize<-c(clusterSize, length(member))
}     

names(cluster)<-clusterId
names(clusterSize)<-clusterId

# sort by cluster size
clusterSize<-clusterSize[clusterSize>=minClusterSize]
cluster<-cluster[names(rev(sort(clusterSize)))]
names(cluster)<-clusterId[1:length(cluster)]
c<-NULL
tmp<-NULL
for(i in 1:length(cluster)){
  tmp<-c(tmp, cluster[[i]])
  c<-c(c, rep(i, length(cluster[[i]])))
}
names(c)<-tmp

if(length(cluster)!=1){
  E<-E[names(c),]
  pdf(paste(outfile, ".pdf", sep=""))
  heatmap(E, scale="none", hclustfun = my.hclust, RowSideColors = mapColor(c), col=my.col(100))
  dev.off()
  E2<-collapse.matrix(E,cluster)
  
  printMatrix(E2, paste(outfile, ".collapsed.tab", sep=""))
  
  pdf(paste(outfile, ".collapsed.pdf", sep=""))
  heatmap(E2, scale="none", hclustfun = my.hclust, col=my.col(100))
  dev.off()
}else{
  pdf(paste(outfile, ".pdf", sep=""))
  heatmap(E, scale="none", hclustfun = my.hclust, col=my.col(100))
  dev.off()
  E2<-collapse.matrix(E,cluster)
  
  printMatrix(E2, paste(outfile, ".collapsed.tab", sep=""))
}

printGmt(cluster, paste(outfile, ".gmt", sep=""))
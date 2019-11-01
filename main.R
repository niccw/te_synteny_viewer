library(gtools)
library(magrittr)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
# source te_matrix.R if use seperatly

plot_param_env <- new.env(parent = emptyenv())
plot_param_env$lim <- 4
plot_param_env$nmin <- 10
plot_param_env$shared <- TRUE
plot_param_env$par1 <- FALSE 
plot_param_env$par2 <- FALSE 
plot_param_env$PARTITION_ngene <- 50
plot_param_env$showlines <- TRUE 

# checked 15oct, this work :)
msynt_matrix_plot <- function(d1,d2,speciesA,speciesB,reactive_env,d_env){ # inherite param from plot_param_env
  # backup d1,d2, not overwrite when reloading parm for same pair of species
  d_env$d1_raw <- d1
  d_env$d2_raw <- d2
  
  # Prepare the data for basic msynt plot
  # The only input is .msynt from both species
  # Filtering: msynt by number of paralogs, scaffold with > n msynt
  
  # Filter msynt by size, filter scaffold by # of msynt ---------------------
  allsame = unique(as.character(d1[d1[,4] %in% d2[,4],4])) # intersect(d1$V4,d2$V4)
  if (plot_param_env$shared) {
    d1=d1[d1[,4] %in% allsame,]
    d2=d2[d2[,4] %in% allsame,]
  }
  
  d1=d1[d1[,5]<=plot_param_env$lim,]
  d2=d2[d2[,5]<=plot_param_env$lim,]
  
  # scaffold ID
  d1scf=as.character(unique(d1[,1])) # unique(d1$V1)
  d2scf=as.character(unique(d2[,1])) # unique(d2$V2)
  
  # mark down number of uniqe msynt on each scaffold
  d1len=NULL
  d2len=NULL
  
  for (i in d1scf) {
    d=d1[d1[,1]==i,]
    d1len=c(d1len,length(unique(as.character(d[,4]))))
  }
  for (i in d2scf) {
    d=d2[d2[,1]==i,]
    d2len=c(d2len,length(unique(as.character(d[,4]))))
  }
  
  # filter by number of msynt on scaffold (only consider scaffold with many msynt)
  d1scf=d1scf[d1len>=plot_param_env$nmin]
  d2scf=d2scf[d2len>=plot_param_env$nmin]
  d1len=d1len[d1len>=plot_param_env$nmin]
  d2len=d2len[d2len>=plot_param_env$nmin]
  
  d1=d1[d1[,1] %in% d1scf,]
  d2=d2[d2[,1] %in% d2scf,]
  
  labelx=d1scf
  labely=d2scf
  
  # Partioning --------------------------------------------------------------
  # The partition can be done by:
  # as.integer((seq_along(1:nrow(sel)) - 1) / ngene)
  
  if (plot_param_env$par1) {
    
    print("Partitioning 1 ! ... ")
    ngene=plot_param_env$PARTITION_ngene
    
    d1scf=as.character(unique(d1[,1]))
    
    newd<-NULL
    for (i in d1scf) {
      count=1
      sel=d1[d1[,1]==i,] # subset scaffold
      n<-NULL
      for (j in c(1:length(sel[,1]))) {
        if ((j/ngene)==round(j/ngene)) { count=count+1 }
        n=c(n,paste(sel[j,1],'.part',count,sep=''))
      }
      sel[,1]=n
      newd=rbind(newd,sel)
    }
    
    d1=newd
    
  }
  if (plot_param_env$par2) {
    print("Partitioning 2 ! ... ")
    ngene=PARTITION_ngene
    
    d2scf=as.character(unique(d2[,1]))
    
    newd<-NULL
    for (i in d2scf) {
      count=1
      sel=d2[d2[,1]==i,]
      n<-NULL
      for (j in c(1:length(sel[,1]))) {
        if ((j/ngene)==round(j/ngene)) { count=count+1 }
        n=c(n,paste(sel[j,1],'.part',count,sep=''))
      }
      sel[,1]=n
      newd=rbind(newd,sel)
    }
    
    d2=newd
  }
  
  # after partition, only consider scaffold.partN shared between both species
  allsame = unique(as.character(d1[d1[,4] %in% d2[,4],4]))
  if (plot_param_env$shared) {
    d1=d1[d1[,4] %in% allsame,]
    d2=d2[d2[,4] %in% allsame,]
  }
  
  
  # Build msynt matrix (to cluster the scaffold) ------------------------------------------------------------
  # convert genomic position to axis (staring from 1, each msynt/partition ++) 
  # more msynt/parition(which is also based on #msynt), wider celler
  # check_abx <- d1 %>% group_by(V1) %>% summarise(n=n())
  
  # mark down the corresponding scaffold, start,end coord of each "base"
  # when plot TE, group them using base, so we can have a plot using same coord system
  # if partition is on, V1 ~ <scaffold_name>.part<n>
  
  base=0
  abx=1
  n=as.character(d1[1,1])
  for (i in c(1:length(d1[,1]))) {
    base=base+1
    if (n==as.character(d1[i,1])) { } else { n=as.character(d1[i,1]); abx=c(abx,base) }
    d1[i,6] = d1[i,3] # back-up the coord
    d1[i,3]=base
  }
  
  
  aby=1
  base=0
  n=as.character(d2[1,1])
  for (i in c(1:length(d2[,1]))) {
    base=base+1
    if (n==as.character(d2[i,1])) { } else { n=as.character(d2[i,1]); aby=c(aby,base) }
    d2[i,6] = d2[i,3]
    d2[i,3]=base
  }
  
  fam = intersect(d1$V4,d2$v4)
  #fam=unique(c(as.character(d1[,4])[as.character(d1[,4]) %in% as.character(d2[,4]) ] ) )
  abx0=abx
  
  # (if using partition) update scaffold vector, use scaffold_partN instead of scaffold
  d1scf=as.character(unique(d1[,1]))
  d2scf=as.character(unique(d2[,1]))
  m <- matrix(0,nrow=length(d1scf),ncol=length(d2scf))
  m_abs <- matrix(0,nrow=length(d1scf),ncol=length(d2scf))
  
  # all to all comparison of #msynt in scaffold/scaffold.part 
  for (i in c(1:length(d1scf))) {
    for (j in c(1:length(d2scf))) {
      fam1=unique(as.character(d1[d1[,1]==d1scf[i],4]))
      fam2=unique(as.character(d2[d2[,1]==d2scf[j],4]))
      allfam=min(length(fam1),length(fam2)) # ? why take the min of family number?
      same=sum(fam1 %in% fam2)
      m[i,j]=same/allfam
      m_abs[i,j]=same
    }
  }
  
  # use heapmap for clustering--------------------------------------------------------------------
  hm=heatmap(m,distfun=function(x) { dist(x,method='euclidean') }, hclustfun=function(x) { hclust(x,method='ward.D2') })
  m_sort=m[hm$rowInd,hm$colInd]
  rownames(m_abs)=d1scf
  
  # Build msynt matrix (for scaled plot) --------------------------------------------------------------
  # use the rowInd and colInd from heatmap
  d0=d1
  d1<-NULL
  ii=0
  
  # reorder based on heatmap clustering, should be better to sort it using rowInd.... for species A
  for (i in d1scf[hm$rowInd]) {
    #ii=ii+1 # what does ii for??
    d1<-rbind(d1,d0[d0[,1]==i,])
  }
  
  base=0
  abx=1
  abx2=1
  n=as.character(d1[1,1])
  
  # mark down the number of msynt per scaffold again; V6 is the coord.... for species A
  for (i in c(1:length(d1[,1]))) {
    base=base+1
    scfname=as.character(d1[i,1])
    if (n==as.character(d1[i,1])) { } else { n=as.character(d1[i,1]); abx=c(abx,base) }
    d1[i,3]=base
  }
  
  
  # reorder based on heatmap clustering.... for species B
  d0=d2
  d2<-NULL
  ii=0
  for (i in d2scf[hm$colInd]) {
    ii=ii+1
    d2<-rbind(d2,d0[d0[,1]==i,])
  }
  
  base=0
  aby=1
  aby2=1
  nc=1
  n=as.character(d2[1,1])
  
  # mark down the number of msynt per scaffold again; V6 is the coord.... for species B
  for (i in c(1:length(d2[,1]))) {
    base=base+1
    scfname=as.character(d2[i,1])
    if (n==as.character(d2[i,1])) { } else { n=as.character(d2[i,1]); aby=c(aby,base) }
    
    d2[i,3]=base
  }
  
  # shared orthogroup between two species
  fam=unique(c(as.character(d1[,4])[as.character(d1[,4]) %in% as.character(d2[,4]) ] ) )
  
  # res is sorted by fam, but since base is calculated based on the cluster order, so it is fine
  res<-NULL
  #res_col<-NULL
  for (i in fam) {
    d1hit=d1[d1[,4] %in% i,]
    d2hit=d2[d2[,4] %in% i,]
    if ((length(d1hit[,1])>0)&(length(d2hit[,1])>0)) { # safe check
      for (x in c(1:length(d1hit[,3]))) { # each of the scaffold in species A with the orthogroup
        #d2c=0;
        for (y0 in c(1:length(d2hit[,3]))) { # each of the scaffold in species B with the orthogroup
          y=d2hit[y0,3]
          #d2c=d2c+1;
          # format of the df:
          # speciesA_base,speciesB_base,speciesA_scaffold,speciesB_scaffold,speciesA_baseName,speciesB_baseName
          res<-rbind(res,c(d1hit[x,3],y ,as.character(d1hit[x,1]),as.character(d2hit[y0,1]),as.character(d1hit[x,2]),as.character(d2hit[y0,2]))) 
        }
      }
    }
  }
  
  # store the res dt to reactiveValues, update d1 and d2 as well
  # reactivevalues (rv) as a environment in global scope
  res_dt <- data.table(res)
  colnames(res_dt) <- c("baseA","baseB","scaffoldA","scaffoldB",
                           "baseNameA","baseNameB")
  res_dt$baseA <- as.numeric(res_dt$baseA)
  res_dt$baseB <- as.numeric(res_dt$baseB)
  #debug 
  print(res_dt)
  reactive_env$res_dt <- res_dt
  d_env$d1 <- d1
  d_env$d2 <- d2
  d_env$abx <- abx
  d_env$aby <- aby
  
  # plot
  par(mar=c(6,6,6,6))
  plot(as.numeric(res[,1]),as.numeric(res[,2]),pch=16,cex=0.1,xlab=speciesA,ylab=speciesB)
  
  if (plot_param_env$showlines) {
    abline(v=abx,lty=3,lwd=0.5)
    abline(h=aby,lty=3,lwd=0.5)
  }
  
  title(sub = paste0("lim=",plot_param_env$lim,";nmin=",plot_param_env$nmin))
  
  axis(3,at=abx,labels=as.character(d1[d1[,3] %in% abx,1]),las=2)
  axis(4,at=aby,labels=as.character(d2[d2[,3] %in% aby,1]),las=2)
  
  p <- recordPlot()
  return(p)
}


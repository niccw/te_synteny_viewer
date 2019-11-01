# libraries
library(gtools)
library(magrittr)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
source("te_matrix.R")


# File and setting --------------------------------------------------------
# "Desktop/vlogin3/hydractinia_2019/sep2019/synteny_break/"

# files
#FOLDER="./"
FOLDER="./data"

FILE0='hech.msynt'
FILE='hsym.msynt'

SCAFFOLD_LEN0 = "Hech_primary_v1.0.fa.genome"
SCAFFOLD_LEN = "Hsym_primary_v1.0.fa.genome"

# TE matrix
# reformated gff3, $3=TE general class
TE_GFF0 = "Hech_primary_v1.0.fa.out.reformat.gff3"
TE_GFF = "Hsym_primary_v1.0.fa.out_reformat.gff3"

# gene gff3
GENE_GFF0 = "Hech_primary_v1.0.gff3"
GENE_GFF = "Hsym_primary_v1.0.gff3"

TE_CLASS = c("DNA")
UPSTREAM_FLANK = TRUE # if FALSE, plot downstream
FLANKING_SIZE = 2500
PLOT_TE_SP0 = TRUE
PLOT_TE_SP1 = TRUE
PLOT_TE_BOTH = TRUE
SAVE_PLOTLY = TRUE
  

#FILTERING
#use families only with lim paralogs or less
lim=4
#include scaffolds with minimum nmin genes
nmin=10  #does not seem to work properly!
#show only shared families
shared=1

#PROCESSING
#partition genome1 and/or genome2?
par1=0
par2=0
#partition scaffolds/chromosomes into chunks with PARTITION_ngene genes
PARTITION_ngene=50

#PLOTTING
#order columns in the plot?
colorder=1
#show separators for chromosomes or scaffolds on the plot?
showlines=1
#prefix for output
PREFIX=''

#SCRIPT

print(paste(FILE0,'vs',FILE))

d1=read.table(paste(FOLDER,FILE0,sep=''),header=F,sep="\t")
d2=read.table(paste(FOLDER,FILE,sep=''),header=F,sep="\t")

# filter scaffold by # of msynt ---------------------
allsame = unique(as.character(d1[d1[,4] %in% d2[,4],4])) # intersect(d1$V4,d2$V4)
if (shared==1) {
d1=d1[d1[,4] %in% allsame,]
d2=d2[d2[,4] %in% allsame,]
}

lim1=lim # filter msynt by # paralogs

d1=d1[d1[,5]<=lim1,]
d2=d2[d2[,5]<=lim,]

# scaffold ID
d1scf=as.character(unique(d1[,1])) # unique(d1$V1)
d2scf=as.character(unique(d2[,1])) # unique(d2$V2)

# number of uniqe msynt on each scaffold ?
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
d1scf=d1scf[d1len>=nmin]
d2scf=d2scf[d2len>=nmin]
d1len=d1len[d1len>=nmin]
d2len=d2len[d2len>=nmin]

d1=d1[d1[,1] %in% d1scf,]
d2=d2[d2[,1] %in% d2scf,]

#FIXME: redundant 
d1scf=as.character(unique(d1[,1]))
d2scf=as.character(unique(d2[,1]))

labelx=d1scf
labely=d2scf


# Partioning --------------------------------------------------------------
# The partition can be done by:
# as.integer((seq_along(1:nrow(sel)) - 1) / ngene)

if (par1==1) {

 print("Partitioning 1 ! ... ")
 ngene=PARTITION_ngene

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
if (par2==1) {

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
if (shared==1) {
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

# fam = intersect(d1$V4,d2$v4)
fam=unique(c(as.character(d1[,4])[as.character(d1[,4]) %in% as.character(d2[,4]) ] ) )

abx0=abx


# (if using partition) update scaffold vector, use scaffold_partN instead of scaffold
d1scf=as.character(unique(d1[,1]))
d2scf=as.character(unique(d2[,1]))
m<-matrix(0,nrow=length(d1scf),ncol=length(d2scf))
m_abs<-matrix(0,nrow=length(d1scf),ncol=length(d2scf))


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


# Plot heatmap and use heapmap for clustering--------------------------------------------------------------------

# (optinal plot)
#pdf(file=paste(FOLDER,PREFIX,FILE0,'-',FILE,'.heatmap.pdf',sep=''))
if (colorder==1) {
hm=heatmap(m,distfun=function(x) { dist(x,method='euclidean') }, hclustfun=function(x) { hclust(x,method='ward.D2') })
} else {
 hm=heatmap(m,Rowv=NA,Colv=NA,distfun=function(x) { dist(x,method='euclidean') }, hclustfun=function(x) { hclust(x,method='ward.D2') })
}
m_sort=m[hm$rowInd,hm$colInd]
#dev.off()

rownames(m_abs)=d1scf
colnames(m_abs)=d2scf
#rite.table(m_abs,file=paste(FOLDER,PREFIX,FILE0,'-',FILE,'.m_abs',sep=''))



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

# build the file name
if ((par2==1)||(par1==1)) {
 fadd2='part'
} else { fadd2='non-part' }
if (colorder==1) {
fadd2=paste(fadd2,'.yordered',sep='') } else {
}

#write.table(res,paste(FOLDER,FILE0,'-',FILE,'.',fadd2,'.res',sep=''),sep='\t',quote=F,row.names=F,col.names=F)

#pdf(file=paste(FOLDER,PREFIX,FILE0,'-',FILE,'.',fadd2,'.pdf',sep=''))

# msynt dot plot (on scale)
par(mar=c(6,6,6,6))

plot(as.numeric(res[,1]),as.numeric(res[,2]),pch=16,cex=0.1,xlab=FILE0,ylab=FILE)

if (showlines==1) {
abline(v=abx,lty=3,lwd=0.5)
abline(h=aby,lty=3,lwd=0.5)
}

axis(3,at=abx,labels=as.character(d1[d1[,3] %in% abx,1]),las=2)
axis(4,at=aby,labels=as.character(d2[d2[,3] %in% aby,1]),las=2)

#dev.off() 

# msynt dot plot (on scale), ggplot ver.
# library(data.table)
# library(ggplot2)
# library(grid)
# 
# # this is important for ggplot!! 
# res_dt$baseA <- as.numeric(res_dt$baseA)
# res_dt$baseB <- as.numeric(res_dt$baseB)
# 
# xlabel <- unique(d1$V1)
# ylabel <- unique(d2$V1)

# ggplot
# g <- ggplot(res_dt,aes(x=baseA,y=baseB))
# p <- g + geom_point(shape = ".") + theme(panel.background=element_blank(),
#                                     #axis.text = element_blank(),
#                                     #axis.ticks = element_blank(),
#                                     plot.margin = unit(c(1,1,3,3),"line")) +
#   geom_vline(xintercept = as.numeric(abx),color="black",linetype="dotted",alpha=0.5) +
#   geom_hline(yintercept = as.numeric(aby),color="black",linetype="dotted",alpha=0.5) 
# 
# 
# 
# # add axis (scaffold name), but if scaffold is short probably will all be overlapped
# for (i in 1:length(abx)) {
#   p <- p + annotation_custom(grob=textGrob(xlabel[i], rot=60, gp = gpar(cex = 0.5)),
#                         xmin=abx[i],xmax=abx[i],
#                         ymin=15,ymax=15)
# }
# 
# gt <- ggplot_gtable(ggplot_build(p))
# gt$layout$clip[gt$layout$name == "panel"] <- "off"
# grid.draw(gt)


# TE cnt ---------------------------------------------------------------
## Use .genome for the scaffold len, make sure end pos not overflow
scaffold_len1 <- read_scaffold_len(SCAFFOLD_LEN0)
scaffold_len2 <- read_scaffold_len(SCAFFOLD_LEN)

# create og to transcript map (only used when considering downstream region as well)
# d1,d2 scaffold already sorted
d1 %<>% separate(V2,c("SP","ID"),sep="\\|")
d2 %<>% separate(V2,c("SP","ID"),sep="\\|")

gene_gff1 <- rtracklayer::import.gff3(GENE_GFF0) # seperate sp name from ID (which match to gene gff)
gene_gff2 <- rtracklayer::import.gff3(GENE_GFF)

# extract the end position of gene
# TODO: if use partition, have to extract the end positionof last OG
gene_gff1 <- gene_gff1[gene_gff1$type == "gene"]
gene_gff2 <- gene_gff2[gene_gff2$type == "gene"]
gene_end1 <- GenomicRanges::end(gene_gff1) # IRange start always smaller than end, don't need to worry the order
gene_end2 <- GenomicRanges::end(gene_gff2)
gene_end1 <- data.table(ID=gene_gff1$Name, gene_end = gene_end1)
gene_end2 <- data.table(ID=gene_gff2$Name, gene_end = gene_end2)
d1 %<>% dplyr::full_join(gene_end1)
d2 %<>% dplyr::full_join(gene_end2)
# filter gene not in msynt
d1 <- d1[!is.na(d1$V1),]
d2 <- d2[!is.na(d2$V1),]
# <...trying

# A map of Base number to a range of genomic Coord
bc_map1 <- base_coord_map(d1,scaffold_len1,upstream = UPSTREAM_FLANK,frank_size = FLANKING_SIZE)
bc_map2 <- base_coord_map(d2,scaffold_len2,upstream = UPSTREAM_FLANK,frank_size = FLANKING_SIZE)

# read repeat gff3
te_gff1 <- rtracklayer::import.gff3(TE_GFF0)
te_gff2 <- rtracklayer::import.gff3(TE_GFF)

# TE count matix, row=base+-window size
# take some time, for 10000 base, ~3 mins
te_base_cnt1 <- te_base_cnt(bc_map1,te_gff1)
te_base_cnt2 <- te_base_cnt(bc_map2,te_gff2)

# te_plot_matrix, filtering TE class of interest
te_df_ls <- te_plot_matrix(tecnt1 = te_base_cnt1, tecnt2 = te_base_cnt2,teClasses = TE_CLASS)
te_df1 <- te_df_ls[1][[1]]
te_df2 <- te_df_ls[2][[1]]
rm(te_df_ls)

res_dt <- data.table(res) # use dt not df, memory efficient!
rm(res)
colnames(res_dt) <- c("baseA","baseB","scaffoldA","scaffoldB",
                      "baseNameA","baseNameB")
res_dt$baseA <- as.numeric(res_dt$baseA)
res_dt$baseB <- as.numeric(res_dt$baseB)

res_dt <- full_join(res_dt,te_df1)
res_dt <- full_join(res_dt,te_df2)
rm(te_df1,te_df2)
res_dt$baseA_n[res_dt$baseA_n==0] <- NA
res_dt$baseB_n[res_dt$baseB_n==0] <- NA
res_dt$baseA_bp[res_dt$baseA_bp==0] <- NA
res_dt$baseB_bp[res_dt$baseB_bp==0] <- NA

# msynt + te plot ---------------------------------------------------------
PLOT_TYPE="n" # c("n","bp")

baseA_plottype <- paste0("baseA_",PLOT_TYPE)
baseB_plottype <- paste0("baseB_",PLOT_TYPE)


if(PLOT_TE_SP0){
  ga <- ggplot(res_dt,aes(x=baseA,
                          text = paste("scaffoldA:",scaffoldA,
                                       "scaffoldB:",scaffoldB,
                                       "baseNameA:",baseNameA,
                                       "baseNameB:",baseNameB
                          )))
  pa <- ga + geom_vline(aes_string(xintercept="baseA",color = baseA_plottype)) + 
    geom_point(aes(y=baseB),shape = ".",size=0.3) + 
    scale_color_gradient(low = "#ffffff", high = "#ff0000",na.value="white") +
    theme(panel.background=element_blank(),
          plot.margin = unit(c(1,1,3,3),"line")) +
    geom_vline(xintercept = as.numeric(abx),color="black",linetype="dotted",alpha=0.5, size=0.3) +
    geom_hline(yintercept = as.numeric(aby),color="black",linetype="dotted",alpha=0.5, size=0.3) +
    labs(colour = paste("TEcnt_ws:", WINDOW_SIZE),
         x = FILE0, y = FILE) 
  
  ggsave(gg_test,file=paste(FOLDER,PREFIX,FILE0,'-',FILE,'.',fadd2,'.ggplot','te_spA','.pdf',sep=''),
         width = 350, height = 350, units = "mm")
  if(SAVE_PLOTLY){
    gpa<- ggplotly(pa,tooltip = c("text", "x", "y"))
    htmlwidgets::saveWidget(gpa, file=paste(FOLDER,PREFIX,FILE0,'-',FILE,'.',fadd2,'.ggplot','te_spA','.html',sep=''))
  }
}

if(PLOT_TE_SP1){
  gb <- ggplot(res_dt,aes(x=baseB,
                          text = paste("scaffoldA:",scaffoldA,
                                       "scaffoldB:",scaffoldB,
                                       "baseNameA:",baseNameA,
                                       "baseNameB:",baseNameB
                          )))
  pb <- gb + geom_hline(aes(yintercept=baseB,color=baseB_tecnt)) + 
    geom_point(aes(y=baseA),shape = ".",size=0.3) + 
    scale_color_gradient(low = "#ffffff", high = "#2053fa",na.value="white") +
    theme(panel.background=element_blank(),
          plot.margin = unit(c(1,1,3,3),"line")) +
    geom_vline(xintercept = as.numeric(abx),color="black",linetype="dotted",alpha=0.5, size=0.3) +
    geom_hline(yintercept = as.numeric(aby),color="black",linetype="dotted",alpha=0.5, size=0.3) +
    labs(colour = paste("TEcnt_ws:", WINDOW_SIZE),
         x = FILE0, y = FILE)
  
  ggsave(gg_test,file=paste(FOLDER,PREFIX,FILE0,'-',FILE,'.',fadd2,'.ggplot','te_spB','.pdf',sep=''),
         width = 350, height = 350, units = "mm")
  if(SAVE_PLOTLY){
    gpa<- ggplotly(gg_test,tooltip = c("text", "x", "y"))
    htmlwidgets::saveWidget(ggp, file=paste(FOLDER,PREFIX,FILE0,'-',FILE,'.',fadd2,'.ggplot','te_spB','.html',sep=''))
  }
}


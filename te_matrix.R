# count TE in range
#library(hash)
library(dplyr)
library(data.table)
library(magrittr)
library(rtracklayer)
library(GenomicRanges)

# read genome.tsv (generated from faidx), return a named vector v[scaffold] = len
read_scaffold_len <- function(genome_tsv_path){
  tb<- read.table(genome_tsv_path,header = F)
  v <- tb[,2]
  names(v) <- tb[,1]
  return(v)
}

# map the base(coord used in the plot, which corresponding to the #msynt) to genomic position (window size)
base_coord_map <- function(df,scaffold_len,upstream=TRUE,frank_size=2500){
  # df: msynt df
  # scaffold_len: named vector of scaffold len
  # window size: window sieze
  start_v <- c()
  end_v <- c()
  # df col: scaffold/scaffold_part, msyntID,base,gene_orthogroupID, #paralog in this species, coord
  for (i in 1:nrow(df)) {
    # start coord
    if(upstream){
      if((df$V6[i] - frank_size) > 1){
        start_v <- c(start_v, df$V6[i] - frank_size)
      }else{
        start_v <- c(start_v, 1)
      }
    }
    else{
      # end coord
      if ((df$gene_end[i] + frank_size) < scaffold_len[df$V1[i]]){ # last record on this scaffold
        end_v <- c(end_v, df$gene_end[i] + frank_size)
      }else{
        end_v <- c(end_v, scaffold_len[df$V1[i]])
      }
    }
  }
  if(upstream){
    return(data.frame(scaffold= df$V1,base_name = df$V3, base_coord = df$V6, window_start = start_v, window_end = df$V6))
  }else{
    return(data.frame(scaffold= df$V1,base_name = df$V3, base_coord = df$V6, window_start = df$gene_end, window_end = end_v))
  }
}

# use sorted gff; count the number of TE per each region
# .....each cell in the heatmap is a scaffold....
# TE count table for all base matix
# 10000 base takes ~3 mins 
te_base_cnt <- function(m,gr,p){
  # p: shiny::Progress
  # create two empty matrix to hold the result
  cnt_m <- matrix(0,nrow=nrow(m),ncol=length(levels(gr$type)))
  bp_m <- matrix(0,nrow=nrow(m),ncol=length(levels(gr$type)))
  
  # use named vector as dictionary for te_class index in the output matrix
  class_nv <- c(1:length(levels(gr$type)))
  names(class_nv) <- levels(gr$type)
  # subset class
  n <- nrow(m) %/% 1000
  for(i in 1:nrow(m)){
    if (i%%1000 == 0){
      p$inc(1/n, detail = paste("Doing part", (i/1000)))
    }
    # scaffold,window_start,window_end
    sub_gr <- gr[seqnames(gr) == as.character(m$scaffold[i]) &
                        start(gr) >= as.numeric(m$window_start[i])&
                        end(gr) <= as.numeric(m$window_end[i])]
    if(length(sub_gr) == 0){
      next
    }
    # te_class, n, bp_sum
    gr_df <- data.table(te_class = sub_gr$type, te_len = width(sub_gr),te_cnt=1) %>% group_by(te_class) %>% summarise(n=n(),bp_sum=sum(te_len))
    
    # FIXME: dt (nrow=#base,col = TE class, value=n/bp_cnt)
    # write cnt_df
    for(j in 1:nrow(gr_df)){
      # col index = class_nv[["te_class"]], row index is i (same as base number)
      # write cnt_m 
      cnt_m[i,class_nv[[as.character(gr_df[[j,1]])]]] <- gr_df[[j,2]]
      # write bp_m
      bp_m[i,class_nv[[as.character(gr_df[[j,1]])]]] <- gr_df[[j,3]]
    }
  }
  colnames(cnt_m) <- levels(gr$type)
  colnames(bp_m) <- levels(gr$type)
  return(list(cnt_m,bp_m))
}

# play with hash (not used in this script)
add_df_to_hash <- function(df,h){
  # df col1:key col2:value (int)
  for(i in 1:nrow(df)){
    # add new key
    te_class <- as.character(df[[i,1]])
    if(is.null(h[[te_class]])){
      h[[te_class]] <- df[[i,2]]
    }else{
      h[[te_class]] <- h[[te_class]] + df[[i,2]]
    }
  }
  return(h)
}

# create TE plot matrix, dimention=baseA x baseB
te_plot_matrix <- function(tecnt1,tecnt2,teClasses){
  
  # subset cnt matix (allow multiple teClass)
  if (length(teClasses) > 1){
    search_phrase = paste(teClasses,collapse="|")
    n1 <- tecnt1[[1]][,grep(search_phrase,colnames(tecnt1[[1]]))] %>% rowSums()
    n2 <- tecnt2[[1]][,grep(search_phrase,colnames(tecnt2[[1]]))] %>% rowSums()
    bp1 <- tecnt1[[2]][,grep(search_phrase,colnames(tecnt1[[2]]))] %>% rowSums()
    bp2 <- tecnt2[[2]][,grep(search_phrase,colnames(tecnt2[[2]]))] %>% rowSums()
  }else{
    search_phrase = as.character(teClasses)
    n1 <- tecnt1[[1]][,grep(search_phrase,colnames(tecnt1[[1]]))] %>% base::as.vector()
    n2 <- tecnt2[[1]][,grep(search_phrase,colnames(tecnt2[[1]]))] %>% base::as.vector()
    bp1 <- tecnt1[[2]][,grep(search_phrase,colnames(tecnt1[[2]]))] %>% base::as.vector()
    bp2 <- tecnt2[[2]][,grep(search_phrase,colnames(tecnt2[[2]]))] %>% base::as.vector()
  }

  # it doesn't make sense to compute all cell, since the unit is msynt for each species, not comparable at all!
  # should be one matrix for one species, showing the number of species class TE regarding 1 species
  # if te class in 1 species only
  if (length(n1) == 0){
    n1 <- NA
    bp1 <- NA
  }
  if(length(n2) == 0){
    n2 <- NA
    bp2 <- NA
  }
  d1 <- data.table(baseA=1:nrow(tecnt1[[1]]), baseA_n=n1, baseA_bp=bp1)
  d2 <- data.table(baseB=1:nrow(tecnt2[[1]]), baseB_n=n2, baseB_bp=bp2)
  
  return(list(d1,d2))
}


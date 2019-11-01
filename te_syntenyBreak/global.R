# 
source("../te_matrix.R")
source("../main.R")
library(magrittr)

# data
data_path="/Users/niccw/Desktop/vlogin3/hydractinia_2019/sep2019/synteny_break/data"
setwd(data_path)

# find the species name at ~/msynt/
species_name <- gsub("\\.msynt","",list.files("./msynt/"))

# find correspoding gene_gff3, te_gff3 and .genome
speciesA_files <- NULL
speciesB_files <- NULL
input_files_by_name <- function(n){
  gene_gff3 <- tolower(list.files("./gene_gff3/")) %>% grep(n, ., value=TRUE) %>% paste0("./gene_gff3/", .)
  te_gff3 <- tolower(list.files("./te_gff3/")) %>% grep(n, ., value=TRUE) %>%  paste0("./te_gff3/", .)
  genome_file <- tolower(list.files("./genome_file/")) %>% grep(n, ., value=TRUE) %>% paste0("./genome_file/", .)
  # If file not exist
  if (gene_gff3 == "./gene_gff3/"){
    gene_gff3 <- NULL
  }
  if(te_gff3 == "./te_gff3/"){
    te_gff3 <- NULL
  }
  if(genome_file == "./genome_file/"){
    genome_file <- NULL
  }
  v <- c(msynt = paste0("./msynt/",n,".msynt"),gene_gff3 = gene_gff3, te_gff3 = te_gff3, genome_file = genome_file,sp_name=n)
  return(v)
}

# load files
data_env <- new.env(parent = emptyenv())
data_env$sp1 <- NULL
data_env$sp2 <- NULL
data_env$d1 <- NULL
data_env$d2 <- NULL
data_env$d1_raw <- NULL
data_env$d2_raw <- NULL
data_env$scaffold_len1 <- NULL
data_env$scaffold_len2 <- NULL
data_env$gene_gff1 <- NULL
data_env$gene_gff2 <- NULL
data_env$te_gff1 <- NULL
data_env$te_gff2 <- NULL
data_env$plot_te <- TRUE
data_env$abx <- NULL
data_env$aby <- NULL
data_env$basic_plot_cp <- NULL

load_files <- function(filelsA,filelsB){
  # mark down species name
  data_env$sp1 <- filelsA[["sp_name"]]
  data_env$sp2 <- filelsB[["sp_name"]]
  # read msynt
  data_env$d1 <- read.table(filelsA[["msynt"]],header=F,sep="\t")
  data_env$d2 <- read.table(filelsB[["msynt"]],header=F,sep="\t")
  # files for te plot
  # read gene gff3
  tryCatch({
    data_env$gene_gff1 <- rtracklayer::import.gff3(filelsA[["gene_gff3"]])
    data_env$gene_gff2 <- rtracklayer::import.gff3(filelsB[["gene_gff3"]])
    # read te gff3
    data_env$te_gff1 <- rtracklayer::import.gff3(filelsA[["te_gff3"]])
    data_env$te_gff2 <- rtracklayer::import.gff3(filelsB[["te_gff3"]])
    # read .genome
    data_env$scaffold_len1 <- read_scaffold_len(filelsA[["genome_file"]])
    data_env$scaffold_len2 <- read_scaffold_len(filelsB[["genome_file"]])
  },
  error=function(e){
    warning_id <<- showNotification("Warning: fail to read gene/te files. Only basic msynt available.", duration=NULL,type = "warning")
    data_env$plot_te <- FALSE
    })
}


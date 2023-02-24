#' ---
#' title: "Initialize"
#' author: "malare"
#' date: "07/03/2022"
#' ---

#' Load packages (make sure you installed all packages!)
library(xlsx)
library(downloader)
library(httr)
library(dplyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(nnet)
library(splines)
library(effects)
library(deSolve)
library(bbmle)
library(SimDesign)
library(MCMCglmm)
library(dplyr)
library(MASS)
library(RDS)
library(doParallel)
library(foreach)


# Set seed
set.seed(60321)


#' Paths
path_function = "R/"
path_saveoutput = "output/"


#' Add custom functions to env
files_functions = list.files(file.path(path_function),full.names = TRUE, include.dirs = TRUE)
files_functions = files_functions[!grepl(pattern = "000|backup|old",files_functions)]
sapply(files_functions[1:5], source)
rm(files_functions)

#' Graphics
theme_set(theme_bw())
col_9 <- (brewer.pal(9,"Set1"))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
yscaling <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("[+]", "", l)
  parse(text=l)
}


yscaling <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "e", l)
  l <- gsub("\\+", "", l)
  l <- gsub("e", "10^", l)
  parse(text=l)
}
g_legend <- function(a.gplot,num){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]] 
  if(num==1){
    legend$grobs[[1]]$grobs[[1]] <-  editGrob(legend$grobs[[1]]$grobs[[1]], gp=gpar(fill ="transparent",col="transparent"))
  }
  if(num==2){
    legend$grobs[[1]]$grobs[[1]] <-  editGrob(legend$grobs[[1]]$grobs[[1]], gp=gpar(fill ="transparent",shape="transparent"))
    legend$grobs[[2]]$grobs[[1]] <-  editGrob(legend$grobs[[2]]$grobs[[1]], gp=gpar(fill ="transparent",shape="transparent"))
  }
  legend
} 


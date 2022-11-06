#!/usr/bin/env Rscript
#Rscript to calculate normalized RF Distance between two trees using TreeDist
list.of.packages <- c("TreeDist", "ape","codetools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos='http://cran.us.r-project.org')
args = commandArgs (trailingOnly=TRUE)
library("TreeDist")
library("ape")
t1 <- ape::read.tree(args[1])
t2 <- ape::read.tree(args[2])
td = TreeDistance(t1,t2)
cat(td)


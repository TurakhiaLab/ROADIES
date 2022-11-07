list.of.packages <- c("TreeDist", "ape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos='http://cran.us.r-project.org')
args = commandArgs (trailingOnly=TRUE)

library("ape")
library("TreeDist")

treePath <- args[1]
filePattern <- args[2]

files <- list.files(path = treePath, pattern=filePattern)
files

trees <- c()
for(file in files){
    filePath <- paste(treePath, '/', file, sep='')
    tree <- ape::read.tree(filePath)
    trees <- append(trees, tree)
}

rfMat <- RobinsonFoulds(trees)
rfMat

cat("mean =", as.character(mean(rfMat)), "\n")
cat("min =", as.character(min(rfMat)), "\n")
cat("max =", as.character(max(rfMat)), "\n")
cat("\n")

# distMat


# refDist <- list()
# for(tree in list(tree_L500_1, tree_L500_2, tree_L500_3)){
#     append(refDist, RobinsonFoulds(refTree, tree))
# }
# refDist
# RobinsonFoulds(list(tree_L200_1, tree_L200_2, tree_L200_3))
# RobinsonFoulds(list(tree_L1000_1, tree_L1000_2, tree_L1000_3))
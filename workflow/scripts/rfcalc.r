library("ape")
library("TreeDist")

param <- "K"
treePath <- "./treedist/trees/"

lenDirs <- list.files(path = treePath, pattern=param)
lens <- c()
for (dir in lenDirs){
    lens <- append(lens, as.integer(substring(dir, 2, 1000)))
}

distMat <- list()
for(len in sort(lens)){
    files <- list.files(path = paste(treePath, param, as.character(len), sep = ""))
    trees <- c()
    for(file in files){
        filePath <- paste(treePath, param, as.character(len), "/", as.character(file), sep = "")
        tree <- ape::read.tree(filePath)
        trees <- append(trees, tree)
    }

    rfMat <- RobinsonFoulds(trees)
    distMat <- append(distMat, rfMat)

    cat("=========", param, "=", as.character(len), "=========\n")
    print(rfMat)
    cat("\n")
    cat("mean =", as.character(mean(rfMat)), "\n")
    cat("min =", as.character(min(rfMat)), "\n")
    cat("max =", as.character(max(rfMat)), "\n")
    cat("\n")
}

# distMat


# refDist <- list()
# for(tree in list(tree_L500_1, tree_L500_2, tree_L500_3)){
#     append(refDist, RobinsonFoulds(refTree, tree))
# }
# refDist
# RobinsonFoulds(list(tree_L200_1, tree_L200_2, tree_L200_3))
# RobinsonFoulds(list(tree_L1000_1, tree_L1000_2, tree_L1000_3))
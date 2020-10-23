require(Seurat)
require(qlcMatrix)

# Function to calculated truncated PC scores, where the  
RunTruncatedPCA <- function(object, n.genes.pc=40, genescale.method="scale", ...){
  
  object <- RunPCA(object, ...) # since Seurat seems to bug out if PCA hasn't initially been run
  old.loadings <- object@reductions$pca@feature.loadings
#  old.loadings <- GetDimReduction(object, reduction.type = "pca", slot="gene.loadings")
  new.pc.scores <- reCalculatePCScores(object, pc.loadings = old.loadings, n.genes.pc = n.genes.pc, use.binary.weights = F, genescale.method = genescale.method)
#  object <- SetDimReduction(object, reduction.type = "pca",slot = "cell.embeddings",new.data = new.pc.scores)
  object@reductions$pca@cell.embeddings<-new.pc.scores
  object
}

reCalculatePCScores <- function(object, pc.loadings, n.genes.pc=30, use.binary.weights=F, genescale.method=c("div.by.max","mean.center","scale","seurat.scale")){
 
#  data <- object@data[rownames(pc.loadings), ] # retain only genes present in the loadings
  data <- object@assays$RNA@data[rownames(pc.loadings), ] # retain only genes present in the loadings
  
    # gene-wise scaling by one of several methods
  if(genescale.method=="div.by.max") data <- data / qlcMatrix::rowMax(data)
  else if(genescale.method=="mean.center") data <- t(scale(t(data), center = T, scale = F))
  else if(genescale.method=="scale") data <- t(scale(t(data)))
#  else data <- object@scale.data[rownames(pc.loadings), ] # use Seurat scaled data
  else data <- object@assays$RNA@scale.data[rownames(pc.loadings), ] # use Seurat scaled data
    
  # if n.genes.pc == 0, then no cutoff is applied and loading is simply recalculated as input loadings * scaled data
  if(n.genes.pc!=0){
    # apply loading cutoffs so that only genes in top n positive or top n negative loadings for each PC get non-zero loadings
    pc.loadings.reweighted <- apply(pc.loadings, 2, cutoff.by.order, n.cutoff=n.genes.pc)
    
    if(use.binary.weights){
      pc.loadings.reweighted[pc.loadings.reweighted > 0] <- 1
    }
    
    pc.scores <- as.matrix( t(data) %*% pc.loadings.reweighted )
  }
  else{
     pc.scores <- as.matrix( t(data) %*% pc.loadings )
  }
  
  
  # recalculate PC scores (cell-wise scores) based on new gene weights
  return(pc.scores)
}

cutoff.by.order <- function(x, n.cutoff=30){
  x.rank <- rank(x)
  len.x <- length(x)
  x[which( (x.rank > n.cutoff) & (x.rank < (len.x-n.cutoff + 1)) )] <- 0
  
  return(x)
}

save.seurat.calc <- function(object, dims.use, file.ext=""){
  if(!dir.exists("seurat_results")) system("mkdir seurat_results")
  
  r_pca <- GetCellEmbeddings(object, reduction.type = 'pca')[, dims.use]
  r_tsne <- GetCellEmbeddings(object, reduction.type = 'tsne')
  r_clusters <- as.character(object@ident)
  cells_use <- as.character(object@cell.names)
  
  write.table(cells_use, file=paste0("seurat_results/cells_use",file.ext, ".txt"), row.names = F, col.names = F, quote = F)
  write.table(r_tsne, file=paste0("seurat_results/tsne", file.ext, ".csv"), row.names = T, col.names = T, quote = F, sep=",")
  write.table(r_pca, file=paste0("seurat_results/pca",file.ext, ".csv"), row.names = T, col.names = T, quote = F, sep=",")
  write.table(r_clusters, file=paste0("seurat_results/clusters", file.ext, ".txt"), row.names = F, col.names = F, quote = F)
}

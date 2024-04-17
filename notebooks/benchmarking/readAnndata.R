require(Seurat)
require(SeuratDisk)

# mb231 = read.table('../../data/toSeurat/mdamb231Jostner.csv', 
# 		   header = TRUE, 
# 		   sep = ',',
# 		   row.names = 1)
# seu <- CreateSeuratObject(mb231)
# seu = NormalizeData(object = seu)
# seu = FindVariableFeatures(object = seu)
# seu = ScaleData(object = seu)
# seu = RunPCA(object = seu)
# seu = FindNeighbors(object = seu, dims = 1:30)
# seu = FindClusters(object = seu)
# seu = RunUMAP(object = seu, dims = 1:30)
# 
# saveRDS(seu, '../../data/toSeurat/231jostnerProcessed.rds')

mb436 = read.table('../../data/toSeurat/mdamb436.csv', 
                   header = TRUE, 
                   sep = ',',
                   row.names = 1)
seu <- CreateSeuratObject(mb436)
seu = NormalizeData(object = seu)
seu = FindVariableFeatures(object = seu)
seu = ScaleData(object = seu)
seu = RunPCA(object = seu)
seu = FindNeighbors(object = seu, dims = 1:30)
seu = FindClusters(object = seu)
seu = RunUMAP(object = seu, dims = 1:30)

saveRDS(seu, '../../data/toSeurat/436jostnerProcessed.rds')
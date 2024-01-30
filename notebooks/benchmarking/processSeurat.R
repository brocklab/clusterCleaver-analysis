require(Seurat)
require(SeuratDisk)


mb231 = LoadH5Seurat('../../data/toSeurat/mdamb231Jostner.h5Seurat')

mb231 = NormalizeData(object = mb231)
mb231 = FindVariableFeatures(object = mb231)
mb231 = ScaleData(object = mb231)
mb231 = RunPCA(object = mb231)
mb231 = FindNeighbors(object = mb231, dims = 1:30)
mb231 = FindClusters(object = mb231)
mb231 = RunUMAP(object = mb231, dims = 1:30)

SaveH5Seurat(mb231, '../../data/toSeurat/mdamb231JostnerProcessed.h5Seurat', overwrite = TRUE)

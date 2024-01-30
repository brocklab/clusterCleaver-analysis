require(Seurat)
require(SeuratDisk)

mb231 = read.table('../../data/toSeurat/mdamb231Jostner.csv', 
		   header = TRUE, 
		   sep = ',',
		   row.names = 1)
seu <- CreateSeuratObject(mb231)
seu = FindVariableFeatures(object = seu)

SaveH5Seurat(seu, '../../data/toSeurat/mdamb231Jostner.h5Seurat', overwrite = TRUE)

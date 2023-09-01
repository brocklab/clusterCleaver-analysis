find_de_MAST_RE <- function(adata_){

    # create a MAST object
    sca <- SceToSingleCellAssay(adata_, class = "SingleCellAssay")
    print("Dimensions before subsetting:")
    print(dim(sca))
    print("")
    # keep genes that are expressed in more than 10% of all cells

    sca <- sca[freq(sca)>0.1,]
    print("Dimensions after subsetting:")
    print(dim(sca))
    print("")
    # add a column to the data which contains scaled number of genes that are expressed in each cell
    cdr2 <- colSums(assay(sca)>0)
    colData(sca)$ngeneson <- scale(cdr2)

    # store the columns that we are interested in as factors
    leiden <- factor(colData(sca)$leiden)

    # set the reference level
    leiden <- relevel(leiden,"0")

    colData(sca)$leiden <- leiden

    # define and fit the model
    zlmCond <- zlm(formula = ~leiden + ngeneson,
                   sca=sca,
                  ) # to speed up calculations

    summary(zlmCond)

    # perform likelihood-ratio test for the condition that we are interested in
    summaryCond <- summary(zlmCond, doLRT='leiden1')

    # get the table with log-fold changes and p-values
    summaryDt <- summaryCond$datatable


    result <- merge(summaryDt[contrast=='leiden1' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                    summaryDt[contrast=='leiden1' & component=='logFC', .(primerid, coef)],
                    by='primerid') # logFC coefficients

    # MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
    result[,coef:=result[,coef]/log(2)]

    # do multiple testing correction
    result[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]

    # result = result[result$FDR<0.01,, drop=F]
    result <- stats::na.omit(as.data.frame(result))

    return(result)
}
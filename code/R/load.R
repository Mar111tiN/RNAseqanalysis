library("pasilla")


load.pasilla <- function() {
    pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package = "pasilla", mustWork = T
                      )
    pasAnno <- system.file("extdata",
                        "pasilla_sample_annotation.csv",
                        package = "pasilla", mustWork = T
                        )

    # load the data
    cts <- as.matrix(read.csv(pasCts, sep="\t", row.names="gene_id"))
    # load the coldata
    # wrangle
    cols <- c("condition", "type")
    coldata <- read.csv(pasAnno, row.names = 1)[,cols]
    for (col in cols) {
        coldata[[col]] <- factor(coldata[[col]])
    }
    # !!!
    # all(rownames(coldata) == colnames(cts))
    # coldata has to be in exact order and have same naming as cols in count matrix !!
    rownames(coldata) <- sub("fb", "", rownames(coldata))
    # adjust order of matrix
    cts <- cts[, rownames(coldata)]

    # all(rownames(coldata) == colnames(cts))
    # check! 
    # create the DESEQ-object
    dds <- DESeqDataSetFromMatrix(
        countData = cts,
        colData = coldata,
        design = ~condition
    )
    return(dds)
}
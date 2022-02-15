#' Covert the ouptut of Seurat's FindMarkers to a table that hciR's fgsea_all can use
#'
#' This function assumes that a user has run Seurat's FindMarkers to
#'      assess differential expression between clusters and would like to
#'      next run GSEA with hciR's fgsea_all and gene to pathway data from
#'      hciRdata.
#'
#' @param markers Table from Seurat's FindMarkers
#' @param contrast The name of the contrast
#' @param database Table of genes to use when gathering annotations. Default is hciRdata::human100
#' @param human_homolog Convert the gene names to human homologs. Default is FALSE
#'
#' @examples
#' \dontrun{
#' sc_results = prep_sc_DE_for_fgsea(markers = markers, contrast = "cluster 4 vs. cluster 0")
#' fgsea_all(res = sc_results, gsets = msig_pathways$KEGG, FDR = 0.1, nperm = 10000)
#' }
#'
#' @export

prep_sc_DE_for_fgsea = function(markers = markers, contrast = NULL, database = hciRdata::human100, human_homolog = FALSE){
    # check that contrast is named
    if(is.null(contrast)){
        stop("contrast must be name: e.g., cluster_1 vs cluster_2")
    }

    # set row names from markers (gene names) to a column
    markers = data.frame(gene_name = row.names(markers), markers)

    # reformat for fgsea_all and merge in the annotations
    if(human_homolog == FALSE){
        res = merge(markers[, c(1:3,6)], database[, c(1:4,8)], by = "gene_name")
    } else {
        res = merge(markers[, c(1:3,6)], database[, c(1:4,8,10)], by = "gene_name")
        res$gene_name = res$human_homolog
        res = res[!is.na(res$gene_name) ,]
    }

    # set column name of the log2fold change to match expected by fgsea_all
    colnames(res) = gsub("avg_log2FC", "log2FoldChange", colnames(res))

    # convert to list as expected by fgsea_all and set contrast name
    res = list(res)
    names(res) = contrast
    return(res)
}

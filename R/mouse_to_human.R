#' Convert Loupe differential gene expression results from mouse to human homologs for GSEA
#'
#' After running differential gene expression in Loupe with mouse data,
#' gene names need to be converted to human homologs for GSEA. Returns
#' list where each element is a cluster from the Loupe export
#'
#' @param loupe_export The CSV table exported from Loupe
#' @param mouse_db The mouse database with human homologs. Default is hciRdata::mouse102
#' @param only_named Return only genes with human homologs. Defalt is TRUE
#'
#' @example
#' human_orthologs = loupe_mouse_to_human(loupe_export = "my_markers.csv")
#'
#' @export

loupe_mouse_to_human = function(loupe_export = NULL, mouse_db = hciRdata::mouse102, only_named = TRUE){
  if(is.null(loupe)){
    stop("loupe export not defined")
  }

  # change gene name column in the loupe export
  colnames(loupe)[2] = "gene_name"

  # merge in annotation table
  named_markers = merge(loupe, mouse_db, by = "gene_name")

  # drop the genes that do not have human homologs
  if(only_named == "TRUE"){
    human_ortholog_markers = human_ortholog_markers[!is.na(human_ortholog_markers$human_homolog) ,]
  }

  # get number of columns in loupe for parsing
  n = ncol(loupe)
  n_clusters = (ncol(loupe)-2)/3

  # return tables for each cluster in the loupe export
  l = vector(mode = "list", length = n_clusters)
  k = 0
  for(i in seq(3,n,3)){
    k = l + 1
    l[[k]] = human_ortholog_markers[,c(1,2,i,i+1,i+2,which(colnames(human_ortholog_markers)=="human_homolog"))]
  }

  return(l)
}

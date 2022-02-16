#' Convert a counts table to fraction of total observations
#'
#' Use either rowSums or colSums to convert each cell to a fraction of total observations by row or columns
#'
#' @param df Either a data frame or a table produced by base::table()
#' @param by Input vectors are either "Row" or "Column"
#' @param transpose Transpose resulting table back to original orientation
#'
#' @examples
#' \dontrun{
#' t = percent_table(df = df, by = "Row", transpose = TRUE)
#' }
#'
#' @export
#'
percent_table = function(df = df, by = "Row", transpose = TRUE) {
  if(by != "Row" & by != "Column"){
    stop("Invalid by agument")
  }
  if(class(df) == "Table"){
    df = apply(as.matrix.noquote(df),2,as.numeric)
  } else{
    if(by == "Row"){
      df = apply(df, 1, FUN = function(x) {x/colSums(df)*100})
    } else {
      df = apply(df, 2, FUN = function(x) {x/rowSums(df)*100})
      }
    }
  if(transpose == TRUE){
    df = t(df)
  }
  return(df)
}


#' Export UMAP embeddings from a Suerat object
#'
#' Extract UMAP embedding from a Suerat object and write to file as csv
#' By defult, writes two files, umap_embeddings.csv and clusters.csv for use in Loupe browser
#' By default, drops a cell ID that is prefixed the cell barcode
#'
#' @param srt A Seurat object
#' @param umap_embeddings_file Name of csv file that contains UMAP embeddings
#' @param cluster_id_file Name of csv file that containst he cluster IDs
#' @param Idents Set the active Ident from Seurat metadata before exporting. Becomes grouping factor in clusters.csv. Default is active Ident.
#' @param drop_cell_ids Drop the cell ID prefixed the the cell barcode by Seurat
#'
#' @examples
#' \dontrun{
#' export_umap(srt = srt, embeddings_file = "umap_embeddings.csv", cluster_id_file = "clusters.csv",
#'     Idents = NULL, drop_cell_ids = TRUE)
#' }
#'
#' @export
#'
export_umap = function(srt = srt, embeddings_file = "umap_embeddings.csv", cluster_id_file = "clusters.csv", Idents = NULL, drop_cell_ids = TRUE){

  # extrace cell embeddings from UMAP dimension reduction slot
  embeddings = as.matrix(srt[["umap"]]@cell.embeddings)

  # reassign Idents if need be
  if(is.null(Idents) == FALSE){
    Idents(srt) = Idents
  }

  # get Idents of each cell
  cluster_id = as.matrix(Idents(srt))

  # drop the cell id added during merge
  if(drop_cell_ids == TRUE){
    row.names(embeddings) = gsub(".*_", "", row.names(embeddings))
    row.names(cluster_id) = gsub(".*_", "", row.names(cluster_id))
  }

  # write to file
  write.csv(embeddings, file = embeddings_file, quote = FALSE)
  write.csv(cluster_id, file = cluster_id_file, quote = FALSE)

}


#' Save the results of FindAllMarkers to Excel
#'
#' The large table produced by FindAllMarkers is split by cluster and saved to Excel
#' When file name is provided the named output is Excel
#' Returns list when file name not provided
#'
#' @param markers The product of Seurat's FindAllMarkers function
#' @param file Name of Excel file. If absent, the default, returns list
#'
#' @examples
#' \dontrun{
#' write_FindAllMarkers(markers = markers, file = "markers.xlsx")
#' }
#'
#' @export
#'
write_FindAllMarkers = function(markers = markers, file = NA) {
  # split the large markers data frame by cluster into list elements
  markers_list = vector(mode = "list", length = length(unique(markers$cluster)))
  for(i in 0:(length(markers_list)-1)) {
    markers_list[[i+1]] = markers[markers$cluster == i ,]
  }

  # if file name not provided then return list else write to file
  if(is.na(file)){
    names(markers_list) = paste("cluster ", seq(0,length(unique(markers$cluster))-1,1), sep = "_")
    return(markers_list)
  } else {
    # write to Excel, one sheet per cluster
    message(paste("saving", file, sep = " "))
    openxlsx::write.xlsx(x = markers_list, file = file, sheetName = paste("cluster ", seq(0,length(unique(markers$cluster))-1,1), sep = ""), row.names = TRUE)
  }
}

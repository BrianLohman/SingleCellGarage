#' Build file manifest from 10X cellranger-atac
#'
#' Assumes that each sample is in its own folder. Pattern match to get multiple samples
#' Returns named list where each sample has paths to peaks.bed, singlecell.csv, and fragments file
#' Also possible to build list or make corrections manually
#'
#' @param pattern Common string among all samples
#'
#' @examples
#' \dontrun{
#' m = build_manifest(pattern = "18655")
#' }
#'
#' @export

build_manifest <- function(pattern = NULL){
  # sample ID: one directory for each sample
  s = list.files(pattern = pattern)
  l = vector(mode = "list", length = length(s))
  names(l) = s

  # peaks.bed
  for(i in 1:length(l)){
    l[[i]]["peaks"] = list.files(path = paste("./", names(l)[i],  sep = ""), pattern = "peaks.bed", full.names = TRUE)
  }

  # singlecell.csv, the metadata
  for(i in 1:length(l)){
    l[[i]]["metadata"] = list.files(path = paste("./", names(l)[i],  sep = ""), pattern = "singlecell.csv", full.names = TRUE)
  }

  # fragments file, get only first match (second is tabix index)
  for(i in 1:length(l)){
    l[[i]]["fragments"] = list.files(path = paste("./", names(l)[i],  sep = ""), pattern = "fragments.tsv.gz", full.names = TRUE)[1]
  }

  return(l)
}

#' Read in multiple 10X ATAC samples with Signac and merge with common peak set
#'
#' Takes file manifest from build_manifest
#' following preffered merge procedure in Signac: https://satijalab.org/signac/articles/merging.html#merge-fragment-files-1
#'
#' @param manifest The file manifest from build_manifest or a list of named lists where each entry is a sample with peaks.bed, singlecell.csv and fragments
#' @param project Name to give merged sample set
#'
#' @examples
#' \dontrun{
#' srt = build_10X_atac(manifest = m, project = "A6440")
#' }
#'
#' @export

build_10x_atac <- function(manifest = NULL, project = NULL){
  if(is.null(manifest)){
    stop("file manifest undefined")
  }

  if(is.null(project)){
    stop("project undefined")
  }

  # assign file manifest
  m = manifest

  # read in peaks files to list
  print("building peaks list")
  peaks_list = vector(mode = "list", length = length(m))
  names(peaks_list) = names(m)

  for(i in names(peaks_list)){
    peaks_list[[i]] = read.table(
      file = as.character(m[[i]]["peaks"]),
      col.names = c("chr", "start", "end")
    )
  }

  # convert to genomic ranges
  print("converting to genomic ranges and building unifed set")
  gr_list = vector(mode = "list", length = length(m))
  names(gr_list) = names(m)

  for(i in names(gr_list)){
    gr_list[[i]] <- makeGRangesFromDataFrame(peaks_list[[i]])
  }

  # Create a unified set of peaks to quantify in each dataset
  combined.peaks <- reduce(GRanges(Reduce("c", gr_list)))

  # Filter out bad peaks based on length
  print("prefiltering")
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

  # load metadata
  print("building metadata list")
  metadata_list = vector(mode = "list", length = length(m))
  names(metadata_list) = names(m)

  for(i in names(metadata_list)){
    metadata_list[[i]] = read.table(
      file = as.character(m[[i]]["metadata"]),
      stringsAsFactors = FALSE,
      sep = ",",
      header = TRUE,
      row.names = 1
    )[-1, ] # remove the first row
  }

  # perform an initial filtering of low count cells
  print("filtering low count cells")
  for(i in names(metadata_list)){
    metadata_list[[i]] = metadata_list[[i]][metadata_list[[i]]$passed_filters > 200, ]
  }

  # create fragment objects
  print("building fragment list")
  fragments_list = vector(mode = "list", length = length(m))
  names(fragments_list) = names(m)

  for(i in names(fragments_list)){
    fragments_list[[i]] <- CreateFragmentObject(
      path = as.character(m[[i]]["fragments"]),
      cells = rownames(metadata_list[[i]])
    )
  }

  # quantify peaks in each data set
  print("quantifying peaks in parallel")
  build_feature_matrix <- function(sample){
    FeatureMatrix(
      fragments = fragments_list[[sample]],
      features = combined.peaks,
      cells = rownames(metadata_list[[sample]])
    )
  }
  counts_list = mclapply(names(m), build_feature_matrix)
  names(counts_list) = names(m)

  # Create Seurat objects
  print("building Seurat objects")
  seurat_list = vector(mode = "list", length = length(m))
  names(seurat_list) = names(m)

  for(i in names(seurat_list)){
    tmp = CreateChromatinAssay(counts_list[[i]], fragments = fragments_list[[i]])
    seurat_list[[i]] = CreateSeuratObject(tmp, assay = "ATAC", project = i, meta.data = metadata_list[[i]])
  }

  # merge all datasets, adding a cell ID to make sure cell names are unique
  print("merging Seurat objects")
  srt <- merge(
    x = seurat_list[[1]],
    y = seurat_list[2:length(m)],
    add.cell.ids = names(m),
    project = project
  )

  return(srt)
}

#' Read in multiple 10X scRNAsea samples
#'
#' Searches for pattern in current working directory and returns list of Seurat objects
#'
#' @param pattern The pattern to search for in current working dicrectory, e.g, run number
#' @param min_cell Seurat min.cells
#' @param min_features Seurat min.features
#' @param outs Include the "outs" in the full path to filtered barcode files. Default is FALSE
#'
#' @examples
#' \dontrun{
#' srt_list = read10X_multiple(pattern = "18560X)
#' }
#'
#' @export
#'
#'
read10X_multiple = function(pattern = NULL, min_cells = 0, min_features = 0, outs = FALSE){
  if(is.null(pattern)){
    stop("sample search pattern not defined")
  }
  # list directories matching pattern
  dirs <- list.files(pattern = pattern)

  # get sample names
  srt_list = vector(length = length(dirs), mode = "list")
  names(srt_list) = dirs

  # load each sample
  for(i in 1:length(dirs)){
    print(paste("loading", names(samples)[i]))
    if(outs == TRUE){
      path = Read10X(data.dir = paste(dirs[i], "outs","filtered_feature_bc_matrix", sep = "/"))
    } else {
      path = Read10X(data.dir = paste(dirs[i], "filtered_feature_bc_matrix", sep = "/"))
    }
    srt_list[[i]] = CreateSeuratObject(counts = path, project = dirs[i], min.cells = min_cells, min.features = min_features)
  }

  return(srt_list)
}


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

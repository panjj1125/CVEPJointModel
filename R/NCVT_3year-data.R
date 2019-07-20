#' 3-year NCVT dataset
#'
#' 3-year data from year 1964 to 1966 for
#' Joint Modeling in Analyzing Highly Unbalanced Multi-Environment Trial Data
#'
#' @docType data
#'
#' @usage data(NCVT_3year)
#'
#'
#' @keywords datasets
#'
#'
#' @source \href{https://phenome.jax.org/projects/Moore1b}{QTL Archive}
#'
#' @examples
#' data(NCVT_3year)
#' re <- CVEP_JM(dat, factors = c("Year","Loc","Rep", "Variety"),
#' TT_mm=list(fixed=c("Year","Loc"),random=c("Variety","Variety:Year","Variety:Loc")),
#' DS_mm=list(fixed=c("Year","Loc"),random=c("Variety")),
#' converg_control=list(nsamp=5,max.iter=5,err=10^(-5),err1=10^(-5)))
#'
"NCVT_3year"

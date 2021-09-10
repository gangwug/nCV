#' Evaluate robustness of clock genes in population scale data
#'
#' This function calculate nCV values of targeted clock genes. The nCV values can indicate
#' the robustness of clock genes. Higher nCV value means more robust of oscillation.
#' Before applying this function, please do the clock functionality test using \code{nCVnet}.
#' If the \code{nCVnet} output significant p-value, it is good to run \code{nCVgene}.
#'
#' @param inputD a data frame. The first column is the gene symbol, and other columns are samples.
#' One row per gene. Do not accept data frame with duplicate gene symbols.
#' @param cgenes a character vector. The targeted clock genes for returning the nCV values.
#' @return A data frame with calculated nCV values of targeted clock genes.
#' @examples
#' ## calculate the nCV values of targeted clock genes
#' cgenes = c("ARNTL", "CLOCK", "NPAS2", "CRY1", "NR1D1", "CIART", "DBP", "PER1", "CRY2", "PER2")
#' ncvD = nCVgene(inputD = nCVegD, cgenes)
#' ncvD
#' @export

nCVgene <- function(inputD, cgenes) {
  ##set row and column names
  inputD = as.data.frame(inputD)
  colnames(inputD)[1] = "geneSym"
  if ( nrow(inputD) == length(unique(inputD$geneSym)) ) {
    rownames(inputD) = inputD$geneSym
    ##get the reformatted data frame
    ncvD = inputD %>%
      tidyr::gather(-geneSym, key = "sampleID", value = "expv") %>%
      split(.$geneSym) %>%
      purrr::map( function(zD) {
        x = zD$expv
        xmean = mean(x)
        xsd = sd(x)
        xcv = xsd/xmean
        return(data.frame( CV = xcv, geneSym = unique(zD$geneSym)) )
      }) %>% bind_rows() %>%
      dplyr::mutate( nCV = CV / mean(CV) ) %>%
      #dplyr::select(geneSym, CV, nCV) %>%
      dplyr::select(geneSym, nCV) %>%
      dplyr::filter(geneSym %in% cgenes)
    return(ncvD)
  }  else  {
    cat("Duplicate gene ids are detected. Please merge rows with duplicate gene ids.\n")
  }
}


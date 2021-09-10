#' Evaluate functionality of circadian clock network in population scale data
#'
#' This function can test whether there is a functional clock network in population scale data.
#' If the p-value of this test is significant, the next step is to get the nCV values of target genes.
#'
#' @param inputD a data frame. The first column is the gene symbol, and other columns are samples.
#' One row per gene. Do not accept data frame with duplicate gene symbols.
#' @param benchD a data frame. The expression correlation values of paired clock genes. The
#' @param hs a logical value. If the samples are from human, set \code{"hs"} as \code{TRUE}.
#' If the samples are from mouse, set \code{hs} as \code{FALSE}.
#' @param nperm a numeric value. The number of permutations used to calculate the p-value. Default is 1000.
#' @param seedN a numeric value. Specify the seeds for random sampling.
#' @return A list containing the following components:
#' \tabular{rl}{
#'              zstat    \tab  the permutation test results of \code{nCVnet}\cr
#'              npermV   \tab  a vector of Z-statistic (Mantel test) of permutations\cr
#'              cmatrix  \tab  a matrix of correlation values of paired clock genes using the input data
#'              }
#' @examples
#' ## test the whether there is a functional clock network running in the example data
#' testD = nCVnet(inputD = nCVegD, benchD = mClockD)
#' testD$zstat
#' @export
#' @importFrom ape mantel.test
#' @importFrom dplyr inner_join


nCVnet <- function(inputD, benchD, hs = TRUE, nperm = 1000, seedN = 10) {
  ##convert mouse gene symbol to human gene symbol
  if (hs) {
    ##since only test the clock genes, it is fine to just use upper case function
    ##for other genes, need to use the homolog gene function for finding the human homolog genes of a mouse gene list
    benchD[,1] <- toupper(benchD[,1])
    colnames(benchD) <- toupper(colnames(benchD))
  }
  ##set row and column names
  inputD <- as.data.frame(inputD)
  benchD <- as.data.frame(benchD)
  colnames(inputD)[1] <- colnames(benchD)[1] <- "geneSym"
  if ( nrow(inputD) == length(unique(inputD$geneSym)) ) {
    rownames(inputD) <- inputD$geneSym
    rownames(benchD) <- benchD$geneSym
    ##get the overlapped genes
    bothID <- inner_join(data.frame(id = benchD$geneSym), data.frame(id = inputD$geneSym), by = "id")
    if ( nrow(bothID) >= 3 )  {
      corR <- as.matrix(benchD[bothID$id,bothID$id])
      corD <- cor(t(inputD[bothID$id,-1]), method = "spearman")
      rownames(corD) = colnames(corD) = bothID$id
      ##calculate the similarity between matrix with 'ape' package
      simD <- mantel.test(corD, corR, nperm = 1)
      ##start the permutation step
      set.seed(seedN)
      indexM <- NULL
      for (np in 1:nperm) {
        indexM <- rbind(indexM, sample(1:nrow(inputD), nrow(bothID)))
      }
      pstat <- apply(indexM, 1, function(indexv) {
        corP <- cor(t(inputD[indexv,-1]), method = "spearman")
        tepD <-  mantel.test(corP, corR, nperm = 1)
        return(tepD$z.stat)
      })
      pva <- sum(pstat > simD$z.stat)/nperm
      return(list(zstat = data.frame(tag = "nCVnet",
                                     geneNumber = nrow(bothID),
                                     zstat = simD$z.stat,
                                     pvalue = pva),
                  npermV = pstat,
                  cmatrix = corD))
    } else {
      cat("Less than 3 genes are overlapped with the bench gene list. Please check the input data.\n")
    }
  } else {
    cat("Duplicate gene ids are detected. Please merge rows with duplicate gene ids.\n")
  }
}

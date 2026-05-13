#' Toy dataset generated through [`pdRCON.simulate`] function
#'
#' Data simulated by the function [`pdRCON.simulate`] with parameters:
#' * `p=20`
#' * `type=c("v", "i","a")`
#' * `force.symm=NULL`
#' * `dens=0.5`
#' * `Sigma` a \eqn{p \times p} matrix with unitary diagonal and 0.5 on off-diagonal elements.
#'
#' @format ## `toy_data`
#' A list with the following components:
#'
#' * `pdColG`, a matrix of class `pdColG`; see [`pdglasso-package`] for details.
#' * `K`, a concentration matrix adapted to `pdColG`.
#' * `sample.data`, a data frame with 250 rows and 20 columns from a multivariate normal distribution
#' with zero mean vector and concentration matrix `K`.
#'
"toy_data"


#' fMRI Dataset
#'
#' This dataset contains residuals obtained through score-driven models on the time-series of brain activity at rest for an individual,
#' measured on 10 regions of interest of the brain in the parietal area (five on the left hemisphere and five paired regions on the right hemisphere).
#' Please check Sections 2 and 7 of Ranciati and Roverato (2024) for more information.
#'
#' @format ## `fMRI_parietal`
#' A \eqn{404 \times 10} dataframe with 404 observations referring to residuals from \eqn{p=5\times 2=10} time-series
#' after filtering through a score-driven model.
#' The residuals are associated to 10 different region of interest of the parietal area of the brain, five on the left empisphere e the paired 
#' homologous regions on the right hemisphere.
#'
#'
#'
"fMRI_parietal"

#' Breast Cancer Dataset
#'
#' The paired samples in this dataset refer to individuals with both tumor and healthy adjacent tissue measurements, in the form of a set of genes of the Hedgehog Pathway
#' and their gene-level transcription estimates, i.e. transformed normalized counts. 
#' Please check Section 8 of Ranciati and Roverato (2026) for more information.
#'
#' @format ## `bcdata`
#' A named \eqn{114 \times 178} matrix, where:
#' * each row refers to one of the 114 individuals in the dataset;
#' * each column is one of the \eqn{p=89 \times 2=178} genes of the Hedgehog Pathway, whose expressions were measured on healthy (first half) and adjacent tumor affected samples (second half)
#' 
"bc_data"
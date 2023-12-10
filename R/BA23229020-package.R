#' @name BA23229020-package
#' @docType package
#' @keywords package
#' @title statistical calculation
#' @description The BA23229020 package implements methods for efficient product quantization 
#'   of high-dimensional data. It is designed to facilitate the compression and retrieval 
#'   of high-dimensional vectors, which is particularly useful in machine learning and 
#'   nearest neighbor search in large-scale databases.
#'
#'   The package includes key functions such as `pqCompress` for compressing data into 
#'   product quantization format and `pqSearch` for performing efficient nearest neighbor 
#'   searches within the compressed space. These functions allow for a significant 
#'   reduction in storage space and search time compared to dealing with raw high-dimensional data.
#'
#' @details The main functions of the package include:
#'   - `pqCompress`: Compresses high-dimensional vectors into a reduced product quantization format.
#'   - `pqSearch`: Performs a nearest neighbor search using the compressed format.
#'   - other functions for HW
#'
#'   Along with these, the package also provides utility functions for analyzing and 
#'   manipulating the compressed data representations.
#'
#' @section Installation: 
#' The package can be installed from a source package or directly from the repository using:
#' \preformatted{
#' devtools::install_github("jzhang-0/zhangjin-BA23229020")
#' }
#'
#' @section Usage:
#' After installation, load the package using:
#' \preformatted{
#' library(BA23229020)
#' }
#' Refer to the documentation for individual functions to learn more about their usage.
#'
#' @references (Jegou, Herve, Matthijs Douze, and Cordelia Schmid. "Product quantization for nearest neighbor search." IEEE transactions on pattern analysis and machine intelligence 33, no. 1 (2010): 117-128.)
#'
#'
#' @import ggplot2
#' @import boot
#' @import bootstrap
#' @import coda
#' @import stats4
#' @import lpSolve
#' @import Rcpp
#' @import microbenchmark
#' @useDynLib BA23229020
NULL

#' Product Quantization Search
#'
#' 这个函数执行基于乘积量化的搜索。
#'
#' @param query 查询向量。
#' @param compressed_data 压缩后的数据。
#' @param centroids 质心列表。
#' @return 返回搜索结果的索引。
#' @export
#' @examples
#' data <- matrix(rnorm(1000), ncol = 100)
#' compressed <- pqCompress(data, 10, 5)
#' query <- rnorm(100)
#' nearest_neightbors <- pqSearchCpp(query, compressed$compressed_data, compressed$centroids)
pqSearchCpp <- function(query, compressed_data, centroids) {
  .Call('_BA23229020_pqSearchCpp', PACKAGE = 'BA23229020', query, compressed_data, centroids)
}

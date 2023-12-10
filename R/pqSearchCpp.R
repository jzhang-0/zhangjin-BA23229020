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
#' # 示例代码...
pqSearchCpp <- function(query, compressed_data, centroids) {
  .Call('_BA23229020_pqSearchCpp', PACKAGE = 'BA23229020', query, compressed_data, centroids)
}

# Help documentation for pqCompress
#' @title Product Quantization Compression
#' @description Compresses high-dimensional data using product quantization.
#' @param data A matrix of high-dimensional data points.
#' @param num_subspaces The number of subspaces.
#' @param num_clusters The number of clusters in each subspace.
#' @return A list containing the centroids and the compressed data representation.
#' @importFrom stats kmeans
#' @export
#' @examples
#' data <- matrix(rnorm(1000), ncol = 100)
#' compressed <- pqCompress(data, 10, 5)
pqCompress <- function(data, num_subspaces, num_clusters) {
  # Validate inputs
  if (!is.matrix(data)) {
    stop("Data must be a matrix.")
  }
  if (num_subspaces <= 0 || num_clusters <= 0) {
    stop("Number of subspaces and clusters must be positive.")
  }

  # Split the data into subspaces
  dim_per_subspace <- ncol(data) / num_subspaces
  if (dim_per_subspace != as.integer(dim_per_subspace)) {
    stop("Number of columns in data must be divisible by num_subspaces.")
  }

  # Placeholder for centroids and compressed data
  centroids <- list()
  compressed_data <- matrix(0, nrow(data), num_subspaces)

  # Apply k-means in each subspace
  for (i in 1:num_subspaces) {
    start_col <- (i - 1) * dim_per_subspace + 1
    end_col <- i * dim_per_subspace
    subspace_data <- data[, start_col:end_col]

    # Apply k-means clustering
    kmeans_result <- kmeans(subspace_data, centers = num_clusters, iter.max = 20)
    centroids[[i]] <- kmeans_result$centers
    for (j in 1:nrow(data)) {
    # 对于每个数据点，计算其与所有质心的距离
        distances <- apply(kmeans_result$centers, 1, function(center) {
        sum((subspace_data[j, ] - center)^2)
    })
    # 找到距离最小的质心
    compressed_data[j, i] <- which.min(distances)
    }

  }

  return(list(centroids = centroids, compressed_data = compressed_data))
}



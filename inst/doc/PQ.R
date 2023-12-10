## -----------------------------------------------------------------------------
library(BA23229020)

## -----------------------------------------------------------------------------
set.seed(123)
d <- 64
data <- matrix(rnorm(1000 * d), ncol = d)

## -----------------------------------------------------------------------------
num_subspaces <- 8
num_clusters <- 256
compressed <- pqCompress(data, num_subspaces, num_clusters)

## -----------------------------------------------------------------------------
head(compressed$compressed_data)

## -----------------------------------------------------------------------------
# Initialize the matrix for reconstructed data
reconstructed_data <- matrix(0, nrow = nrow(compressed$compressed_data), ncol = d)

# Reconstruct the original data for each subspace
for (i in 1:num_subspaces) {
  # Get the compressed indices for the current subspace
  subspace_indices <- compressed$compressed_data[, i]
  # Create a matrix of centroids for the current subspace
  subspace_centroids <- matrix(unlist(compressed$centroids[[i]]), byrow = FALSE, ncol = d / num_subspaces)
  # Ensure that the indices are within the bounds of the centroids matrix
  if (max(subspace_indices) > nrow(subspace_centroids)) {
    stop("One of the indices is out of bounds of the centroids matrix.")
  }
  # Add the centroids of each subspace to the corresponding position in the reconstructed data matrix
  for (j in 1:nrow(compressed$compressed_data)) {
    reconstructed_data[j, ((i - 1) * (d / num_subspaces) + 1):(i * (d / num_subspaces))] <-
      subspace_centroids[subspace_indices[j], ]
  }
}

## -----------------------------------------------------------------------------
mse <- sum((data - reconstructed_data)^2) / sum((data)^2)

print(paste("Relative Mean Squared Error:", mse))

## -----------------------------------------------------------------------------
set.seed(123)
queries <- matrix(rnorm(100 * d), ncol = d)

nearest_neighbors_list <- lapply(1:100, function(i) pqSearchCpp(queries[i, ], compressed$compressed_data, compressed$centroids))

## -----------------------------------------------------------------------------
# Function to calculate recall at k for a single relevant item
calculate_recall_at_k <- function(query, data, compressed_data, true_rankings, retrieved_neighbors, k) {
  # Calculate the true k nearest neighbors based on the original data
  true_neighbors <- head(true_rankings, 1)

  # Calculate the intersection of true and retrieved neighbors
  common_neighbors <- length(intersect(true_neighbors, retrieved_neighbors[1:k]))

  # Recall@k is the proportion of true neighbors that were retrieved
  recall_at_k <- common_neighbors
  return(recall_at_k)
}

## -----------------------------------------------------------------------------
distance <- function(x, y) {
  sqrt(sum((x - y)^2))
}

# k values
k_values <- c(1, 4, 8, 16, 32, 64)

average_recalls <- numeric(length(k_values))

# compute recall
for (idx in seq_along(k_values)) {
  k <- k_values[idx]
  recalls_at_k <- sapply(1:100, function(i) {
    # Compute distances from the query to all points in the original data
    distances <- apply(data, 1, distance, queries[i, ])

    # Get the ranking of the original points based on their distance to the query
    true_rankings <- order(distances)

    # Get the top k retrieved neighbors from the compressed data search
    retrieved_neighbors <- head(nearest_neighbors_list[[i]], k)

    # Calculate recall at k for the current query
    calculate_recall_at_k(queries[i, ], data, compressed$compressed_data, true_rankings, retrieved_neighbors, k)
  })

  # Store the average recall at the current k
  average_recalls[idx] <- mean(recalls_at_k)
}

# Create a data frame for plotting
recall_data <- data.frame(k = k_values, Recall = average_recalls)


library(ggplot2)
ggplot(recall_data, aes(x = k, y = Recall)) +
  geom_point() +
  geom_line(group = 1, colour = "blue") + 
  scale_x_continuous(breaks = k_values) + 
  labs(title = "Recall Curve", x = "k", y = "Average Recall@k") +
  theme(
    panel.grid.major = element_blank(), 
  )


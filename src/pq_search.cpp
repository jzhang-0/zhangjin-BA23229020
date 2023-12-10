#include <Rcpp.h>
using namespace Rcpp;

// 辅助函数用于比较两个值对
bool cmp(const std::pair<double, int> &a, const std::pair<double, int> &b)
{
    return a.first < b.first;
}

// [[Rcpp::export]]
NumericVector pqSearchCpp(NumericVector query, NumericMatrix compressed_data, List centroids)
{
    int num_subspaces = centroids.size();
    int query_dim = query.size() / num_subspaces;
    std::vector<double> distances(compressed_data.nrow(), 0.0);

    //逐个空间建表，计算距离
    for (int i = 0; i < num_subspaces; ++i)
    {
        NumericVector subspace_query = query[Range(i * query_dim, (i + 1) * query_dim - 1)];
        NumericMatrix subspace_centroids = as<NumericMatrix>(centroids[i]); // K * d_sub
        NumericVector lookup_table(subspace_centroids.nrow());

        for (int j = 0; j < subspace_centroids.nrow(); ++j)
        {
            NumericVector centroid = subspace_centroids(j, _);
            lookup_table[j] = sum(pow(subspace_query - centroid, 2.0));
        }

        for (int k = 0; k < compressed_data.nrow(); ++k)
        {
            distances[k] += lookup_table[compressed_data(k, i) - 1];
        }
    }

    // 创建一个索引向量和对应的距离
    std::vector<std::pair<double, int> > indexed_distances;
    for (int i = 0; i < distances.size(); ++i)
    {
        indexed_distances.push_back(std::make_pair(distances[i], i + 1)); // R中的索引从1开始
    }

    // 排序索引向量
    std::sort(indexed_distances.begin(), indexed_distances.end(), cmp);

    // 提取排序后的索引
    NumericVector ordered_indices(distances.size());
    for (int i = 0; i < distances.size(); ++i)
    {
        ordered_indices[i] = indexed_distances[i].second;
    }

    return ordered_indices;
    // return List::create(Named("distances") = distances, Named("ordered_indices") = ordered_indices);
}

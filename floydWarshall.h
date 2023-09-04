#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <climits>
#include <cassert>

// Authored by Xuhui Li, https://github.com/Emerlee

// Floyd-Warshall algorithm helps finding the path from a node in the graph to another with the minimal cost (distance).
// The basic idea is when traveling from point A to point C, if point B is on the optimal path, then we have /distance(A, C) = distance(A, B) + distance(B, C)/.
// Unlike Dijkstra algorithm, the cost (distance) between two nodes can be negative. Also, this algorithm can also get the minimal distance between any two nodes.


// Two tables are needed:
//   -- distance table (D[i][j] stores the minimal distance from node [i] to node [j], initialized as inf)
//   -- sequence table (S[i][j] stores the next node from node [i] to node [j], e.g. if S[i][j] == k (k != j), then the path is i->k->j. initialized as [j])

template <typename T>
void floydWarshall(
    const std::vector<std::vector<T>>& adjacency_matrix,
    std::vector<std::vector<T>>& distance,
    std::vector<std::vector<int>>& sequence
    ) {
    // If the graph is directional, both adjacency list and adjacency matrix are fine. But if the graph is non-directional, we need to transfer the adjacency list into an adjacency table;
    
    // initialization
    int n = distance.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            distance[i][j] = adjacency_matrix[i][j];
            sequence[i][j] = j;
        }
    }
        
    // dynamic programming
    for (int k = 0; k < n; k++) { // loop over all possible middle nodes
        for (int i = 0; i < n; i++) { // loop over all start nodes
            for (int j = 0; j < n; j++) { // loop over all destination nodes
                if (distance[i][k] != INT_MAX && distance[k][j] != INT_MAX && distance[i][k] + distance[k][j] < distance[i][j]) {
                    distance[i][j] = distance[i][k] + distance[k][j];
                    sequence[i][j] = sequence[i][k];
                }
            }
        }
    }
}


template <typename T>
std::vector<std::vector<T>> adjacencyList2Matrix(const std::vector<std::vector<T>>& adjacency_list, int n, bool directional=false) {
    // This function transforms a bi-directional adjacency list into an adjacency matrix
    assert(adjacency_list[0].size() == 3);
    
    std::vector<std::vector<T>> adjacency_matrix(n, std::vector<T>(n, INT_MAX));
    for (int i = 0; i < n; i++) {
        adjacency_matrix[i][i] = 0;
    }
    
    if (directional) {
        for (int i = 0; i < adjacency_list.size(); i++) {
            adjacency_matrix[adjacency_list[i][0]][adjacency_list[i][1]] = adjacency_list[i][2];
        }
    } else {
        for (int i = 0; i < adjacency_list.size(); i++) {
            adjacency_matrix[adjacency_list[i][0]][adjacency_list[i][1]] = adjacency_list[i][2];
            adjacency_matrix[adjacency_list[i][1]][adjacency_list[i][0]] = adjacency_list[i][2];
        }
    }
    
    return adjacency_matrix;
}

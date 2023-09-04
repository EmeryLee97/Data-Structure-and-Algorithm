#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>
#include <cassert>

// Authored by Xuhui Li, https://github.com/Emerlee

// Dijkstra algorithm is based on greedy strategy. It helps finding the smallest cost from start node to any destination node, and also records the shortest path.
//            8       7
//        1⃣️-----2⃣️------3⃣️        -- The dis from 0 to 1 = 4. 0->1
//       /|       | \     | \      -- The min dis from 0 to 2 = 12. 0->1->2
//    4 / |     2 |  \    |  \9    -- The min dis from 0 to 3 = 19. 0->1->2->3
//     /  |       |   \ 4 |   \    -- The min dis from 0 to 4 = 21. 0->7->6->5->4
//    0⃣️  |11   /8⃣️    \  |14 4⃣️   -- The min dis from 0 to 5 = 11. 0->7->6->5
//     \  |    /  |     \ |   /    -- The min dis from 0 to 6 = 9. 0->7->6
//    8 \ |  / 7  | 6    \|  /10   -- The min dis from 0 to 7 = 8. 0->7
//       7⃣️/-----6⃣️------5⃣️/       -- The min dis from 0 to 8 = 14. 0->1->2->8
//            1       2
// https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-greedy-algo-7/

// Three lists are needed to store nodes information
//   -- visited  (One node will not be revisited, initialized as false)
//   -- distance (The distance from the start node to any node, initialized as INF)
//   -- parent (From which node we arrive the current node, initialized as -1)

// Steps:
// 1. Each time a new node is visited, update the lowest cost to it's adjacent nodes, also change their parent as the current node.
// 2. Find the node which is unvisited and has the lowest cost, repeat step 1.

template <typename T>
T dijkstra(const std::vector<std::vector<T>>& adjacency_matrix, int n, int node) {
    // If the graph is directional, both adjacency list and adjacency matrix are fine. But if the graph is non-directional, we need to transfer the adjacency list into an adjacency table;
    
    assert(("Please input a valid node number!", node < n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            assert("The cost between two nodes can not be negative!", adjacency_matrix[i][j] >= 0;)
        }
    }
    
    std::vector<bool> visited(n, false);
    std::vector<T> distance(n, INT_MAX);
    std::vector<int> parent(n, -1);
    
    // Initialization
    distance[0] = 0;
    parent[0] = 0;
    
    for (int i = 0; i < n; i++) {
        // Find the unvisited node with minimun distance
        int curr_node = 0;
        T curr_distance = INT_MAX;
        for (int j = 0; j < n; j++) {
            if (!visited[j] && distance[j] < curr_distance) {
                curr_node = j;
                curr_distance = distance[j];
            }
        }
        visited[curr_node] = true;
        
        // Update it's adjacent nodes
        for (int j = 0; j < n; j++) {
            if (visited[j] || adjacency_matrix[curr_node][j] == INT_MAX) continue;
            int next_distance = distance[curr_node] + adjacency_matrix[curr_node][j];
            if (next_distance < distance[j]) {
                distance[j] = next_distance; // update distance
                parent[j] = curr_node; // update parent
            }
        }
    }
    // Returns the minimal distance from node0 to any other node;
    return distance[node];
}


//   adjacency list:              adjacency matrix:
//     0----1----4             0  1  2  3  4  5  6  7  8
//     0----7----8          0  0  4  N  N  N  N  N  8  N
//     1----2----8          1  4  0  8  N  N  N  N  11 N
//     1----7----11         2  N  8  0  7  N  4  N  N  2
//     2----3----7          3  N  N  7  0  9  14 N  N  N
//     2----5----4          4  N  N  N  9  0  10 N  N  N
//     2----8----2          5  N  N  4  14 10 0  2  N  N
//     3----4----9          6  N  N  N  N  N  2  0  1  6
//     3----5----14         7  8  11 N  N  N  N  1  0  7
//     4----5----10         8  N  N  2  N  N  N  6  7  0
//     5----6----2
//     6----7----1                  (N = INF)
//     6----8----6
//     7----8----7

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

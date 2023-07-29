# include <iostream>
# include <vector>
# include <algorithm>
# include "unionFind.h"

// minimun spanning tree, Kruskal's algorithm
// https://en.wikipedia.org/wiki/Kruskal%27s_algorithm

// This function takes a (n, 3) vector as input, where vec[i][0] and vec[i][1]
// means the two nodes to be connected together, while vec[i][2] means the cost
// to connect them.
// It returns the minimun cost to connect all nodes together.

int minimunSpanningTree_K(n, std::vector<std::vector<int>>& vec) {
    // 1. Sort the given vector according to the 2nd colomn (costs)
    std::sort(
        vec.begin(),
        vec.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) {
            return a.at(2) < b.at(2);}
    );
    
    // 2. iterate throuth the sorted vector and union the two nodes together.
    //    if the connection doesn't form a circular, add the cost to sum,
    //    otherwise continue;
    int sum_cost{0};
    int count{0};
    UnionFind uf(n);
    for (const auto& elem : vec) {
        if (uf.find(elem[0]-1) != uf.find(elem[1]-1)) {
            uf.Union(elem[0]-1, elem[1]-1);
            sum_cost += elem[2];
            count++;
        } else {
            continue;
        }
    }
    if (count == n-1) return sum_cost;
    else return -1;
}


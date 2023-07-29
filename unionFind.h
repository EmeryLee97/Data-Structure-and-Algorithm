#include <iostream>
#include <vector>

// Union Find Set including quick find
// (union according to their ranks)

class UnionFind {
private:
    std::vector<int> _parents;
    std::vector<int> _rank;
    
public:
    UnionFind(int size) {
        for (int i = 0; i < size; i++) {
            // initialize the index as itselfs parent
            _parents.push_back(i);
            // initialize all nodes' ranks as 1
            // (The tree onlt contains itself)
            _rank.push_back(1);
        }
    }
    
    // This function returns the root(not the parent) of a node
    int find(int node) {
        if (_parents[node] == node) {
            return node;
        } else {
            return _parents[node] = find(_parents[node]);
        }
    }
    
    // This function unions two trees together
    // "union is already used in c++"
    void Union(int x, int y) {
        int rootx = find(x);
        int rooty = find(y);
        if (rootx != rooty) {
            if (_rank[rootx] > _rank[rooty]) {
                // The tree x belongs to is higher than the tree y belongs to
                // Put y's root under x's root
                _parents[rooty] = rootx;
            } else if (_rank[rootx] < _rank[rooty]) {
                // The tree x belongs to is shorter than the tree y belongs to
                // Put x's root under y's root
                _parents[rootx] = rooty;
            } else {
                // Two trees have the same height
                // It doesn't matter to choose who as the root
                _parents[rooty] = rootx;
                _rank[rootx]++;
            }
        }
    }
};

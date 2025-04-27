//
// Created by jawad on 26/04/2025.
//

#include "preprocessing.h"
#include <iostream>
#include <vector>
#include <utility> // for pair

using namespace std;

void preprocessInsertedEdges(const vector<pair<int, int>>& insertedEdges, int numVertices, vector<vector<pair<int, int>>>& I) {
    I.clear();
    I.resize(numVertices);  // Initialize I[v] lists

    for (const auto& edge : insertedEdges) {
        int u = edge.first;
        int v = edge.second;

        I[v].push_back({u, v});  // Correct: add into I[v], not I[u]
    }
}

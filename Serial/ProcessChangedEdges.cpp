//
// Created by jawad on 26/04/2025.
//

#include "ProcessChangedEdges.h"
#include <iostream>
#include <vector>
#include <utility> // for pair
#include <limits>  // for infinity
#include <algorithm> // for min
#include <map>

using namespace std;

void processChangedEdges(
    const vector<vector<pair<int, int>>>& I,
    vector<float>& delta,
    vector<int>& parent,
    vector<int>& marked,
    vector<int>& Aff,
    const map<pair<int, int>, float>& weightMap
) {
    int V = delta.size();
    Aff.clear();
    marked.assign(V, 0);  // Step 5: initialize marked with 0

    // Step 6: for each vertex v in parallel (simulate sequentially for now)
    for (int v = 0; v < V; ++v) {
        // Step 7: for each incoming edge (u, v) in I[v]
        for (const auto& edge : I[v]) {
            int u = edge.first; // source node
            int v = edge.second; // destination node (same as outer v)

            auto weightIt = weightMap.find({u, v});
            if (weightIt == weightMap.end()) {
                continue; // If no weight, skip
            }

            float weight = weightIt->second;

            // Step 8: Check if shorter path found
            if (delta[v] > delta[u] + weight) {
                // Step 9: Add v into Aff
                Aff.push_back(v);

                // Step 10: Update distance
                delta[v] = delta[u] + weight;

                // Step 11: Update parent
                parent[v] = u;

                // Step 12: Mark as updated
                marked[v] = 1;
            }
        }
    }
}

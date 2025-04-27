//
// Created by jawad on 26/04/2025.
//

#include "PropagateUpdate.h"
#include <vector>
#include <map>
#include <iostream>
#include <limits>

using namespace std;

void propagateUpdate(
    const vector<vector<pair<int, int>>>& adjacencyList, // I[v] (incoming edges)
    vector<float>& delta,         // Distance array
    vector<int>& parent,          // Parent array
    vector<int>& marked,          // Marked nodes
    vector<int>& Aff,             // Current affected nodes
    const map<pair<int, int>, float>& weightMap // Edge weight lookup
) {
    while (!Aff.empty()) { // Step 14: While Aff is not empty
        vector<int> N, AffPrime; // Step 15: Initialize empty N and Aff'

        // Step 16-17: Find neighbors of affected nodes
        for (int v : Aff) {
            for (const auto& edge : adjacencyList[v]) {
                int u = edge.first; // Predecessor node
                N.push_back(u);  // Add to neighbors list
            }
        }

        // Step 18-26: Process all neighbors in parallel
        for (int v : N) {
            for (const auto& edge : adjacencyList[v]) {
                int u = edge.first; // Predecessor neighbor

                if (marked[u] != 1) continue; // Step 20-21: Skip if not marked

                auto weightIt = weightMap.find({u, v});
                if (weightIt == weightMap.end()) continue; // If no weight, skip

                float weight = weightIt->second;

                // Step 22: Relax edge (u,v) if a shorter path is found
                if (delta[v] > delta[u] + weight) {
                    AffPrime.push_back(v); // Step 23: Add v to Aff'
                    delta[v] = delta[u] + weight; // Step 24: Update distance
                    parent[v] = u;  // Step 25: Update parent
                    marked[v] = 1;  // Step 26: Mark as updated
                }
            }
        }

        Aff = AffPrime; // Step 27: Update Aff with Aff'
        AffPrime.clear(); // Step 28: Reset Aff'
    }
}

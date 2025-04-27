#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility> // for pair
#include <limits>  // for infinity
#include <string>
#include <unordered_set>
#include "preprocessing.h"
#include "File_Reading.h"
#include "ProcessChangedEdges.h"

using namespace std;

// Function to build adjacency lists
void buildAdjacencyLists(
    int numVertices,
    const vector<pair<int, int>>& edges,
    vector<vector<int>>& adjList,
    vector<vector<int>>& reverseAdjList
) {
    adjList.assign(numVertices, {});
    reverseAdjList.assign(numVertices, {});

    for (const auto& edge : edges) {
        int u = edge.first;
        int v = edge.second;
        adjList[u].push_back(v);
        reverseAdjList[v].push_back(u);
    }
}

// Function to propagate the updates (Step 2)
void propagateUpdates(
    vector<int>& Aff,
    vector<float>& delta,
    vector<int>& parent,
    vector<int>& marked,
    const map<pair<int, int>, float>& weightMap,
    const vector<vector<int>>& adjList,
    const vector<vector<int>>& reverseAdjList
) {
    while (!Aff.empty()) {
        vector<int> AffPrime;
        unordered_set<int> neighborsSet;

        // Step 16-17: Gather neighbors of affected nodes
        for (int v : Aff) {
            for (int neighbor : adjList[v]) {
                neighborsSet.insert(neighbor);
            }
        }

        // Convert set to vector
        vector<int> N(neighborsSet.begin(), neighborsSet.end());

        // Step 18-26: Check updates for neighbors
        for (int v : N) {
            for (int u : reverseAdjList[v]) { // Predecessor neighbors
                if (marked[u] != 1)
                    continue;

                auto weightIt = weightMap.find({u, v});
                if (weightIt == weightMap.end())
                    continue;

                float weight = weightIt->second;

                if (delta[v] > delta[u] + weight) {
                    delta[v] = delta[u] + weight;
                    parent[v] = u;

                    if (!marked[v]) {
                        marked[v] = 1;
                        AffPrime.push_back(v);
                    }
                }
            }
        }

        // Step 27-28: Update Aff for next round
        Aff = AffPrime;
    }
}

int main() {
    string filename = "edges.txt"; // Your file with src,dest, etc.

    vector<pair<int, int>> insertedEdges;
    map<pair<int, int>, float> weightMap;

    file_reading(filename, insertedEdges, weightMap);

    int numVertices = 5; // (adjust this based on your actual graph)

    vector<vector<pair<int, int>>> I;
    preprocessInsertedEdges(insertedEdges, numVertices, I);

    // Display the I array
    cout << "--- Inserted Edges ---" << endl;
    for (int v = 0; v < numVertices; ++v) {
        cout << "I[" << v << "]: ";
        for (const auto& e : I[v]) {
            cout << "(" << e.first << " -> " << e.second << ") ";
        }
        cout << endl;
    }

    // Initialize
    vector<float> delta(numVertices, numeric_limits<float>::infinity());
    delta[1] = 0.0;  // Assume node 1 as source

    vector<int> parent(numVertices, -1);
    vector<int> marked;
    vector<int> Aff;

    processChangedEdges(I, delta, parent, marked, Aff, weightMap);

    cout << "\n--- After ProcessChangedEdges ---" << endl;
    cout << "Affected Nodes: ";
    for (int node : Aff) cout << node << " ";
    cout << endl;

    cout << "Distances: ";
    for (float d : delta) cout << d << " ";
    cout << endl;

    cout << "Parents: ";
    for (int p : parent) cout << p << " ";
    cout << endl;

    // Build adjacency lists for propagation
    vector<vector<int>> adjList, reverseAdjList;
    buildAdjacencyLists(numVertices, insertedEdges, adjList, reverseAdjList);

    // Step 2: Propagate the update
    propagateUpdates(Aff, delta, parent, marked, weightMap, adjList, reverseAdjList);

    // Final output after propagation
    cout << "\n--- After PropagateUpdates ---" << endl;
    cout << "Distances: ";
    for (float d : delta) cout << d << " ";
    cout << endl;

    cout << "Parents: ";
    for (int p : parent) cout << p << " ";
    cout << endl;

    return 0;
}

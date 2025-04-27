//
// Created by jawad on 26/04/2025.
//

#ifndef PROPAGATEUPDATE_H
#define PROPAGATEUPDATE_H
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
);
#endif //PROPAGATEUPDATE_H

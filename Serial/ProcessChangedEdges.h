//
// Created by jawad on 26/04/2025.
//

#ifndef PROCESSCHANGEDEDGES_H
#define PROCESSCHANGEDEDGES_H
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
);
#endif //PROCESSCHANGEDEDGES_H
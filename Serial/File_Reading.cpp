//
// Created by jawad on 25/04/2025.
//

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <string>

using namespace std;

void file_reading(const string& filename, vector<pair<int, int>>& insertedEdges, map<pair<int, int>, float>& weightMap) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;
    bool firstLine = true;

    while (getline(file, line)) {
        if (firstLine) {
            firstLine = false; // Skip header line
            continue;
        }

        stringstream ss(line);
        string token;
        vector<string> tokens;

        while (getline(ss, token, ',')) {
            tokens.push_back(token);
        }

        if (tokens.size() < 3) {  // We need at least src, dest, and distance
            cerr << "Invalid line: " << line << endl;
            continue;
        }

        int src = stoi(tokens[0]);
        int dest = stoi(tokens[1]);
        float distance = stof(tokens[2]);  // Read distance as weight

        insertedEdges.push_back({src, dest});
        weightMap[{src, dest}] = distance;  // Save weight into map
    }

    file.close();
}

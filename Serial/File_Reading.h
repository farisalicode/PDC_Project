//
// Created by jawad on 25/04/2025.
//

#ifndef FILE_READING_H
#define FILE_READING_H
#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
void file_reading(const string& filename, vector<pair<int, int>>& insertedEdges, map<pair<int, int>, float>& weightMap);

#endif //FILE_READING_H

#ifndef DATAFUNCTIONS_H
#define DATAFUNCTIONS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "Graph.h"

using namespace std;

// Function to load the graph data from a file
// Returns true if the graph was loaded successfully, false otherwise
bool loadGraphFromFile(const string& filename, Graph& G, int a);

#endif // DATAFUNCTIONS_H
// In both dataFunctions.h and sequentialAlgorithm.h
#ifndef GRAPH_H
#define GRAPH_H

#include <vector>

struct Graph { // Graph structure
    int n, a; // Number of vertices and number of vertices to be assigned to X
    std::vector<std::vector<int>> adjMatrix; // Adjacency matrix
};

#endif // GRAPH_H
#ifndef MWCSOLVER_H
#define MWCSOLVER_H

#include <vector>
#include <limits>
#include <iostream>
#include "Graph.h"

using namespace std;

struct state { // State of the search
    vector<bool> nodes; // Nodes assigned to X (0) or Y (1) or not assigned (false)
    int weight; // Weight of the cut between X and Y
    int amountNodes; // Number of nodes assigned to X
    int depth; // Depth of the search tree
    state(vector<bool> nodes, int weight, int amountNodes, int depth) : nodes(nodes), weight(weight), amountNodes(amountNodes), depth(depth) {}
};

class MWCSolver {
    struct solution {
        vector<bool> nodes; // Best partition of the graph
        int weight; // Weight of the best cut
        int recursionCalls; // Number of recursion calls
    } bestSolution;

    public:
        MWCSolver(const Graph& graph, int numThreads); // Constructor for the MWCSolver class

        void solve();  // Solve the Minimum Weighted Cut problem

    private:
        const Graph& G;  // Reference to the graph object
        int numThreads;  // Number of threads to use
        int maxDepth; // max depth
        vector<state> states;

        void dfs(state currentState);  // Depth-first search for DFS
        void parallelDFS();
        void dfsAlmostSeq(state currentState);  // Depth-first search for DFS with almost sequential execution
        int computeLowerBound(const vector<bool> &nodes, int size);  // Compute lower bound
        int getEdgeCutWeight(const vector<bool> &nodes, int index, bool value);  // Get edge cut weight for a given vector of nodes
        void printSolution();  // Print the solution (best partition and cut weight)
};

#endif // MWCSOLVER_H
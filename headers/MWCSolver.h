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

struct MPIState {
    bool nodes[100];
    int size;
    int weight;
    int amountNodes;
    int depth;
};

class MWCSolver {
    struct solution {
        vector<bool> nodes; // Best partition of the graph
        int weight; // Weight of the best cut
        int recursionCalls; // Number of recursion calls
    } bestSolution;

    struct MPIsolution {
        bool nodes[100];
        int size;
        int weight;
        int recursionCalls;
    };

    public:
        MWCSolver(const Graph& graph, int numThreads); // Constructor for the MWCSolver class

        void masterSolve(int numProcesses); // Master process solve
        void slaveSolve(int rank); // Slave process solve

    private:
        const Graph& G;  // Reference to the graph object
        int numThreads;  // Number of threads to use
        int maxDepth; // max depth
        vector<state> states;

        MPIState normalToMPIState(state s); // Convert state to MPIState
        state MPItoNormalState(MPIState s); // Convert MPIState to state
        MPIsolution normalToMPIsolution(solution s); // Convert solution to MPIsolution
        solution MPItoNormalsolution(MPIsolution s); // Convert MPIsolution to solution

        void dfs(state currentState);  // Depth-first search for DFS
        void dfsAlmostSeq(state currentState);  // Depth-first search for DFS with almost sequential execution
        int computeLowerBound(const vector<bool> &nodes, int size);  // Compute lower bound
        int getEdgeCutWeight(const vector<bool> &nodes, int index, bool value);  // Get edge cut weight for a given vector of nodes
        void printSolution();  // Print the solution (best partition and cut weight)
};

#endif // MWCSOLVER_H
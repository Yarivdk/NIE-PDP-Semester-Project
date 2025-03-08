#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include "Graph.h"
#include "MWCSolver.h"  // Include the header file

using namespace std;


// Constructor and method implementations
MWCSolver::MWCSolver(const Graph& graph, int numThreads) : G(graph), numThreads(numThreads) { // Constructor for the MWCSolver class
    bestSolution.nodes.resize(G.n, 0); // Initialize the best partition vector
    bestSolution.weight = numeric_limits<int>::max(); // Initialize the best cut weight to maximum value
    bestSolution.recursionCalls = 0; // Initialize the number of recursion calls

    if (G.n > 9) { // If the graph has more than 4 vertices set the max depth to 4
        maxDepth = 9;
    }
    else {
        maxDepth = G.n-1;
    }
}

void MWCSolver::solve() {
    state s = state(vector<bool>(G.n, false), 0, 0, 0); // Initialize the initial state
    dfs(s); // Start the search
    if (states.size() > 0) {
        parallelDFS(); // Start the parallel search
    }
    printSolution(); // Print the solution
}

void MWCSolver::dfs(state currentState) {
    // Assign vertex to Y (0) and search next state
    currentState.nodes[currentState.depth] = false; // Assign vertex to Y
    int newWeight = currentState.weight + getEdgeCutWeight(currentState.nodes, currentState.depth, false); // Update the cut weight
    state leftState = state(currentState.nodes, newWeight, currentState.amountNodes, currentState.depth + 1); // Create the next state

    if ((leftState.amountNodes <= G.a) && (G.a - leftState.amountNodes) <= (G.n - leftState.depth)) { // if the amount of nodes assigned to X is less than a and the remaining nodes can be assigned to X
        int lb = newWeight + computeLowerBound(leftState.nodes, leftState.depth); // Compute the lower bound
        if (newWeight < bestSolution.weight && lb < bestSolution.weight) { // If the lower bound is better than the best we have so far
            if (currentState.depth < maxDepth) { // If not all vertices are assigned
                dfs(leftState); // Continue the search
            } else {
                states.push_back(leftState); // Continue the search
            }
        }
    }
    
    
    // Assign vertex to X (1) and update cut weight
    currentState.nodes[currentState.depth] = true; // Assign vertex to X
    newWeight = currentState.weight + getEdgeCutWeight(currentState.nodes, currentState.depth, true); // Update the cut weight
    state rightState = state(currentState.nodes, newWeight, currentState.amountNodes + 1, currentState.depth + 1); // Create the next state

    if ((rightState.amountNodes <= G.a) && (G.a - rightState.amountNodes) <= (G.n - rightState.depth))  { // if the amount of nodes assigned to X is less than a and the remaining nodes can be assigned to X
        int lb = newWeight + computeLowerBound(rightState.nodes, rightState.depth); // Compute the lower bound
        if (newWeight < bestSolution.weight && lb < bestSolution.weight) { // If the lower bound is better than the best we have so far
            if (currentState.depth < maxDepth) { // If not all vertices are assigned
                dfs(rightState); // Continue the search
            } else {
                states.push_back(rightState); // Continue the search
            }
        }
    }

}

void MWCSolver::parallelDFS() {
    #pragma omp parallel for schedule(dynamic) num_threads(numThreads) // Start parallel region
    for (size_t i = 0; i < states.size(); i++) {
        dfsAlmostSeq(states[i]); // Continue the search
    }
}

void MWCSolver::dfsAlmostSeq(state currentState) {
    #pragma omp atomic // Start atomic region (ensures that the operation is executed atomically)
    bestSolution.recursionCalls++; // Increment the number of recursion calls
    #pragma omp critical // Start critical region (only one thread can execute this block at a time)
    { 
        if (currentState.depth == G.n) { // If all vertices are assigned
            if (currentState.amountNodes == G.a && currentState.weight < bestSolution.weight) { // If the partition is valid and the cut weight is better than the best we have so far 
                bestSolution.nodes = currentState.nodes; // Update the best partition
                bestSolution.weight = currentState.weight; // Update the best cut weight
            }
            return;
        }
    }
    
    // cout << "Depth: " << currentState.depth << " AmountNodes: " << currentState.amountNodes << " Weight: " << currentState.weight << endl;
    
    // Assign vertex to Y (0) and search next state
    currentState.nodes[currentState.depth] = false; // Assign vertex to Y
    int newWeight = currentState.weight + getEdgeCutWeight(currentState.nodes, currentState.depth, false); // Update the cut weight
    state leftState = state(currentState.nodes, newWeight, currentState.amountNodes, currentState.depth + 1); // Create the next state

    if ((leftState.amountNodes <= G.a) && (G.a - leftState.amountNodes) <= (G.n - leftState.depth)) { // if the amount of nodes assigned to X is less than a and the remaining nodes can be assigned to X
        int lb = newWeight + computeLowerBound(leftState.nodes, leftState.depth); // Compute the lower bound
        if (newWeight < bestSolution.weight && lb < bestSolution.weight) { // If the lower bound is better than the best we have so far
            dfsAlmostSeq(leftState); // Continue the search
        }
    }
    
    
    // Assign vertex to X (1) and update cut weight
    currentState.nodes[currentState.depth] = true; // Assign vertex to X
    newWeight = currentState.weight + getEdgeCutWeight(currentState.nodes, currentState.depth, true); // Update the cut weight
    state rightState = state(currentState.nodes, newWeight, currentState.amountNodes + 1, currentState.depth + 1); // Create the next state

    if ((rightState.amountNodes <= G.a) && (G.a - rightState.amountNodes) <= (G.n - rightState.depth))  { // if the amount of nodes assigned to X is less than a and the remaining nodes can be assigned to X
        int lb = newWeight + computeLowerBound(rightState.nodes, rightState.depth); // Compute the lower bound
        if (newWeight < bestSolution.weight && lb < bestSolution.weight) { // If the lower bound is better than the best we have so far
            dfsAlmostSeq(rightState); // Continue the search
        }
    }
}

int MWCSolver::getEdgeCutWeight(const vector<bool> &nodes, int index, bool value) { // Get edge cut weight for a given vector of nodes
    int weight = 0; // Initialize the weight
    for (int u = 0; u < index; ++u) {  // only consider the vertices before index
        weight += (nodes[u] != value) * G.adjMatrix[u][index]; // if vertex is not in partition, then add the weight
    }
    return weight;
}


int MWCSolver::computeLowerBound(const vector<bool> &nodes, int size) { // Compute lower bound
    int bound = 0; // Initialize the bound 
    for (int v = size; v < G.n; v++) { // for each vertex after size
        int cutWeightFalse = 0, cutWeightTrue = 0; // Initialize the cut weights
        for (int u = 0; u < size; ++u) { // for each vertex before size
            int w = G.adjMatrix[u][v]; // get the weight of the edge
            cutWeightFalse += (nodes[u]) ? w : 0; // if vertex is in partition, then add the weight
            cutWeightTrue  += (!nodes[u]) ? w : 0; // if vertex is not in partition, then add the weight
        }
        bound += std::min(cutWeightFalse, cutWeightTrue); // add the minimum of the two cuts
    }
    return bound;
}


void MWCSolver::printSolution() { // Print the solution (best partition, cut weight and the amount of recursion calls)
    cout << "Minimum Edge Cut: " << bestSolution.weight << endl;
    cout << "Partition X: ";
    for (int i = 0; i < G.n; ++i) {
        if (bestSolution.nodes[i] == 0) cout << i << " ";
    }
    cout << "\nPartition Y: ";
    for (int i = 0; i < G.n; ++i) {
        if (bestSolution.nodes[i] == 1) cout << i << " ";
    }
    cout << endl;
    cout << "Recursion calls: " << bestSolution.recursionCalls << endl;
}
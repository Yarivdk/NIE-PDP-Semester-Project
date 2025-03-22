#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include "Graph.h"
#include "MWCSolver.h"  // Include the header file
#include <mpi.h>
#include <omp.h>

using namespace std;


// Constructor and method implementations
MWCSolver::MWCSolver(const Graph& graph, int numThreads) : G(graph), numThreads(numThreads) { // Constructor for the MWCSolver class
    bestSolution.nodes.resize(G.n, 0); // Initialize the best partition vector
    bestSolution.weight = numeric_limits<int>::max(); // Initialize the best cut weight to maximum value
    bestSolution.recursionCalls = 0; // Initialize the number of recursion calls

    if (G.n > 9) { // If the graph has more than 9 vertices set the max depth to 9
        maxDepth = 9;
    }
    else {
        maxDepth = G.n-1;
    }
}


MPIState MWCSolver::normalToMPIState(state s) { // Convert state to MPIState
    MPIState mpiState;
    for (int i = 0; i < s.nodes.size(); i++) { // Convert from vector to bool array, because MPI does not support vectors
        mpiState.nodes[i] = s.nodes[i];
    }
    mpiState.size = s.nodes.size();
    mpiState.weight = s.weight;
    mpiState.amountNodes = s.amountNodes;
    mpiState.depth = s.depth;

    return mpiState;
}

state MWCSolver::MPItoNormalState(MPIState s) { // Convert MPIState to state
    vector<bool> nodes;
    for (int i = 0; i < s.size; i++) {
        nodes.push_back(s.nodes[i]);
    }
    state normalState = state(nodes, s.weight, s.amountNodes, s.depth);

    return normalState;
}

MWCSolver::MPIsolution MWCSolver::normalToMPIsolution(solution s) { // Convert solution to MPIsolution
    MPIsolution mpiSolution;
    for (int i = 0; i < s.nodes.size(); i++) { // Convert from vector to bool array, because MPI does not support vectors
        mpiSolution.nodes[i] = s.nodes[i];
    }
    mpiSolution.size = s.nodes.size();
    mpiSolution.weight = s.weight;
    mpiSolution.recursionCalls = s.recursionCalls;

    return mpiSolution;
}

MWCSolver::solution MWCSolver::MPItoNormalsolution(MPIsolution s) { // Convert MPIsolution to solution
    vector<bool> nodes;
    for (int i = 0; i < s.size; i++) {
        nodes.push_back(s.nodes[i]);
    }
    solution normalSolution = {nodes, s.weight, s.recursionCalls};

    return normalSolution;
}

void MWCSolver::masterSolve(int numProcesses) {
    state s = state(vector<bool>(G.n, false), 0, 0, 0); // Initialize the initial state

    dfs(s) ; // Run the DFS algorithm

    int i = 0;
    int amountItems = states.size();

    for (int rank = 1; rank < numProcesses; rank++) { // Send states to slaves
        if (i < amountItems) {
            MPIState mpiState = normalToMPIState(states[i]);
            MPI_Send(&mpiState, sizeof(MPIState), MPI_BYTE, rank, 1, MPI_COMM_WORLD);
            i++;
        }
    }

    MPI_Status status;
    int working = numProcesses-1;

    while (working > 0) { // Receive results from slaves
        MPIsolution mpiSolution;
        MPI_Recv(&mpiSolution, sizeof(MPIsolution), MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

        solution localSolution = MPItoNormalsolution(mpiSolution);
        if (localSolution.weight < bestSolution.weight) { // Update the best solution
            bestSolution.weight = localSolution.weight;
            bestSolution.nodes = localSolution.nodes;
            bestSolution.recursionCalls += localSolution.recursionCalls; // Update the number of recursion calls
        }   

        if (i < amountItems) { // Send new state to the slave
            MPIState mpiState = normalToMPIState(states[i]);
            MPI_Send(&mpiState, sizeof(MPIState), MPI_BYTE, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
            i++;
        } else {
            MPI_Send(NULL, 0, MPI_BYTE, status.MPI_SOURCE, 2, MPI_COMM_WORLD); // Send termination signal
            working--;
        }
    }

    printSolution(); // Print the solution
}

void MWCSolver::slaveSolve(int rank) {
    MPI_Status status;
    while (true) {
        MPIState mpiState;
        MPI_Recv(&mpiState, sizeof(MPIState), MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // Receive state from master

        if (status.MPI_TAG == 2) { // If termination signal is received
            break;
        }

        state s = MPItoNormalState(mpiState);

        #pragma cmp parallel num_threads(numThreads) // start parallel region to run DFS 
        {
            int numThreadsInRegion;
            #pragma omp single
            {
                numThreadsInRegion = omp_get_num_threads();
                cout << "Number of threads in the region: " << numThreadsInRegion << endl;
                dfsAlmostSeq(s); // Run the DFS algorithm
            }
            
        }

        MPIsolution mpiSolution = normalToMPIsolution(bestSolution);
        MPI_Send(&mpiSolution, sizeof(MPIsolution), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }
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

void MWCSolver::dfsAlmostSeq(state currentState) { 
    #pragma omp atomic // Start atomic region (ensures that the operation is executed atomically)
    bestSolution.recursionCalls++; // Increment the number of recursion calls
    if (currentState.depth == G.n) { // If all vertices are assigned
        #pragma omp critical // Start critical region (only one thread can execute this block at a time)
        {
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
            if (currentState.depth < 0.75*G.n) { // If not all vertices are assigned
                #pragma omp task // Start a task region
                dfsAlmostSeq(leftState); // Continue the search
            } else {
                dfsAlmostSeq(leftState); // Continue the search
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
            if (currentState.depth < 0.75*G.n) { // If not all vertices are assigned
                #pragma omp task // Start a task region
                dfsAlmostSeq(rightState); // Continue the search
            } else {
                dfsAlmostSeq(rightState); // Continue the search
            }
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
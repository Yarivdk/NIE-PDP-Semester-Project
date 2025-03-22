#include <iostream>
#include <vector>
#include <chrono>
#include "dataFunctions.h"
#include "MWCSolver.h"
#include "Graph.h"
#include <omp.h>
#include <mpi.h>    

// USAGE: ./output.exe <filename> <a>

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <filename> <a> <number of threads>" << std::endl;
        return 1;
    }

    std::string filename = "input/" + std::string(argv[1]);
    int a = std::stoi(argv[2]);  // Convert string argument to integer
    int numThreads = std::stoi(argv[3]);  // Convert string argument to integer
    Graph G;

    if (loadGraphFromFile(filename, G, a)) {
        // Successfully loaded the graph
        cout << "Graph loaded successfully!" << endl;
        cout << "Number of vertices: " << G.n << endl;
    } else {
        cerr << "Failed to load the graph." << endl;
    }


    MWCSolver solver(G, numThreads);  // Create a solver object

    int rank, size;
    MPI_Init(NULL, NULL); // Initialize MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of the process
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the total number of processes

    if (rank == 0) {
        auto start = std::chrono::high_resolution_clock::now(); // Start measuring time

        solver.masterSolve(size);  // Solve the problem

        auto stop = std::chrono::high_resolution_clock::now(); // Stop measuring time

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); // Calculate duration

        std::cout << "Execution Time: " << duration.count() << " ms" << std::endl;
    }
    else {
        solver.slaveSolve(rank);  // Solve the problem
    }

    MPI_Finalize(); // Finalize MPI
    
    return 0;
}
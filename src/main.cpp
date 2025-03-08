#include <iostream>
#include <vector>
#include <chrono>
#include "dataFunctions.h"
#include "MWCSolver.h"
#include "Graph.h"
#include <omp.h>

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
    
    auto start = std::chrono::high_resolution_clock::now(); // Start measuring time

    solver.solve();  // Solve the problem

    auto stop = std::chrono::high_resolution_clock::now(); // Stop measuring time

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); // Calculate duration

    std::cout << "Execution Time: " << duration.count() << " ms" << std::endl;
    
    return 0;
}
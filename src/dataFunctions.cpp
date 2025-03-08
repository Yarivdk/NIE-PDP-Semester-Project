#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "dataFunctions.h"
#include "Graph.h"

using namespace std;

// Function to load graph data from a file
bool loadGraphFromFile(const string& filename, Graph& G, int a) { // Load the graph data from a file
    ifstream file(filename);  // Open the file

    if (!file.is_open()) { // Check if the file was opened successfully
        cerr << "Error opening file!" << endl;
        return false;  // Return false if the file could not be opened
    }

    string line; // String to store each line of the file
    int row = 0;

    if (getline(file, line)) { // Read the first line of the file
        istringstream(line) >> G.n;  // Read the number of vertices (n)
        G.a = a;  // Set the value of a
        G.adjMatrix.resize(G.n, vector<int>(G.n, 0));  // Resize adjacency matrix based on the number of vertices
    } else {
        cerr << "Error reading the number of vertices!" << endl;
        return false;
    }

    while (getline(file, line) && row < G.n) { // Read the remaining lines of the file
        istringstream ss(line); // Create a string stream from the line
        int value, col = 0; 

        while (ss >> value && col < G.n) { // Read the values for the current row of the adjacency matrix
            G.adjMatrix[row][col] = value; // Store the value in the adjacency matrix
            col++; // Move to the next column
        }
        row++; // Move to the next row
    }

    if (row != G.n) { // Check if the number of rows in the adjacency matrix matches the number of vertices
        cerr << "Error: Number of rows in the adjacency matrix does not match the number of vertices!" << endl;
        return false;
    }

    file.close();  // Close the file
    return true;  // Return true if the file was loaded successfully
}


#ifndef GRAPH
#define GRAPH

#include <iostream>
#include <list>
#include <stack>
#include <queue>
#include <algorithm>
#include <numeric> 
#include <tuple> 

// === Disjoint Set Data Structure ===
// This class is used for union-find operations that are necessary for Kruskal's algorithm
class DisjointSet {
public:
    // Constructor initializes the parent and rank vectors
    DisjointSet(int n) {
        parent.resize(n);          // Resizing the parent vector to hold n elements
        rank.resize(n, 0);         // Initializing all ranks to 0
        for (int i = 0; i < n; ++i) {
            parent[i] = i;         // Initially, every vertex is its own parent
        }
    }

    // --- Find operation with Path Compression ---
    // The function returns the root of the set to which u belongs, 
    // compressing the path for efficiency.
    int find(int u) {
        if (parent[u] != u) {
            parent[u] = find(parent[u]);  // Path compression
        }
        return parent[u];
    }

    // --- Union operation by Rank ---
    // This function unites two sets (u and v) by comparing their ranks
    void unite(int u, int v) {
        int rootU = find(u);
        int rootV = find(v);
        
        if (rootU != rootV) { // If they belong to different sets
            if (rank[rootU] > rank[rootV]) {
                parent[rootV] = rootU;
            } else if (rank[rootU] < rank[rootV]) {
                parent[rootU] = rootV;
            } else { 
                parent[rootV] = rootU;
                ++rank[rootU];  // Increment the rank when both have the same rank
            }
        }
    }

private:
    std::vector<int> parent, rank;  // Parent and rank vectors
};

// === Graph Class ===
// This class represents a Graph using an adjacency list
class Graph {
public:
    // Constructor initializes the graph with the specified number of vertices
    Graph(int);

    // --- Add an edge between two vertices ---
    void addEdge(int, int);
    
    // --- Remove an edge between two vertices ---
    void removeEdge(int, int);
    
    // --- Check if an edge exists between two vertices ---
    bool hasEdge(int, int) const;
    
    // --- Get the neighbors of a vertex ---
    std::vector<int> getNeighbors(int) const;
    
    // --- Print the graph in adjacency list format ---
    void print() const;

    // --- Get the number of vertices in the graph ---
    int getVertexCount() const;

    // --- Get the degree (number of neighbors) of a vertex ---
    int getDegree(int) const;

    // --- Depth First Search (Recursive version) ---
    void depthFirstSearchRecursive(int);

    // --- Depth First Search (Iterative version using stack) ---
    void depthFirstSearchIterative(int) const;

    // --- Breadth First Search ---
    void breadthFirstSearch(int);
    
    // --- Find the number of connected components in the graph ---
    int connectedComponentsCount();
    
    // --- Count the number of vertices at a given level ---
    int vertexCountAtLevel(int, int);
    
    // --- Find the shortest path between two vertices ---
    std::vector<int> shortestPath(int, int);

    // --- Count the number of distinct paths between two vertices ---
    int pathCount(int, int);
    
    // --- Find all paths between two vertices ---
    std::vector<std::vector<int>> allPaths(int, int);
    
    // --- Check if the graph has a cycle ---
    bool isCyclic();

    // --- Kruskal's Minimum Spanning Tree algorithm ---
    void KruskalMST() const;
    
    // --- Prim's Minimum Spanning Tree algorithm ---
    void PrimMST() const;

private:
    int countOfVertices;            // Number of vertices in the graph
    std::vector<std::list<int>> adjList;  // Adjacency list representing the graph

    // --- Helper functions for DFS ---
    void dfsRecursiveHelper(int, std::vector<bool>&);

    // --- Helper function for counting paths ---
    void pathCountHelper(int, int, int&, std::vector<bool>&);
    
    // --- Helper function for finding all paths ---
    void allPathsHelper(int, int, std::vector<std::vector<int>>&, std::vector<int>, std::vector<bool>&);
    
    // --- Helper function for cycle detection ---
    bool cycleDetectionHelper(std::vector<bool>&, int, int);
    
    // --- Validate if a vertex is within bounds ---
    void validateVertex(int) const;
};

#endif
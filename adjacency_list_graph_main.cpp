#include <iostream>
#include "adjacency_list_graph.h"

int main() {
    // Create a graph with 6 vertices
    Graph graph(6);

    // Add edges
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);
    graph.addEdge(1, 3);
    graph.addEdge(1, 4);
    graph.addEdge(2, 5);
    graph.addEdge(3, 4);

    // Print the graph's adjacency list
    std::cout << "Graph Adjacency List:\n";
    graph.print();

    // Check if there is an edge between two vertices
    std::cout << "Has edge between 0 and 1? " << (graph.hasEdge(0, 1) ? "Yes" : "No") << std::endl;

    // Get neighbors of a vertex
    std::vector<int> neighbors = graph.getNeighbors(1);
    std::cout << "Neighbors of vertex 1: ";
    for (int neighbor : neighbors) {
        std::cout << neighbor << " ";
    }
    std::cout << std::endl;

    // Get the degree of a vertex
    std::cout << "Degree of vertex 1: " << graph.getDegree(1) << std::endl;

    // Perform DFS (Recursive)
    std::cout << "DFS (Recursive) starting from vertex 0: ";
    graph.depthFirstSearchRecursive(0);

    // Perform DFS (Iterative)
    std::cout << "DFS (Iterative) starting from vertex 0: ";
    graph.depthFirstSearchIterative(0);

    // Perform BFS
    std::cout << "BFS starting from vertex 0: ";
    graph.breadthFirstSearch(0);

    // Find connected components count
    std::cout << "Number of connected components: " << graph.connectedComponentsCount() << std::endl;

    // Count vertices at a certain level in BFS
    std::cout << "Vertices at level 2 from vertex 0: " << graph.vertexCountAtLevel(0, 2) << std::endl;

    // Find shortest path between two vertices
    std::vector<int> shortestPath = graph.shortestPath(0, 4);
    std::cout << "Shortest path from vertex 0 to vertex 4: ";
    for (int vertex : shortestPath) {
        std::cout << vertex << " ";
    }
    std::cout << std::endl;

    // Count all paths between two vertices
    std::cout << "Number of paths between vertex 0 and vertex 4: " << graph.pathCount(0, 4) << std::endl;

    // Get all paths between two vertices
    std::vector<std::vector<int>> allPaths = graph.allPaths(0, 4);
    std::cout << "All paths from vertex 0 to vertex 4:\n";
    for (const auto& path : allPaths) {
        for (int vertex : path) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }

    // Check if the graph is cyclic
    std::cout << "Is the graph cyclic? " << (graph.isCyclic() ? "Yes" : "No") << std::endl;

    // Kruskal's MST
    std::cout << "Kruskal's MST:\n";
    graph.KruskalMST();

    // Prim's MST
    std::cout << "Prim's MST:\n";
    graph.PrimMST();

    return 0;
}

#include "adjacency_list_graph.h"

// === Constructor: Initialize the graph ===
// The constructor takes the number of vertices and initializes the adjacency list
Graph::Graph(int countOfVertices) : countOfVertices(countOfVertices), adjList(countOfVertices) { 
    if (countOfVertices <= 0) {
        throw std::invalid_argument("Number of vertices must be positive.");
    }
}

// === Add Edge ===
// This function adds an edge between two vertices (undirected graph)
void Graph::addEdge(int source, int destination) {
    try {
        validateVertex(source);  // Validate source vertex
        validateVertex(destination);  // Validate destination vertex
        
        adjList[source].push_back(destination);  // Add destination to the adjacency list of source
        adjList[destination].push_back(source);  // Add source to the adjacency list of destination
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Error in addEdge: " << e.what() << std::endl;
    }
}

// === Remove Edge ===
// This function removes the edge between two vertices (undirected graph)
void Graph::removeEdge(int source, int destination) {
    try {
        validateVertex(source);  // Validate source vertex
        validateVertex(destination);  // Validate destination vertex
        
        adjList[source].remove(destination);  // Remove destination from the adjacency list of source
        adjList[destination].remove(source);  // Remove source from the adjacency list of destination
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Error in removeEdge: " << e.what() << std::endl;
    }
}

// === Has Edge? ===
// This function checks if an edge exists between two vertices
bool Graph::hasEdge(int source, int destination) const {
    try {
        validateVertex(source);  // Validate source vertex
        validateVertex(destination);  // Validate destination vertex
        
        for (int neighbor : adjList[source]) {
            if (neighbor == destination) {
                return true;  // Edge exists
            }
        }
    } 
    catch (const std::out_of_range&) {
        // In case the vertex is invalid, we return false.
    }
    return false;  // Edge doesn't exist
}

// === Get Neighbors ===
// This function returns a list of neighbors for a given vertex
std::vector<int> Graph::getNeighbors(int vertex) const {
    try {
        validateVertex(vertex);  // Validate the vertex
        return {adjList[vertex].begin(), adjList[vertex].end()};  // Return the neighbors as a vector
    } catch (const std::out_of_range& e) {
        std::cerr << "Error in getNeighbors: " << e.what() << std::endl;
        return {};  // Return empty vector in case of error
    }
}

// === Print Graph ===
// This function prints the adjacency list of the graph
void Graph::print() const {
    for (int i{}; i < countOfVertices; ++i) {
        std::cout << i << ": ";  // Print the vertex
        for (int neighbor : adjList[i]) {
            std::cout << neighbor << " ";  // Print each neighbor
        }
        std::cout << std::endl;
    }
}

// === Get Vertex Count ===
// This function returns the total number of vertices in the graph
int Graph::getVertexCount() const {
    return countOfVertices;
}

// === Get Degree ===
// This function returns the degree (number of edges) of a vertex
int Graph::getDegree(int vertex) const {
    try {
        validateVertex(vertex);  // Validate the vertex
        return adjList[vertex].size();  // Return the degree of the vertex
    } 
    catch (const std::out_of_range& e) {
        std::cerr << "Error in getDegree: " << e.what() << std::endl;
        return -1;  // Return -1 in case of error
    }
}

// === Depth First Search (Recursive) ===
// This function performs DFS recursively starting from the given vertex
void Graph::depthFirstSearchRecursive(int startVertex) {
    try {
        validateVertex(startVertex);  // Validate the starting vertex
        
        std::vector<bool> visited(countOfVertices, false);  // Track visited vertices
        dfsRecursiveHelper(startVertex, visited);  // Call the helper function for DFS
        std::cout << std::endl;  // Print new line after DFS completion
    } catch (const std::out_of_range& e) {
        std::cerr << "Error in depthFirstSearchRecursive: " << e.what() << std::endl;
    }
}

// === Depth First Search (Iterative) ===
// This function performs DFS iteratively using a stack
void Graph::depthFirstSearchIterative(int startVertex) const {
    try {
        validateVertex(startVertex);  // Validate the starting vertex
        
        std::vector<bool> visited(countOfVertices, false);  // Track visited vertices
        std::stack<int> stack;  // Stack to simulate recursive calls
        
        stack.push(startVertex);  // Push the starting vertex to the stack
        
        while (!stack.empty()) {
            int vertex = stack.top();
            stack.pop();
            
            if (!visited[vertex]) {
                visited[vertex] = true;  // Mark the vertex as visited
                std::cout << vertex << " ";  // Print the vertex
            }
            
            // Push all unvisited neighbors of the current vertex to the stack
            for (int neighbor : adjList[vertex]) {
                if (!visited[neighbor]) {
                    stack.push(neighbor);
                }
            }
        }
        
        std::cout << std::endl;  // Print new line after DFS completion
    } catch (const std::out_of_range& e) {
        std::cerr << "Error in depthFirstSearchIterative: " << e.what() << std::endl;
    }
}

// === Breadth First Search ===
// This function performs BFS starting from the given vertex
void Graph::breadthFirstSearch(int startVertex) {
    try {
        validateVertex(startVertex);  // Validate the starting vertex
        
        std::vector<bool> visited(countOfVertices, false);  // Track visited vertices
        std::queue<int> q;  // Queue to hold vertices for BFS
        
        visited[startVertex] = true;  // Mark the start vertex as visited
        q.push(startVertex);  // Push the start vertex into the queue
        
        while (!q.empty()) {
            int vertex = q.front(); 
            q.pop();  // Remove the front element from the queue
            
            std::cout << vertex << " ";  // Print the vertex
            
            // Add unvisited neighbors to the queue
            for (int neighbor : adjList[vertex]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;  // Mark neighbor as visited
                    q.push(neighbor);  // Add to queue for further exploration
                }
            }
        }
        std::cout << std::endl;  // Print new line after BFS completion
    } catch (const std::out_of_range& e) {
        std::cerr << "Error in breadthFirstSearch: " << e.what() << std::endl;
    }
}

// === Connected Components Count ===
// Function to count the number of connected components in the graph
int Graph::connectedComponentsCount() {
    // Vector to keep track of visited vertices
    std::vector<bool> visited(countOfVertices, false);
    
    // Variable to store the count of connected components
    int componentCount{};
    
    // Iterate through all vertices
    for (int i{}; i < countOfVertices; ++i) {
        // If the vertex hasn't been visited, it means it's the start of a new component
        if (!visited[i]) {
            // Perform DFS starting from this vertex to visit all vertices in this component
            dfsRecursiveHelper(i, visited);
            // Increment the component count after exploring all vertices in the component
            ++componentCount;
        }
    }
    // Return the total number of connected components
    return componentCount;
}

// === Vertex CountAt A Given Level ===
// Function to count the number of vertices at a specific level starting from a given vertex
int Graph::vertexCountAtLevel(int startVertex, int level) {
    try {
        // Validate the starting vertex
        validateVertex(startVertex);

        // Check if the level is non-negative
        if (level < 0) {
            throw std::invalid_argument("Level cannot be negative.");
        }

        // Vector to keep track of visited vertices
        std::vector<bool> visited(countOfVertices, false);
        
        // Queue to perform BFS
        std::queue<int> q;
        
        // Variable to keep track of the current level in BFS
        int currentLevel{};

        // Push the start vertex into the queue and mark it as visited
        q.push(startVertex);
        visited[startVertex] = true;

        // Perform BFS
        while (!q.empty()) {
            int levelSize = q.size(); // Get the number of vertices at the current level

            // If the current level matches the target level, return the count of vertices at this level
            if (currentLevel == level) {
                return levelSize;
            }

            // Process all vertices at the current level
            while (levelSize--) {
                int u = q.front();
                q.pop();

                // Explore the neighbors of the current vertex
                for (int neighbor : adjList[u]) {
                    // If a neighbor has not been visited, mark it and add it to the queue
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        q.push(neighbor);
                    }
                }
            }
            // Increment the level after processing the current level's vertices
            ++currentLevel;
        }

        // If the level does not exist in the graph, return 0
        return 0; 
    } catch (const std::exception& e) {
        // Catch and handle any exceptions (like invalid vertex or negative level)
        std::cerr << "Error in vertexCountAtLevel: " << e.what() << std::endl;
        return -1; // Return -1 to indicate an error
    }
}

// === Shortest Path Source To Destination ===
// Function to find the shortest path between two vertices in the graph
std::vector<int> Graph::shortestPath(int source, int destination) {
    try {
        // Validate the source and destination vertices
        validateVertex(source);
        validateVertex(destination);

        // If the source and destination are the same, return a path containing only the source
        if (source == destination) {
            return {source}; 
        }

        // Vector to track visited vertices
        std::vector<bool> visited(countOfVertices, false);
        
        // Vector to store the parent of each vertex (for path reconstruction)
        std::vector<int> parents(countOfVertices, -1);
        
        // Queue for BFS traversal
        std::queue<int> q;

        // Push the source vertex into the queue and mark it as visited
        q.push(source);
        visited[source] = true;

        // Perform BFS
        while (!q.empty()) {
            int u = q.front();
            q.pop();

            // Explore the neighbors of the current vertex
            for (int neighbor : adjList[u]) {
                // If a neighbor has not been visited, mark it as visited and set its parent
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    parents[neighbor] = u;
                    q.push(neighbor);

                    // If the destination vertex is reached, reconstruct and return the path
                    if (neighbor == destination) {
                        std::vector<int> path;
                        for (int v = destination; v != -1; v = parents[v]) {
                            path.insert(path.begin(), v); // Insert the vertex at the beginning of the path
                        }
                        return path; // Return the reconstructed path
                    }
                }
            }
        }

        // If no path exists from source to destination, return an empty vector
        return {}; 
    } catch (const std::exception& e) {
        // Catch and handle any exceptions (like invalid vertices)
        std::cerr << "Error in shortestPath: " << e.what() << std::endl;
        return {}; // Return an empty vector to indicate failure
    }
}

// === Number Of Paths Form the Source Vertex to the Destination Vertex ===
// Function to count the number of paths from the source vertex to the destination vertex
int Graph::pathCount(int source, int destination) {
    try {
        // Validate that the source and destination vertices exist in the graph
        validateVertex(source);
        validateVertex(destination);

        // Variable to store the number of paths found
        int count{};

        // Vector to track visited vertices during the path search
        std::vector<bool> visited(countOfVertices, false);

        // Call the helper function to recursively find all paths and update the count
        pathCountHelper(source, destination, count, visited);

        // Return the total number of paths found
        return count;
    } 
    catch (const std::exception& e) {
        // Catch and handle any exceptions (e.g., invalid vertex or graph errors)
        std::cerr << "Error in pathCount: " << e.what() << std::endl;
        return -1; // Return -1 to indicate an error
    }
}

// === All Possible Paths Source to Destination ===
// Function to find and return all paths from the source vertex to the destination vertex
std::vector<std::vector<int>> Graph::allPaths(int source, int destination) {
    validateVertex(source);
    validateVertex(destination);

    // Vector to store all paths found from source to destination
    std::vector<std::vector<int>> paths;

    // Vector to store the current path as we explore it
    std::vector<int> path;

    // Vector to track visited vertices during the path search
    std::vector<bool> visited(countOfVertices, false);

    // Call the helper function to recursively find all paths and store them in 'paths'
    allPathsHelper(source, destination, paths, path, visited);

    // Return the list of all paths found
    return paths;
}

// === Check If the Graph contains Any Cycle ===
// Function to check if the graph contains any cycles (i.e., if the graph is cyclic)
bool Graph::isCyclic() {
    try {
        // Vector to track visited vertices during the cycle detection
        std::vector<bool> visited(countOfVertices, false);

        // Iterate through all vertices in the graph
        for (int i{}; i < countOfVertices; ++i) {
            // If the vertex has not been visited, start DFS from it
            if (!visited[i]) {
                // If a cycle is detected starting from the current vertex, return true
                if (cycleDetectionHelper(visited, i, -1)) {
                    return true;
                }
            }
        }

        // Return false if no cycle is detected
        return false;
    } catch (const std::exception& e) {
        // Catch and handle any exceptions (e.g., invalid vertex or graph errors)
        std::cerr << "Error in isCyclic: " << e.what() << std::endl;
        return false;  // Return false in case of an error
    }
}

// === Krsukal's MST Algorithm ===
// Function to find and print the Minimum Spanning Tree (MST) using Kruskal's Algorithm with edge weights
void Graph::KruskalMST() const {
    // Vector to store edges, now with weights (weight, u, v)
    std::vector<std::tuple<int, int, int>> edges;

    // Collect all edges from the adjacency list, adding weight (set to 1 for simplicity)
    for (int u{}; u < countOfVertices; ++u) {
        for (int v : adjList[u]) {
            if (u < v) {  // Avoid duplicate edges
                int weight{1};  // Assign weight to the edge (this can be modified if needed)
                edges.push_back({weight, u, v});  // Add edge to the list with weight
            }
        }
    }

    // Sort the edges by weight
    std::sort(edges.begin(), edges.end()); 

    // Initialize DSU structure for cycle detection
    DisjointSet dsu(countOfVertices);
    std::vector<std::tuple<int, int>> mst;  // Vector to store MST edges

    // Process edges to build the MST, avoiding cycles
    for (const auto& edge : edges) {
        int u{std::get<1>(edge)};
        int v{std::get<2>(edge)};
        
        // If u and v belong to different sets, add the edge to the MST
        if (dsu.find(u) != dsu.find(v)) {
            mst.push_back({u, v});  // Add edge to the MST
            dsu.unite(u, v);  // Union the sets of u and v
        }
    }

    // Print the edges in the MST
    std::cout << "Kruskal's MST Edges: \n";
    for (const auto& edge : mst) {
        std::cout << std::get<0>(edge) << " - " << std::get<1>(edge) << std::endl;
    }
}

// === Prim's MST Algorithm ===
// Function to find and print the Minimum Spanning Tree (MST) using Prim's Algorithm
void Graph::PrimMST() const {
    // Vector to track if a vertex is included in the MST
    std::vector<bool> inMST(countOfVertices, false);

    // Vector to store the parent of each vertex (for MST construction)
    std::vector<int> parent(countOfVertices, -1);

    // Vector to store the minimum edge weight to each vertex (key)
    std::vector<int> key(countOfVertices, INT_MAX);

    // Start from vertex 0
    key[0] = 0;

    // Min-heap (priority queue) to store vertices with the minimum key values
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;

    // Push the starting vertex (vertex 0) into the priority queue with key 0
    pq.push({0, 0});

    // While there are vertices to process
    while (!pq.empty()) {
        // Get the vertex with the minimum key value
        int u{pq.top().second};
        pq.pop();

        // If u is already in the MST, skip it
        if (inMST[u]) continue;

        // Mark u as included in the MST
        inMST[u] = true;

        // Explore all adjacent vertices of u
        for (int v : adjList[u]) {
            int weight{1};  // Set the weight of each edge (this can be adjusted)
            
            // If v is not yet in the MST and the edge weight is less than the current key for v
            if (!inMST[v] && weight < key[v]) {
                // Update the key for v and set u as its parent
                key[v] = weight;
                pq.push({key[v], v});  // Push v into the priority queue
                parent[v] = u;  // Set u as the parent of v
            }
        }
    }

    // Print the edges in the MST formed by Prim's Algorithm
    std::cout << "Prim's MST Edges:\n";
    for (int i{1}; i < countOfVertices; ++i) {
        std::cout << parent[i] << " - " << i << std::endl;  // Print parent-child edge in MST
    }
}

// Helper function for DFS (Depth First Search) to traverse the graph recursively
void Graph::dfsRecursiveHelper(int source, std::vector<bool>& visited) {
    // Mark the source vertex as visited
    visited[source] = true;
    
    // Print the source vertex (part of the DFS traversal)
    std::cout << source << ' ';

    // Recursively visit all unvisited adjacent vertices
    for (int v : adjList[source]) {
        if (!visited[v]) {
            dfsRecursiveHelper(v, visited);  // Recurse for the adjacent vertex
        }
    }
}

// Helper function for counting the number of paths from 'source' to 'destination' in the graph.
void Graph::pathCountHelper(int source, int destination, int& count, std::vector<bool>& visited) {
    // Mark the current vertex as visited
    visited[source] = true;

    // If the source is equal to destination, increment the path count
    if (source == destination) {
        ++count;
    }
    else {
        // Recursively visit all unvisited adjacent vertices
        for (int v : adjList[source]) {
            if (!visited[v]) {
                pathCountHelper(v, destination, count, visited);  // Recurse with the next vertex
            }
        }
    }
    // Mark the current vertex as unvisited for other paths
    visited[source] = false;
}

// Helper function to find all paths from 'source' to 'destination' and store them.
void Graph::allPathsHelper(
    int source,
    int destination,
    std::vector<std::vector<int>>& paths,
    std::vector<int> currentPath,
    std::vector<bool>& visited)
{
    // Mark the current vertex as visited
    visited[source] = true;

    // Add the current vertex to the current path
    currentPath.push_back(source);

    // If the destination is reached, add the current path to the list of paths
    if (source == destination) {
        paths.push_back(currentPath);
    }
    else {
        // Recursively visit all unvisited adjacent vertices
        for (int v : adjList[source]) {
            if (!visited[v]) {
                // Recurse with the next vertex
                allPathsHelper(v, destination, paths, currentPath, visited);
            }
        }
    }

    // Backtrack: remove the current vertex from the path and mark it as unvisited
    currentPath.pop_back();
    visited[source] = false;
}

// Helper function to detect if there is a cycle in the graph using DFS.
bool Graph::cycleDetectionHelper(std::vector<bool>& visited, int vertex, int parent) {
    // Mark the current vertex as visited
    visited[vertex] = true;
    
    // Recursively visit all adjacent vertices
    for (int neighbor : adjList[vertex]) {
        // If the neighbor is not visited, recurse with the neighbor
        if (!visited[neighbor]) {
            if (cycleDetectionHelper(visited, neighbor, vertex)) {
                return true;  // If cycle is detected, return true
            }
        }
        // If the neighbor is visited and not the parent, a cycle is detected
        else if (neighbor != parent) {
            return true;  // Cycle detected
        }
    }
    return false;  // No cycle detected from the current vertex
}

// Function to validate if a vertex is within the valid range (0 to countOfVertices - 1).
void Graph::validateVertex(int vertex) const {
    // Check if the vertex is out of range
    if (vertex < 0 || vertex >= countOfVertices) {
        throw std::out_of_range(
            "Vertex " + std::to_string(vertex) + " is out of valid range [0, " + 
            std::to_string(countOfVertices - 1) + "]."
        );
    }
}

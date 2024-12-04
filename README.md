# undirected_graph_implementation

# Graph Implementation in C++

This project provides a **C++ implementation of a graph** data structure. It includes functionalities for representing a graph, adding/removing vertices and edges, performing graph traversal algorithms (Depth-First Search, Breadth-First Search), and finding all paths between two vertices. The graph is represented using an adjacency list, and this project is designed for educational purposes to demonstrate common graph algorithms and operations.

## Table of Contents
- [Graph Data Structure Overview](#graph-data-structure-overview)
- [Key Concepts](#key-concepts)
- [Graph Class Implementation](#graph-class-implementation)
  - [Constructor & Destructor](#constructor--destructor)
  - [Adding/Removing Vertices and Edges](#addingremoving-vertices-and-edges)
  - [Graph Traversal](#graph-traversal)
  - [Find All Paths](#find-all-paths)
- [Usage](#usage)
- [Example](#example)
- [Compiling the Code](#compiling-the-code)
- [License](#license)

## Graph Data Structure Overview

A **graph** is a non-linear data structure consisting of a set of **vertices (nodes)** connected by **edges (arcs)**. Graphs are used to represent relationships, networks, and connections between data points. For example, graphs can be used to represent social networks, transportation systems, and even the structure of the internet.

Graphs can be classified into different types:
- **Directed Graphs**: In these graphs, edges have a direction (i.e., the edge points from one vertex to another).
- **Undirected Graphs**: In these graphs, edges are bidirectional (i.e., the edge connects two vertices in both directions).
- **Weighted Graphs**: Each edge in the graph has an associated weight (usually representing distance, cost, or other metrics).
- **Unweighted Graphs**: In these graphs, edges have no weight and represent simple connections between vertices.

This implementation uses an **adjacency list** to store the graph, where each vertex points to a list of adjacent vertices. This representation is efficient for sparse graphs and provides easy access to a vertex’s neighbors.

## Key Concepts

- **Vertex**: A node in the graph. It can represent any entity, such as a city, person, or computer.
- **Edge**: A connection between two vertices. It represents a relationship or interaction between two entities.
- **Adjacency List**: A data structure that stores the list of adjacent vertices for each vertex.
- **Path**: A sequence of vertices connected by edges.
- **DFS (Depth-First Search)**: A traversal algorithm that explores as far as possible along each branch before backtracking.
- **BFS (Breadth-First Search)**: A traversal algorithm that explores all neighbors of a vertex before moving to the next level of vertices.

## Graph Class Implementation

The `Graph` class is designed to perform various operations on the graph, such as adding vertices and edges, performing DFS and BFS, and finding all paths between two vertices.

### Constructor & Destructor

```cpp
Graph(int n); // Constructor: Initializes the graph with `n` vertices
~Graph();     // Destructor: Cleans up memory when the graph object is destroyed

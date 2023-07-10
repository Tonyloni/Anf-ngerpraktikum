#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <iostream>

using namespace std;

struct Node {
    int index; //index for adjacency array
    int adj_nodes; //number of adjacent nodes
    double coordinate_x;
    double coordinate_y;

    Node();

    void setCoords(double x, double y);
};

struct Edge {
    int target;
    double weight;

    Edge();
    Edge(int target);
    Edge(int target, double weight);
};

class Graph {
    
    protected:
    int numNodes;
    int numEdges;

    //adjacency array
    std::vector<Node> Nodes;
    std::vector<Edge> Edges;

    public:
 
    Graph(int n, int m);
    Graph() {}

    int getNumNodes();
    int getNumEdges();
    Node getNode(int i);
    Edge getEdge(int i);
    void addEdge(int sourcenode, int targetnode);
    void addCoords(int node, double x, double y);
    void printGraph();
};

void readCoordinates(string filename, Graph* g);
Graph* readGraph(string filename_graph, string filename_coords);

#endif //GRAPH_H
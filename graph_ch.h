#ifndef GRAPH_CH_H
#define GRAPH_CH_H

#include "graph.h"
#include <iostream>
#include <vector>
using namespace std;

extern int num_er; //number edge reductions

struct Node_ch : Node {
    int id; //safe id of the node to access node array directly
    int cnc; //contracted neighbour counter
    int importance; //importance for hierarchy
    int level; //level of hierarchy
    int end_index; //index to the end of the edge group, following free space for additional edges
    int adj_nodes_ch; //number of adjacent nodes that are not contracted
    //bool contracted;
    int end; //end of edge group (after additional space)
    int second_index;

    Node_ch();

    bool operator==(const Node_ch& other) const;
};

struct Edge_ch : Edge {
    Edge_ch(int target);
    Edge_ch();
    Edge_ch(int target, double weight);
};

class Graph_ch {
    private:
    int level;
    int numNodes;
    int numEdges;

    //adjacency array
    vector<Node_ch> Nodes;
    vector<Edge_ch> Edges;

    public:
 
    Graph_ch(int n, int m);
    Graph_ch() {}

    int getLevel();
    void increaseLevel();
    int getNumNodes();
    int getNumEdges();
    Node_ch* getNode(int i);
    void setNode(int index, Node_ch* u);
    Edge_ch* getEdge(int i);
    void setEdge(int index, Edge_ch* e);
    void swapEdges(int edge1, int edge2);
    void addEdge(int sourcenode, int targetnode);
    void removeEdge(int edge_index, int node_index);
    void addCoords(int node, double x, double y);
    void addShortcut(int node_u, int node_w, double d_sc);
    void printGraph();
    void printGraph2();
    void printEdges();
    void writeGraphIntoFile();
};

Graph_ch* readContractionHierarchieFromFile(string nodes_file, string edges_file);
Graph_ch* readGraph_ch(string filename_graph, string filename_coords);
void readCoordinates(string filename, Graph_ch* g);

#endif // GRAPH_CH_H
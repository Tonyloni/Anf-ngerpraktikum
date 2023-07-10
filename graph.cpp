#include <iostream>
#include <list> 
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include "graph.h"

using namespace std;

Node::Node() {
        coordinate_x = 0;
        coordinate_y = 0;
        adj_nodes = 0;
        index = 0;
    }

void Node::setCoords(double x, double y) {
        this->coordinate_x = x;
        this->coordinate_y = y;
    }

Edge::Edge() {
        weight = INFINITY;
}

Edge::Edge(int target) {
        weight = INFINITY;
        this->target = target;
    }

Edge::Edge(int target, double weight) {
        this->target = target;
        this->weight = weight;
    }
 
Graph::Graph(int n, int m) {
        Nodes = vector<Node>(n);
        Edges = vector<Edge>(m);
        this->numNodes = n;
        this->numEdges = m;
    }

int Graph::getNumNodes() {
        return this->numNodes;
    }

int Graph::getNumEdges() {
        return this->numEdges;
    }

Node Graph::getNode(int i) {
        return this->Nodes[i];
    }

Edge Graph::getEdge(int i) {
        return this->Edges[i];
    }

void Graph::addEdge(int sourcenode, int targetnode) {

        if (sourcenode == 0) {
            Nodes[sourcenode].index = 0;
        }
        else {
            Nodes[sourcenode].index = Nodes[sourcenode-1].index + Nodes[sourcenode-1].adj_nodes;
        }

        //create new edge, calculate weight and add it to the adjacency array
        Edge e = Edge(targetnode);

        //euclidian norm? idk dude, but should be
        
        e.weight = sqrt(    pow(Nodes[targetnode-1].coordinate_x - Nodes[sourcenode].coordinate_x, 2.0)
                        +   pow(Nodes[targetnode-1].coordinate_y - Nodes[sourcenode].coordinate_y, 2.0));

        //fill new edge in adjacency array 
        Edges[Nodes[sourcenode].index + Nodes[sourcenode].adj_nodes] = e;
        Nodes[sourcenode].adj_nodes++;
    }

void Graph::addCoords(int node, double x, double y) {
        Nodes[node].setCoords(x, y);
    }

void Graph::printGraph() {
        for (int i = 0; i < numNodes; i++) {
            cout << i+1 << " -> ";
            for (int j = 0; j < Nodes[i].adj_nodes; j++) {
                cout << Edges[Nodes[i].index + j].target << " /w:" << Edges[Nodes[i].index + j].weight << " ";
            }
            cout << ", ";
            cout << Nodes[i].coordinate_x << "/" << Nodes[i].coordinate_y;
            cout << endl;
        }
    }

void readCoordinates(string filename, Graph* g) {

    //read coordinate file and add the coordinates to the Node array of the graph

    ifstream file(filename);

    if (!file.is_open()) {
        cout << "error: file not open" << endl;
    }

    string line;

    int i = 0;

    while (getline(file, line)) {

        if (line[0] == '%') {
            continue;
        }

        stringstream ss(line);

        long double coord_x;
        long double coord_y;

        ss >> coord_x;
        ss >> coord_y;

        g->addCoords(i, coord_x, coord_y);

        i++;
    }

    file.close();

}

Graph* readGraph(string filename_graph, string filename_coords) {

    int nodes;
    int edges;
    
    ifstream file(filename_graph);

    if (!file.is_open()) {
        cout << "error: file not open" << endl;
    }

    //read first line of the graph file to initialize the graph with nodes and edges
    string line;

    getline(file, line);

    while (line[0] == '%') {
        getline(file, line);
    }

    stringstream ss(line);
    ss >> nodes;
    ss >> edges;

    //initialize graph
    Graph* g = new Graph(nodes, 2*edges);

    //read coordinate file and initialize coordinates of the nodes in the graph
    readCoordinates(filename_coords, g);

    int i = 0;

    //read adjacencies of the nodes in the graph file
    while (getline(file, line)) {

        if (line[0] == '%') {
            continue;
        }

        stringstream ss(line);

        int target;

        while ( ss >> target) {
            g->addEdge(i, target);
        }

        i++;
    }

    file.close();

    //g.printGraph();

    return g;

}




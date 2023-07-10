#ifndef GRAPH_CH2_H
#define GRAPH_CH2_H

#include "graph_ch.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

struct Node_ch2 {
    int id;
    double coord_x;
    double coord_y;
    int second_index;
    int importance;
    int cnc; //contracted neighbour counter
    int level;
    std::vector<Edge_ch> edges;
};

template <class Iterator> struct Range {
Iterator _b, _e;
auto begin() { return _b; }
auto end() { return _e; }
const auto begin() const { return _b; }
const auto end() const { return _e; }
using value_type = typename Iterator::value_type;
Range(Iterator begin, Iterator end) : _b(begin), _e(end) {}
};

class Graph_ch2 {
    private: 
    std::vector<Node_ch2> Nodes;
    int level;

    public:
    Graph_ch2() {}
    Graph_ch2(int n);
    int getLevel();
    void increaseLevel();
    int getNumNodes();
    Node_ch2* getNode(int i);
    void addEdge(int source, int target);
    void addCoords(int node, double x, double y);
    void addShortcut(int node_u, int node_w, double d_sc);
    std::vector<Node_ch2>::iterator getBegin();
    std::vector<Node_ch2>::iterator getEnd();
    void printGraph();
    void printGraphNotContracted();
    void writeGraphIntoFile2(string nodefile, string edegfile);
};

void readCoordinates(string filename, Graph_ch2* g);
Graph_ch2* readGraph_ch2(string filename_graph, string filename_coords);
Graph_ch* readContractionHierarchieFromFile2(string nodes_file, string edges_file);

#endif //GRAPH_CH2_H
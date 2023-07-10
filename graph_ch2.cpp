#include "graph_ch.h"
#include "graph_ch2.h"
#include <iostream>
#include <iomanip>

Graph_ch2::Graph_ch2(int n) {
    Nodes = std::vector<Node_ch2>(n);
    level = 0;
}

int Graph_ch2::getLevel() {
    return this->level;
}

void Graph_ch2::increaseLevel() {
    this->level++;
}

int Graph_ch2::getNumNodes() {
    return Nodes.size();
}

Node_ch2* Graph_ch2::getNode(int i) {
    Node_ch2* node = &(this->Nodes[i]);
    return node;
}

std::vector<Node_ch2>::iterator Graph_ch2::getBegin() {
    return Nodes.begin();
}

std::vector<Node_ch2>::iterator Graph_ch2::getEnd() {
    return Nodes.end();
}

void Graph_ch2::addEdge(int source, int target) {
    Nodes[source].id = source;
    Nodes[source].second_index = 0;
    Nodes[source].level = -1;

    Edge_ch e = Edge_ch(target);
    //calculate euclidian norm //MULTIPLICATION BETTER
    e.weight = sqrt(    pow(Nodes[target].coord_x - Nodes[source].coord_x, 2.0)
                    +   pow(Nodes[target].coord_y - Nodes[source].coord_y, 2.0));
   
    Nodes[source].edges.push_back(e);
}

void Graph_ch2::addCoords(int node, double x, double y) {
    Nodes[node].coord_x = x;
    Nodes[node].coord_y = y;
}

void Graph_ch2::addShortcut(int node_u, int node_w, double d_sc) {
    Edge_ch e_w = Edge_ch(node_w, d_sc);
    Edge_ch e_u = Edge_ch(node_u, d_sc);

    Nodes[node_u].edges.push_back(e_w);
    Nodes[node_w].edges.push_back(e_u);
}

void Graph_ch2::printGraph() {
    for (auto& node : Nodes) {
        std::cout << node.id << "/i:" << node.importance << "/seci:" << node.second_index << "/l:" << node.level << " -> ";
        for (auto& edge : node.edges) {
            std::cout << edge.target <<  "/l:" << Nodes[edge.target].level << "/w:" << edge.weight << ", ";
        } 
        std::cout << std::endl;
    }
}

void Graph_ch2::printGraphNotContracted() {
    for (auto& node : Nodes) {
        std::cout << node.id << "/i:" << node.importance << "/level:" << node.level << " -> ";
        Range<std::vector<Edge_ch>::iterator> range_nc_edges(node.edges.begin() + node.second_index, node.edges.end());
        for (auto& edge : range_nc_edges) {
            std::cout << edge.target << "/w:" << edge.weight << ", ";
        }
        std::cout << std::endl;
    }
}

void Graph_ch2::writeGraphIntoFile2(string nodefile, string edgefile) {   
    int numEdges = 0;
    for (const auto& node : Nodes) {
        numEdges += node.edges.size();
    }

    std::ofstream NodeFile(nodefile);
    int edge_index = 0;
    if (NodeFile.is_open()) {
        NodeFile << Nodes.size() << " " << numEdges << " " << "\n";
        for (const auto& node : Nodes) {
            NodeFile << edge_index << " " << edge_index + node.second_index << " " << edge_index + node.edges.size()-1 << " " << node.level << "\n";
            edge_index += node.edges.size();
        }
        /*for (int i = 0; i < Nodes.size(); i++) {
            NodeFile << Nodes[i].index << " " << Nodes[i].second_index << " " << Nodes[i].end_index << "\n";
        }*/

        NodeFile.close();
    } else {
        std::cout << "unable to open file ch_nodes.txt" << std::endl;
    }

    std::ofstream EdgeFile(edgefile);

    if (EdgeFile.is_open()) {
        for (const auto& node : Nodes) {
            for (const auto& edge : node.edges) {
                EdgeFile << edge.target << " " << std::setprecision(15) << edge.weight << "\n";
            }
        }
        EdgeFile.close();
    } else {
        std::cout << "unable to open file ch_edges.txt" << std::endl;
    }
    std::cout << "graph written into file" << std::endl;
}

void readCoordinates(string filename, Graph_ch2* g) {

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

        double coord_x;
        double coord_y;

        ss >> coord_x;
        ss >> coord_y;

        g->addCoords(i, coord_x, coord_y);

        i++;
    }
    file.close();
}

Graph_ch2* readGraph_ch2(string filename_graph, string filename_coords) {
    int nodes;
    ifstream file(filename_graph);

    if (!file.is_open()) {
        cout << "error: file not open" << endl;
    }
    //read first line of the graph file to initialize the graph with nodes only -> edges are added dynamically via pushback
    string line;

    getline(file, line);

    while (line[0] == '%') {
        getline(file, line);
    }

    stringstream ss(line);
    ss >> nodes;

    //initialize graph
    Graph_ch2* g = new Graph_ch2(nodes);

    //read coordinate file and initialize coordinates of the nodes in the graph
    std::cout << "read Coordinates..." << std::endl;
    readCoordinates(filename_coords, g);

    int i = 0;
    std::cout << "read graph... " << std::endl;
    //read adjacencies of the nodes in the graph file
    while (getline(file, line)) {

        if (line[0] == '%') {
            continue;
        }

        stringstream ss(line);
        int target;
        while ( ss >> target) {
            g->addEdge(i, target-1); 
        }
        i++;
    }
    file.close();

    return g;
}

Graph_ch* readContractionHierarchieFromFile2(string nodes_file, string edges_file) {
    long int nodes;
    long int edges;
    
    ifstream nodeFile(nodes_file);

    if (!nodeFile.is_open()) {
        cout << "error: file not open" << endl;
    }

    //read first line of the nodefile to initialize the graph with nodes and edges
    string line;

    getline(nodeFile, line);

    while (line[0] == '%') {
        getline(nodeFile, line);
    }

    stringstream ss(line);
    ss >> nodes;
    ss >> edges;

    //initialize graph
    Graph_ch* g = new Graph_ch(nodes, edges);

    int i = 0;

    //read indices of the nodes in the node file
    while (getline(nodeFile, line)) {

        if (line[0] == '%') {
            continue;
        }

        stringstream ss(line);

        int index_;
        int second_index_;
        int end_index_;
        int level_;

        ss >> index_;
        ss >> second_index_;
        ss >> end_index_; 
        ss >> level_;   

    	g->getNode(i)->index = index_;
        g->getNode(i)->second_index = second_index_;
        g->getNode(i)->end_index = end_index_;
        g->getNode(i)->level = level_;

        i++;
    }

    nodeFile.close();

    int j = 0;

    ifstream edgeFile(edges_file);

    while (getline(edgeFile, line)) {

        if (line[0] == '%') {
            continue;
        }

        stringstream ss(line);

        int target_;
        double weight_;

        ss >> target_;
        ss >> weight_;

        g->getEdge(j)->target = target_;
        g->getEdge(j)->weight = weight_;
        j++;
    }

    edgeFile.close();

    return g;
} 



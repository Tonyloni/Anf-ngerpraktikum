#include <iostream>
#include "argtable3.h"
#include <fstream>
#include <chrono>
#include "graph.h"
#include "graph_ch.h"
#include "graph_ch2.h"
//#include "contraction_h.cpp"
#include "contraction_h2.cpp"
#include <cmath>
#include <numeric>
void printPath(std::vector<int> path) {
    for (int i = 0; i < path.size(); i++) {
        std::cout << path[i] << " -> ";
    }
    std::cout << std::endl;
}



double getWeight(Graph_ch* g, int source, int target) {
    for (int i = g->getNode(source)->index; i <= g->getNode(source)->end_index; i++) {
        if (g->getEdge(i)->target == target) {
            return g->getEdge(i)->weight;
        }
    }
    return -1;
}

double getWeight2(Graph_ch2* g, int source, int target) {
    std::pair<std::vector<int>, double> weight;
    weight.second = 0.0;
    weight = bidirec_dijkstra2(g, source, target);
    return weight.second;
}



int main(int argc, char* argv[]) {

    
    int source_node = -1;
    int target_node = -1;
    const char* graph_filename = NULL;
    const char* coordinates_filename = NULL;
    const char* nodes_ch_filename = NULL;
    const char* edges_ch_filename = NULL;

    // Define the argument structures
    struct arg_int* arg_source = arg_int0("s", "source", "<s>", "Source node");
    struct arg_int* arg_target = arg_int0("t", "target", "<t>", "Target node");
    struct arg_str* arg_graph = arg_str0("g", "graph", "<file>", "Graph filename");
    struct arg_str* arg_coordinates = arg_str0("c", "coordinates", "<file>", "Coordinates filename");
    struct arg_str* arg_ch_nodes = arg_str0("n", "nodes", "<file>", "ch nodes filename");
    struct arg_str* arg_ch_edges = arg_str0("e", "edegs", "<file>", "ch edges filename");
    struct arg_end* end = arg_end(20);

    // Create an array of pointers to the argument structures
    void* argtable[] = { arg_source, arg_target, arg_graph, arg_coordinates, arg_ch_nodes, arg_ch_edges, end }; //, 

    // Parse the command-line arguments
    int nerrors = arg_parse(argc, argv, argtable);

    // Check for parsing errors
    if (nerrors > 0) {
        arg_print_errors(stderr, end, "my_program");
        arg_print_syntax(stderr, argtable, "\n");
        arg_print_glossary(stderr, argtable, "  %-25s %s\n");
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        exit(1);
    }

    // Get the values of the parsed arguments
    source_node = arg_source->count > 0 ? arg_source->ival[0] : -1;
    target_node = arg_target->count > 0 ? arg_target->ival[0] : -1;
    graph_filename = arg_graph->count > 0 ? arg_graph->sval[0] : NULL;
    coordinates_filename = arg_coordinates->count > 0 ? arg_coordinates->sval[0] : NULL;
    nodes_ch_filename = arg_ch_nodes->count > 0 ? arg_ch_nodes->sval[0] : NULL;
    edges_ch_filename = arg_ch_edges->count > 0 ? arg_ch_edges->sval[0] : NULL;

    // Free the memory allocated by Argtable
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

    // Check if all required arguments are provided
    if (source_node == -1 || target_node == -1 || graph_filename == NULL || coordinates_filename == NULL || nodes_ch_filename == NULL || edges_ch_filename == NULL) { //
        printf("Error: Missing required arguments\n");
        exit(1);
    }

    printf("Source node: %d\n", source_node);
    printf("Target node: %d\n", target_node);
    printf("Graph filename: %s\n", graph_filename);
    printf("Coordinates filename: %s\n", coordinates_filename);
    printf("ch nodes filename: %s\n", nodes_ch_filename);
    printf("ch edges filename: %s\n", edges_ch_filename);

    
    Graph_ch2* g1 = new Graph_ch2();
    g1 = readGraph_ch2(graph_filename, coordinates_filename);
        Graph_ch* g2 = new Graph_ch();
    g2 = readContractionHierarchieFromFile2(nodes_ch_filename, edges_ch_filename);
    //contractionHierarchy(g1);
    std::vector<double> times_ch;
    std::vector<double> times_normal;

    for(int tries = 0;tries<100;tries++){
        if(tries %10==0){
            std::cout<<"Versuch:"<< tries<<std::endl;
        }
source_node = rand()%g1->getNumNodes()+1;
target_node = rand()%g1->getNumNodes()+1;

  
    auto start_time_query_normal = std::chrono::high_resolution_clock::now();
    std::pair<std::vector<int>, double> path_normal = dijkstra(g1, source_node-1, target_node-1);
    auto end_time_query_normal = std::chrono::high_resolution_clock::now();
    auto duration_query_normal = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_query_normal - start_time_query_normal).count();
    times_normal.push_back(duration_query_normal/1000.0);
   // g1->writeGraphIntoFile2(nodes_ch_filename, edges_ch_filename);
   
    //std::pair<std::vector<int>, double> path_normal = bidirec_dijkstra2(g1, source_node-1, target_node-1);

    //double distance_ch = modifiedBidirectional(g2, source_node-1, target_node-1);
    auto start_time_query_ch = std::chrono::high_resolution_clock::now();
    std::pair<std::vector<int>, double> path_ch = bidirec_dijkstra(g2, source_node-1, target_node-1);
    auto end_time_query_ch = std::chrono::high_resolution_clock::now();
    auto duration_query_ch = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_query_ch - start_time_query_ch).count();
    times_ch.push_back(duration_query_ch/1000.0);
if(std::abs( path_normal.second- path_ch.second)>1e-5&&path_normal.second!=0.0){
    std::cerr<<source_node<<" to "<<target_node<<" is off by "<<std::abs( path_normal.second- path_ch.second)<<std::endl;
      printPath(path_normal.first);
    std::cout << "pathsize_normal: " << path_normal.first.size() << std::endl;
    std::cout << "distance_normal: " << path_normal.second << std::endl;
    std::cout << "time_normal: " << duration_query_normal << std::endl;


    //std::cout << "distance ch: " << distance_ch << std::endl;
    printPath(path_ch.first);
    std::cout << "pathsize_ch: " << path_ch.first.size() << std::endl;
    std::cout << "distance_ch: " << path_ch.second << std::endl;
    std::cout << "time_ch: " << duration_query_ch << std::endl;


    std::cout << "weights path_ch:" << std::endl;
    for (int i = 0; i < path_ch.first.size()-1; i++) {
        std::cout << getWeight(g2, path_ch.first[i], path_ch.first[i+1]) << ", ";
    }
    std::cout << std::endl;
}
      }
      //calculate means
      std::cout<<"AVG normal "<< std::accumulate(times_normal.begin(),times_normal.end(),0.0)/times_normal.size()<<std::endl;
      std::cout<<"AVG ch "<< std::accumulate(times_ch.begin(),times_ch.end(),0.0)/times_ch.size()<<std::endl;
        std::cout<<"GEO normal "<< exp(std::transform_reduce(times_normal.begin(),times_normal.end(),0.0,[](auto a,auto b){return a+b;},[](auto a){return log(a+0.00001);})/times_normal.size())<<std::endl;
      std::cout<<"GEO ch "<< exp(std::transform_reduce(times_ch.begin(),times_ch.end(),0.0,[](auto a,auto b){return a+b;},[](auto a){return log(a+0.00001);})/times_ch.size())<<std::endl;
    //std::cout << "distance: " << path_ch.second << std::endl;
    //double path_ch = bidirec_dijkstra2(g1, source_node-1, target_node-1);

  

    /*std::cout << "weights path_normal:" << std::endl;
    for (int i = 0; i < path_ch.first.size()-1; i++) {
        std::cout << getWeight2(g1, path_ch.first[i], path_ch.first[i+1]) << ", ";
    }
    std::cout << std::endl;

    bool contains = false;
    for (int i = g2->getNode(28072)->index; i <= g2->getNode(28072)->end_index; i++) {
        //std::cout << "HERE" << std::endl;
        if (g2->getEdge(i)->target == 27889) {
            contains = true;
        }
    }
    if (contains) {
        std::cout << "contained" << std::endl;
    } else {
        std::cout << "NOT CONTAINED" << std::endl;
    }
*/
    
/*
double epsilon = 1e-6;
    for (int i = 1; i < g2->getNumNodes(); i++) {
        std::pair<std::vector<int>, double> path_ch = bidirec_dijkstra(g2, 0, i);
        std::pair<std::vector<int>, double> path_normal = bidirec_dijkstra2(g1, 0, i);
        if  (std::abs(path_ch.second - path_normal.second) >= epsilon) {
            std::cout << "not same distance to target: " << i << std::endl;
            std::cout << "path_ch.second: " << path_ch.second << std::endl;
            std::cout << "path_normal.second: " << path_normal.second << std::endl;
            break;
        }
    }
    */

    //contractionHierarchy2(g1, source_node, target_node);
    //g->printGraph();
    //std::cout << bidirec_dijkstra(g, source_node-1, target_node-1) << std::endl;
    //g->printGraph();
    //g1->writeGraphIntoFile2();

    //g->printGraph();
     
    
    auto start_time_query = std::chrono::high_resolution_clock::now();

    auto end_time_query = std::chrono::high_resolution_clock::now();
    auto duration_query = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_query - start_time_query).count();

    //std::cout << "Time taken for query: " << duration_query << " milliseconds" << std::endl;

    return 0;
}

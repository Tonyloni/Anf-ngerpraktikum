#include <iostream>
#include "graph_ch.h"
#include "graph_ch2.h"
#include "dijkstra.cpp"
#include <fstream>
#include <sstream>
#include <chrono>

int main(int argc, char* argv[]) {
    std::vector<std::pair<int, int>> TestNodes;
    std::ifstream testNodes("../random10000bel");

    if (!testNodes.is_open()) {
        cout << "error: file not open" << endl;
    }

    //read first line of the nodefile to initialize the graph with nodes and edges
    string line;

    while (getline(testNodes, line)) {

        stringstream ss(line);

        int source;
        int target;

        ss >> source;
        ss >> target; 

        TestNodes.push_back({source, target}); 
    }
    std::cout << TestNodes.size() << std::endl;

    std::vector<int> v;
    

    Graph_ch2* g_normal = new Graph_ch2();
    g_normal = readGraph_ch2("../bel.graph", "../bel.graph.xyz");
    Graph_ch* g_ch = new Graph_ch();
    g_ch = readContractionHierarchieFromFile2("../ch_nodes_bel2.txt", "../ch_edges_bel2.txt");

    std::ofstream TestResults("../testresults10000bel.txt");

    if (TestResults.is_open()) {
        for (int i = 0; i < TestNodes.size(); i++) {
            auto start_time_query_normal = std::chrono::high_resolution_clock::now();
            std::pair<std::vector<int>, double> path_normal = dijkstra(g_normal, TestNodes[i].first-1, TestNodes[i].second-1);
            auto end_time_query_normal = std::chrono::high_resolution_clock::now();
            auto duration_query_normal = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_query_normal - start_time_query_normal).count();

            auto start_time_query_ch = std::chrono::high_resolution_clock::now();
            std::pair<std::vector<int>, double> path_ch = bidirec_dijkstra(g_ch, TestNodes[i].first-1, TestNodes[i].second-1);
            auto end_time_query_ch = std::chrono::high_resolution_clock::now();
            auto duration_query_ch = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_query_ch - start_time_query_ch).count();


            if (std::abs( path_normal.second- path_ch.second)>1e-5&&path_normal.second!=0.0) {
                std::cout << "DISTANCE NOT THE SAME" << std::endl;
                break;
            }
            TestResults << duration_query_normal << " " << duration_query_ch << " " << path_normal.first.size() << "/n";
            if (i % 100 == 0) {
                std::cout << "test run: " << i << std::endl;
            }
        }
    }

    TestResults.close();
    std::cout << "results written into file" << std::endl;

    return 0;
}
#include <iostream>
#include <algorithm>
#include <limits>
#include "priorityqueue.cpp"
#include "adressable_pq.cpp"
#include <utility>
#include <cmath>
#include "graph_ch.h"
#include "graph.h"
#include "graph_ch2.h"
#include "minheap.h"

int num_er = 0;

std::vector<int> reconstructPathDijkstra(const std::vector<int>& pred, int source, int target) {
    std::vector<int> path;
    int currentNode = target;

    // Traverse the predecessor array from the target node to the source node
    while (currentNode != source) {
        path.push_back(currentNode);
        currentNode = pred[currentNode];
    }

    // Add the source node to the path
    path.push_back(source);

    // Reverse the path vector to get the correct order
    std::reverse(path.begin(), path.end());

    //print path
    /*for (int i = 0; i < path.size(); i++) {
        std::cout << path[i] << " -> ";
    }
    std::cout << "\n";*/

    return path;
}



std::pair<std::vector<int>, double> dijkstra(Graph_ch2* g, int source, int target) {

    //when accessing the d array, there has to be -1 to match the index of the graph file (there is no node 0)

    int n = g->getNumNodes();

    //create a priority queue where the element type is the node id (int) and the key is the distance, d[i] (double)
    //Adr_PriorityQueue<int, double>* p = new Adr_PriorityQueue<int, double>();
    MinHeap* p = new MinHeap();

    std::vector<double> d(n, INFINITY);
    std::vector<int> parent(n, -1);
    //double* d = new double[n];
    //parent array not really used here but is in the pseudocode
    //int parent[n];


    //parent[source] = source;
    d[source] = 0;

    p->insert(source, d[source]);

    while (!p->isEmpty()) {

        int u = p->deleteMin(); 

        if (u == target) {
            break;
        }
        Node_ch2* node_u = g->getNode(u);

        Range<std::vector<Edge_ch>::iterator> range_edges_u(node_u->edges.begin() + node_u->second_index, node_u->edges.end());
        for (auto& v : range_edges_u) { 
            if ((d[u] + v.weight) < d[v.target]) {
                d[v.target] = d[u] + v.weight;
                parent[v.target] = u;

                if (p->isElementOf(v.target)) {
                    p->updateKey(v.target, d[v.target]);
                }
                else {
                    p->insert(v.target, d[v.target]);
                }

            }
        }
    }
    if (d[target] == INFINITY) {
        return {std::vector<int>(), d[target]};
    }
    std::vector<int> shortest_path = reconstructPathDijkstra(parent, source, target);

    //delete[] d;

    //std::cout <<  "distance to target: " << target_dist << std::endl;

    //return shortest_path;
    return {shortest_path, d[target]};
}
/*

double dijkstra_ch(Graph_ch* g, int source, int target, double d_sc) { //modified dijkstra for local search

    //when accessing the d array, there has to be -1 to match the index of the graph file (there is no node 0)

    int n = g->getNumNodes();

    //create a priority queue where the element type is the node id (int) and the key is the distance, d[i] (double)
    Adr_PriorityQueue<int, double>* p = new Adr_PriorityQueue<int, double>();

    std::vector<double> d(n, INFINITY);
    //double* d = new double[n];
    //parent array not really used here but is in the pseudocode
    //int parent[n];


    //parent[source] = source;
    d[source] = 0;

    p->insert(source, d[source]);

    while (!p->isEmpty()) {

        int u = p->deleteMin(); 

        if (u == target) {
            break;
        }

        //relaxation CHANGE ADJ_NODES TO INDEX/ENDINDEX
        for (int i = g->getNode(u)->index; i <=  g->getNode(u)->end_index; i++) {
            //std::cout << "u: " << u << ", g->getNode(u)->id: " << g->getNode(u)->id << std::endl;

            if (d[u] > d_sc) { //break dijkstra when shortcut distance is smaller
                Edge_Reduction(g, d, g->getNode(source));
                    //std::cout << "d[u]: " << d[u] << std::endl;
                    return -1;
            }

            //std::cout << "g->getNode(g->getEdge(i)->target): " << g->getNode(g->getEdge(i)->target)->id << ", contracted: " << g->getNode(g->getEdge(i)->target)->contracted << std::endl;
            if (!g->getNode(g->getEdge(i)->target)->contracted) {
                //std::cout << "adjacent node: " << g->getNode(g->getEdge(i)->target)->id << std::endl;
                Edge_ch* v = g->getEdge(i);

                if ((d[u] + v->weight) < d[v->target]) {
                    
                    d[v->target] = d[u] + v->weight;
                    //std::cout << "d[v->target]: " << d[v->target] << std::endl;
                    //parent[v.target-1] = u;

                    if (v->target == target && d[v->target] <= d_sc) {
                        Edge_Reduction(g, d, g->getNode(source));
                            //std::cout << "RETURN d[v->target]: " << d[v->target] << std::endl;
                        return d[v->target];
                    }

                    if (p->isElementOf(v->target)) {
                        p->updateKey(v->target, d[v->target]);
                    }
                    else {
                        p->insert(v->target, d[v->target]);
                        //std::cout << "insert: " << v->target << ", d[v->target]: " << d[v->target] <<std::endl;
                    }

                }

            }
            

        }
    }

    delete p;
    Edge_Reduction(g, d, g->getNode(source));

    if (d[target] > d_sc) {
        return -1;
    }
    //delete[] d;

    //std::cout <<  "distance to target: " << target_dist << std::endl;

    return d[target];
}
*/

std::vector<int> reconstructPathBidirecDijkstra(const std::vector<int>& pred_f, const std::vector<int>& pred_b, int source, int target, int meetingNode) {
    std::vector<int> path;
    int currentNode;

    if (meetingNode == -1) {
        return path;
    }

    // Reconstruct path from source to meetingNode
    currentNode = meetingNode;
    while (currentNode != source) {
        path.push_back(currentNode);
        currentNode = pred_f[currentNode];
    }
    path.push_back(source);

    // Reverse the path vector
    std::reverse(path.begin(), path.end());

    // Reconstruct path from target to meetingNode
    currentNode = meetingNode;
    while (currentNode != target) {
        currentNode = pred_b[currentNode];
        path.push_back(currentNode);
    }

    //print path
    /*for (int i = 0; i < path.size(); i++) {
        std::cout << path[i] << " -> ";
    }
    std::cout << "\n";*/
    return path;
}


std::pair<std::vector<int>, double> bidirec_dijkstra2(Graph_ch2* g, int source, int target) {

    int n = g->getNumNodes();

    //create two pq, one for forward search, one for backward search
    /*Adr_PriorityQueue<int, double>* p_f = new Adr_PriorityQueue<int, double>();
    Adr_PriorityQueue<int, double>* p_b = new Adr_PriorityQueue<int, double>();*/

    MinHeap* p_f = new MinHeap();
    MinHeap* p_b = new MinHeap();

    //two distance arrays
    std::vector<double> d_f(n, INFINITY);    //double* d = new double[n];
    std::vector<double> d_b(n, INFINITY);

    //parent array 
    std::vector<int> parent_f(n, -1);
    std::vector<int> parent_b(n, -1);

    //parent[source] = source;
    d_f[source] = 0;
    d_b[target] = 0;

    p_f->insert(source, d_f[source]);
    p_b->insert(target, d_b[target]);

    double ten_sc = INFINITY; //tentative shortest path
    int meetingNode = -1;

    while (!p_f->isEmpty() && !p_b->isEmpty()) {

        int u = p_f->deleteMin();
        int v = p_b->deleteMin();
        Node_ch2* node_u = g->getNode(u);
        Node_ch2* node_v = g->getNode(v);

        //relaxation for all (u, x) - only relax edges with higher level
        Range<std::vector<Edge_ch>::iterator> range_edges_u(node_u->edges.begin() + node_u->second_index, node_u->edges.end());
        for (auto& x : range_edges_u) { 

            //relax (u, x)
            if ((d_f[u] + x.weight) < d_f[x.target]) {
                d_f[x.target] = d_f[u] + x.weight;
                parent_f[x.target] = u;

                if (p_f->isElementOf(x.target)) {
                    p_f->updateKey(x.target, d_f[x.target]);
                } else {
                    p_f->insert(x.target, d_f[x.target]);
                }

            }

            if (d_b[x.target] < INFINITY && (d_f[u] + x.weight + d_b[x.target]) < ten_sc) {
                ten_sc = d_f[u] + x.weight + d_b[x.target];
                meetingNode = x.target;
            }             
        }

        Range<std::vector<Edge_ch>::iterator> range_edges_v(node_v->edges.begin() + node_v->second_index, node_v->edges.end());
        //relaxation for all (v, x) - only relax edges with higher level
        for (auto& x : range_edges_v) {

            //relax (v, x)
            if ((d_b[v] + x.weight) < d_b[x.target]) {
                d_b[x.target] = d_b[v] + x.weight;
                parent_b[x.target] = v;
                //parent[v.target-1] = u;
                if (p_b->isElementOf(x.target)) {
                    p_b->updateKey(x.target, d_b[x.target]);
                } else {
                    p_b->insert(x.target, d_b[x.target]);
                }
            }
            if (d_f[x.target] < INFINITY && (d_b[v] + x.weight + d_f[x.target]) < ten_sc) {
                ten_sc = d_b[v] + x.weight + d_f[x.target];
                meetingNode = x.target;
            } 
        }
        //ten_sc is true shortest path distance    
        if (d_f[u] >= ten_sc && d_b[v] >= ten_sc) {
            //std::cout << "MEETING NODE: " << meetingNode << std::endl;
            //std::cout << "ten_sc: " << ten_sc << std::endl;
            break;
        }    
    }

    std::vector<int> path = reconstructPathBidirecDijkstra(parent_f, parent_b, source, target, meetingNode);
    return {path, ten_sc};
    //return ten_sc;
}

std::pair<std::vector<int>, double> bidirec_dijkstra(Graph_ch* g, int source, int target) {

    int n = g->getNumNodes();

    //create two pq, one for forward search, one for backward search
    Adr_PriorityQueue<int, double>* p_f = new Adr_PriorityQueue<int, double>();
    Adr_PriorityQueue<int, double>* p_b = new Adr_PriorityQueue<int, double>();
    
    //MinHeap* p_f = new MinHeap();
    //MinHeap* p_b = new MinHeap();

    //two distance arrays
    std::vector<double> d_f(n, INFINITY);    //double* d = new double[n];
    std::vector<double> d_b(n, INFINITY);

    //parent array 
    std::vector<int> parent_f(n, -1);
    std::vector<int> parent_b(n, -1);

    //parent[source] = source;
    d_f[source] = 0;
    d_b[target] = 0;

    p_f->insert(source, d_f[source]);
    p_b->insert(target, d_b[target]);

    double ten_sc = INFINITY; //tentative shortest path
    int meetingNode = -1;

    int u;
    int v;

    while (!(p_f->isEmpty() && p_b->isEmpty())) {

        if (!p_f->isEmpty()) {
            u = p_f->deleteMin();
        }
        if (!p_b->isEmpty()) {
            v = p_b->deleteMin();
        } 


        //relaxation for all (u, x) - only relax edges with higher level
        for (int i = g->getNode(u)->second_index; i <= g->getNode(u)->end_index; i++) {

            Edge_ch* x = g->getEdge(i);

            //relax (u, x)
            if ((d_f[u] + x->weight) < d_f[x->target]) {
                d_f[x->target] = d_f[u] + x->weight;
                parent_f[x->target] = u;

                if (p_f->isElementOf(x->target)) {
                    p_f->updateKey(x->target, d_f[x->target]);
                } else {
                    p_f->insert(x->target, d_f[x->target]);
                }

            }

            if (d_b[x->target] < INFINITY && (d_f[u] + x->weight + d_b[x->target]) < ten_sc) {
                ten_sc = d_f[u] + x->weight + d_b[x->target];
                meetingNode = x->target;
            }             
        }
        //relaxation for all (v, x) - only relax edges with higher level
        for (int i = g->getNode(v)->second_index; i <= g->getNode(v)->end_index; i++) {

            Edge_ch* x = g->getEdge(i);

            //relax (v, x)
            if ((d_b[v] + x->weight) < d_b[x->target]) {
                d_b[x->target] = d_b[v] + x->weight;
                parent_b[x->target] = v;
                //parent[v.target-1] = u;
                if (p_b->isElementOf(x->target)) {
                    p_b->updateKey(x->target, d_b[x->target]);
                } else {
                    p_b->insert(x->target, d_b[x->target]);
                }
            }
            if (d_f[x->target] < INFINITY && (d_b[v] + x->weight + d_f[x->target]) < ten_sc) {
                ten_sc = d_b[v] + x->weight + d_f[x->target];
                meetingNode = x->target;
            } 
        }   
        //std::cout << "d_f[u]: " << d_f[u] << ", d_b[u]: " << d_b[v] << std::endl;
        //std::cout << "p_f.size: " << p_f->size() << ", p_b.size: " << p_b->size() << std::endl;
        //ten_sc is true shortest path distance  
        if (d_f[u] > ten_sc && d_b[v] > ten_sc) {
            //std::cout << "MEETING NODE: " << meetingNode << std::endl;
            //std::cout << "ten_sc2: " << ten_sc << std::endl;
            break;
        }  
        /*if (d_f[u] + d_b[v] >= ten_sc) {
            break;
        }*/    
    }

    std::vector<int> path = reconstructPathBidirecDijkstra(parent_f, parent_b, source, target, meetingNode);
    //std::cout << "ten_sc: " << ten_sc << std::endl;
    return {path, ten_sc};
    //return ten_sc;
}

std::pair<std::vector<int>, double> modifiedBidirectional(Graph_ch* g, int source, int target) {
    int n = g->getNumNodes();

    //create two pq, one for forward search, one for backward search
    Adr_PriorityQueue<int, double>* p_f = new Adr_PriorityQueue<int, double>();
    Adr_PriorityQueue<int, double>* p_b = new Adr_PriorityQueue<int, double>();

    //two distance arrays
    std::vector<double> d_f(n, INFINITY);    //double* d = new double[n];
    std::vector<double> d_b(n, INFINITY);

    //parent array 
    std::vector<int> parent_f(n, -1);
    std::vector<int> parent_b(n, -1);

    std::vector<int> settledNodes;

    //parent[source] = source;
    d_f[source] = 0;
    d_b[target] = 0;

    p_f->insert(source, d_f[source]);
    p_b->insert(target, d_b[target]);

    double ten_sc = INFINITY; //tentative shortest path
    int meetingNode = -1;

    while (!p_f->isEmpty()) {
        int u = p_f->deleteMin();
        for (int i = g->getNode(u)->second_index; i <= g->getNode(u)->end_index; i++) {

            Edge_ch* x = g->getEdge(i);

            //relax (u, x)
            if ((d_f[u] + x->weight) < d_f[x->target]) {
                d_f[x->target] = d_f[u] + x->weight;
                parent_f[x->target] = u;

                if (p_f->isElementOf(x->target)) {
                    p_f->updateKey(x->target, d_f[x->target]);
                } else {
                    p_f->insert(x->target, d_f[x->target]);
                }

            }
        }
    }

    while (!p_b->isEmpty()) {
        int v = p_b->deleteMin();
        for (int i = g->getNode(v)->second_index; i <= g->getNode(v)->end_index; i++) {

            Edge_ch* x = g->getEdge(i);

            //relax (v, x)
            if ((d_b[v] + x->weight) < d_b[x->target]) {
                d_b[x->target] = d_b[v] + x->weight;
                parent_b[x->target] = v;
                //parent[v.target-1] = u;
                if (p_b->isElementOf(x->target)) {
                    p_b->updateKey(x->target, d_b[x->target]);
                } else {
                    p_b->insert(x->target, d_b[x->target]);
                }
            }
        }
        if (d_f[v] < INFINITY) {
            settledNodes.push_back(v);
        }  
    }
    std::cout << "settledNodes.size(): " << settledNodes.size() << std::endl;

    double min = d_f[settledNodes[0]] + d_b[settledNodes[0]];
    for (int i = 1; i < settledNodes.size(); i++) {
        //std::cout << "d_f[settledNodes[i]] + d_b[settledNodes[i]]: " << d_f[settledNodes[i]] + d_b[settledNodes[i]] << std::endl;
        if (d_f[settledNodes[i]] + d_b[settledNodes[i]] < min) {
            min = d_f[settledNodes[i]] + d_b[settledNodes[i]];
            meetingNode = settledNodes[i];
        }
    }

    std::vector<int> path = reconstructPathBidirecDijkstra(parent_f, parent_b, source, target, meetingNode);

    return {path, min};
}

#include "graph_ch2.h"
#include <iostream>
#include <cmath>
#include <iterator>
#include "dijkstra.cpp"
//#include "priorityqueue.cpp"

std::vector<double> d;
int num_sc = 0;
double time_ls = 0;
double time1 = 0;
double time2 = 0;
double time3 = 0;
double time4 = 0;
double time_btw = 0;

void calculateImportance(Graph_ch2* g, Node_ch2* node) {
    
    int adj_nodes_nc = node->edges.size() - node->second_index;
    //number of potential shortcuts
    int num_sc = (adj_nodes_nc * (adj_nodes_nc - 1)) / 2;
    node->importance = (num_sc - adj_nodes_nc) + 3 * node->cnc; //contracted neghbour counter
}

void localsearch(Graph_ch2* g, Node_ch2* node_u, Node_ch2* c_node, int edge_cnode_u) { //edge_cnode_u is edge c_node->u, c_node is the node that gets contracted
    //calculate maxPath := w(u, v) + max{w(v, w)}
    int c_node_sec_i = c_node->second_index; //contracted node second index
    double maxPath = c_node->edges[c_node_sec_i].weight;
    Range<std::vector<Edge_ch>::iterator> c_node_range(c_node->edges.begin() + c_node_sec_i+1, c_node->edges.end());

    for (auto& c_node_edge : c_node_range) {
        maxPath = std::max(maxPath, c_node_edge.weight);;
    }
    maxPath += c_node->edges[edge_cnode_u].weight;

    Adr_PriorityQueue<int, double>* p = new Adr_PriorityQueue<int, double>();
    std::fill(d.begin(), d.end(), INFINITY);
    d[node_u->id] = 0;
    p->insert(node_u->id, d[node_u->id]);

    //iterate through different targets
    int hop_limit = 0;
    Range<std::vector<Edge_ch>::iterator> ls_range(c_node->edges.begin() + edge_cnode_u+1, c_node->edges.end());
    for (auto& edge : ls_range) {
        double d_sc = c_node->edges[edge_cnode_u].weight + edge.weight;
        int target = edge.target;

        //check if edge already exists and only weight needs to be changed
        bool contained = false;
        Range<std::vector<Edge_ch>::iterator> node_u_range(node_u->edges.begin() + node_u->second_index, node_u->edges.end());

        for (auto& node_u_edge : node_u_range) {
            if (node_u_edge.target == target) {
                contained = true;
                if (d_sc < node_u_edge.weight) {
                    node_u_edge.weight = d_sc;

                    Node_ch2* target_node = g->getNode(target);
                    Range<std::vector<Edge_ch>::iterator> target_node_range(target_node->edges.begin() + target_node->second_index, target_node->edges.end());
                    for (auto& target_node_edge : target_node_range) {
                        if (target_node_edge.target == node_u->id) {
                            target_node_edge.weight = d_sc;
                            break;
                        }
                    }
                }
                break;
            }
        }
        if (contained) { //no dijkstra search if edge already exists
            continue;
        }

        //bool no_sc = false;
        bool sc = false;
        while (!p->isEmpty()) {

            /*if (no_sc) {
                no_sc = false;
                break;
            }*/

            //check hop-limit to reduce local search
            /*if (hop_limit == 500) {
                sc = true;
                hop_limit = 0;
                break;
            } */
                
            int x = p->deleteMin(); 

            if (x == target) { // || d[target] <= d_sc
                break;
            }

            if (d[x] > maxPath) { //break dijkstra when Pmax is smaller than current distance
                sc = true;
                break;
            }  

            //hop_limit++;
            Node_ch2* node_x = g->getNode(x);
            Range<std::vector<Edge_ch>::iterator> range_x(node_x->edges.begin() + node_x->second_index, node_x->edges.end());

            //realx all edges (x, v) (only not contracted nodes)
            for (auto& v : range_x) {
                if ((d[x] + v.weight) < d[v.target]) {
                    d[v.target] = d[x] + v.weight;
                            
                    /*if (v.target == target && d[v.target] <= d_sc) {
                        no_sc = true;
                        break;
                    }*/

                    if (p->isElementOf(v.target)) {
                        p->updateKey(v.target, d[v.target]);
                    } else {
                            p->insert(v.target, d[v.target]);
                    }

                }
            }
        }
        if ((d[target] > d_sc && d[target] < INFINITY) || sc) {
                //std::cout << "add shortcut:" << node_u->id << " -> " << target << std::endl;
                g->addShortcut(node_u->id, target, d_sc); 
                num_sc++;
        }
    }
    delete p;   
}

double localDijkstra(Graph_ch2* g, int target, double dist_sc, double maxPath, Adr_PriorityQueue<int, double>* pq) {
    int hoplimit = 0;
    while (!pq->isEmpty()) {

        int u = pq->deleteMin(); 

        if (u == target || d[u] > dist_sc) {
            break;
        }

        if (hoplimit == 500) {
            return dist_sc+1;
        }
        //relaxation
        Node_ch2* node_u = g->getNode(u);
        Range<std::vector<Edge_ch>::iterator> range_u(node_u->edges.begin() + node_u->second_index, node_u->edges.end());
        hoplimit++;
        for (auto& v : range_u) {
            if ((d[u] + v.weight) < d[v.target]) {
                d[v.target] = d[u] + v.weight;
                //parent[v.target] = u;

                if (v.target == target && d[v.target] <= dist_sc) {
                        return -1;
                    }

                if (pq->isElementOf(v.target)) {
                    pq->updateKey(v.target, d[v.target]);
                }
                else {
                    pq->insert(v.target, d[v.target]);
                }

            }
        }
    }
    return d[target];
}

void localsearch(Graph_ch2* g, Node_ch2* contracted_node) {

    //calculate Pmax
    double Pmax = 0;
    /*auto start_time1 = std::chrono::high_resolution_clock::now();
    for (int i = contracted_node->second_index; i < contracted_node->edges.size()-1; i++) {
        for (int j = i+1; j < contracted_node->edges.size(); j++) {
            if (contracted_node->edges[i].weight + contracted_node->edges[j].weight > Pmax) {
                Pmax = contracted_node->edges[i].weight + contracted_node->edges[j].weight;
            }
        }
    }
    auto end_time1 = std::chrono::high_resolution_clock::now();
    time1 += (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_time1- start_time1).count();
*/
    for (int i = contracted_node->second_index; i < contracted_node->edges.size(); i++) {
        auto start_timebtw = std::chrono::high_resolution_clock::now();
        double d_u_v = contracted_node->edges[i].weight;
        Node_ch2* u = g->getNode(contracted_node->edges[i].target);

        Adr_PriorityQueue<int, double>* p = new Adr_PriorityQueue<int, double>();
        std::fill(d.begin(), d.end(), INFINITY);
        d[u->id] = 0;
        p->insert(u->id, d[u->id]);
        auto end_timebtw = std::chrono::high_resolution_clock::now();
        time_btw += (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_timebtw - start_timebtw).count();

        //localSearch from u to different neighbours of v - add range index so targets are only searched in one direction and not two times
        Range<std::vector<Edge_ch>::iterator> targets(contracted_node->edges.begin() + i+1, contracted_node->edges.end());
        for (auto& target_w : targets) {
            double d_sc = d_u_v + target_w.weight; ///shortcut distance
            Node_ch2* w = g->getNode(target_w.target);

            auto start_time2 = std::chrono::high_resolution_clock::now();
            //check if edge already exists
            bool exists = false;
            Range<std::vector<Edge_ch>::iterator> range_u(u->edges.begin() + u->second_index, u->edges.end());
            for (auto& neighbour_u : range_u) {
                if (neighbour_u.target == w->id) {
                    exists = true;
                    if (d_sc < neighbour_u.weight) {
                        neighbour_u.weight = d_sc;
                        
                        Range<std::vector<Edge_ch>::iterator> range_w(w->edges.begin() + w->second_index, w->edges.end());
                        for (auto& neighbour_w : range_w) {
                            if (neighbour_w.target == u->id) {
                                neighbour_w.weight = d_sc;
                                break;
                            }
                        } 
                        break;
                    }
                    
                }
            }
            auto end_time2 = std::chrono::high_resolution_clock::now();
            time2 += (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_time2- start_time2).count();
            if (!exists) {
                auto start_time3 = std::chrono::high_resolution_clock::now();
                double d_dijkstra = localDijkstra(g, w->id, d_sc, Pmax, p);
                auto end_time3 = std::chrono::high_resolution_clock::now();
                time3 += (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_time3- start_time3).count();
                if (d_dijkstra > d_sc && d_dijkstra > 0) {
                    
                    g->addShortcut(u->id, w->id, d_sc);
                    num_sc++;
                    //std::cout << "add shortcut " << std::endl;
                }
            }
        }
        delete p;
    }   
}


/*void contractNode(Graph_ch2* g, Node_ch2* node) { //node is v, the node that gets contracted
    node->level = g->getLevel();
    g->increaseLevel();
    Range<std::vector<Edge_ch>::iterator> range_edges(node->edges.begin() + node->second_index, node->edges.end());

    //change second index of all adjacent nodes with higher level
    for (auto& edge : range_edges) {

        Node_ch2* u = g->getNode(edge.target);
        Range<std::vector<Edge_ch>::iterator> range_targets(u->edges.begin() + u->second_index, u->edges.end());

        for (auto& targets : range_targets) {
            if (targets.target == node->id) {
                std::swap(targets, u->edges[u->second_index]);
                u->second_index++;
                break;
            }
        }
    }

    for (std::vector<Edge_ch>::iterator it = node->edges.begin() + node->second_index; it != node->edges.end(); it++) {
        Node_ch2* u = g->getNode(it->target);
        //increase contracted neighbour counter and adjacencies counter 
        u->cnc++;
        localsearch(g, u, node, std::distance(node->edges.begin(), it)); //std::distance(node->edges.begin(), it) := edge index u->nodebt
    }
}*/


void contractNode(Graph_ch2* g, Node_ch2* node, Adr_PriorityQueue<Node_ch2*, int>* p) {
    node->level = g->getLevel();
    g->increaseLevel(); 

    //update importance and second index of all neighbours
    auto start_time4 = std::chrono::high_resolution_clock::now();
    Range<std::vector<Edge_ch>::iterator> neighbours(node->edges.begin() + node->second_index, node->edges.end());
    for (auto& edge : neighbours) {
        Node_ch2* neighbour = g->getNode(edge.target);
        neighbour->cnc++;

        Range<std::vector<Edge_ch>::iterator> neighbours_edges(neighbour->edges.begin() + neighbour->second_index, neighbour->edges.end());
        for (auto& n_edge : neighbours_edges) {
            if (n_edge.target == node->id) {
                std::swap(n_edge, neighbour->edges[neighbour->second_index]);
                neighbour->second_index++;
                break;
            }
        }
    }
    auto end_time4 = std::chrono::high_resolution_clock::now();
    time4 += (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_time4- start_time4).count();
}



/*void contractionHierarchy2(Graph_ch2* g, int source, int target) {
    std::cout << "contraction hierarchy" << std::endl;
    d.resize(g->getNumNodes());

    //create pq and insert elements
    Adr_PriorityQueue<Node_ch2*, int> p;
    Range<std::vector<Node_ch2>::iterator> nodes(g->getBegin(), g->getEnd());
    for (auto& node : nodes) {
        calculateImportance(g, &node);
        p.insert(&node, node.importance);
    }

    Node_ch2* u = nullptr;
    int i = 0;
    while (!p.isEmpty()) {

        u = p.top();
        Node_ch2 u_ = *u; //copy

        calculateImportance(g, u);
        //lazy update
        int trigger = 0;
        while(u_.importance != u->importance) {
            //if too many successfull lazy updates, update all remaining nodes
            if (trigger == 10) {
                for (int k = 0; k < p.size(); k++) {
                    calculateImportance(g, g->getNode(k));
                    p.updateKey(g->getNode(k), g->getNode(k)->importance);
                    trigger = 0;
                }
                break;
            }

            p.updateKey(p.top(), u->importance);
            trigger++;
            //j++;
            u = p.top();
            u_ = *u;
            calculateImportance(g, u);
        }

        u = p.deleteMin();
        //std::cout << "contract node: " << u->id << std::endl;
        contractNode(g, u);

        //update all direct neighbours
        Range<std::vector<Edge_ch>::iterator> range_edges(u->edges.begin() + u->second_index, u->edges.end());
        for (auto& edge : range_edges) {
            Node_ch2* node = g->getNode(edge.target);
            calculateImportance(g, node);
            p.updateKey(node, node->importance);
        }
        i++;
        if (i % 1000 == 0) {
            std::cout << i << std::endl;
            std::cout << "num_sc: " << num_sc << std::endl;
            num_sc = 0;
        }
    } 
}*/

void contractionHierarchy(Graph_ch2* g) {
    std::cout << "start contraction hierarchy" << std::endl;
    d.resize(g->getNumNodes());
    Adr_PriorityQueue<Node_ch2*, int>* p = new Adr_PriorityQueue<Node_ch2*, int>;
    Range<std::vector<Node_ch2>::iterator> nodes(g->getBegin(), g->getEnd());
    for (auto& node : nodes) {
        calculateImportance(g, &node);
        p->insert(&node, node.importance);
    }
    int i = 0;
    while (!p->isEmpty()) {
        Node_ch2* u = p->deleteMin();
        contractNode(g, u, p);
        i++;
        auto start_time_localsearch = std::chrono::high_resolution_clock::now();
        localsearch(g, u);
        auto end_time_localsearch = std::chrono::high_resolution_clock::now();
        time_ls += (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_time_localsearch - start_time_localsearch).count();

        Range<std::vector<Edge_ch>::iterator> neighbours(u->edges.begin() + u->second_index, u->edges.end());
        for (auto& edge : neighbours) {
            Node_ch2* neighbour = g->getNode(edge.target);
            neighbour->cnc++;
            calculateImportance(g, neighbour);
            p->updateKey(neighbour, neighbour->importance);
        }

        if (i % 1000 == 0) {
        std::cout << i << std::endl;
        std::cout << "num_sc: " << num_sc << std::endl;
        std::cout << "time_ls: " << time_ls << std::endl;
        std::cout << "time1: " << time1 << std::endl;
        std::cout << "time2: " << time2 << std::endl;
        std::cout << "time3: " << time3 << std::endl;
        std::cout << "time4: " << time4 << std::endl;
        std::cout << "timebtw: " << time_btw << std::endl;

        time_ls = 0;
        time3 = 0;
        time2 = 0;
        time4 = 0;
        time_btw = 0;
        num_sc = 0;
        }
    }
    
}
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "graphs.h"
#include "minheap.h"
#include <algorithm> //for swap of vector elements
#include <cstdlib>
#include <chrono> //for timing



using namespace std;
using namespace std::chrono;
using namespace graphs;
using namespace heap; //no overlap




void heapTests();
void graphTests();
void vectorTests();


vector<edge> prim(Graph og) //Prim's Algorithm implementation
{
    int n = og.numvertices;
    double * distances;
    int * S;
    distances = new double[n]; //despite fixed size, arrays must be dynamically allocated on heap to prevent stack overflow for large n
    S = new int[n];
    //S[v] = 1 means v is in S; 0 means otherwise
    for (int i = 0; i < n; i++){
        distances[i] = INT_MAX;
        S[i] = 0;
    }
    vector<edge> prev;
    prev.reserve(n);
    for (int k = 0; k < n; k++){
        prev.push_back({0,0,0.0}); //dummy edge, to be replaced for all v > 0 by Prim
    }
    MinHeap heap(n);

    heap.insrt({0, 0.0}); //vertices as integers, 0-indexed
    distances[0] = 0;

    while (heap.H.size() > 0){
        heapNode v = heap.extract_min();
        S[v.v] = 1;
        for (unsigned int w = 0; w < og.adj[v.v].size(); w++){


            edge ew = og.adj[v.v][w]; //ew is an edge in adj[v], which is a vector
            if (S[ew.u] == 0){ //u \in V\S
                if (distances[ew.u] > ew.elen){
                    distances[ew.u] = ew.elen;
                    prev[ew.u] = ew;
                    heap.insrt({ew.u, distances[ew.u]});
                }
            }
        }
    }
    delete distances;
    delete S;
    //cout << "finished prim, about to return \n";
 return prev;
}

int main(int argc, char* argv[])
{
    srand (time (0)); //seeding RAND for random graphs

//    for (int i = 0; i < argc; i++){
//        cout << "\n" << argv[i];
//    }


    int n = atoi(argv[2]);
    int dim = atoi(argv[4]);
    int trials = atoi(argv[3]);

    double avgTotalWeight = 0;
    double avgMaxEdge = 0;
    double avgTime = 0;


    Graph graph(n);

    for (int i = 0; i < trials; i++){
        graph.genGraph(dim); //see graph.cpp - new graph generated each time

        auto start = high_resolution_clock::now();
        vector<edge> ms = prim(graph);
        auto stop = high_resolution_clock::now(); //measures time elapsed for Prim's algorithm; timing function from https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
        auto duration = duration_cast<microseconds>(stop - start);

        avgTime += duration.count() / trials;


        double maxWeight = 0.0;
        double totalWeight = 0;
        for (unsigned int i = 1; i < ms.size(); i++){

            //std::cout << "\n Edge from " << ms[i].v << " to " << ms[i].u << " of length " << ms[i].elen << "\n";
            totalWeight += ms[i].elen;
            if (ms[i].elen > maxWeight){
                maxWeight = ms[i].elen;
            }
        }
        avgTotalWeight += totalWeight / double(trials);
        avgMaxEdge += maxWeight / double(trials);
    }

    //cout << "\n average total weight: " << avgTotalWeight << ", Average Max Weight: " << avgMaxEdge << " Prim execution time: " << avgTime << "\n";

    cout << avgTotalWeight << " " << n << " " << trials << " " << dim;
    return 0;
}

void heapTests()
{
        cout << "\n";

    MinHeap heap(10); //ANY VALUE GREATER THAN 6 MAKES THIS BREAK...WHY???
    heap.insrt({0, .2});
    cout << "\n insert {0, .2}: ";
    //heap.printHeap();
    cout << "\n insert {1, .42}: ";
    heap.insrt({1, .42});
    //heap.printHeap();
    cout << "\n insert {2, .1}: ";
    heap.insrt({2, .1});
    //heap.printHeap();
//    heap.insrt({3, .01});
//    heap.insrt({4, .6});
//    heap.insrt({5, .0001});
    heap.insrt({6, .2});
//    heap.insrt({7, .42});
//    heap.insrt({8, .1});
//    heap.insrt({9, .51});
//    heap.insrt({10, .12});
//    heap.insrt({11, .01});
//    heap.insrt({12, .6});
//    heap.insrt({13, .0001});
cout<<"FINAL: \t";
    //heap.printHeap();


//    heap.printHeap();
//
//    heap.insrt({2, 0.001});
//    heap.printHeap();
//
//    heapNode n;
//
//    for (int i = 0; i < heap.H.size(); i++){
//        n = heap.extract_min();
//        cout << n.v << ", " << n.dist << "\n";
//    }
}

void graphTests()
{


 //push_back doubles vector size if not enough room - same as realloc * 2

  Graph graph1(7);
  graph1.adj[0].push_back({0,1,0.5});
  graph1.printGraph();

    Graph g1(4);
    g1.genGraph(4);
    g1.printGraph();
    std::cout << "\n";
    //g1.getWeights();
  //delete(g1);
}

void vectorTests()
{
    vector<int> test;
    test.resize(10);
    for (int i = 0; i < 10; i++){
        test[i] = i+1;
    }
    cout << "initial: \n";
    for (int j = 0; j < 10; j++){
    cout << test[j] << "\t";
    }
    std::iter_swap(test.begin() + 3, test.begin() + 8);
    cout << "\n Final: \n";
    for (int j = 0; j < 10; j++){
    cout << test[j] << "\t";
    }

    test.pop_back();
    test.shrink_to_fit();
    cout << "\n Final: \n";
    for (int j = 0; j < 16; j++){
    cout << test[j] << "\t";
    }
}

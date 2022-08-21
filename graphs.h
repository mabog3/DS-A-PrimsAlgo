#ifndef GRAPHS_H_INCLUDED
#define GRAPHS_H_INCLUDED
#include <vector>
#include <string>
#include <utility>


using namespace std;
namespace graphs
{
    struct edge
    {
        int v;
        int u;
        double elen;
    };

    class Graph //implementations in graphs.cpp
    {
    public:
        vector<vector<edge>> adj;  //adjacency list representation; outside vector holds all vertices, and adj[v] is a vector of all edges from v
        int numvertices; // number of vertices

        Graph(int numvert); //constructor - implemented in graphs.cpp

        vector<double> getWeights(); //get weight of all edges in graph and weight of largest edge in graph
        void genGraph(int dimension); //changes the empty graph object of size n into a populated complete graph (minus too-large edges)
        void printGraph(); //testing purposes
    };
}

#endif // GRAPHS_H_INCLUDED

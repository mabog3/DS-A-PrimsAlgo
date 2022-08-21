#include <vector>
#include <iostream>
#include <utility>
#include "graphs.h"
#include <cmath>
#include <math.h>       /* exp */


using namespace graphs;

Graph::Graph(int numvert)
{
    numvertices = numvert;
    adj.resize(numvert); //pre-allocate memory & empty edge vectors for all of the vertices
}

void Graph::printGraph()
{
    for (int i = 0; i < numvertices; i++){
        std::cout << "\n Edge from " << i << " to: \n";
        for (unsigned int v = 0; v < adj[i].size(); v++) {
            cout << "\t " << adj[i][v].u << " has length " << adj[i][v].elen;
        }
    }
}

void Graph::genGraph(int dimension) //assumes dimension is 0,2,3, or 4
{
    if (adj[0].size() > 0){//if a graph was already generated
        for (int i = 0; i < numvertices; i++){
            adj[i].clear();
            //adj[i].shrink_to_fit() -- PERFORMANCE???
        }
    }
    int dim = dimension;
    if (dim == 0){
        dim = 1; //idk guys, a line from 0 to 1 seems 1 dimensional to me, not 0 dimensional
    }

    double *x1 = new double[numvertices](); //these arrays represent the coordinates (x1, x2, x3, x4) in R^{dim}; value-initialized to 0
    double *x2 = new double[numvertices]();
    double *x3 = new double[numvertices]();
    double *x4 = new double[numvertices]();

    for (int i = 0; i < numvertices; i++){ //give each vertex a valid location - at first wanted to just generate edges, but then could have had
            // impossible graph, like (a,b) = .9, (b,c) = .01, (c, a) = .01
            x1[i] = double(std::rand()) / double(RAND_MAX);
            if (dim > 1){
                x2[i] = double(std::rand()) / double(RAND_MAX); //R^k is isomorphic to the subspace of R^n, n - k = d > 0, where d coordinates of R^n are 0
            }
            if (dim > 2){
                x3[i] = double(std::rand()) / double(RAND_MAX);
            }
            if (dim > 3){
                x4[i] = double(std::rand()) / double(RAND_MAX);
            }

    }

    double maxLen = 1;


    if (numvertices > 5000){
        if (dim == 1){
            maxLen = exp(-0.000006*numvertices - 5);
        }

        else if (dim == 2){
            maxLen = exp(-0.0000060965*numvertices - 3.7);
        }
        else if (dim == 3){
            maxLen = exp(-0.0000041965*numvertices - 2.4);
        }
        else{
            maxLen = exp(-0.00000592*numvertices*(.5) - 1.71);
        }
    }



    //cout << "maxlen: " << maxLen;
    for (int v = 0; v < numvertices; v++){
        for (int u = 0; u < v; u++) {
            double templen = sqrt(pow(x1[v] - x1[u], 2) + pow(x2[v] - x2[u], 2) + pow(x3[v] - x3[u], 2) + pow(x4[v] - x4[u], 2) );
            if (templen < maxLen){
                adj[v].push_back({v, u, templen});
                adj[u].push_back({u, v, templen});
            }

        }
    }


    delete x1; //deallocate location array memory
    delete x2;
    delete x3;
    delete x4;
}



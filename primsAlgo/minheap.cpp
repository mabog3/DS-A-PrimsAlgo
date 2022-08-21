#include "minheap.h"
#include <vector>
#include <algorithm> //for swap of vector elements
#include <iostream>


using namespace heap;

//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_ALWAYS_DF );


MinHeap::MinHeap(int numv)
{
    H.reserve(numv/4); //for performance reasons - reserve enough space for a large amount of vertices ahead of time
    posit = new int[numv];
    numvert = numv;
    for (int i = 0; i < numvert; i++){
        posit[i] = -1;
    }
}

int MinHeap::left(int i)
{
    return 2*i + 1; // because actual arrays are 0-indexed, not 1-indexed
}

int MinHeap::right(int i)
{
    return 2*i + 2;
}

int MinHeap::parent(int i)
{
    return (i - 1)/2; //1-indexed would be i/2
}

heapNode MinHeap::peek()
{
    return H[0]; //top element of the actual heap
}

void MinHeap::swapN(int a, int b)
{
    posit[H[a].v] = b;
    posit[H[b].v] = a; //swap positions in posit
    std::iter_swap(H.begin() + b, H.begin() + a);
}

void MinHeap::min_heapify(int N)
{
    unsigned int l = left(N);
    unsigned int r = right(N);

    int smallest = N;

    if (l < H.size()){
        if (H[l].dist < H[N].dist){
              smallest = l;
            }
        } //if Left exists & is smaller than N
    if (r < H.size()){
            if (H[r].dist < H[smallest].dist){
                smallest = r;
            }

    }
    if (smallest != N){
        swapN(smallest, N);
        min_heapify(smallest);
    }

}

heapNode MinHeap::extract_min()
{
    heapNode m = H[0];
    H[0] = H[H.size() - 1]; //last element moves to first
    posit[H[0].v] = 0; //updte position of last node to 0
    H.resize(H.size() - 1); //deletes last element in vector - DOES NOT ACTUALLY RELEASE VALUE FROM MEMORY, SO BEWARE MEM LEAK!!!
//    H.shrink_to_fit();
    posit[m.v] = -1; //node no longer in heap
    min_heapify(0);
    return m;
}

void MinHeap::insrt(heapNode node)
{
    if (posit[node.v] == -1) {//node not already in heap
        H.push_back(node); //adds new node to back; O(1) unless reallocation, so equivalent to checking & doubling array size via malloc
        int N = H.size() - 1;
        posit[node.v] = N;
        while (N != 0 && H[parent(N)].dist > H[N].dist){
            swapN(parent(N), N);
            N = parent(N);
        }
    }
    else{ //node already in heap - update its dist value
//        std::cout << "\n Found a repeat: v = " << node.v << "currently at H[" << << "\n";
        int N = posit[node.v];
        if (node.dist < H[N].dist){
            H[posit[node.v]].dist = node.dist;
            while (N && H[parent(N)].dist > H[N].dist){
                swapN(parent(N), N);
                N = parent(N);
            }

        }
    }

}

void MinHeap::printHeap()
{
    std::cout << "\n";
    for (unsigned int i = 0; i < H.size(); i++)
    {
        std::cout << i << ": {" << H[i].v << ", " << H[i].dist << " }" << "\t";
    }
    for (int v = 0; v < numvert; v++)
    {
        std::cout << "v" << v << " at H[" << posit[v] << "]; " ;
    }
}


MinHeap::~MinHeap() //destructor
{
    H.clear();
    H.shrink_to_fit();
    std::vector<heapNode>().swap(H);
    //std::cout << H.size();
    delete posit;
}



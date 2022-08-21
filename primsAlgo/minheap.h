#ifndef MINHEAP_H_INCLUDED
#define MINHEAP_H_INCLUDED
#include <vector>

using namespace std;
namespace heap
{
    struct heapNode
    {
        int v; //vertex
        double dist; //'distance,' for use in Prim's algo
    };

    class MinHeap //implementation in minheap.cpp
    {
    public:

        vector<heapNode> H;
//        int len;
//        int capacity;
        int* posit;
        int numvert;

        MinHeap(int numv);


        int left(int i);
        int right(int i);
        int parent(int i);
        void swapN(int a, int b); //swaps node at heap index a with node at heap index b

        void min_heapify(int N); //index of "Node N" from algorithm
        heapNode peek();
        heapNode extract_min();
        void insrt(heapNode node);

        ~MinHeap();
        void printHeap(); //testing purposes

    };
}

#endif // MINHEAP_H_INCLUDED

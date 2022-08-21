#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono> //for timing

#define lint long long int
#define max_iter 25000

using namespace std;
using namespace std::chrono;

struct heapNode
    {
        int v; //vertex
        lint dist; //'distance,' for use in Prim's algo
    };

class MaxHeap
    {
    public:

        vector<heapNode> H;
//        int len;
//        int capacity;
        int* posit;
        int numvert;

        MaxHeap(int numv);


        int left(int i);
        int right(int i);
        int parent(int i);
        void swapN(int a, int b); //swaps node at heap index a with node at heap index b

        void max_heapify(int N); //index of "Node N" from algorithm
        heapNode peek();
        heapNode extract_max();
        void insrt(heapNode node);

        ~MaxHeap();
        void printHeap(); //testing purposes

    };



MaxHeap::MaxHeap(int numv)
{
    H.reserve(numv);
    posit = new int[numv];
    numvert = numv;
    for (int i = 0; i < numvert; i++){
        posit[i] = -1;
    }
}

int MaxHeap::left(int i)
{
    return 2*i + 1; // because actual arrays are 0-indexed, not 1-indexed
}

int MaxHeap::right(int i)
{
    return 2*i + 2;
}

int MaxHeap::parent(int i)
{
    return (i - 1)/2; //1-indexed would be i/2
}

heapNode MaxHeap::peek()
{
    return H[0]; //top element of the actual heap
}

void MaxHeap::swapN(int a, int b)
{
    posit[H[a].v] = b;
    posit[H[b].v] = a; //swap positions in posit
    std::iter_swap(H.begin() + b, H.begin() + a);
}

void MaxHeap::max_heapify(int N)
{
    unsigned int l = left(N);
    unsigned int r = right(N);

    int largest = N;

    if (l < H.size()){
        if (H[l].dist > H[N].dist){
              largest = l;
            }
        } //if Left exists & is smaller than N
    if (r < H.size()){
            if (H[r].dist > H[largest].dist){
                largest = r;
            }

    }
    if (largest != N){
        swapN(largest, N);
        max_heapify(largest);
    }

}

heapNode MaxHeap::extract_max()
{
    heapNode m = H[0];
    H[0] = H[H.size() - 1]; //last element moves to first
    posit[H[0].v] = 0; //updte position of last node to 0
    H.resize(H.size() - 1); //deletes last element in vector - DOES NOT ACTUALLY RELEASE VALUE FROM MEMORY, SO BEWARE MEM LEAK!!!
//    H.shrink_to_fit();
    posit[m.v] = -1; //node no longer in heap
    max_heapify(0);
    return m;
}

void MaxHeap::insrt(heapNode node)
{
    if (true){//(posit[node.v] == -1) {//node not already in heap
        H.push_back(node); //adds new node to back; O(1) unless reallocation, so equivalent to checking & doubling array size via malloc
        int N = H.size() - 1;
        posit[node.v] = N;
        while (N != 0 && H[parent(N)].dist < H[N].dist){
            swapN(parent(N), N);
            N = parent(N);
        }
}
}


MaxHeap::~MaxHeap() //destructor
{
    H.clear();
    H.shrink_to_fit();
    std::vector<heapNode>().swap(H);
    //std::cout << H.size();
    delete posit;
}


lint kkarp(int n, const vector<lint> & A){//returns residue

    MaxHeap heap(n);

    for (int i = 0; i < n; i++){
        heap.insrt({i, A[i]}); //first value is the original position in A, 2nd value is the actual num in A
    }

//    cout << "testing heap: \n";
//    for (int i = 0; i < n; i++){
//        heapNode currMax = heap.extract_max();
//        cout << currMax.dist << "\n";
//    }
//
    while(true){
        heapNode max1 = heap.extract_max();
        heapNode max2 = heap.extract_max();

        if (max2.dist == 0){
            //cout << "residue: " << max1.dist << "\n";
            return max1.dist;
        }

        lint diff = llabs(max1.dist - max2.dist);
        heap.insrt({max1.v, diff});
        heap.insrt({max2.v, 0});
    }

    return 0;
}


lint solResidue(int n, vector<int> & sol, const vector<lint> & A, bool part = false){
    lint residue = 0;

    if (!part){
        for (int i = 0; i < n; i++){
            residue = residue + (A[i] * sol[i]);
        }
        return llabs(residue); //absolute value of residue
    }

    else{
        vector<lint> Aprime(n);

        for (int j = 0; j < n; j++){
            Aprime[sol[j]] = Aprime[sol[j]] + A[j];
        }

        return kkarp(n, Aprime);
    }


}

void randomSolution(int n, vector<int> & sol, bool part = false) //changes input array to a random solution of size n
{
    if (!part){
        double pr;
        for (int i = 0; i < n; i++){
            pr = double(rand()) / double(RAND_MAX);
            if (pr < 0.5){
                sol[i] = 1;
            }
            else{
                sol[i] = -1;
            }
        //cout << sol[i] << "\n";
        }
    }
    else {
        for (int i = 0; i < n; i++){
            sol[i] = rand() % n;
        }
    }
}


double cool(int iter) {
    return pow(10, 10) * pow(0.8, iter / 300); //cooling function for simulated annealing
}

void arrcpy(int n, vector<int> & source, vector<int> & target) //copies source into target
{
    for (int i = 0; i < n; i++){
        target[i] = source[i];
    }
}

void genNeighbor(int n, vector<int> & sol, vector<int> & neighbor, bool part = false) //generates a neighbor of sol
{
    arrcpy(n, sol, neighbor); //copies sol into neighbor

    int idx1 = rand() % n;  //get first index
    int idx2 = idx1;

    if (!part){
        neighbor[idx1] = -1 * sol[idx1];

        if (double(rand()) / double(RAND_MAX) > 0.5){
            idx2 = idx1;
            while (idx1 == idx2){
                idx2 = rand() % n; //ensures new idx2 is picked if its the same as idx1
            }
            neighbor[idx2] = -1* sol[idx2];
        }
    }


    else{
        while (idx1 == idx2){
            idx2 = rand() % n; //ensures new idx2 is picked if its the same as idx1
        }
        neighbor[idx1] = idx2;
    }
}


lint rrandom(int n, const vector<lint> & A, bool part = false){
    vector<int> sol(n);
    //cout << "part: " << part << "\n";

    lint currRes;
    lint newRes;

    randomSolution(n, sol, part);

    currRes = solResidue(n, sol, A, part);
    //cout << "1st res: " << currRes << "\n";

    for (int i = 1; i < max_iter; i++){
        randomSolution(n, sol, part);
        newRes = solResidue(n, sol, A, part);
        if (newRes < currRes){
            currRes = newRes;
            //cout << "New res: " << currRes << "\n";
        }
    }
    return currRes;
}

lint hill(int n, const vector<lint> &A, bool part = false){
    vector<int> sol(n);
    vector<int> neighbor(n);

    lint currRes;
    lint newRes;


    randomSolution(n, sol, part);
    currRes = solResidue(n, sol, A, part);

    for (int i = 0; i < max_iter; i++){
        genNeighbor(n, sol, neighbor, part);
        newRes = solResidue(n, neighbor, A, part);

        if (newRes < currRes){
            arrcpy(n, neighbor, sol);
            currRes = newRes;
            //cout << "New res: " << currRes << "\n";
        }
    }


    return currRes;

}

lint simanneal(int n, const vector<lint> & A, bool part = false)
{
    vector<int> sol(n);
    vector<int> bestSol(n);
    vector<int> neighbor(n);

    lint bestRes;
    lint currRes;
    lint newRes;

    randomSolution(n, sol, part);
    arrcpy(n, sol, bestSol);


    currRes = solResidue(n, sol, A, part);
    bestRes = solResidue(n, bestSol, A, part);

    for (int i = 0; i < max_iter; i++){

        genNeighbor(n, sol, neighbor, part); // sets neighbor ot a random neighbor of current sol

        currRes = solResidue(n, sol, A, part);
        newRes = solResidue(n, neighbor, A, part);


        if (newRes < currRes){
            arrcpy(n, neighbor, sol);
            currRes = newRes;
            //cout << "New res: " << currRes << "\n";
        }
        else {
            long double pr = exp(-1 * (newRes - currRes) / cool(i));
            if (double(rand()) / double(RAND_MAX) < pr){
               arrcpy(n, neighbor, sol);
               currRes = newRes;
            }
        }

        if (currRes < bestRes){
            arrcpy(n, sol, bestSol);
            bestRes = currRes;
        }
    }

    return currRes;
}


void genRandomInstance(int n, vector<lint> &A)
{
    lint x = 10e12; //code for generating random large ints adapted from https://stackoverflow.com/questions/37396278/how-to-generate-very-large-random-number-in-c

    random_device dev;
    default_random_engine generator(dev());
    uniform_int_distribution<long long unsigned> distribution(0,0xFFFFFFFFFFFFFFFF);

    for (int i; i < n; i++){
        A[i] = distribution(generator) % (x + 1);
    }
}

int main(int argc, char *argv[])
{
	if (argc != 4) {
        //cout << "/partition flag algorithm inputfile\n";
        return 1;
	}
    srand (time (0)); //seeding RAND

	int algo = stoi(argv[2]); //stoi may cause problems with older c++
	int test = stoi(argv[1]); //optional toggle for tests

	//cout << "algo: " << algo << ", test: " << test << "\n";

	int n = 0; //size of input list

	string inpt = argv[3]; //input file name
	ifstream input;
	input.open(inpt);
	string inputString;
	vector<lint> A; //padded with 0s to nearest power of 2;

	if ( input.is_open() ) { // reads file, adds to vector, adds up dimension
        while (getline(input, inputString)) {
            A.push_back(stoll(inputString)); //adds line as a long long int to A
            n = n + 1;
        }
    }

    //cout << "n: " << n << "size: " << A.size() << "\n";

    if (test == 1){

        vector<double> times(7); //karp, rand, hill, sa, prand, phill, psa
        vector<vector<lint>> avgResidues(7, vector<lint> (51, 0));;//last spot holds avg

        n = 100;
        vector<lint> B(n);
        genRandomInstance(n, B);

        lint res = 0;

        for (int i = 0; i < 50; i++)
        {

            genRandomInstance(n, B);

            auto start = high_resolution_clock::now();
            res = kkarp(n, B);
            auto stop = high_resolution_clock::now(); //measures time elapsed for Strassen's algorithm; timing function from https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            auto duration = duration_cast<milliseconds>(stop - start);
            times[0] += duration.count() / 50;
            avgResidues[0][i] = res;
            avgResidues[0][50] += res;

            auto start1 = high_resolution_clock::now();
            res = rrandom(n, B);
            auto stop1 = high_resolution_clock::now(); //measures time elapsed for Strassen's algorithm; timing function from https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            auto duration1 = duration_cast<milliseconds>(stop1 - start1);
            times[1] += duration1.count() / 50;
            avgResidues[1][i] = res;
            avgResidues[1][50] += res;

            auto start2 = high_resolution_clock::now();
            res = hill(n, B);
            auto stop2 = high_resolution_clock::now(); //measures time elapsed for Strassen's algorithm; timing function from https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            auto duration2 = duration_cast<milliseconds>(stop2 - start2);
            times[2] += duration2.count() / 50;
            avgResidues[2][i] = res;
            avgResidues[2][50] += res;

            auto start3 = high_resolution_clock::now();
            res = simanneal(n, B);
            auto stop3 = high_resolution_clock::now(); //measures time elapsed for Strassen's algorithm; timing function from https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            auto duration3 = duration_cast<milliseconds>(stop3 - start3);
            times[3] += duration3.count() / 50;
            avgResidues[3][i] = res;
            avgResidues[3][50] += res;

            auto start4 = high_resolution_clock::now();
            res = rrandom(n, B, true);
            auto stop4 = high_resolution_clock::now(); //measures time elapsed for Strassen's algorithm; timing function from https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            auto duration4 = duration_cast<milliseconds>(stop4 - start4);
            times[4] += duration4.count() / 50;
            avgResidues[4][i] = res;
            avgResidues[4][50] += res;

            auto start5 = high_resolution_clock::now();
            res = hill(n, B, true);
            auto stop5 = high_resolution_clock::now(); //measures time elapsed for Strassen's algorithm; timing function from https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            auto duration5 = duration_cast<milliseconds>(stop5 - start5);
            times[5] += duration5.count() / 50;
            avgResidues[5][i] = res;
            avgResidues[5][50] += res;

            auto start6 = high_resolution_clock::now();
            res = simanneal(n, B, true);
            auto stop6 = high_resolution_clock::now(); //measures time elapsed for Strassen's algorithm; timing function from https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            auto duration6 = duration_cast<milliseconds>(stop6 - start6);
            times[6] += duration6.count() / 50;
            avgResidues[6][i] = res;
            avgResidues[6][50] += res;



        }

        cout << "Average Times: ";
            for (int j = 0; j < 7; j++){
                cout << times[j] << "\n";
            }

            for (int j = 0; j < 7; j++){
                cout << "Algo " << j << "\n";
                for (int k = 0; k < 50; k++){
                    cout << avgResidues[j][k] << ", ";
                }
                cout << "\n Avg: " << avgResidues[j][50] / 50 << "\n\n";
            }
        //lint x = 10e12;

//        for (int i = 0; i < n; i++){
//            lint diff = x - B[i];
//            cout << "B[" << i << "] = " << B[i] << ", diff = " << diff << "\n";
//        }
    }

    if (test == 2){
            vector<lint> avg(50);
            //lint res = 0;
            vector<lint> B(n);
            genRandomInstance(n, B);

            for (int i = 0; i < 50; i++){
                genRandomInstance(n, B);
                cout << simanneal(n, B, true) << ", ";
            }
        }

    if (algo == 0){
        int karp = kkarp(n, A);
        cout << karp;

    }

    if (algo == 1){
        lint rrand = rrandom(n, A);
        cout << rrand;
    }

    if (algo == 2){
        lint hll = hill(n, A);
        cout << hll;
    }

    if (algo == 3){
        lint sa = simanneal(n, A);
        cout << sa;
    }

    if (algo == 11){
        cout << rrandom(n, A, true);
    }

    if (algo == 12){
        cout << hill(n, A, true);
    }
    if (algo == 13){
        cout << simanneal(n, A, true);
    }

	return 0;
}

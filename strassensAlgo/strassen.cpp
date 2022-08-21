
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <chrono> //for timing

using namespace std;
using namespace std::chrono;

void printMatrix(const vector<vector<int>> & A, int dim = 0){ //for testing purposes
    if (dim == 0){
        dim = A.size();
    }
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << "\n";
    }

}

vector<vector<int>> randMatrix(int dim, int range){ //testing purposes; creates a random square matrix with entires 0-range &
    vector<vector<int>> reslt( dim , vector<int> (dim)); //empty vector is initialized to 0's by itself
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            reslt[i][j] = rand() % range;
        }
    }
    return reslt;
}



void printDiagonals(const vector<vector<int>>  &A, int dim = 0){ //for returning final answer
    if (dim == 0){
        dim = A.size();
    }
    for (int i = 0; i < dim; i++){
            cout << A[i][i] << "\n";
    }
}


void basicMult(const vector<vector<int>> & A, const vector<vector<int>> & B, vector<vector<int>> &reslt, int dim){
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            for (int k = 0; k < dim; k++){
                reslt[i][j] = reslt[i][j] + (A[i][k] * B[k][j]);
            }
        }
    }

    return;
}

void addMat(const vector<vector<int>> & Mat1, const vector<vector<int>> & Mat2, vector<vector<int>> & reslt, int dim){

    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            reslt[i][j] = Mat1[i][j] + Mat2[i][j];
        }
    }
    return;
}

void subtractMat(const vector<vector<int>> & Mat1, const vector<vector<int>> & Mat2, vector<vector<int>> & reslt, int dim){

    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            reslt[i][j] = Mat1[i][j] - Mat2[i][j];
        }
    }
    return;
}


void strassen(const vector<vector<int>> MAT1, const vector<vector<int>> & MAT2, vector<vector<int>> & reslt, int dim, int padStart, int currRow1 = 0, int currCol1 = 0, int currRow2 = 0, int currCol2 = 0){

    if (currCol1 > padStart && currRow1 > padStart){//one of the input vectors is all 0s - no further operations are needed
        for (int i = 0; i < dim; i++){
            for (int j = 0; j < dim; j++){
                reslt[i][j] = 0;
            }
        }
    }
    else if (currCol2 > padStart && currRow2 > padStart){//one of the input vectors is all 0s - no further operations are needed
        for (int i = 0; i < dim; i++){
            for (int j = 0; j < dim; j++){
                reslt[i][j] = 0;
            }
        }
    }


    else if (dim <= 32){ //threshold
        basicMult(MAT1, MAT2, reslt, dim);
    }
else{
    int dim2 = dim / 2; //new dim; half of old dim

    vector<vector<int>> A( dim2 , vector<int> (dim2)); //MAT1 upper left
    vector<vector<int>> B( dim2 , vector<int> (dim2)); //MAT1 upper right
    vector<vector<int>> C( dim2 , vector<int> (dim2)); //MAT1 lower left
    vector<vector<int>> D( dim2 , vector<int> (dim2)); //MAT1 lower right
    vector<vector<int>> E( dim2 , vector<int> (dim2)); //MAT2 upper left
    vector<vector<int>> F( dim2 , vector<int> (dim2)); //MAT2 upper right
    vector<vector<int>> G( dim2 , vector<int> (dim2)); //MAT2 lower left
    vector<vector<int>> H( dim2 , vector<int> (dim2)); //MAT2 lower right

    for (int i = 0; i< dim2; i++){
        for (int j = 0; j < dim2; j++){
            A[i][j] = MAT1[i][j];
            E[i][j] = MAT2[i][j];

            B[i][j] = MAT1[i][j + dim2];
            F[i][j] = MAT2[i][j + dim2];

            C[i][j] = MAT1[i + dim2][j];
            G[i][j] = MAT2[i + dim2][j];

            D[i][j] = MAT1[i + dim2][j + dim2];
            H[i][j] = MAT2[i + dim2][j + dim2];
        }
    }

    vector<vector<int>> P1( dim2 , vector<int> (dim2));
    vector<vector<int>> P2( dim2 , vector<int> (dim2));
    vector<vector<int>> P3( dim2 , vector<int> (dim2));
    vector<vector<int>> P4( dim2 , vector<int> (dim2));
    vector<vector<int>> P5( dim2 , vector<int> (dim2));
    vector<vector<int>> P6( dim2 , vector<int> (dim2));
    vector<vector<int>> P7( dim2 , vector<int> (dim2));

    vector<vector<int>> temp1( dim2 , vector<int> (dim2)); //vectors of integers are intialized to all 0s automatically
    vector<vector<int>> temp2( dim2 , vector<int> (dim2));

  //P1 = A(F - H)
  subtractMat(F, H, temp1, dim2);
  strassen(A, temp1, P1, dim2, padStart, currRow1, currCol1, currRow2, currCol2 + dim2); //F is "less far out" than H

  //P2 = (A + B)H
  addMat(A, B, temp1, dim2);
  strassen(temp1, H, P2, dim2, padStart, currRow1, currCol1, currRow2 + dim2, currCol2 + dim2);

  //P3 = (C + D)E
  addMat(C, D, temp1, dim2);
  strassen(temp1, E, P3, dim2, padStart, currRow1 + dim2, currCol1, currRow2, currCol2);

  //P4 = D(G - E)
  subtractMat(G, E, temp1, dim2);
  strassen(D, temp1, P4, dim2, padStart, currRow1 + dim2, currCol1 + dim2, currRow2, currCol2);

  //P5 = (A + D)(E + H)
  addMat(A, D, temp1, dim2);
  addMat(E, H, temp2, dim2);
  strassen(temp1, temp2, P5, dim2, padStart, currRow1, currCol1, currRow2, currCol2);

  //P6 = (B - D)(G + H)
  subtractMat(B, D, temp1, dim2);
  addMat(G, H, temp2, dim2);
  strassen(temp1, temp2, P6, dim2, padStart, currRow1, currCol2 + dim2, currRow2 + dim2, currCol2);

  //P7 = (C - A)(E + F)
  subtractMat(C, A, temp1, dim2);
  addMat(E, F, temp2, dim2);
  strassen(temp1, temp2, P7, dim2, padStart, currRow1, currCol1, currRow2, currCol2);

  for (int i = 0; i < dim2; i++){
    for (int j = 0; j < dim2; j++){
        reslt[i][j] = P5[i][j] + P4[i][j] + P6[i][j] - P2[i][j]; // P5 + (P4 - P2) + P6
        reslt[i][j + dim2] = P1[i][j] + P2[i][j]; //Q2 = P1 + P2
        reslt[i + dim2][j] = P3[i][j] + P4[i][j]; //P3 + P4
        reslt[i + dim2][j + dim2] = P1[i][j] - P3[i][j] + P5[i][j] + P7[i][j]; //P1 - P3 + P5 + P7
    }
  }
}
}

void multTests(int t, int range, int trials) {

    int dimp = t; //only works for powers of 2
    dimp = int(pow(2, ceil(log2(double(t))))); // gets smallest the power of 2 greater than or equal to the dimension
	vector<vector<int>> mat_A = randMatrix(t, range);
	vector<vector<int>> mat_B = randMatrix(t, range);

	vector<vector<int>> A(dimp, vector<int> (dimp, 0));
	vector<vector<int>> B(dimp, vector<int> (dimp, 0));


	for (int i = 0; i < t; i++){
	    for (int j = 0; j < t; j++){
            A[i][j] = mat_A[i][j];
            B[i][j] = mat_B[i][j];
	    }
	}
	vector<vector<int>> str( dimp , vector<int> (dimp));
	vector<vector<int>> basic( t , vector<int> (t)); //**** DIM

	if (t < 50){
	    strassen(A, B, str, dimp, t);
	    basicMult(A, B, basic, t);

        cout << "A: \n";
        printMatrix(A);
        cout << "B: \n";
        printMatrix(B);


        cout << "\n basic mult: \n";
        printDiagonals(basic);

        cout << "STRASSEN MULT:\n";

        cout << "\n strassen mult: \n";
        printDiagonals(str, t);
	}

    double avgStrasTime = 0;
    double avgBasicTime = 0;

    for (int i = 0; i < trials; i++){ //timing
        mat_A = randMatrix(t, range);
        mat_B = randMatrix(t, range);
        for (int i = 0; i < t; i++){
            for (int j = 0; j < t; j++){
                A[i][j] = mat_A[i][j];
                B[i][j] = mat_B[i][j];
            }
	    }

        auto start = high_resolution_clock::now();
        strassen(A, B, str, dimp, t);
        auto stop = high_resolution_clock::now(); //measures time elapsed for Strassen's algorithm; timing function from https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
        auto duration = duration_cast<microseconds>(stop - start);

        avgStrasTime += duration.count() / trials; //stras time

        cout << "did stras \n";

        start = high_resolution_clock::now();
//        if (i == 0){
//            basicMult(A, B, basic, t);
//        }
        //
        stop = high_resolution_clock::now(); //measures time elapsed for Strassen's algorithm; timing function from https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
        duration = duration_cast<microseconds>(stop - start);

        cout << "did basic \n";
        avgBasicTime += duration.count();// / trials; //basic time

    }
    cout << "\n dim: " << t;
    cout << "\n stras time (microsec): " << avgStrasTime;
    cout << "\n basic time (microsec): " << avgBasicTime;
}

double graphCount(double p, int dim){
    vector<vector<int>> reslt( dim , vector<int> (dim));
    double pr;

    for (int i = 0; i < dim; i++){ //generating random adjacency matrix
        for (int j = 0; j < dim; j++){
            pr = double(std::rand()) / double(RAND_MAX);
            if (pr <= p){
                reslt[i][j] = 1;
            }
            else{
                reslt[i][j] = 0;
            }
        }
    }

    //getting A^3

    vector<vector<int>> temp1( dim , vector<int> (dim));
    vector<vector<int>> temp2( dim , vector<int> (dim));

    strassen(reslt, reslt, temp1, dim, dim);
    strassen(temp1, reslt, temp2, dim, dim);

    double total = 0;

    //getting sum of diagonal entries

    for (int i = 0; i < dim; i++){
        total += temp2[i][i];
    }

    return (total / 6);

}


int main(int argc, char *argv[])
{
	if (argc != 4) {
        cout << "./strassen 0 dimension inputfile\n";
        return 1;
	}
    srand (time (0)); //seeding RAND

	int dim = stoi(argv[2]); //stoi may cause problems with older c++
	int test = stoi(argv[1]); //optional toggle for tests
	int dimp = dim;


    dimp = int(pow(2, ceil(log2(double(dim))))); // gets smallest the power of 2 greater than or equal to the dimension
        //cout << "dimp: " << dimp << "\n";

	string inpt = argv[3];

	ifstream input;
	input.open(inpt);

	vector<vector<int>> reslt(dimp, vector<int> (dimp, 0)); //padded with 0s to nearest power of 2;

	string inputString;

	vector<vector<int>> mat_A(dimp, vector<int> (dimp, 0)); //padded with 0s to nearest power of 2;
	vector<vector<int>> mat_B(dimp, vector<int> (dimp, 0)); //padded with 0s to nearest power of 2;

	int n = 0;
	int r = 0;

	if ( input.is_open() ) { // reads file and turns it into 2 matrices
        while (getline(input, inputString)) {
            if (n < dim * dim){
                mat_A[r][n - r * dim] = stoi(inputString);
                if ((n + 1) % dim == 0){
                    r = r + 1;
                }
            } //first matrix
            else {
                if (n == dim * dim){
                    r = 0;
                    //cout << "reset \n";
                }
                mat_B[r][n - (r + dim)*dim] = stoi(inputString);
                if ((n + 1) % dim == 0){
                        r = r + 1;
                }
            }
            n = n + 1;
        }
    }
    if (test == 1){
        cout << "\n A: \n";
        printMatrix(mat_A);
        cout << "\n B: \n";
        printMatrix(mat_B);
    }


    //vector<vector<int>> result(dim, vector<int> (dim, 0));
    vector<vector<int>> results(dimp, vector<int> (dimp, 0));
//    basicMult(mat_A, mat_B, result, dim); //basic uses normal dim, not padded
//    cout << "Result: \n";
//    printMatrix(result);
//
//    cout << "strass: \n";
//
//    vector<vector<int>> results(dimp, vector<int> (dimp, 0));
//    printMatrix(results);

    strassen(mat_A, mat_B, results, dimp, dim + 1); //padding starts at row/column right after dim
    printDiagonals(results, dim);







    if (test == 1){
        printMatrix(results);
        cout << "\n Tests: \n";

        vector<vector<int>> addTest(dimp, vector<int> (dimp, 0));
        cout << "\n AD TEST \n";
        addMat(mat_A, mat_B, addTest, dimp);
        printMatrix(addTest, dimp);

        cout << "\n SUB TEST \n";
        subtractMat(mat_A, mat_B, addTest, dimp);
        printMatrix(addTest, dimp);
        multTests(4096, 3, 5); //dim, range, trials
    }

    if (test == 2){
        double p01 = 0;
        double p02 = 0;
        double p03 = 0;
        double p04 = 0;
        double p05 = 0;
        for (int i = 0; i < 3; i++){
            p01 += graphCount(0.01, 1024) / 3;
            p02 += graphCount(0.02, 1024) / 3;
            p03 += graphCount(0.03, 1024) / 3;
            p04 += graphCount(0.04, 1024) / 3;
            p05 += graphCount(0.05, 1024) / 3;
        }

        cout << "p01: " << p01 << "\n";
        cout << "p02: " << p02 << "\n";
        cout << "p03: " << p03 << "\n";
        cout << "p04: " << p04 << "\n";
        cout << "p05: " << p05 << "\n";
    }




	return 0;
}

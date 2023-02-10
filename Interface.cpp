#include <iostream>
#include<fstream>
#include<cstdlib>
#include<stdio.h>
#include <cstring>
#include <ctime>
#include <vector>
#include "Interface.h"
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include "Solver.h"
#include "Solver.cpp"

using namespace std;

Interface::Interface() {}

Interface::~Interface() {}

void Interface::interfaceStart() {
    cout << endl;
    cout << " -------------------------------------------" << endl;
    cout << "|  MATRIX  LINEAR  SOLVER  COMMANDLINE  TOOL |" << endl;
    cout << " -------------------------------------------" << endl;
    cout << endl;
    cout << " ------------------------------------------------------" << endl;
    cout << "|                     Matrix Data                      |" << endl;
    cout << " ------------------------------------------------------" << endl;
    cout << "| Data inside the matrix can be randomly generated or  |" << endl;
    cout << "| entered by user through keyboard and text file.  And |" << endl;
    cout << "| matrix data range can be decided by user.            |" << endl;
    cout << " ------------------------------------------------------" << endl;
    cout << endl;
    cout << " ------------------------------------------------------" << endl;
    cout << "|                     Matrix Size                      |" << endl;
    cout << " ------------------------------------------------------" << endl;
    cout << "| Matrix size are set to be larger than 10 X 10        |" << endl;
    cout << "| entered by user through keyboard and text file.      |" << endl;
    cout << "| Two types of matrices are provided                   |" << endl;
    cout << " ------------------------------------------------------" << endl;
    cout << endl;
    cout << " ---------------------------------- " << endl;
    cout << "|          Matrix Type             |" << endl;
    cout << " ---------------------------------- " << endl;
    cout << "| 1. Dense matrix                  |" << endl;
    cout << "| 2. Sparse matrix                 |" << endl;
    cout << " ----------------------------------" << endl;
    cout << endl;
    cout << " ----------------------------------------------------- " << endl;
    cout << "| Do you want to continue and select the matrix (y/n) |" << endl;
    cout << " ----------------------------------------------------- " << endl;
    cout << "| Enter y to continue                                 |" << endl;
    cout << "| Enter n to exit                                     |" << endl;
    cout << " ---------------------------------------------------- " << endl;
    cout << ">> ";
    char option;
    cin >> option;
    if (cin.fail()){
        // if input is not a char, remind user the type of input
        cerr << "Please enter a char.";
    } else if (option == 'y'){
        interfacePartOne();
    } else if (option == 'n'){
        exit(0);
    } else{
        // if input is not y or n
        // re-ask the user
        interfaceStart();
    }
}

void Interface::interfacePartOne() {
    // clears the error flag on cin
    // so that future I/O operations will work correctly
    cin.clear();
    cout << endl;
    cout << " ----------------------------------------------------------" << endl;
    cout << "|               Select matrix storage type                  |" << endl;
    cout << " ----------------------------------------------------------" << endl;
    cout << "| 1: Dense Matrix                                           |" << endl;
    cout << "| 2: Sparse Matrix                                          |" << endl;
    cout << "| 3: Test Matrix (Run testing on each solver automatically) |" << endl;
    cout << "| b: Back                                                   |" << endl;
    cout << "| x: Exit                                                   |" << endl;
    cout << " ----------------------------------------------------------" << endl;
    cout << ">> ";
    char option;
    cin >> option;
    if (cin.fail()){
        // if input is not a char, remind user the type of input
        cerr << "Please enter a char.";
    }
    switch (option) {
        case '1':{
            interfacePartTwo(option);
        }
        case '2':{
            interfacePartTwo(option);
        }
        case '3':{
            cout << endl;
            cout << " ----------------------------------------------------------" << endl;
            interfaceMatrixTest();
            cout << " ----------------------------------------------------------" << endl;
            cout << endl;
            exit(0);
        }
        case 'b':{
            // back to the start
            interfaceStart();
        }
        case 'x':{
            exit(0);
        }
        default:{
            cerr << "Invalid input";
        }
    }
}

void Interface::interfacePartTwo(char matrixType) {
    if (matrixType == '1'){
        // Dense matrix
        cout << endl;
        cout << " -----------------------------------------------" << endl;
        cout << "|           Set up the matrix size              |:" << endl;
        cout << " -----------------------------------------------" << endl;
        cout << "| 1: 10 x 10                                    |" << endl;
        cout << "| 2: 15 x 15                                    |" << endl;
        cout << "| 3: Enter the size (greater than 10 X 10)      |" << endl;
        cout << "| b: Back                                       |" << endl;
        cout << "| x: Exit                                       |" << endl;
        cout << " -----------------------------------------------" << endl;
        cout << ">> ";

        char option;
        cin >> option;
        if (cin.fail()){
            // if input is not a char, remind user the type of input
            cerr << "Please enter a char.";
        }

        // matrix size variables
        int rows, cols;
        switch (option) {
            case '1':{
                rows = 10;
                cols = 10;
                interfacePartThree(rows, cols);
            }
            case '2':{
                rows = 15;
                cols = 15;
                interfacePartThree(rows, cols);
            }
            case '3':{
                cout << " ------------------------------------------------------------------ " << endl;
                cout << "| Enter the row number or column number since it is a square matrix |" << endl;
                cout << " ------------------------------------------------------------------ " << endl;
                cout << ">> ";
                // read the row number
                cin >> rows;
                // column number equals row number
                cols = rows;
                interfacePartThree(rows, cols);
            }
            case 'b':{
                // re-select the martix type
                interfacePartOne();
            }
            case 'x':{
                exit(0);
            }
            default:{
                cerr << "Invalid input";
            }
        }
    } else if (matrixType == '2'){
        // Sparse matrix

        cout << endl;
        cout << " -----------------------------------------------" << endl;
        cout << "|           Set up the matrix size              |:" << endl;
        cout << " -----------------------------------------------" << endl;
        cout << "| 1: 10 x 10                                    |" << endl;
        cout << "| 2: 15 x 15                                    |" << endl;
        cout << "| 3: Enter the size (greater than 10 X 10)      |" << endl;
        cout << "| b: Back                                       |" << endl;
        cout << "| x: Exit                                       |" << endl;
        cout << " -----------------------------------------------" << endl;
        cout << ">> ";

        char option;
        cin >> option;
        if (cin.fail()){
            // if input is not a char, remind user the type of input
            cerr << "Please enter a char.";
        }

        // matrix size variables
        int rows, cols;
        switch (option) {
            case '1':{
                rows = 10;
                cols = 10;
                interfacePartThree_sparse(rows, cols);
            }
            case '2':{
                rows = 15;
                cols = 15;
                interfacePartThree_sparse(rows, cols);
            }
            case '3':{
                cout << " ------------------------------------------------------------------ " << endl;
                cout << "| Enter the row number or column number since it is a square matrix |" << endl;
                cout << " ------------------------------------------------------------------ " << endl;
                cout << ">> ";
                // read the row number
                cin >> rows;
                // column number equals row number
                cols = rows;
                interfacePartThree_sparse(rows, cols);
            }
            case 'b':{
                // re-select the martix type
                interfacePartOne();
            }
            case 'x':{
                exit(0);
            }
            default:{
                cerr << "Invalid input";
            }
        }
    }
}

void Interface::interfacePartThree_sparse(int rows, int cols){
    // fill the sparse matrix data

    cout << endl;
    cout << " -----------------------------------------------" << endl;
    cout << "|        Fill up the sparse matrix data         |:" << endl;
    cout << " -----------------------------------------------" << endl;
    cout << "| 1: Randomly generated data                    |" << endl;
    cout << "| 2: Enter the data                             |" << endl;
    cout << "| 3: Text file                                  |" << endl;
    cout << "| b: Back                                       |" << endl;
    cout << "| x: Exit                                       |" << endl;
    cout << " -----------------------------------------------" << endl;
    cout << ">> ";

    char option;
    cin >> option;
    if (cin.fail()){
        // if input is not a char, remind user the type of input
        cerr << "Please enter a char.";
    }

    switch (option) {
        case '1':{
            // random generation of sparse matrix function to be written...、
            double left_boundary = 0;
            double right_boundary = 0;
            cout << endl;
            cout << " ------------------------------------------------------------------------- " << endl;
            cout << "| Generate randomly generated data between left_boundary and right_boundary |" << endl;
            cout << " ------------------------------------------------------------------------- " << endl;
            cout << "  Please enter the left_boundary first" << endl;
            cout << ">> ";
            cin >> left_boundary;
            cout << "  Then please enter the right_boundary first" << endl;
            cout << ">> ";
            cin >> right_boundary;

            // generate random csr matrix
            CSRMatrix<double> mat;
            CSRMatrix<double> randomMat =
                    mat.randomCSRMatrix(rows, cols, left_boundary, right_boundary);

            // print out the random csr matrix
            cout << endl;
            cout << "------------------------------------------------" << endl;
            cout << "Print entered matrix A";
            randomMat.printMatrix();
            cout << "------------------------------------------------" << endl;

            cout << "------------------------------------------------" << endl;
            cout << "Print entered matrix A in CSR format:" << endl;
            randomMat.printMatrixCSR();
            cout << "------------------------------------------------" << endl;
            cout << endl;

            // generate and print the random vector
            cout << "------------------------------------------------" << endl;
            double a[rows];
            randomVectb(a, rows, left_boundary, right_boundary);
            cout << "Print randomly generated vector b:" << endl;
            cout << "[";
            for (int i = 0; i < randomMat.cols; ++i) {
                if (i == randomMat.cols-1){
                    cout << a[i];
                } else{
                    cout << a[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "------------------------------------------------" << endl;

            interfacePartFour_sparse(randomMat, a);
        }
        case '2':{

            cout << endl;
            cout << " -----------------------------------------------------" << endl;
            cout << "|     The way to fill up sparse matrix data           |:" << endl;
            cout << " -----------------------------------------------------" << endl;
            cout << "| 1: Enter the values, row position, col index arrays |" << endl;
            cout << "| 2: Enter Row-major ordering matrix                  |" << endl;
            cout << "| b: Back                                             |" << endl;
            cout << "| x: Exit                                             |" << endl;
            cout << " -----------------------------------------------------" << endl;
            cout << ">> ";

            char selection;
            cin >> selection;
            if (cin.fail()){
                // if input is not a char, remind user the type of input
                cerr << "Please enter a char.";
            }

            switch (selection) {
                case '1':{
                    cout << endl;
                    cout << " ---------------------------------------------------------------------------- " << endl;
                    cout << "|     Follow instructions below to enter the matrix data through keyboard     |" << endl;
                    cout << " ---------------------------------------------------------------------------- " << endl;
                    cout << "| User is expected to ensure that the inputs given to the jacobi method       |" << endl;
                    cout << "| and the gauss seidel method are diagonally dominant to ensure convergence.  |" << endl;
                    cout << " ---------------------------------------------------------------------------- " << endl;

                    int nnzs;
                    cout << "  Please enter the value of non-zero entries in the sparse matrix" << endl;
                    cout << ">> ";
                    cin >> nnzs;

                    // enter the sparse matrix data through keyboard
                    auto *cvalues = new double[rows * cols];
                    for (int i = 0; i < rows * cols; i++) {
                        cout << "  Please enter the value stored in value array[" << i << "]" << endl;
                        cout << ">> ";
                        cin >> cvalues[i];
                    }

                    int *row_position = new int[nnzs];
                    for (int i = 0; i < rows + 1; i++) {
                        cout << "  Please enter the value stored in row position array[" << i << "]" << endl;
                        cout << ">> ";
                        cin >> row_position[i];
                    }

                    int *col_index = new int[rows * cols];
                    for (int i = 0; i < rows * cols; i++) {
                        cout << "  Please enter the value stored in column index array[" << i << "]" << endl;
                        cout << ">> ";
                        cin >> col_index[i];
                    }

                    // generate rows X cols sparse matrix with nnzs
                    CSRMatrix<double> mat(rows, cols, nnzs, cvalues, row_position, col_index);


                    // print out the sparse matrix
                    cout << endl;
                    cout << "------------------------------------------------" << endl;
                    cout << "Print entered matrix A";
                    mat.printMatrix();
                    cout << "------------------------------------------------" << endl;

                    cout << "------------------------------------------------" << endl;
                    cout << "Print entered matrix A in CSR format:" << endl;
                    mat.printMatrixCSR();
                    cout << "------------------------------------------------" << endl;
                    cout << endl;

                    // generate a vector with zeros
                    double a[rows];
                    // looping to enter values in vector b
                    for (int i = 0; i < rows; i++) {
                        cout << "  Please enter the value stored in vector[" << i << "]" << endl;
                        cout << ">> ";
                        cin >> a[i];
                    }

                    // generate and print the entered vector
                    cout << "------------------------------------------------" << endl;
                    cout << "Print entered vector b:" << endl;
                    cout << "[";
                    for (int i = 0; i < mat.cols; ++i) {
                        if (i == mat.cols-1){
                            cout << a[i];
                        } else{
                            cout << a[i] << ",";
                        }
                    }
                    cout << "]" << endl;
                    cout << "------------------------------------------------" << endl;

                    interfacePartFour_sparse(mat, a);

                }

                case '2':{
                    cout << endl;
                    cout << " ---------------------------------------------------------------------------- " << endl;
                    cout << "|     Follow instructions below to enter the matrix data through keyboard     |" << endl;
                    cout << " ---------------------------------------------------------------------------- " << endl;
                    cout << "| User is expected to ensure that the inputs given to the jacobi method       |" << endl;
                    cout << "| and the gauss seidel method are diagonally dominant to ensure convergence.  |" << endl;
                    cout << " ---------------------------------------------------------------------------- " << endl;

                    auto* values = new double[rows * cols];
                    for (int i=0; i < rows; i++){
                        for (int j=0; j < cols; j++){
                            cout << "  Please enter the value stored in matrix[" << i << "," << j << "]" << endl;
                            cout << ">> ";
                            cin >> values[i * cols + j];
                        }
                    }

                    CSRMatrix<double> mat(rows, cols, values);

                    // print out the sparse matrix
                    cout << endl;
                    cout << "------------------------------------------------" << endl;
                    cout << "Print entered matrix A";
                    mat.printMatrix();
                    cout << "------------------------------------------------" << endl;

                    cout << "------------------------------------------------" << endl;
                    cout << "Print entered matrix A in CSR format:" << endl;
                    mat.printMatrixCSR();
                    cout << "------------------------------------------------" << endl;
                    cout << endl;

                    // generate a vector with zeros
                    double a[rows];
                    // looping to enter values in vector b
                    for (int i = 0; i < rows; i++) {
                        cout << "  Please enter the value stored in vector[" << i << "]" << endl;
                        cout << ">> ";
                        cin >> a[i];
                    }

                    // generate and print the entered vector
                    cout << "------------------------------------------------" << endl;
                    cout << "Print entered vector b:" << endl;
                    cout << "[";
                    for (int i = 0; i < mat.cols; ++i) {
                        if (i == mat.cols-1){
                            cout << a[i];
                        } else{
                            cout << a[i] << ",";
                        }
                    }
                    cout << "]" << endl;
                    cout << "------------------------------------------------" << endl;

                    interfacePartFour_sparse(mat, a);

                }

                case 'b':{
                    interfacePartOne();
                }

                case 'x':{
                    exit(0);
                }

                default:{
                    cerr << "Invalid input";
                }
            }

//            interfacePartFour(mat, a);
        }
        case '3':{
            // read sparse matrix data from text file function to be written...

            // read matrix data from txt file
            ifstream ifs;
            // edit the run configuration
            // set the working directory into file directory
            ifs.open("matrix_data.txt",ios::in);

            if (!ifs.is_open()){
                cerr << "Open matrix data file failed." << endl;
                exit(0);
            }

            string buf;
            while (getline(ifs, buf)){
            }

            // convert string into char array
            char arr[buf.size() + 1];
            strcpy(arr, buf.c_str());

            // creating a vector for storing each matrix element
            vector<int> matrixData;

            // segment string with delimiter ","
            char *p;
            p = strtok(arr, ",");
            while(p != nullptr){
                if (*p == '-'){
                    // if the element is -5, the split char is stored as ['-','5']
                    char tmp = *(p+1);
                    // it is a negative number
                    int val = 0 - (tmp - '0');
                    matrixData.push_back(val);
                } else{
                    int val = *p - '0';
                    matrixData.push_back(val);
                }
                p = strtok(nullptr, ",");
            }

            int *values = new int[matrixData.size()];
            if (!matrixData.empty()){
                memcpy(values, &matrixData[0], matrixData.size()*sizeof(int));
            }

            // enter the matrix data through txt file
            CSRMatrix<int> mat(rows, cols, values);

            // print out the matrix
            cout << endl;
            cout << "------------------------------------------------" << endl;
            cout << "Print entered matrix A";
            mat.printMatrix();
            cout << "------------------------------------------------" << endl;

            cout << "------------------------------------------------" << endl;
            cout << "Print entered matrix A in CSR format:" << endl;
            mat.printMatrixCSR();
            cout << "------------------------------------------------" << endl;
            cout << endl;

            // close the first stream
            ifs.close();

            // read vector data file
            ifstream ifs1;
            ifs1.open("vector_data.txt",ios::in);

            if (!ifs1.is_open()){
                cerr << "Open vector data file failed." << endl;
                exit(0);
            }


            string buf1;
            while (getline(ifs1, buf1)){
            }

            // convert string into char array
            char arr1[buf1.size() + 1];
            strcpy(arr1, buf1.c_str());

            // creating a vector for storing each matrix element
            vector<int> vectorData;

            // segment string with delimiter ","
            char *p1;
            p1 = strtok(arr1, ",");
            while(p1 != nullptr){
                if (*p1 == '-'){
                    // if the element is -5, the split char is stored as ['-','5']
                    char tmp = *(p1+1);
                    // it is a negative number
                    int val = 0 - (tmp - '0');
                    vectorData.push_back(val);
                } else{
                    int val = *p1 - '0';
                    vectorData.push_back(val);
                }
                p1 = strtok(nullptr, ",");
            }

            // generate a vector with zeros
            double a[rows];
            for (int i = 0; i < rows; i++) {
                a[i] = vectorData[i];
            }

            // generate and print the random vector
            cout << "------------------------------------------------" << endl;
            cout << "Print vector b from vector_data.txt file:" << endl;
            cout << "[";
            for (int i = 0; i < mat.cols; ++i) {
                if (i == mat.cols-1){
                    cout << a[i];
                } else{
                    cout << a[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "------------------------------------------------" << endl;

            ifs1.close();

//            interfacePartFour_sparse(mat, a);

            exit(0);
        }
        case 'b':{
            interfacePartOne();
        }
        case 'x':{
            exit(0);
        }
        default:{
            cerr << "Invalid input";
        }
    }
}

void Interface::interfacePartThree(int rows, int cols) {
    cout << endl;
    cout << " -----------------------------------------------" << endl;
    cout << "|        Fill up the dense matrix data          |:" << endl;
    cout << " -----------------------------------------------" << endl;
    cout << "| 1: Randomly generated data                    |" << endl;
    cout << "| 2: Enter the data                             |" << endl;
    cout << "| 3: Text file                                  |" << endl;
    cout << "| b: Back                                       |" << endl;
    cout << "| x: Exit                                       |" << endl;
    cout << " -----------------------------------------------" << endl;
    cout << ">> ";

    char option;
    cin >> option;
    if (cin.fail()){
        // if input is not a char, remind user the type of input
        cerr << "Please enter a char.";
    }

    switch (option) {
        case '1':{
            int left_boundary = 0;
            int right_boundary = 0;
            cout << endl;
            cout << " ------------------------------------------------------------------------- " << endl;
            cout << "| Generate randomly generated data between left_boundary and right_boundary |" << endl;
            cout << " ------------------------------------------------------------------------- " << endl;
            cout << "  Please enter the left_boundary first" << endl;
            cout << ">> ";
            cin >> left_boundary;
            cout << "  Then please enter the right_boundary first" << endl;
            cout << ">> ";
            cin >> right_boundary;

            // generate random matrix
            Matrix<double> mat(rows, cols);
            mat = mat + mat.randomMatrix(rows, cols, left_boundary, right_boundary);

            // print out the random matrix
            cout << endl;
            cout << "------------------------------------------------" << endl;
            cout << "Print randomly generated matrix A";
            mat.printMatrix();
            cout << "------------------------------------------------" << endl;
            cout << endl;

            // generate and print the random vector
            cout << "------------------------------------------------" << endl;
            double a[rows];
            randomVectb(a, rows, left_boundary, right_boundary);
            cout << "Print randomly generated vector b:" << endl;
            cout << "[";
            for (int i = 0; i < mat.cols; ++i) {
                if (i == mat.cols-1){
                    cout << a[i];
                } else{
                    cout << a[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "------------------------------------------------" << endl;

            interfacePartFour(mat, a);
        }
        case '2':{

            cout << endl;
            cout << " ---------------------------------------------------------------------------- " << endl;
            cout << "|     Follow instructions below to enter the matrix data through keyboard     |" << endl;
            cout << " ---------------------------------------------------------------------------- " << endl;
            cout << "| User is expected to ensure that the inputs given to the jacobi method       |" << endl;
            cout << "| and the gauss seidel method are diagonally dominant to ensure convergence.  |" << endl;
            cout << " ---------------------------------------------------------------------------- " << endl;

            // enter the matrix data through keyboard
            // generate rows X cols matrix with zeros
            Matrix<double> mat(rows, cols);
            // looping to enter the value
            for (int i=0; i < rows; i++){
                for (int j=0; j < cols; j++){
                    cout << "  Please enter the value stored in matrix[" << i << "," << j << "]" << endl;
                    cout << ">> ";
                    cin >> mat.values[i * cols + j];
                }
            }

            // print out the matrix
            cout << endl;
            cout << "------------------------------------------------" << endl;
            cout << "Print entered matrix A";
            mat.printMatrix();
            cout << "------------------------------------------------" << endl;
            cout << endl;

            // generate a vector with zeros
            double a[rows];
            for (int i = 0; i < rows; i++) {
                cout << "  Please enter the value stored in vector[" << i << "]" << endl;
                cout << ">> ";
                cin >> a[i];
            }

            // generate and print the random vector
            cout << "------------------------------------------------" << endl;
            cout << "Print entered vector b:" << endl;
            cout << "[";
            for (int i = 0; i < mat.cols; ++i) {
                if (i == mat.cols-1){
                    cout << a[i];
                } else{
                    cout << a[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "------------------------------------------------" << endl;

            interfacePartFour(mat, a);
        }
        case '3':{
            // read matrix data from txt file
            ifstream ifs;
            // edit the run configuration
            // set the working directory into file directory
            ifs.open("matrix_data.txt",ios::in);

            if (!ifs.is_open()){
                cerr << "Open matrix data file failed." << endl;
                exit(0);
            }

            string buf;
            while (getline(ifs, buf)){
            }

            // convert string into char array
            char arr[buf.size() + 1];
            strcpy(arr, buf.c_str());

            // creating a vector for storing each matrix element
            vector<int> matrixData;

            // segment string with delimiter ","
            char *p;
            p = strtok(arr, ",");
            while(p != nullptr){
                if (*p == '-'){
                    // if the element is -5, the split char is stored as ['-','5']
                    char tmp = *(p+1);
                    // it is a negative number
                    int val = 0 - (tmp - '0');
                    matrixData.push_back(val);
                } else{
                    int val = *p - '0';
                    matrixData.push_back(val);
                }
                p = strtok(nullptr, ",");
            }

            // enter the matrix data through txt file
            Matrix<double> mat(rows, cols);
            // looping to enter the value
            for (int i=0; i < rows; i++){
                for (int j=0; j < cols; j++){
                    mat.values[i * cols + j] = matrixData[i * cols + j];
                }
            }

            // print out the matrix
            cout << endl;
            cout << "------------------------------------------------" << endl;
            cout << "Print matrix A from matrix_data.txt";
            mat.printMatrix();
            cout << "------------------------------------------------" << endl;
            cout << endl;

            // close the first stream
            ifs.close();

            // read vector data file
            ifstream ifs1;
            ifs1.open("vector_data.txt",ios::in);

            if (!ifs1.is_open()){
                cerr << "Open vector data file failed." << endl;
                exit(0);
            }


            string buf1;
            while (getline(ifs1, buf1)){
            }

            // convert string into char array
            char arr1[buf1.size() + 1];
            strcpy(arr1, buf1.c_str());

            // creating a vector for storing each matrix element
            vector<int> vectorData;

            // segment string with delimiter ","
            char *p1;
            p1 = strtok(arr1, ",");
            while(p1 != nullptr){
                if (*p1 == '-'){
                    // if the element is -5, the split char is stored as ['-','5']
                    char tmp = *(p1+1);
                    // it is a negative number
                    int val = 0 - (tmp - '0');
                    vectorData.push_back(val);
                } else{
                    int val = *p1 - '0';
                    vectorData.push_back(val);
                }
                p1 = strtok(nullptr, ",");
            }

            // generate a vector with zeros
            double a[rows];
            for (int i = 0; i < rows; i++) {
                a[i] = vectorData[i];
            }

            // generate and print the random vector
            cout << "------------------------------------------------" << endl;
            cout << "Print vector b from vector_data.txt file:" << endl;
            cout << "[";
            for (int i = 0; i < mat.cols; ++i) {
                if (i == mat.cols-1){
                    cout << a[i];
                } else{
                    cout << a[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "------------------------------------------------" << endl;

            ifs1.close();

            interfacePartFour(mat, a);
        }
        case 'b':{
            interfacePartOne();
        }
        case 'x':{
            exit(0);
        }
        default:{
            cerr << "Invalid input";
        }
    }
}

template<class U>
void Interface::interfacePartFour_sparse(CSRMatrix<U> &matA, U *b) {
    cout << endl;
    cout << " ----------------------------------------------" << endl;
    cout << "|  Select the sparse matrix solver              |:" << endl;
    cout << " ----------------------------------------------" << endl;
    cout << "| 1: Gaussian Elimination                       |" << endl;
    cout << "| 2: Gaussian Jordan                            |" << endl;
    cout << "| 3: Gaussian Seidel                            |" << endl;
    cout << "| 4: LU factorisation                           |" << endl;
    cout << "| 5: Jacobi                                     |" << endl;
    cout << "| b: Back                                       |" << endl;
    cout << "| x: Exit                                       |" << endl;
    cout << " -----------------------------------------------" << endl;
    cout << ">> ";

    char option;
    cin >> option;
    if (cin.fail()){
        // if input is not a char, remind user the type of input
        cerr << "Please enter a char.";
    }

    switch (option) {
        case '1':{
            auto* output = new double[matA.rows];

            Solver s1;
            clock_t start = clock();
            s1.Gaussian_Elimination_CSRMatrix(matA, b, output);
            clock_t end = clock();

            cout << "The output solved by Gaussian Elimination CSRMatrix: [";
            for (int i = 0; i < matA.cols; ++i) {
                if (i == matA.cols-1){
                    cout << output[i];
                } else{
                    cout << output[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
            cout << "------------------------------------------------" << endl;

            delete[] output;

            exit(0);
        }

        case '2':{
            auto* output = new double[matA.rows];

            Solver s1;
            clock_t start = clock();
            s1.Gauss_Jordan_CSRMatrix(matA, b, output);
            clock_t end = clock();

            cout << "The output solved by Gaussian Jordan CSRMatrix: [";
            for (int i = 0; i < matA.cols; ++i) {
                if (i == matA.cols-1){
                    cout << output[i];
                } else{
                    cout << output[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
            cout << "------------------------------------------------" << endl;

            delete[] output;

            exit(0);
        }

        case '3':{
            auto* output = new double[matA.rows];

            Solver s1;
            clock_t start = clock();
            s1.Gauss_Seidel_CSRMatrix(matA, b, output);
            clock_t end = clock();

            cout << "The output solved by Gaussian Seidel CSRMatrix: [";
            for (int i = 0; i < matA.cols; ++i) {
                if (i == matA.cols-1){
                    cout << output[i];
                } else{
                    cout << output[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
            cout << "------------------------------------------------" << endl;

            delete[] output;

            exit(0);
        }

        case '4':{
            auto* output = new double[matA.rows];

            Solver s1;
            clock_t start = clock();
            s1.LU_solver_CSRMatrix( matA,  b, output);
            clock_t end = clock();

            cout << "The output solved by LU Solver CSRMatrix: [";
            for (int i = 0; i < matA.cols; ++i) {
                if (i == matA.cols-1){
                    cout << output[i];
                } else{
                    cout << output[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
            cout << "------------------------------------------------" << endl;

            delete[] output;

            exit(0);
        }

        case '5':{
            auto* output = new double[matA.rows];

            Solver s1;
            int it_max = 1000;
            double it_tol = 1.e-6;

            clock_t start = clock();
            s1.Jacobi_Method_CSRMatrix(matA, b, it_max, it_tol, output);
            clock_t end = clock();

            cout << "The output solved by Jacobi Method CSRMatrix: [";
            for (int i = 0; i < matA.cols; ++i) {
                if (i == matA.cols-1){
                    cout << output[i];
                } else{
                    cout << output[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
            cout << "------------------------------------------------" << endl;

            delete[] output;

            exit(0);
        }

        case 'b':{
            interfacePartOne();
        }

        case 'x':{
            exit(0);
        }

        default:{
            cerr << "Invalid input";
        }
    }
}

template <class U>
void Interface::interfacePartFour(Matrix<U>& matA, U* b) {
    cout << endl;
    cout << " ----------------------------------------------" << endl;
    cout << "|  Select the dense matrix solver               |:" << endl;
    cout << " ----------------------------------------------" << endl;
    cout << "| 1: Gaussian Elimination                       |" << endl;
    cout << "| 2: Gaussian Jordan                            |" << endl;
    cout << "| 3: Gaussian Seidel                            |" << endl;
    cout << "| 4: LU factorisation                           |" << endl;
    cout << "| 5: LU factorisation with partial pivoting     |" << endl;
    cout << "| 6: Jacobi                                     |" << endl;
    cout << "| b: Back                                       |" << endl;
    cout << "| x: Exit                                       |" << endl;
    cout << " -----------------------------------------------" << endl;
    cout << ">> ";
    char option;
    cin >> option;
    if (cin.fail()){
        // if input is not a char, remind user the type of input
        cerr << "Please enter a char.";
    }
    switch (option) {
        case '1':{
            auto* output = new double[matA.rows];
            Solver s1;
            clock_t start = clock();
            s1.Gaussian_Elimination(matA, b, output);
            clock_t end = clock();
            cout << endl;
            cout << "------------------------------------------------" << endl;
            cout << "The output solved by Gaussian Elimination:" << endl;
            cout << "[ ";
            for (int i = 0; i < matA.rows; i++) {
                if (i == matA.rows - 1){
                    cout << output[i];
                } else{
                    cout << output[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
            cout << "------------------------------------------------" << endl;
            delete[] output;
            exit(0);
//            interfacePartFour(matA, b);
        }
        case '2':{
            auto* output = new double[matA.rows];
            Solver s1;
            clock_t start = clock();
            s1.Gauss_Jordan(matA, b, output);
            clock_t end = clock();
            cout << endl;
            cout << "------------------------------------------------" << endl;
            cout << "The output solved by Gaussian Jordan:" << endl;
            cout << "[ ";
            for (int i = 0; i < matA.rows; i++) {
                if (i == matA.rows - 1){
                    cout << output[i];
                } else{
                    cout << output[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
            cout << "------------------------------------------------" << endl;
            delete[] output;
            exit(0);
//            interfacePartFour(matA, b);
        }
        case '3':{
            auto* output = new double[matA.rows];
            Solver s1;
            clock_t start = clock();
            s1.gauss_seidel(matA, b, output);
            clock_t end = clock();
            cout << endl;
            cout << "------------------------------------------------" << endl;
            cout << "The output solved by Gaussian Seidel:" << endl;
            cout << "[ ";
            for (int i = 0; i < matA.rows; i++) {
                if (i == matA.rows - 1){
                    cout << output[i];
                } else{
                    cout << output[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
            cout << "------------------------------------------------" << endl;
            delete[] output;
            exit(0);
//            interfacePartFour(matA, b);
        }
        case '4':{
            auto* output = new double[matA.rows];
            Solver s1;
            clock_t start = clock();
//            s1.LU_solver(mat_test, b_test, output);
            s1.LU_solver(matA, b, output);
            clock_t end = clock();
            cout << endl;
            cout << "------------------------------------------------" << endl;
            cout << "The output solved by LU Solver:" << endl;
            cout << "[ ";
            for (int i = 0; i < matA.rows; i++) {
                if (i == matA.rows - 1){
                    cout << output[i];
                } else{
                    cout << output[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
            cout << "------------------------------------------------" << endl;
            delete[] output;
            exit(0);
//            interfacePartFour(matA, b);
        }
        case '5':{
            auto* output = new double[matA.rows];
            Solver s1;
            clock_t start = clock();
            //            s1.LU_solver(mat_test, b_test, output);
            s1.LU_solver_partial_pivoting(matA, b, output);
            clock_t end = clock();
            cout << endl;
            cout << "------------------------------------------------" << endl;
            cout << "The output solved by LU Partial Pivoting Solver:" << endl;
            cout << "[ ";
            for (int i = 0; i < matA.rows; i++) {
                if (i == matA.rows - 1){
                    cout << output[i];
                } else{
                    cout << output[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
            cout << "------------------------------------------------" << endl;
            delete[] output;
            exit(0);
//            interfacePartFour(matA, b);
        }
        case '6':{
            int it_max = 1000;
            double it_tol = 1.e-6 ;
            auto* output = new double[matA.rows];
            Solver s1;
            clock_t start = clock();
            s1.template Jacobi_Method(matA, b, it_max, it_tol,output);
            clock_t end = clock();
            cout << endl;
            cout << "------------------------------------------------" << endl;
            cout << "The output solved by Jacobi Solver:" << endl;
            cout << "[ ";
            for (int i = 0; i < matA.rows; i++) {
                if (i == matA.rows - 1){
                    cout << output[i];
                } else{
                    cout << output[i] << ",";
                }
            }
            cout << "]" << endl;
            cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
            cout << "------------------------------------------------" << endl;
            delete[] output;
            exit(0);
//            interfacePartFour(matA, b);
        }
        case 'b':{
            // back to the start
            interfaceStart();
        }
        case 'x':{
            exit(0);
        }
        default:{
            cerr << "Invalid input";
        }
    }
}

void randomVectb(double *a, int n, double left, double right)//生成范围在l~r的随机数
{
    srand(time(0));
    for(int i=0;i<n;i++){
        a[i] = left + 1.0 * ( rand() % RAND_MAX ) / RAND_MAX * (right - left);
    }
}

void Interface::test_Gaussian_Elimination() {
    // correct output: [-6,5,-0.5]
    double values[9] = {
            2,3,-4,
            6,8,2,
            4,8,-6
    };
    Matrix<double> A(3,3,values);

    auto* b = new double[3] {5,3,19};

    auto* output = new double[3] {0,0,0};

    Solver s1;
    clock_t start = clock();
    s1.Gaussian_Elimination(A, b, output);
    clock_t end = clock();

    cout << "Input matrix A: ";
    A.printMatrix();

    cout << "Input vector b: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << b[i];
        } else{
            cout << b[i] << ",";
        }
    }
    cout << "]" << endl;

    cout << "The output solved by Gaussian Elimination: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << output[i];
        } else{
            cout << output[i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
    cout << "------------------------------------------------" << endl;

    delete[] b;
    delete[] output;
}

void Interface::test_Gaussian_Jordan() {
    // correct output: [1.86957,1.04348,-0.782609]
    double values[9] = {
            2,3,-4,
            3,-1,2,
            4,2,2
    };
    Matrix<double> A(3,3,values);

    auto* b = new double[3] {10,3,8};

    auto* output = new double[3] {0,0,0};

    Solver s1;
    clock_t start = clock();
    s1.Gauss_Jordan(A, b, output);
    clock_t end = clock();

    cout << "Input matrix A: ";
    A.printMatrix();

    cout << "Input vector b: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << b[i];
        } else{
            cout << b[i] << ",";
        }
    }
    cout << "]" << endl;

    cout << "The output solved by Gaussian Jordan: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << output[i];
        } else{
            cout << output[i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
    cout << "------------------------------------------------" << endl;

    delete[] b;
    delete[] output;
}

void Interface::test_LU_solver() {
    // correct output: [-6,5,-0.5]
    double values[9] = {
        2,3,-4,
        6,8,2,
        4,8,-6
    };
    Matrix<double> A(3,3,values);

    auto* b = new double[3] {5,3,19};

    auto* output = new double[3] {0,0,0};

    cout << "Input matrix A: ";
    A.printMatrix();

    Solver s1;
    clock_t start = clock();
    s1.LU_solver(A, b, output);
    clock_t end = clock();

    cout << "Input vector b: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << b[i];
        } else{
            cout << b[i] << ",";
        }
    }
    cout << "]" << endl;

    cout << "The output solved by LU Solver: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << output[i];
        } else{
            cout << output[i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
    cout << "------------------------------------------------" << endl;

    delete[] b;
    delete[] output;
}

void Interface::test_LU_partial_pivoting_solver() {
    // correct output: [1.86957,1.04348,-0.782609]
    double values[9] = {
            2,3,-4,
            6,8,2,
            4,8,-6
    };
    Matrix<double> A(3,3,values);

    auto* b = new double[3] {5,3,19};

    auto* output = new double[3] {0,0,0};

    Solver s1;
    clock_t start = clock();
    s1.LU_solver_partial_pivoting(A, b, output);
    clock_t end = clock();

    cout << "Input matrix A: ";
    A.printMatrix();

    cout << "Input vector b: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << b[i];
        } else{
            cout << b[i] << ",";
        }
    }
    cout << "]" << endl;

    cout << "The output solved by LU Partial Pivoting Solver: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << output[i];
        } else{
            cout << output[i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
    cout << "------------------------------------------------" << endl;

    delete[] b;
    delete[] output;
}

void Interface::test_Gaussian_Seidel() {
    // correct output: [-0.163408,-0.0153271,0.273353,0.368936]
    double values[16] = {
            10,2,3,5,
            1,14,6,2,
            -1,4,16,-4,
            5,4,3,11
    };
    Matrix<double> A(4,4,values);

    double* b = new double[4] {1,2,3,4};

    double* x = new double[4] {0,0,0,0};

    Solver s1;
    clock_t start = clock();
    s1.gauss_seidel<double>(A, b, x);
    clock_t end = clock();

    cout << "Input matrix A: ";
    A.printMatrix();

    cout << "Input vector b: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << b[i];
        } else{
            cout << b[i] << ",";
        }
    }
    cout << "]" << endl;

    cout << "The output solved by Gauss Seidel: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << x[i];
        } else{
            cout << x[i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
    cout << "------------------------------------------------" << endl;

    delete[] b;
    delete[] x;
}

void Interface::test_Gaussian_Elimination_CSRMatrix() {
    double values[9] = {
            2,3,-4,
            6,8,2,
            4,8,-6
    };
    CSRMatrix<double> A(3,3,values);

    auto* b = new double[3] {5,3,19};

    auto* output = new double[3] {0,0,0};

    Solver s1;
    clock_t start = clock();
    s1.Gaussian_Elimination_CSRMatrix(A, b, output);
    clock_t end = clock();

    cout << "Input matrix A in CSRMatrix format:" << endl;
    A.printMatrixCSR();

    cout << "Input vector b: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << b[i];
        } else{
            cout << b[i] << ",";
        }
    }
    cout << "]" << endl;

    cout << "The output solved by Gaussian Elimination CSRMatrix: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << output[i];
        } else{
            cout << output[i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
    cout << "------------------------------------------------" << endl;

    delete[] b;
    delete[] output;
}

void Interface::test_Gaussian_Seidel_CSRMatrix() {
    double values[16] = {
            10,2,3,5,
            1,14,6,2,
            -1,4,16,-4,
            5,4,3,11
    };
    CSRMatrix<double> A(4,4,values);

    double* b = new double[4] {1,2,3,4};

    double* x = new double[4] {0,0,0,0};

    Solver s1;
    clock_t start = clock();
    s1.Gauss_Seidel_CSRMatrix(A, b, x);
    clock_t end = clock();

    cout << "Input matrix A in CSRMatrix format:" << endl;
    A.printMatrixCSR();

    cout << "Input vector b: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << b[i];
        } else{
            cout << b[i] << ",";
        }
    }
    cout << "]" << endl;

    cout << "The output solved by Gauss Seidel CSRMatrix: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << x[i];
        } else{
            cout << x[i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
    cout << "------------------------------------------------" << endl;

    delete[] b;
    delete[] x;
}
void Interface::test_Jacobi_Method() {
    double values[16] = {
            10,2,3, 5,
            1,14,6, 2,
            -1, 4, 16, -4,
            5,4,3, 11
    };
    Matrix<double> A(4,4,values);

    double* b = new double[4] {1,2,3, 4};
    double* x = new double[4] {0,0,0,0};

    Solver s1;
    int it_max = 1000;
    double it_tol = 1.e-6 ;
    clock_t start = clock();
    s1.Jacobi_Method(A, b, it_max, it_tol, x);
    clock_t end = clock();

    cout << "Input matrix A in CSRMatrix format:";
    A.printMatrix();

    cout << "Input vector b: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << b[i];
        } else{
            cout << b[i] << ",";
        }
    }
    cout << "]" << endl;

    cout << "The output solved by Jacobi Method: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << x[i];
        } else{
            cout << x[i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
    cout << "------------------------------------------------" << endl;

    delete[] b;
    delete[] x;
}

void Interface::test_Jacobi_Method_CSRMatrix() {
    double *cvalues = new double[16]{10,2,3,5, 1, 14, 6, 2,
                                     -1, 4, 16, -4, 5, 4, 3, 11};

    CSRMatrix<double> A(4, 4, cvalues);

    double* b = new double[4] {1,2,3, 4};
    double* x = new double[4] {0,0,0,0};

    cout << "Input matrix A in CSRMatrix format:" << endl;
    A.printMatrixCSR();

    cout << "Input vector b: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << b[i];
        } else{
            cout << b[i] << ",";
        }
    }
    cout << "]" << endl;

    Solver s1;
    int it_max = 1000;
    double it_tol = 1.e-6;
    clock_t start = clock();
    s1.Jacobi_Method_CSRMatrix(A, b, it_max, it_tol, x);
    clock_t end = clock();

    cout << "The output solved by Jacobi Method CSRMatrix: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << x[i];
        } else{
            cout << x[i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
    cout << "------------------------------------------------" << endl;

    delete[] b;
    delete[] x;
}

void Interface::interfaceMatrixTest() {
    Interface anInterface;
    anInterface.test_Gaussian_Elimination();
    cout << endl;
    cout << endl;
    anInterface.test_Gaussian_Jordan();
    cout << endl;
    cout << endl;
    anInterface.test_LU_solver();
    cout << endl;
    cout << endl;
    anInterface.test_LU_partial_pivoting_solver();
    cout << endl;
    cout << endl;
    anInterface.test_Gaussian_Seidel();
    cout << endl;
    cout << endl;
    anInterface.test_Jacobi_Method();
    cout << endl;
    cout << endl;
    anInterface.test_Gaussian_Elimination_CSRMatrix();
    cout << endl;
    cout << endl;
    anInterface.test_Gaussian_Jordan_CSRMatrix();
    cout << endl;
    cout << endl;
    anInterface.test_Gaussian_Seidel_CSRMatrix();
    cout << endl;
    cout << endl;
    anInterface.test_Jacobi_Method_CSRMatrix();
    cout << endl;
    cout << endl;
    anInterface.test_LU_Solver_CSRMatrix();
}

void Interface::test_Gaussian_Jordan_CSRMatrix() {
    double values[9] = {
            2, 3, -4,
            3, -1, 2,
            4, 2, 2
    };
    CSRMatrix<double> A(3, 3, values);

    auto *b = new double[3]{10, 3, 8};

    auto *output = new double[3]{0, 0, 0};

    Solver s1;
    clock_t start = clock();
    s1.Gauss_Jordan_CSRMatrix(A, b, output);
    clock_t end = clock();

    cout << "Input matrix A in CSRMatrix format:" << endl;
    A.printMatrixCSR();

    cout << "Input vector b: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols - 1) {
            cout << b[i];
        } else {
            cout << b[i] << ",";
        }
    }
    cout << "]" << endl;

    cout << "The output solved by Gaussian Jordan CSRMatrix: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols - 1) {
            cout << output[i];
        } else {
            cout << output[i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;
    cout << "------------------------------------------------" << endl;

    delete[] b;
    delete[] output;
}

void Interface::test_LU_Solver_CSRMatrix()
{
    double values[9] = {
            2,3,-4,
            6,8,2,
            4,8,-6
    };
    CSRMatrix<double> A(3,3,values);

    auto* b = new double[3] {5,3,19};

    auto* output = new double[3] {0,0,0};

    cout << "Input matrix A in CSRMatrix format:" << endl;
    A.printMatrixCSR();

    cout << "Input vector b: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << b[i];
        } else{
            cout << b[i] << ",";
        }
    }
    cout << "]" << endl;

    Solver s1;
    clock_t start = clock();
    s1.LU_solver_CSRMatrix( A,  b, output);
    clock_t end = clock();

    cout << "The output solved by LU Solver CSRMatrix: [";
    for (int i = 0; i < A.cols; ++i) {
        if (i == A.cols-1){
            cout << output[i];
        } else{
            cout << output[i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "Time spent: " << (double) (end-start) / (double)(CLOCKS_PER_SEC) * 1000.0 << " milliseconds" << endl;

    delete[] b;
    delete[] output;
}


#include "Matrix.h"

#ifndef CPP_TUTORIAL_CSRMATRIX_H
#define CPP_TUTORIAL_CSRMATRIX_H

#pragma once
/**
    A sparse matrix class that can initialise sparse matrices (stored in CSR format) and 
    perform operations related to these matrices. Our sparse matrix class is templated so 
    that it can accept multiple types of input (eg int, doubles).
*/

template <class T>
class CSRMatrix: public Matrix<T>{
public:
    // ------------------------Constructors------------------------
    // default constructor
    CSRMatrix<T>(){}
    // constructor where we want to preallocate memory ourselves
    CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
    // constructor where already have allocated memory outside
    CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index);
    // constructor that can convert a Dense Matrix (with values given by the values pointer) to 
    // a Sparse (CSR) Matrix automatically
    CSRMatrix(int rows, int cols, T *values_ptr);
    // deep copy constructor
    CSRMatrix(const CSRMatrix<T> &A);
    // destructor
    ~CSRMatrix();
    // ------------------------Printing Functions------------------------
    // print the Matrix in CSR form
    virtual void printMatrixCSR();
    // override the function from parent class
    virtual void printMatrix();

    // ------------------------Functions for matrix operations------------------------
    // get the value of a matrix element given the specified row and column 
    T getValue(int row, int col) const;
    // set the value of a matrix element given the specified row and column
    void setValue(int row, int col, T value);

    // perform matrix vector multiplication
    void matVecMult(T *input, T *output);

    // random generation of CSR matrix
    CSRMatrix<T> randomCSRMatrix(int rowNum, int colNum, T min, T max);

    // creates an identity matrix of a specified size (nxn)
    CSRMatrix<T> identity(int n);

    // ------------------------Class variables------------------------
    int *row_position = nullptr;
    int *col_index = nullptr;
    // nnzs is the number of non-zeros
    int nnzs = -1;

private:

};


#endif //CPP_TUTORIAL_CSRMATRIX_H

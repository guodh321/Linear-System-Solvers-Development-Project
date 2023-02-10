#include <iostream>
#include<cmath>
#include "Matrix.h"
/**
    This file contains the definitions of the constructors and functions specified in the Matrix.h file.
*/
using namespace std;

// ------------------------Constructors------------------------
// constructor where we want to preallocate memory ourselves
template <class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate): rows(rows), cols(cols), size_of_values(rows * cols), preallocated(preallocate)
{
    if (this->preallocated)
    {
        this->values = new T[size_of_values];
    }
}

// constructor where we already have allocated memory outside
template <class T>
Matrix<T>::Matrix(int rows, int cols, T *values_ptr): rows(rows), cols(cols), size_of_values(rows * cols), values(values_ptr)
{}

// deep copy constructor
template<class T>
Matrix<T>::Matrix(const Matrix<T> &mat) {
    this->cols = mat.cols;
    this->rows = mat.rows;
    this->size_of_values = mat.rows * mat.cols;
    // deep copy
    T *new_values = new T[mat.size_of_values];
    for (int i = 0; i < mat.size_of_values; i++) {
        new_values[i] = mat.values[i];
    }
    this->values = new_values;
}

// creates a zero Matrix with specified dimensions
template<class T>
Matrix<T>::Matrix(int rows, int cols): rows(rows), cols(cols), size_of_values(rows * cols){
    this->values = new T [this->size_of_values];
    this->preallocated = true;
    for (int i=0; i < rows; i++){
        for (int j=0; j < cols; j++){
            this->values[i * cols + j] = 0.0;
        }
    }
}

// destructor
template <class T>
Matrix<T>::~Matrix()
{
    // Delete the values array
    if (this->preallocated){
        delete[] this->values;
    }
}

// print out the values in our matrix
template <class T>
void Matrix<T>::printValues()
{
    std::cout << "Printing values" << std::endl;
    for (int i = 0; i< this->size_of_values; i++)
    {
        std::cout << this->values[i] << " ";
    }
    std::cout << std::endl;
}

// print out the values in our matrix
template <class T>
void Matrix<T>::printMatrix()
{
    for (int j = 0; j< this->rows; j++)
    {
        std::cout << std::endl;
        for (int i = 0; i< this->cols; i++)
        {
            // We have explicitly used a row-major ordering here
            std::cout << this->values[i + j * this->cols] << " ";
        }
    }
    std::cout << std::endl;
}

// creates an identity matrix of a specified size (nxn)
template<class T>
Matrix<T> Matrix<T>::identity(int n) {
    // creating an nxn matrix
    Matrix<T> matE(n,n);
    int c = matE.cols;
    for (int i = 0; i < n; i++) {
        // setting all diagonal elements to 1 (creating identity matrix)
        matE.values[i * c + i] = 1.0;
    }
    return matE;
}

template<class T>
Matrix<T> Matrix<T>::randomMatrix(int rowNum, int colNum, T left, T right) {
    Matrix<T> matE(rowNum, colNum);

    srand(time(0));

    for (int i = 0; i < rowNum; i++) {
        for (int j = 0; j < colNum; j++) {
            matE.values[i * matE.cols + j] = left + 1.0 * ( rand() % RAND_MAX ) / RAND_MAX * (right - left);
        }
    }

    return matE;
}

// get the value of a matrix element given the specified row and column 
template<class T>
T Matrix<T>::getValue(int row, int col) const {
    return this->values[row * this->cols + col];
}

// set the value of a matrix element given the specified row and column
template<class T>
void Matrix<T>::setValue(int row, int col, T value) {
    this->values[row * this->cols + col] = value;
}

// calculate the norm of the Matrix
template<class T>
T Matrix<T>::matNorm() {
    T norm = 0;
    T temp;
    for (int i = 0; i < this->rows * this->cols; i++)
    {   // squaring matrix values and adding to norm variable
        temp = pow(this->values[i], 2);
        norm += temp;
    }
    // calulating the norm by finding the square root
    norm = sqrt(norm);
    return norm;
}

// ------------------------Functions for matrix operations------------------------
// Do matrix matrix multiplication
// output = mat_left * this
// m * k = m * n * n * k
template <class T>
void Matrix<T>::matMatMult(Matrix<T>& mat_left, Matrix<T>& output)
{

    // Check our dimensions match
    if (this->cols != output.cols)
    {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
        return;
    }

    // Check if our output matrix has had space allocated to it
    if (output.values != nullptr)
    {
        // Check our dimensions match
        if (this->rows != mat_left.cols || mat_left.rows != output.rows)
        {
            std::cerr << "Input dimensions for matrices don't match" << std::endl;
            return;
        }
    }
        // The output hasn't been preallocated, so we are going to do that
    else
    {
        output.values = new T[this->cols * mat_left.rows];
        output.preallocated = true;
    }

    // Set values to zero before hand
    for (int i = 0; i < output.size_of_values; i++)
    {
        output.values[i] = 0;
    }

    // Do matrix-matrix multiplication
    // Change this for loop ordering around
    for(int i = 0; i < mat_left.rows; i++)
    {
        for(int j = 0; j < this->cols; j++)
        {
            for(int k = 0; k < mat_left.cols; k++)
            {
                output.values[i * this->cols + j] += mat_left.values[i * mat_left.cols + k] * this->values[k * this->cols + j];
            }
        }
    }
}


// Do matrix vector multiplication
// output = this * vec
// m * 1 = m * n * n * 1
template <class T>
void Matrix<T>::matVecMult(const T* vec, T* output) {
    // set the output vector to 0 first
    for (int i = 0; i < this->cols; i++)
    {
        output[i] = 0;
    }
    // Do matrix-vector multiplication
    for (int i = 0; i < this->rows; i++) {
        for (int k = 0; k < this->cols; k++) {
            output[i] += this->values[i * this->cols + k] * vec[k];
        }
    }
}

// overload = to matrix assigment
template<class T>
void Matrix<T>::operator=(const Matrix<T> &A) {
    // assigning matrix to have values given by the input matrix values 
    for(int i = 0; i < A.rows * A.cols; i++){
        this->values[i] = A.values[i];
    }
}
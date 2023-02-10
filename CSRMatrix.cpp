#include <iostream>
#include "CSRMatrix.h"
/**
    This file contains the definitions of the constructors and functions specified in the CSRMatrix.h file.
*/
using  namespace std;

// ------------------------Constructors------------------------
// constructor where we want to preallocate memory ourselves
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate): Matrix<T>(rows, cols, false), nnzs(nnzs)
{
    this->preallocated = preallocate;
    this->size_of_values = nnzs;

    if (this->preallocated)
    {   // setting the values pointer, row_position and column_index to be related to the number of non zeros in the matrix
        this->values = new T[this->nnzs];
        this->row_position = new int[this->rows +1];
        this->col_index = new int[this->nnzs];
    }
}

// constructor where already have allocated memory outside
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index): Matrix<T>(rows, cols, values_ptr), nnzs(nnzs), row_position(row_position), col_index(col_index)
{
    // this is set to be rows*cols in our parent class Matrix
    // which is not correct for our CSRMatrix
    this->size_of_values = nnzs;
}

// constructor that can convert a Dense Matrix (with values given by the values pointer) to 
// a Sparse (CSR) Matrix automatically
template<class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, T *values_ptr): Matrix<T>(rows, cols, values_ptr)
{
    int i_nn = 0;
    this->nnzs = 0;
    int temp_row, temp_col;
    T *values_nn = new T[this->size_of_values];
    int *col_index_nn = new int[this->size_of_values];
    this->row_position = new int[this->rows + 1];

    // set the row_position to zero vector first
    for ( int i = 0; i < (this->rows + 1); i++)
    {
        row_position[i] = 0;
    }

    for (int i = 0; i < this->size_of_values; i++)
    {   // if statement to specify that we only want to store non zero elements from our dense matrix
        if (this->values[i] != 0)
        {
            // calculate the position in the dense array of the current element ( eg (3,2)))
            temp_col = i%this->cols;
            if ( temp_col == -1)
            {
                temp_col = this->cols - 1;
            }

            temp_row = (i + 1 - (temp_col + 1))/this->cols;
            
            // modifying the values vector, row position vector and column index vector of the sparse matrix 
            // given the value of the dense array
            values_nn[i_nn] = this->values[i];
            col_index_nn[i_nn] = temp_col;
            this->row_position[temp_row + 1]++;
            i_nn++;
            nnzs++;

        }

    }

    for (int i = 1; i < rows + 1; i++)
    {
        row_position[i] += row_position[i-1];
    }

    this->values = new T[this->nnzs];
    this->col_index = new int[this->nnzs];

    for (int i = 0; i < nnzs; i++)
    {
        this->values[i] = values_nn[i];
        this->col_index[i] = col_index_nn[i];
    }

    delete[] values_nn;
    delete[] col_index_nn;
}

// deep copy constructor
template<class T>
CSRMatrix<T>::CSRMatrix(const CSRMatrix<T> &mat) :Matrix<T>(mat.rows, mat.cols, true)
{

    this->cols = mat.cols;
    this->rows = mat.rows;
    this->size_of_values = mat.size_of_values;
    this->nnzs = mat.nnzs;
    // deep copy
    T *new_values = new T[mat.size_of_values];
    int *new_col_index = new int[mat.size_of_values];
    int *new_row_position = new int[mat.rows + 1];
    for (int i = 0; i < mat.size_of_values; i++) {
        new_values[i] = mat.values[i];
        new_col_index[i] = mat.col_index[i];
    }
    for (int i = 0; i < mat.rows + 1; i++){
        new_row_position[i] = mat.row_position[i];
    }

    this->col_index = new int[nnzs];
    this->values = new T[nnzs];
    this->row_position = new int[this->rows+1];

    for (int i = 0; i < mat.size_of_values; i++) {
        this->values[i] = new_values[i];
        this->col_index[i] = new_col_index[i];
    }
    for (int i = 0; i < mat.rows + 1; i++){
        this->row_position[i] = new_row_position[i];
    }

    delete[] new_values;
    delete[] new_col_index;
    delete[] new_row_position;

}

// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix<T>()
{
    if (this->preallocated)
    {
        // Matrix class destructor
        delete[] this->row_position;
        delete[] this->col_index;
    }
}

// ------------------------Printing Functions------------------------
// print the Matrix in CSR form
template <class T>
void CSRMatrix<T>::printMatrixCSR()
{
    std::cout << "Printing matrix";
    std::cout << "Values: ";
    for (int j = 0; j< this->nnzs; j++)
    {
        std::cout << this->values[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "row_position: ";
    for (int j = 0; j< this->rows+1; j++)
    {
        std::cout << this->row_position[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "col_index: ";
    for (int j = 0; j< this->nnzs; j++)
    {
        std::cout << this->col_index[j] << " ";
    }
    std::cout << std::endl;
}

// print the CSRMatrix in dense form
template<class U>
void CSRMatrix<U>::printMatrix()
{
    int row_begin, row_end, flag;
    for (int i = 0; i < this->rows; i++) {
        cout << endl;
        for (int j = 0; j < this->cols; j++)
        {
            flag = 0;
            row_begin = this->row_position[i];
            row_end = this->row_position[i+1];
            for (int k = row_begin; k < row_end; k++)
            {
                if (j == this->col_index[k])
                {
                    cout << this->values[k] << " ";
                    flag = 1;
                }
            }
            if (flag == 0)
            {
                cout << "0" << " ";
            }
        }

    }
    cout << endl;
}

// ------------------------Functions for matrix operations------------------------
// get the value of a matrix element given the specified row and column 
template<class T>
T CSRMatrix<T>::getValue(int row, int col) const
{
    T result;
    int row_begin, row_end, flag = 0;
    row_begin = this->row_position[row];
    row_end = this->row_position[row + 1];
    for (int k = row_begin; k < row_end; k++)
    {
        if ((col == this->col_index[k]) && (flag == 0)) {
            result = this->values[k];
            flag = 1;
        }
    }
    if (flag == 0) {
        result = 0;
    }

    return result;
}

// set the value of a matrix element given the specified row and column
template<class T>
void CSRMatrix<T>::setValue(int row, int col, T value)
{
    int row_begin, row_end, temp_position, change_position, flag = 0, nn;
    row_begin = row_position[row];
    row_end = row_position[row + 1];

    temp_position = row_begin;
    change_position = temp_position;

    for (int k = row_begin; k < row_end; k++)
    {
        if ((col > this->col_index[k]))
        {
            temp_position = row_begin + 1;
            change_position = temp_position;
            nn++;
        }
        if (col == this->col_index[k]) flag = 1;
    }

    if (flag == 1) {
        if (value != 0)
        {
            this->values[change_position]  = value;
        }
        else
        {
            this->nnzs -= 1;
            T *col_index_new = new T[nnzs];
            T *values_new = new T[nnzs];

            // change the row_position vector
            for (int r = (row + 1); r < (this->rows + 1); r++ )
            {
                row_position[r] -= 1;
            }

            // change the values vector and col_index vector
            for (int i = 0; i < change_position; i++)
            {
                values_new[i] = this->values[i];
                col_index_new[i] = this->col_index[i];
            }

            for (int i = change_position; i < this->nnzs; i++)
            {
                values_new[i] = this->values[i+1];
                col_index_new[i] = this->col_index[i+1];
            }

            // reset the new vector pointers
            delete[] this->values;
            delete[] this->col_index;

            this->col_index = new int[nnzs];
            this->values = new T[nnzs];

            for (int i = 0; i < nnzs; i++)
            {
                this->values[i] = values_new[i];
                this->col_index[i] = col_index_new[i];
            }

            delete[] values_new;
            delete[] col_index_new;

        }

    }

    if (flag == 0 && value != 0) {

        this->nnzs += 1;
        T *col_index_new = new T[nnzs];
        T *values_new = new T[nnzs];

        // change the row_position vector
        for (int r = (row + 1); r < (this->rows + 1); r++ )
        {
            row_position[r] += 1;
        }

        // change the values vector and col_index vector
        for (int i = 0; i < change_position; i++)
        {
            values_new[i] = this->values[i];
            col_index_new[i] = this->col_index[i];
        }

        values_new[change_position] = value;
        col_index_new[change_position] = col;

        for (int i = (change_position + 1); i < this->nnzs; i++)
        {
            values_new[i] = this->values[i-1];
            col_index_new[i] = this->col_index[i-1];
        }

        // reset the new vector pointers
        delete[] this->values;
        delete[] this->col_index;

        this->col_index = new int[nnzs];
        this->values = new T[nnzs];

        for (int i = 0; i < nnzs; i++)
        {
            this->values[i] = values_new[i];
            this->col_index[i] = col_index_new[i];
        }

        delete[] values_new;
        delete[] col_index_new;

    }

}

// perform matrix vector multiplication
template<class T>
void CSRMatrix<T>::matVecMult(T *input, T *output)
{
    if (input == nullptr || output == nullptr)
    {
        std::cerr << "Input or output haven't been created" << std::endl;
        return;
    }

    // Set the output to zero
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0.0;
    }

    int val_counter = 0;
    // loop over each row
    for (int i = 0; i < this->rows; i++)
    {
        // loop over all the entries in this col
        for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
        {
            output[i] += this->values[val_index] * input[this->col_index[val_index]];

        }
    }
}

template<class T>
CSRMatrix<T> CSRMatrix<T>::randomCSRMatrix(int rowNum, int colNum, T min, T max) {
    // values array
    int n = rowNum * colNum;
    T a[n];

    // random fill values in the array
    srand(time(0));
    for(int i=0; i < n; i++){
        if (1 == rand() % 10){
            // since a sparse matrix is a matrix that is comprised of mostly zero values
            // there are 10% chance of fill random value between min and max
            a[i] = min + 1.0 * ( rand() % RAND_MAX ) / RAND_MAX * (max - min);
        } else{
            // 90% chance of generating zero
            a[i] = 0;
        }
    }

    CSRMatrix<T> matE(rowNum, colNum, a);

    return matE;
}

// creates an identity matrix of a specified size (nxn)
template<class T>
CSRMatrix<T> CSRMatrix<T>::identity(int n) {
    T *i_values = new T[n*n];
    for (int i = 0; i < n; i++) {
        // setting all diagonal elements to 1 (creating identity matrix)
        i_values[i * n + i] = 1.0;
    }

    // creating an nxn CSRMatrix
    CSRMatrix<T> matE(n,n, i_values);

    delete[] i_values;
    return matE;
}

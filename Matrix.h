

#pragma once
/**
    A matrix class that can initialise matrices and perform operations related to matrices.
    Our matrix class is templated so that it can accept multiple types of input (eg int, doubles).
*/
using namespace std;

template <class T>
class Matrix
{
public:
    // ------------------------Constructors------------------------
    // default constructor
    Matrix(){}
    // constructor where we want to preallocate memory ourselves
    Matrix(int rows, int cols, bool preallocate);
    // constructor where we already have allocated memory outside
    Matrix(int rows, int cols, T *values_ptr);
    // deep copy constructor
    Matrix(const Matrix<T> &A);
    // creates a zero Matrix with specified dimensions
    Matrix(int rows,int cols);
    // destructor
    virtual ~Matrix();

    // ------------------------Printing Functions------------------------
    // print out the values in our matrix
    void printValues();
    virtual void printMatrix();

    // ------------------------Functions for matrix operations------------------------
    // creates an identity matrix of a specified size (nxn)
    Matrix<T> identity(int n);
    // creates a matrix of random values with specified size and the minimum 
    // and maximum size of any element required
    Matrix<T> randomMatrix(int rowNum, int colNum, T min, T max);

    // get the value of a matrix element given the specified row and column 
    T getValue(int row, int col) const;
    // set the value of a matrix element given the specified row and column
    void setValue(int row, int col, T value);
    // calculate the norm of the Matrix
    T matNorm( );

    // perform matrix matrix multiplication
    void matMatMult(Matrix<T>& mat_left, Matrix<T>& output);
    // perform matrix vector multiplication
    void matVecMult(const T* vec, T* ouput);


    // overload = to matrix assigment
    void operator=(const Matrix<T> &A);

    //friend functions
    //R = A + B
    friend Matrix<T> operator+(const Matrix<T> & A, const Matrix<T> & B){

        Matrix<T> R(A.rows, A.cols, true);

        for(int i = 0; i < A.rows * B.cols; i++)
        {
            R.values[i] = A.values[i] + B.values[i];
        }

        return R;
    }

    // R = A - B
    friend Matrix<T> operator-(const Matrix<T> & A, const Matrix<T> & B)
    {

        Matrix<T> R(A.rows, A.cols, true);

        for(int i = 0; i < A.rows * B.cols; i++)
        {
            R.values[i] = A.values[i] - B.values[i];
        }
        return R;
    }

    // overloading the matmatMult operation
    // R = A@B
    friend Matrix<T> operator*(const Matrix<T> & A, const Matrix<T> & B)
    {
        Matrix<T> R(A.rows, B.cols, true);

        // Set values to zero before hand
        for (int i = 0; i < R.size_of_values; i++)
        {
            R.values[i] = 0;
        }

        // Change this for loop ordering around
        for(int i = 0; i < A.rows; i++)
        {
            for(int j = 0; j < B.cols; j++)
            {
                for(int k = 0; k < A.cols; k++)
                {
                    R.values[i * B.cols + j] += A.values[i * A.cols + k] * B.values[k * B.cols + j];
                }
            }
        }
        return R;
    }

    // ------------------------Class variablesie------------------------
    // Initialising the values pointer to be a nullpointer (until changed via values being inputted)
    T *values = nullptr;
    int rows = -1;
    int cols = -1;
    int size_of_values = -1;

protected:
    bool preallocated = false;

};

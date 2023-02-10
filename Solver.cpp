#include <iostream>
#include "Solver.h"
/**
    This file contains the definitions of the constructors and functions specified in the Solver.h file.
*/
using namespace std;

// --------------auxiliary function-----------------
template<class U>
void Solver::back_substitution(Matrix<U> &A, U* b, U* output) {
    // n is the size of vector b which is the same as matrix A's column number
    int n = A.rows;
    for (int k = n-1; k != -1; k--){
        // loop from the end to front (backwards)
        U sum = 0;
        for (int j = k+1; j != n; j++){
            sum = sum + A.values[k * A.cols + j] * output[j];
        }
        output[k] = (b[k] - sum) / A.values[k * A.cols + k];
    }
}

template<class U>
void Solver::back_substitution_CSRMatrix(CSRMatrix<U> &A, U* b, U* output) {
    // n is the size of vector b which is the same as matrix A's column number
    int n = A.rows;
    for (int k = n-1; k != -1; k--){
        // loop from the end to front (backwards)
        U sum = 0;
        for (int j = k+1; j != n; j++){
            sum = sum + A.getValue(k,j) * output[j];
        }
        output[k] = (b[k] - sum) / A.getValue(k,k);
    }
}

template<class U>
void Solver::forward_substitution(Matrix<U> &A, U *b, U *output) {
    // n is the size of vector b which is the same as matrix A's column number
    int n = A.rows;
    for (int k = 0; k != n; k++){
        // loop from the end to front (backwards)
        U sum = 0;
        for (int j = 0; j != k; j++){
            sum = sum + A.values[k * A.cols + j] * output[j];
        }
        output[k] = (b[k] - sum) / A.values[k * A.cols + k];
    }
}

template<class U>
void Solver::upper_triangle(Matrix<U> &A, U* b) {
    // n is the size of vector b which is the same as matrix A's column number
    int n = (int) A.cols;
    // Loop through each pivot row
    // except the last one, which will never be used as a pivot.
    for (int k = 0; k != n-1; k++){
        for ( int i = k+1; i != n; i++){
            // Define the scaling factor for this row outside the innermost loop otherwise
            U sum = A.values[i * A.cols + k] / A.values[k * A.cols + k];
            for ( int j = k; j != n; j++){
                // update the current row of A
                A.values[i * A.cols + j] -= sum * A.values[k * A.cols + j];
            }
            // update entries of vector b
            b[i] -= sum * b[k];
        }
    }
}

template<class U>
void Solver::upper_triangle_CSRMatrix(CSRMatrix<U> &A, U* b) {
    // n is the size of vector b which is the same as matrix A's column number
    int n = (int) A.cols;
    // Loop through each pivot row
    // except the last one, which will never be used as a pivot.
    for (int k = 0; k != n-1; k++){
        for ( int i = k+1; i != n; i++){
            // Define the scaling factor for this row outside the innermost loop otherwise
            U sum = A.getValue(i,k)/ A.getValue(k ,k);
            for ( int j = k; j != n; j++){
                // update the current row of A
                A.setValue(i,j,A.getValue(i,j)-sum * A.getValue(k,j));
            }
            // update entries of vector b
            b[i] -= sum * b[k];
        }
    }
}

template<class U>
void Solver::swap_row(Matrix<U> &mat, U* b,int cur_row, int max_row) {
    U tmp;
    U tmp1;
    // swapping matrix A
    for (int j = 0; j < mat.cols; j++){
        // temporary store the value to avoid overwriting it
        tmp = mat.values[cur_row * mat.cols + j];
        mat.values[cur_row * mat.cols + j] = mat.values[max_row * mat.cols + j];
        mat.values[max_row * mat.cols + j] = tmp;
    }
    // swapping vector b
    tmp1 = b[cur_row];
    b[cur_row] = b[max_row];
    b[max_row] = tmp1;
}

template<class U>
void Solver::swap_row(Matrix<U> &matrix, int current_row, int max_row) {
    U temp;
    for (int j = 0; j < matrix.cols; j++)
    {
        temp = matrix.values[current_row * matrix.cols + j];
        matrix.values[current_row * matrix.cols + j] = matrix.values[max_row * matrix.cols + j];
        matrix.values[max_row * matrix.cols + j] = temp;
    }
}

template<class U>
void Solver::swap_row_CSRMatrix(CSRMatrix<U> &matrix, int current_row, int max_row) {
    U temp;
    for (int j = 0; j < matrix.cols; j++)
    {
        temp = matrix.values[current_row * matrix.cols + j];
        matrix.values[current_row * matrix.cols + j] = matrix.values[max_row * matrix.cols + j];
        matrix.values[max_row * matrix.cols + j] = temp;
    }
}

template<class U>
void Solver::upper_triangle_pp(Matrix<U> &A, U* b) {
    int num = (int) A.cols;
    U s;
    for (int k = 0; k < num - 1; k++) {
        // Swap the rows so that we're always dividing by the highest number.
        int kmax = k; // initiatise kmax with the current pivot row (k)
        for (int i = k + 1; i < num; i++) {
            // Select the k with the highest abs value from all entries below the pivot.
            if (abs(A.values[kmax * A.cols + k]) < abs(A.values[i * A.cols + k])){
                kmax = i;
            }
        }
        //replace the current pivot row (k) with the row with the biggest abs value.
        swap_row(A, b, kmax, k);
        for (int i = k + 1; i < num; i++) {
            s = A.values[i * A.cols + k] / A.values[k * A.cols + k];
            for (int j = k; j < num; j++) {
                A.values[i * A.cols + j] -= A.values[k * A.cols + j] * s;
            }
            b[i] -= s * b[k];
        }
    }
}

template<class U>
void Solver::upper_triangle_pp_CSRMatrix(CSRMatrix<U> &A, U* b) {
    int num = (int) A.cols;
    U s;
    for (int k = 0; k < num - 1; k++) {
        // Swap the rows so that we're always dividing by the highest number.
        int kmax = k; // initiatise kmax with the current pivot row (k)
        for (int i = k + 1; i < num; i++) {
            // Select the k with the highest abs value from all entries below the pivot.
            if (abs(A.getValue(kmax,k)) < abs(A.getValue(i,k))){
                kmax = i;
            }
        }
        //replace the current pivot row (k) with the row with the biggest abs value.
        swap_row(A, b, kmax, k);
        for (int i = k + 1; i < num; i++) {
            s = A.getValue(i,k)/ A.getValue(k,k);
            for (int j = k; j < num; j++) {
                A.setValue(i,j,A.getValue(i,j)-A.getValue(k,j) * s);
            }
            b[i] -= s * b[k];
        }
    }
}

template<class U>
U Solver::gauss_seidel_addition(Matrix<U> &A, U *b,U * x, int i) {
    //initialisng a total value to be zero 
    U total = 0;
    int c = A.cols;
    int r = A.rows;
    if(i > 0){
        //looping over columns less than i and calculating required algorithmic values
        for(int a = 0; a < i; a++){
            total += A.values[i*c+a]* x[a];
        }
    }
    //looping over columns greater than i and calculating required algorithmic values
    for (int b = i + 1; b < r; b++) {
        total += A.values[i*c+b] * x[b];
    }

    return total;
}

template<class U>
U Solver::gauss_seidel_addition_CSRMatrix(CSRMatrix<U> &A, U *b,U *x, int i) {
    U total = 0;
    int c = A.cols;
    int r = A.rows;
    if(i > 0){
        for(int a = 0; a < i; a++){
            total += A.values[i*c+a]* x[a];
        }
    }

    for (int b = i + 1; b < r; b++) {
        total += A.values[i*c+b] * x[b];
    }

    return total;
}

// ---------------------Dense linear solvers---------------------
// Gaussian Elimination Solver
template<class U>
void Solver::Gaussian_Elimination(const Matrix<U> &matA, const U* vecb, U* output) {
    // assume input matrix is a square matrix
    int n = (int) matA.cols;

    // deep copy of matrix A
    Matrix<U> A(matA);

    // deep copy of vector b
    U* b = new U[n];
    for (int i = 0; i < n; i++) {
        b[i] = vecb[i];
    }

    // perform operations on A and b
    // store results in output array
    upper_triangle(A, b);
    back_substitution(A, b, output);
}

// Gauss Jordan Solver
template<class U>
void Solver::Gauss_Jordan(const Matrix<U> &matA, const U* vecb, U* output) {
    // assume input matrix is a square matrix
    int n = (int) matA.cols;

    // copy matrix
    Matrix<U> A(matA);

    // copy vector
    U* b = new U[n];
    for (int i = 0; i < n; i++) {
        b[i] = vecb[i];
    }

    // perform operations on A and b
    // store results in output array
    upper_triangle_pp(A, b);
    back_substitution(A, b, output);

    delete[] b;
}

// LU Decomposition Solver
template<class U>
void Solver::LU_decomposition(Matrix<U>& matA, Matrix<U>& matL) {
    int n = matA.cols;
    int m = matL.cols;
    for (int k = 0; k < n - 1; k++){
        for (int i = k + 1; i < n; i++){
            U sum = matA.values[i * n + k] / matA.values[k * n + k];
            for (int j = k; j != n; j++)
                matA.values[i * n + j] -= sum * matA.values[k * n + j];
            matL.values[i * m + k] = sum;
        }
    }
    matL = matL + matL.identity(n);
}

// LU Solver
template<class U>
void Solver::LU_solver( Matrix<U> &matA,  U* vecb, U* output) {
    // assume input matrix is a square matrix
    int n = (int) matA.rows;

    // copy matrix
    Matrix<U> A(matA);

    // copy vector
    U* b = new U[n];
    for (int i = 0; i < n; i++) {
        b[i] = vecb[i];
    }

    // initialize L into nXn matrix with zero
    Matrix<U> matL(n, n);

    // apply operations on matrix A and matrix L in place
    LU_decomposition(matA, matL);

    // declare an array for storing the output from forward substitution
    U* d = new U[n];
    forward_substitution(matL, b, d);
    back_substitution(matA, d, output);

    delete[] d;
    delete[] b;
}

template<class U>
void Solver::LU_decomposition_partial_pivoting(Matrix<U> &A, Matrix<U> &matP, Matrix<U> &matL, Matrix<U> &matU) {
    // square matrix is cheked by the call function
    int m = A.rows;

    matP = A.identity(m);

    for ( int k = 0; k < m-1; k++){
        vector<U> vec1;
        for (int q = k; q != m; q++){
            vec1.push_back(abs(A.values[q * A.cols + k]));
        }

        // find the max element in vector
        int j = 0; // the index of max element
        int max = 0;
        for (int i = 0; i < vec1.size(); i++) {
            if (vec1[i] > max){
                max = vec1[i];
                j = i;
            }
        }
        j += k;

        // swap operations
        swap_row(A, j, k);
        swap_row(matP, j, k);
        swap_row(matL, j, k);
        for (int i = k+1; i < m; i++){
            U sum = A.values[i * A.cols + k] / A.values[k * A.cols + k];
            for (int q = k; q != m; q++)
                A.values[i * A.cols + q] -= A.values[k * A.cols + q] * sum;
            matL.values[i * matL.cols + k] = sum;
        }
    }

    matU = A;
    matL = matL + A.identity(m);
}

// LU factorisation with partial pivoting
template<class U>
void Solver::LU_solver_partial_pivoting(Matrix<U>& matA, U* vecb, U* output) {
    // assume input matrix is a square matrix
    int n = (int) matA.rows;

    // copy matrix
    Matrix<U> A(matA);

    // copy vector
    U* b = new U[n];
    for (int i = 0; i < n; i++) {
        b[i] = vecb[i];
    }

    // initialize L into nXn matrix with zero
    Matrix<U> matL(n, n);
    Matrix<U> matU(n, n);
    Matrix<U> matP(n, n);

    // apply operations on matrix A and matrix L in place
    LU_decomposition_partial_pivoting(A, matP, matL, matU);

    // declare an array for storing the output from forward substitution
    U* d = new U[n];
    U* c = new U[n];

    // dot product
    U temp = 0;
    for (int i = 0; i < matP.rows; i++)
    {
        temp = 0;
        for (int j = 0; j < matP.cols; j++)
            temp += matP.values[i * matP.cols + j] * vecb[j];
        c[i] = temp;
    }

    forward_substitution(matL, c, d);
    back_substitution(matU, d, output);

    delete[] d;
    delete[] b;
    delete[] c;
}

// Gauss Seidel Solver
template<class U>
void Solver::gauss_seidel(Matrix<U>& A,U *b,U *x) {
    int m=A.rows;
    int n=A.cols;
    int maxit=500;
    if (b == nullptr || x == nullptr)
    {
        std::cout << "Input or output haven't been created" << std::endl;
        return;
    }

    for(int k=0;k<maxit;k++){
        //looping until maxmimum number of iterations reached
        for(int i=0;i<m;i++){
            // note that here we will be re-using the x as we update it each i
            //calculating the next value in the iteration and storing it
            x[i] = ( b[i]- gauss_seidel_addition(A,b,x,i) )/A.values[i*n+i];

        }

    }
}

// Jacobi Method
template<class U>
void Solver::Jacobi_Method(Matrix<U> & matA, const  U* vecb, int it_max, double it_tol, U* output)
{
    int n = (int) matA.rows;
    U temp_value, residual;

    Matrix<U> M_matB = Matrix<U>(matA.rows, matA.cols, true);
    U* vecg = new U [matA.rows];
    U* vecx = new U [matA.rows];
    U* vecx_new = new U [matA.rows];
    U* vec_resultOfMul1 = new U [matA.rows];
    U* vec_resultOfMul2 = new U [matA.rows];
    U* vec_resultOfAdd = new U [matA.rows];

// calculate B
    for(int i = 0;i < n;i++)
    {
        for(int j = 0;j < n;j++)
        {
            if(i != j)
            {
                temp_value = -1*(matA.getValue(i, j))/matA.getValue(i, i);
                M_matB.setValue(i, j, temp_value);
            }
            else
                M_matB.setValue(i, j, 0);
        }
    }

// calculate g
    for(int i = 0; i < n; i++)
    {
        temp_value = vecb[i] / matA.getValue(i, i);
        vecg[i] = temp_value;
    }

// initialize vectors to 0
    for(int i = 0;i < n;i++)
    {
        vecx[i] = 0;
        vecx_new[i] = 0;
    }

// implement the iterative methods
    for (int k = 0; k < it_max; k++)
    {
        // x_new = g + (B @ x)
        M_matB.matVecMult( vecx, vec_resultOfMul1);
        for (int i = 0; i < n; i++)
        {
            vecx_new[i] = vecg[i] + vec_resultOfMul1[i];
        }

        for (int i = 0; i < n; i++)
        {
            vec_resultOfMul1[i] = 0;
        }

        // calculate residual
        // residual = norm((A @ x_new) - b)
        residual = 0;
        matA.matVecMult( vecx_new, vec_resultOfMul2);

        for (int i = 0; i < n; i++)
        {
            vec_resultOfAdd[i] =  vec_resultOfMul2[i] -vecb[i];
        }

        for (int i = 0; i < n; i++)
        {
            temp_value = pow(vec_resultOfAdd[i], 2);
            residual += temp_value;
        }

        for (int i = 0; i < n; i++)
        {
            vec_resultOfMul2[i] = 0;
            vec_resultOfAdd[i] = 0;
        }

        residual = sqrt(residual);

        // if residual < it_tol, end iteration
        if ( residual <= it_tol)
        {
            break;
        }

        for (int i = 0; i < n; i++)
        {
            vecx[i] = vecx_new[i];
        }

    }

    for (int i = 0; i < n; i++)
    {
        output[i] = vecx_new[i];

    }

    delete[] vecg;
    delete[] vecx;
    delete[] vecx_new;
    delete[] vec_resultOfMul1;
    delete[] vec_resultOfMul2;
    delete[] vec_resultOfAdd;

}

// ---------------------Sparse linear solvers---------------------
// Gaussian Elimination Solver for CSRMatrix
template<class U>
void Solver::Gaussian_Elimination_CSRMatrix(CSRMatrix<U>& matA, const  U* vecb, U* output)
{
// assume input matrix is a square matrix
    int n = (int) matA.cols;

// deep copy of matrix A
    CSRMatrix<U> A(matA);

// deep copy of vector b
    U* b = new U[n];
    for (int i = 0; i < n; i++) {
        b[i] = vecb[i];
    }

// perform operations on A and b
// store results in output array
    upper_triangle_CSRMatrix(matA, b);
    back_substitution_CSRMatrix(matA, b, output);

    delete[] b;
}

// Gauss Jordan Solver for CSRMatrix
template<class U>
void Solver::Gauss_Jordan_CSRMatrix(CSRMatrix<U> &matA, const U* vecb, U* output) {
// assume input matrix is a square matrix
    int n = (int) matA.cols;

// copy matrix
    CSRMatrix<U> A(matA);

// copy vector
    U* b = new U[n];
    for (int i = 0; i < n; i++) {
        b[i] = vecb[i];
    }

// perform operations on A and b
// store results in output array
    upper_triangle_pp_CSRMatrix(A, b);
    back_substitution_CSRMatrix(A, b, output);

    delete[] b;
}

// Gauss Seidel Solver for CSRMatrix
template<class U>
void Solver::Gauss_Seidel_CSRMatrix(CSRMatrix<U> &A, U *b, U *x)
{
    int m=A.rows;
    int n=A.cols;
    int maxit=500;
    if (b == nullptr || x == nullptr)
    {
        std::cout << "Input or output haven't been created" << std::endl;
        return;
    }

    for(int k=0;k<maxit;k++){
        //while loop here instead and counter
        for(int i=0;i<m;i++){
            // note that here we will be re-using the x as we update it each i
            x[i] = ( b[i]- gauss_seidel_addition_CSRMatrix(A,b,x,i) )/A.values[i*n+i];
        }

    }
}

// Jacobi Method for CSRMatrix
template<class U>
void Solver::Jacobi_Method_CSRMatrix(CSRMatrix<U> &matA, const U *vecb, int it_max, double it_tol, U *output)
{
    int n = (int) matA.rows;
    U temp_value, residual;

    Matrix<U> M_matB = Matrix<U>(matA.rows, matA.cols, true);
    U* vecg = new U [matA.rows];
    U* vecx = new U [matA.rows];
    U* vecx_new = new U [matA.rows];
    U* vec_resultOfMul1 = new U [matA.rows];
    U* vec_resultOfMul2 = new U [matA.rows];
    U* vec_resultOfAdd = new U [matA.rows];

// calculate B
    for(int i = 0;i < n;i++)
    {
        for(int j = 0;j < n;j++)
        {
            if(i != j)
            {
                temp_value = -1*(matA.getValue(i, j))/matA.getValue(i, i);
                M_matB.setValue(i, j, temp_value);
            }
            else
                M_matB.setValue(i, j, 0);
        }
    }

// calculate g
    for(int i = 0; i < n; i++)
    {
        temp_value = vecb[i] / matA.getValue(i, i);
        vecg[i] = temp_value;
    }


// initialize vectors to 0
    for(int i = 0;i < n;i++)
    {
        vecx[i] = 0;
        vecx_new[i] = 0;
    }

// implement the iterative methods
    for (int k = 0; k < it_max; k++)
    {
        // x_new = g + (B @ x)
        M_matB.matVecMult( vecx, vec_resultOfMul1);
        for (int i = 0; i < n; i++)
        {
            vecx_new[i] = vecg[i] + vec_resultOfMul1[i];
        }

        for (int i = 0; i < n; i++)
        {
            vec_resultOfMul1[i] = 0;
        }

        // calculate residual
        // residual = norm((A @ x_new) - b)
        residual = 0;
        matA.matVecMult(vecx_new, vec_resultOfMul2);

        for (int i = 0; i < n; i++)
        {
            vec_resultOfAdd[i] =  vec_resultOfMul2[i] -vecb[i];
        }

        for (int i = 0; i < n; i++)
        {
            temp_value = pow(vec_resultOfAdd[i], 2);
            residual += temp_value;
        }

        for (int i = 0; i < n; i++)
        {
            vec_resultOfMul2[i] = 0;
            vec_resultOfAdd[i] = 0;
        }

        residual = sqrt(residual);

        // if residual < it_tol, end iteration
        if ( residual <= it_tol)
        {
            break;
        }

        for (int i = 0; i < n; i++)
        {
            vecx[i] = vecx_new[i];
        }

    }

    for (int i = 0; i < n; i++)
    {
        output[i] = vecx_new[i];

    }

    delete[] vecg;
    delete[] vecx;
    delete[] vecx_new;
    delete[] vec_resultOfMul1;
    delete[] vec_resultOfMul2;
    delete[] vec_resultOfAdd;

}

// LU Solver for CSRMatrix
template<class U>
void Solver::LU_solver_CSRMatrix(CSRMatrix<U>& matA,  U* vecb, U* output) {
    // assume input matrix is a square matrix
    int n = (int) matA.rows;

    // copy matrix
//    CSRMatrix<U> A(matA);

    // copy vector
    U* b = new U[n];
    for (int i = 0; i < n; i++) {
        b[i] = vecb[i];
    }

    // initialize L into nXn matrix with zero
    U* zero_values = new U[n * n];
    for (int i = 0; i < n * n; i++) {
        zero_values[i] = 0;
    }
    CSRMatrix<U> matL(n, n, zero_values);

    // apply operations on matrix A and matrix L in place
    LU_decomposition_CSRMatrix(matA, matL);

    // declare an array for storing the output from forward substitution
    U* d = new U[n];
    forward_substitution_CSRMatrix(matL, b, d);
    back_substitution_CSRMatrix(matA, d, output);

    delete[] zero_values;
    delete[] d;
    delete[] b;
}

// LU decomposition CSRMatrix
template<class U>
void Solver::LU_decomposition_CSRMatrix(CSRMatrix<U>& matA, CSRMatrix<U>& matL) {
    int n = matA.cols;
    int m = matL.cols;
    U sum, temp_values1, temp_values2;
    for (int k = 0; k < n - 1; k++){
        for (int i = k + 1; i < n; i++){
            sum = matA.values[i * n + k] / matA.values[k * n + k];
            for (int j = k; j != n; j++)
            {
                matA.values[i * n + j] -= sum * matA.values[k * n + j];
            }
            matL.setValue(i, k, sum);
        }
    }

    for ( int i = 0; i < n; i++)
    {
        temp_values2 = matL.getValue(i, i) + 1;
        matL.setValue(i, i, temp_values2);
        temp_values2 = 0;
    }
}

template<class U>
void Solver::forward_substitution_CSRMatrix(CSRMatrix<U> &A, U *b, U *output) {
    // n is the size of vector b which is the same as matrix A's column number
    int n = A.rows;
    U temp_value;
    for (int i = 0; i< n;i++)
        output[i] = 0;
    for (int k = 0; k != n; k++){
        // loop from the end to front (backwards)
        U sum = 0;
        for (int j = 0; j != k; j++){
            temp_value = A.getValue(k, j) ;
            sum = sum + temp_value * output[j];
        }
        output[k] = (b[k] - sum) / A.getValue(k, k);
    }
}

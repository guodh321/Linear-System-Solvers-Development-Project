#include "Matrix.h"
#include "CSRMatrix.h"
#include<vector>

/**
    A solver class that has functions that correspond to various linear solver algorithms.
*/
using namespace std;

#ifndef GROUP_ASSIGNMENT_LSG_SOLVER_H
#define GROUP_ASSIGNMENT_LSG_SOLVER_H

class Solver {
public:

    // --------------auxiliary function-----------------
    /**
     * @tparam U
     * @param A The input square matrix on upper triangular form
     * @param b The vector b from upper_triangle(Matrix<U> & A, U* b)
     * @param output The output array to store the vector x in linear system Ax=b
     */
    template <class U>
    void back_substitution(Matrix<U>& A, U* b, U* output);

    /**
    * @tparam U
    * @param A The input square matrix on upper triangular form
    * @param b The vector b from upper_triangle(Matrix<U> & A, U* b)
    * @param output The output array to store the vector x in linear system Ax=b
    */
    template <class U>
    void back_substitution_CSRMatrix(CSRMatrix<U>& A, U* b, U* output);

    /**
     * @tparam U
     * @param A The input square matrix on upper triangular form
     * @param b The vector b from upper_triangle(Matrix<U> & A, U* b)
     * @param output The output array to store the vector x in linear system Ax=b
     */
    template <class U>
    void forward_substitution(Matrix<U>& A, U* b, U* output);

    /**
     * @tparam U
     * @param A The input square matrix A that will be converted into upper triangular form through row operations
     * @param b The same row operations are performed on the vector b
     */
    template<class U>
    void upper_triangle(Matrix<U>& A, U* b);

    /**
     * @tparam U
     * @param A The input square matrix A that will be converted into upper triangular form through row operations
     * @param b The same row operations are performed on the vector b
     */
    template<class U>
    void upper_triangle_CSRMatrix(CSRMatrix<U>& A, U* b);

    /**
     * This version takes a matrix and an array (vector)
     * Used in Gaussian Jordan solver
     * @tparam U
     * @param mat Input matrix A
     * @param b Input vector b
     * @param cur_row Swap rows cur_row and max_row
     * @param max_row Swap rows cur_row and max_row
     */
    template<class U>
    void swap_row(Matrix<U>& mat, U* b, int cur_row, int max_row);

    /**
     This version takes three parameters, which is used in LU decomposition pp
     * @tparam U
     * @param matrix Input matrix
     * @param current_row Swap rows cur_row and max_row
     * @param max_row Swap cur_row and max_row
     */
    template<class U>
    void swap_row(Matrix<U>& matrix, int current_row, int max_row);

    template<class U>
    void swap_row_CSRMatrix(CSRMatrix<U>& matrix, int current_row, int max_row);

    /** LU decomposition
    * @tparam U
    * @param matA Input matrix A which we would apply operations to change it into upper triangular matrix U in place
    * @param matL Input matrix L is initially zero matrix and will be changed into lower triangular matrix L in place
    */
    template <class U>
    void LU_decomposition(Matrix<U> & matA, Matrix<U> & matL);

    /**
     * This version uses partial pivoting
     * @tparam U
     * @param A The input square matrix A that will be converted into upper triangular form through row operations
     * @param b The same row operations are performed on the vector b
     */
    template<class U>
    void upper_triangle_pp(Matrix<U> & A, U* b);

    template<class U>
    void upper_triangle_pp_CSRMatrix(CSRMatrix<U> & A, U* b);

    /**
     * @tparam U
     * @param matA
     * @param matP
     * @param matL
     * @param matU
     */
    template <class U>
    void LU_decomposition_partial_pivoting(Matrix<U>& matA, Matrix<U>& matP, Matrix<U>& matL, Matrix<U>& matU);

    template <class U>
    U gauss_seidel_addition(Matrix<U>& A, U *b, U *x, int i);

    template <class U>
    U gauss_seidel_addition_CSRMatrix(CSRMatrix<U> &A, U *b, U *x, int i);

    template <class U>
    void LU_decomposition_CSRMatrix(CSRMatrix<U> & matA, CSRMatrix<U> & matL);

    template <class U>
    void forward_substitution_CSRMatrix(CSRMatrix<U>& A, U* b, U* output);

    // ---------------------Dense linear solvers---------------------
    /** Gaussian Elimination Solver
     * @tparam U
     * @param matA The input square matrix A in the linear system Ax=b
     * @param vecb The vector b in the linear system Ax=b
     * @param output The vector x to be solved in the linear system Ax=b
     */
    template<class U>
    void Gaussian_Elimination(const Matrix<U>& matA, const  U* vecb, U* output);

    /** Gauss Jordan Solver
     * This version uses partial pivoting
     * @tparam U
     * @param matA The input square matrix A in the linear system Ax=b
     * @param vecb The vector b in the linear system Ax=b
     * @param output The vector x to be solved in the linear system Ax=b
     */
    template <class U>
    void Gauss_Jordan(const Matrix<U>& matA, const U* vecb, U* output);

    /** LU Solver
     * @tparam U
     * @param matA The input square matrix A in the linear system Ax=b
     * @param vecb The vector b in the linear system Ax=b
     * @param output The vector x to be solved in the linear system Ax=b
     */
    template <class U>
    void LU_solver( Matrix<U>& matA,  U* vecb, U* output);

    /** LU Solver with partial pivoting
     * This version uses partial pivoting
     * @tparam U
     * @param matA The input square matrix A in the linear system Ax=b
     * @param vecb The vector b in the linear system Ax=b
     * @param output The vector x to be solved in the linear system Ax=b
     */
    template <class U>
    void LU_solver_partial_pivoting(Matrix<U>& matA,  U* vecb, U* output);

    /** Gauss Seidel Solver
     * @tparam U
     * @param matA The input square matrix A in the linear system Ax=b
     * @param vecb The vector b in the linear system Ax=b
     * @param output The vector x to be solved in the linear system Ax=b
     */
    template <class U>
    void gauss_seidel(Matrix<U>& matA, U* vecb, U* output);

    /** Jacobi Method Solver
     * @tparam U
     * @param matA The input square matrix A in the linear system Ax=b
     * @param vecb The vector b in the linear system Ax=b
     * @param it_max Maximum number of iterations
     * @param it_tol Maxium tolerance for residual (Once residual < it_tol, end iteration)
     */
    template <class U>
    void Jacobi_Method(Matrix<U> & matA, const  U* vecb, int it_max, double it_tol, U* output);


    // ---------------------Sparse linear solvers---------------------
    /** Gaussian Elimination Solver for CSRMatrix
     * @tparam U
     * @param matA The input square matrix A in the linear system Ax=b
     * @param vecb The vector b in the linear system Ax=b
     * @param output The vector x to be solved in the linear system Ax=b
     */
    template <class U>
    void Gaussian_Elimination_CSRMatrix(CSRMatrix<U>& matA, const  U* vecb, U* output);

    /** Gauss Jordan Solver for CSRMatrix
     * This version uses partial pivoting
     * @tparam U
     * @param matA The input square matrix A in the linear system Ax=b
     * @param vecb The vector b in the linear system Ax=b
     * @param output The vector x to be solved in the linear system Ax=b
     */
    template <class U>
    void Gauss_Jordan_CSRMatrix(CSRMatrix<U>& matA, const U* vecb, U* output);

    /** Gauss Seidel Solver for CSRMatrix
     * @tparam U
     * @param matA The input square matrix A in the linear system Ax=b
     * @param vecb The vector b in the linear system Ax=b
     * @param output The vector x to be solved in the linear system Ax=b
     */
    template <class U>
    void Gauss_Seidel_CSRMatrix(CSRMatrix<U>& matA, U *vecb, U * output);

    /** Jacobi Method Solver for CSRMatrix
     * @tparam U
     * @param matA The input square matrix A in the linear system Ax=b
     * @param vecb The vector b in the linear system Ax=b
     * @param it_max Maximum number of iterations
     * @param it_tol Maxium tolerance for residual (Once residual < it_tol, end iteration)
     */
    template <class U>
    void Jacobi_Method_CSRMatrix(CSRMatrix<U> & matA, const  U* vecb, int it_max, double it_tol, U* output);

    /** LU Solver for CSRMatrix
    * @tparam U
    * @param matA The input square matrix A in the linear system Ax=b
    * @param vecb The vector b in the linear system Ax=b
    * @param output The vector x to be solved in the linear system Ax=b
    */
    template <class U>
    void LU_solver_CSRMatrix( CSRMatrix<U>& matA,  U* vecb, U* output);

};

//function to check convergence, needed for Gauss-Seidel method
template <class T>
bool check_error(Matrix<T>& mat, Matrix<T>& vect, Matrix<T>& vect_output);

#endif //GROUP_ASSIGNMENT_LSG_SOLVER_H

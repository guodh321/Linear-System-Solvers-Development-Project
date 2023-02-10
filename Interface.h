#include "Matrix.h"
#include "CSRMatrix.h"

#ifndef GROUP_ASSIGNMENT_LSG_INTERFACE_H
#define GROUP_ASSIGNMENT_LSG_INTERFACE_H

/**
 * A interface class that provides all the necessary functions to initialize the user interface
 * Users can select the matrix type, enter matrix data, and use linear solvers in the interface.
 */

class Interface{

public:
    // constructor
    Interface();
    // destructor
    virtual ~Interface();

    // ------display user interface on the terminal------
    void interfaceStart();

    // ---------------Solver test functions---------------
    void test_Gaussian_Elimination();
    void test_Gaussian_Jordan();
    void test_LU_solver();
    void test_LU_partial_pivoting_solver();
    void test_Gaussian_Seidel();
    void test_Jacobi_Method();

    void test_Gaussian_Elimination_CSRMatrix();
    void test_Gaussian_Jordan_CSRMatrix();
    void test_Gaussian_Seidel_CSRMatrix();
    void test_Jacobi_Method_CSRMatrix();
    void test_LU_Solver_CSRMatrix();


private:
    // part 1 is used to select the type of matrix
    void interfacePartOne();

    // part 2 is used to select matrix size
    void interfacePartTwo(char matrixType);

    // part 3 is used to fill data into matrix
    void interfacePartThree(int rows, int cols);

    // part 3 is used to fill data into sparse matrix
    void interfacePartThree_sparse(int rows, int cols);

    // part 4 is used to select linear solvers for dense matrix
    template<class U>
    void interfacePartFour(Matrix<U>& matA, U* b);

    // part 4 is used to select linear solvers for dense matrix
    template<class U>
    void interfacePartFour_sparse(CSRMatrix<U>& matA, U* b);

    // run the testing cases for each solver
    void interfaceMatrixTest();

};

// randomly generate double array
void randomVectb(double *a, int n, double l, double r);

#endif //GROUP_ASSIGNMENT_LSG_INTERFACE_H

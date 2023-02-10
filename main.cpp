//
// Created by Yifan Cheung on 2022/1/20.
//

#include <fstream>
#include "Interface.h"
#include "Solver.h"
#include "Matrix.h"
#include "CSRMatrix.h"
#include "Solver.cpp"
#include "Matrix.cpp"
#include "CSRMatrix.cpp"

using namespace std;

int main(){
    Interface anInterface;
    anInterface.interfaceStart();
//    anInterface.test_Gaussian_Elimination();
//    cout << endl << endl;
//    anInterface.test_Gaussian_Jordan();
//    cout << endl << endl;
//    anInterface.test_LU_solver();
//    cout << endl << endl;
//    anInterface.test_LU_partial_pivoting_solver();
//    cout << endl << endl;
//    anInterface.test_Gaussian_Seidel();
//    cout << endl << endl;
//    anInterface.test_Jacobi_Method();
//    cout << endl;
//
//    cout << endl;
//    anInterface.test_Gaussian_Elimination_CSRMatrix();
//    cout << endl << endl;
    // 有bug，会打印无关的数组
//    anInterface.test_Gaussian_Seidel_CSRMatrix();
//    cout << endl << endl;
//    anInterface.test_Jacobi_Method_CSRMatrix();
}
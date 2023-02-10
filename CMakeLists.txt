cmake_minimum_required(VERSION 3.21)
project(group_assignment_lsg)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(group_assignment_lsg
        main.cpp
        Matrix.cpp
        Matrix.h
        Solver.cpp
        Solver.h
        CSRMatrix.cpp
        CSRMatrix.h
        Interface.cpp
        Interface.h


        )

#add_executable(group_assignment_lsg
#
#        Matrix.cpp
#        Matrix.h
#        Solver.cpp
#        Solver.h
#        CSRMatrix.cpp
#        CSRMatrix.h
##        Interface.cpp
##        Interface.h
#
##        Solver_d.h
##        Solver_d.cpp
#        Test_d.cpp
#        )
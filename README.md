[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-f059dc9a6f8d3a56e377f745f24479a46679e63a5d9fe6f495e02850cd0d8118.svg)](https://classroom.github.com/online_ide?assignment_repo_id=6667355&assignment_repo_type=AssignmentRepo)

# MATRIX  LINEAR  SOLVER  COMMANDLINE  TOOL

The objective of this assignment was to produce an efficient C++ library for solving linear systems of the form Ax = b where A is a matrix and x and b are vectors respectively. We use multiple algorithms for dense and sparse matrices to achieve this objective.

| **Dense solvers implemented**          | **Sparse solver implemented** |
| -------------------------------------- | ----------------------------- |
| Gaussian Elimination                   | Gaussian Elimination          |
| Gaussian Jordan                        | Gaussian Jordan               |
| Gaussian-Seidel                        | Gaussian-Seidel               |
| LU factorisation                       | LU factorisation              |
| LU factorisation with partial pivoting | Jacobi method                 |
| Jacobi method                          |                               |



## Live Demonstration

![LiveDemo](https://media.giphy.com/media/NtEGUePGphJb7YO7BE/giphy-downsized-large.gif)





## Getting Started

- **Prerequisite**

  - Compiler: Clang version 12.0.0 or above

- **Clone the project**

  - ```bash
    cd /path-to-the-cloned-repo
    ```

  - ```bash
    git clone https://github.com/acse-2020/group-assignment-lsg.git
    ```

- **Open the project**

  1. We recommend to open the project with IDE `CLion` since we included the correponding IDE configuration file in the repo which helps you to run the command line interface with one click!
  2. You can also run via the command line through a sequnence of commands below

  ```bash
  g++-11 -o output Matrix.cpp Solver.cpp Interface.cpp main.cpp
  ```

  ```bash
  ./output
  ```





## User Instructions

After running **main.cpp** in the repo, an introduction of this command line tool is shown

```bash
 -------------------------------------------
|  MATRIX  LINEAR  SOLVER  COMMANDLINE  TOOL |
 -------------------------------------------
| Developed by : Team LSG                    |
|       _       _______  _______             |
|      ( \      (  ____ \(  ____ \           |
|      | (      | (    \/| (    \/           |
|      | |      | (_____ | |                 |
|      | |      (_____  )| | ____            |
|      | |            ) || | \_  )           |
|      | (____/\/\____) || (___) |           |
|      (_______/\_______)(_______)           |
 -------------------------------------------

 ------------------------------------------------------
|                     Matrix Data                      |
 ------------------------------------------------------
| Data inside the matrix can be randomly generated or  |
| entered by user through keyboard and text file.  And |
| matrix data range can be decided by user.            |
 ------------------------------------------------------

 ------------------------------------------------------
|                     Matrix Size                      |
 ------------------------------------------------------
| Matrix size are set to be larger than 10 X 10        |
| entered by user through keyboard and text file.      |
| Two types of matrices are provided                   |
 ------------------------------------------------------

 ---------------------------------- 
|          Matrix Type             |
 ---------------------------------- 
| 1. Dense matrix                  |
| 2. Sparse matrix                 |
 ----------------------------------

 ----------------------------------------------------- 
| Do you want to continue and select the matrix (y/n) |
 ----------------------------------------------------- 
| Enter y to continue                                 |
| Enter n to exit                                     |
 ---------------------------------------------------- 
>> y
```

User enter `y` to continue and select the matrix storage type. Users can select dense or sparse matrix. Predefined testing matrix cases for each solver are also provided.

```bash
 ----------------------------------------------------------
|               Select matrix storage type                  |
 ----------------------------------------------------------
| 1: Dense Matrix                                           |
| 2: Sparse Matrix                                          |
| 3: Test Matrix (Run testing on each solver automatically) |
| b: Back                                                   |
| x: Exit                                                   |
 ----------------------------------------------------------
>> 1
```

Next, user can define the matrix size by themselves or use the predefined matrix size. Here we enter the size we want.

```bash
 -----------------------------------------------
|           Set up the matrix size              |:
 -----------------------------------------------
| 1: 10 x 10                                    |
| 2: 15 x 15                                    |
| 3: Enter the size (greater than 10 X 10)      |
| b: Back                                       |
| x: Exit                                       |
 -----------------------------------------------
>> 3
 ------------------------------------------------------------------ 
| Enter the row number or column number since it is a square matrix |
 ------------------------------------------------------------------ 
>> 3
```

There are three ways to fill up the matrix and vector b. The first option is random generation of matrix data, we can define the minimum and maximum values present within the randomly generated matrix.

```bash
 -----------------------------------------------
|        Fill up the dense matrix data          |:
 -----------------------------------------------
| 1: Randomly generated data                    |
| 2: Enter the data                             |
| 3: Text file                                  |
| b: Back                                       |
| x: Exit                                       |
 -----------------------------------------------
>> 1

 ------------------------------------------------------------------------- 
| Generate randomly generated data between left_boundary and right_boundary |
 ------------------------------------------------------------------------- 
  Please enter the left_boundary first
>> 10
  Then please enter the right_boundary first
>> 20

------------------------------------------------
Print randomly generated matrix A
11.3784 17.1028 16.0673 
13.0392 19.2356 12.8096 
10.5268 14.0893 18.0547 
------------------------------------------------

------------------------------------------------
Print randomly generated vector b:
[11.3784,17.1028,16.0673]
------------------------------------------------
```

Finally user can select the solvers to solve the linear system.

```bash
 ----------------------------------------------
|  Select the dense matrix solver               |:
 ----------------------------------------------
| 1: Gaussian Elimination                       |
| 2: Gaussian Jordan                            |
| 3: Gaussian Seidel                            |
| 4: LU factorisation                           |
| 5: LU factorisation with partial pivoting     |
| 6: Jacobi                                     |
| b: Back                                       |
| x: Exit                                       |
 -----------------------------------------------
>> 1

------------------------------------------------
The output solved by Gaussian Elimination:
[ 7.73753,-4.04774,-0.46274]
Time spent: 0.006 milliseconds
------------------------------------------------
```



## Documentation

A well-written 4-page project report that gives a brief description of the project architecture, code structure, input/output, and how to run it can be checked herer [report link](https://github.com/acse-2020/group-assignment-lsg/blob/main/report.pdf)



## Testing

1. **Testing on each solvers**

Each solver testing function is written in **Interface::interfaceMatrixTest( )**. Users can select the **Test matrix** in the interface commnd line tool.

2. **Timing test**

The run time for each solver is timed every time we execute them to solve linear systems. The timing test is integrated in the **Interface::interfaceMatrixTest( )** as well.



## Authors

- **Team LSG**
  - [Yifan Zhang](https://github.com/edsml-yz1921)
  - [Donghu Guo](https://github.com/edsml-dg321)
  - [Sanjeet Chhokar](https://github.com/SC7818-EDSML)



## License

The scripts and documentation in this project are released under the [MIT License](https://github.com/actions/upload-artifact/blob/main/LICENSE).
# A Finite Element Project

## The blueprint
---

<br/>

The main goal of the project is to implement a solver for some simple **Partial Differential Problems in lower dimensions (up to two dimensions)**, specifically **stationary elliptic PDEs involving multivariate scalar-valued functions using the Spectral Element method**, depicted in CHQZ2[^1]. 

[^1]: C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,"Spectral Methods. Fundamentals in Single Domains"Springer Verlag, Berlin Heidelberg New York, 2006.

The rectangular meshes used to analyse the problem are created by modifing the MATLAB script written by prof. Paola Gervasio[^2] and suitably converted into `csv` files form which the code extrapolates the necessary information to compute the approximated solution of the problem considered.

[^2]: [Paola Gervasio](gervasio@ing.unibs.it), Department of Mathematics, University of Brescia, 25133 Brescia (Italy). 

The solution is then written onto a `.vtp` file and can be handled with visualization softwares such as ParaView.

It's worth mentioning the premises of the Spectral Element method, i.e. the compuation of the solution on the mesh is carried out in a **matrix free** way: while considering the diffent elements and computing their specific contribution on the structure of the gloabl system matrix, we avoid explicitly building the local matrix but, instead, we exploit the local stiffness and mass matrix of the reference element $[-1,1]\times[-1,1]$ and then apply the transformation to the specific element of the mesh considered (such trasformation is performed as the algebraic product between matrixes).

<br/>

## Requirements
---
<br/>

Since the code heavily relies on the `Eigen` library, to deal with algrabraic structures, as  well as the `VTK` library, to export the solution to a file (for later visualization and processing), to properly run the program, these libraries must be in scope when compiling using the incorporated CMake file.

<br/>

## Assumptions
---
<br/>

Some assuptions and generality restrictions  where taken along the way, to simplify the coding process:

- First and foremost, the code has been written to work, without any changes, to solve a fixed two-dimensional problem. The solver class, in particular, is endowed with a general structure al long as the different in-class methods are concered, still its lines of code have been written to some specific problems: **this means that changing the problem considered would also require to change the solver class to some (shallow) extent**

- The project was carried out with the intention of preserving the generality in dimensions, degree and number of quadrature nodes per cell considered - this is why most of the classes and objects defined are strongly dependent on some global constant parameters.
Yet along the development of the code, some assumpions where taken:
   - The degree $r$ of the finite element space used is supposed to be the same for both $x$ an $y$ directions: it is possible to change it in the first lines of the `main.cpp` file, inside the `src` directory.

    - The same holds true for the number of quadrature points along the two directions, set to be equal to $r$ (A total of $r+1$ points along each side of each cell).

<br/>

## How to run the code
---

<br/>

The code comes with a CMake file to make the compilation process easier: to run the code, follow these steps form terminal inside the directory of the project

1. `mkdir build`

2. `cd build`

3. `cmake ..`

4. `make`

5. Finally launch the project `./A_Finite_Element_project`

The solution is then stored with the mesh inside a `solution.vtp` file inside the `build` directory and can be easily visualized using ParaView.

>> insert ooption selection from terminal

<br/>

## Test case
---

<br/>


The test case considered is the following:
Let $\Omega = (0, 1) × (0, 1)$, and let us consider the following Poisson problem with homogeneous Dirichlet boundary conditions:

$$ 
\begin{cases}
    \nabla(\mu \nabla u) + \sigma u = f && x\space\epsilon\space\Omega \\ 
     u = 0 && on\space\delta\Omega
\end{cases}
$$

where $x = (x, y)^T$ , $\mu(x) = 1$, $\sigma = 1$ and
$$f(x)=(20π2 +1)sin(2πx)sin(4πy)$$
The exact solution to this problem is
$$u_{ex}(x, y) = sin(2\pi x) sin(4\pi y)$$

>> insert more test cases
<br/>


## The parallel code
---
The parallelisation was carried out using OpenMP:
two main approaches where considered when parallelising the code, both of them centered around the assembly phase of global system matrix.

- The first approach is based on the use of **OpenMP locks**: when the single threads compute a subset of the local matrixes, it may appen that two or more of them were to write on the same entry of the matrix. To avoid data races, locks where used.
Idealy it would be better to have a single lock for each coefficinet of the matrix, but doing so would result in terrible performances due to the use of lock's overheads. To mitigate such phenomena, the idea would be to define a resticted number of locks: specifically, we define a set of lock to partition the global matrix into sub-matrixes. Each time a thread accesses an element of one of these blocks, it acquires the lock for the specifc sub-matrix and then returns it when it has finished updating the coefficient.
To obtain a better performance we also use a **_try lock_ mechanism**.

- The second approach would require using some **local support** varibles, private to each thread: each one of them creates a vector of local matrixes and updates them. Once all the threads have finished, we reduce the matrixes on the global one.

## References
---
<br/>

* _A. Quarteroni, [Numerical Models for Differential Problems, MS&A 16](https://doi.org/10.1007/978-3-319-49316-9)_

* _Paola Gervasio, [Spectral Element Library - CHQZ lib](https://paola-gervasio.unibs.it/software/)_


## Thanks
---

<br/>

Special thanks to PHD student [Matteo Caldana](matteo.caldana@polimi.it) and prof. [Luca Formaggia](luca.formaggia@polimi.it), MOX Politecnico di Milano, Department of Mathematics, Via Bonardi n9, 20133 Milano

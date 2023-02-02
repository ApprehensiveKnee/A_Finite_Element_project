# A Finite Element Project

>## The blueprint

<br/>

The main goal of the project is to implement a solver for some simple **Partial Differential Problems in lower dimensions**, specifically **stationary elliptic PDEs using the Spectral Element method** depicted in CHQZ2[^1]. 

[^1]: C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,"Spectral Methods. Fundamentals in Single Domains"Springer Verlag, Berlin Heidelberg New York, 2006.

The rectangular meshes used to analyse the problem are created by modifing the MATLAB script written by prof. Paola Gervasio[^2] and suitably converted into `csv` files form which the code extrapolates the necessary information to compute the approximated solution of the problem considered.

[^2]: [Paola Gervasio](gervasio@ing.unibs.it), Department of Mathematics, University of Brescia, 25133 Brescia (Italy). 

The solution is then written onto a `.vtp` file and can be handled with visualization sofwares such as ParaView.

It's worth mentioning the premises of the Spectral Element method, i.e. the fact that the compuation of the solution on the mesh is carried out in a **matrix free** way: while considering the diffent elements and computing their specific contribution on the structure of the gloabl system matrix, we avoid explicitly building the local matrix but, instead, we exploit the local stiffness and mass matrix of the reference element $[-1,1]\times[-1,1]$ and then apply the transformation to the specific element of the mesh considered (such trasformation is performed as the algebraic product between matrixes).

<br/>

>## Requirements
<br/>

Since the code heavily relies on the `Eigen` library, to deal with algrabraic structures, as  well as the `VTK` library, to export the solution to a file (for later visualization and processing), to properly run the program, these libraries must be in scope when compiling using the incorporated CMake file.

<br/>

>## Assumptions
<br/>

Some assuptions and generality restrictions  where taken along the way, to simplify the coding process:

- First and foremost, the code has been written to work without any changes to solve a fixed two-dimensional problem. The solver class, in particular, is endowed with a general structure al long as the different in-class methods are concered, still its lines of code have been written to solve a specific problem: **this means that changing the problem considered would also require to change the solver class to some (shallow) extent**

-  The degree $r$ of the finite element space used is supposed to be the same for both $x$ an $y$ directions: it is possible to change it in the first lines of the `main.cpp` file, inside the `src` directory.

- The same holds true for the number of quadrature points along the two directions, set to be equal to $r$ (A total of $r+1$ points along each side of each cell).

<br/>

>## How to run the code

<br/>

The code comes with a CMake file to make the compilation process easier: to run the code, follow these steps form terminal inside the directory of the project

1. `mkdir build`

2. `cd build`

3. `cmake ..`

4. `make`

5. Finally launch the project `./A_Finite_Element_project`

The solution is then stored with the mesh inside a `solution.vtp` file inside the `build` directory and can be easily visualized using ParaView.

<br/>

>## Test case

<br/>
The test case considered is the following:

Let $\Omega = (0, 1) × (0, 1)$, and let us consider the following Poisson problem with homogeneous Dirichlet boundary conditions:

$$ 
\begin{cases}
    \nabla(\mu \nabla u) + \sigma u && x\space\epsilon\space\Omega \\ 
    \left\lvert\ u = 0 && on\space\delta\Omega
\end{cases}
$$

where $x = (x, y)^T$ , $\mu(x) = 1$, $\sigma = 1$ and
$$f(x)=(20π2 +1)sin(2πx)sin(4πy)$$
The exact solution to this problem is
$$u_{ex}(x, y) = sin(2\pi x) sin(4\pi y)$$

<br/>

>## References
<br/>

* _A. Quarteroni, [Numerical Models for Differential Problems, MS&A 16](https://doi.org/10.1007/978-3-319-49316-9)_

* _Paola Gervasio, [Spectral Element Library - CHQZ lib](https://paola-gervasio.unibs.it/software/)_



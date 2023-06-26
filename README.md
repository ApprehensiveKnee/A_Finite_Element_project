# A Finite Element Project

## The blueprint
---

<br/>

The main goal of the project is to implement a solver for some simple <span style="color:darkorange">**Partial Differential Problems in lower dimensions (1 and 2D)**</span>, specifically <span style="color:darkorange">**stationary elliptic PDEs involving multivariate scalar-valued functions using the Spectral Element method**</span>, depicted in CHQZ2[^1]. 

[^1]: C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,[_"Spectral Methods. Fundamentals in Single Domains"_](http://dx.doi.org/10.1007/978-3-540-30726-6) Springer Verlag, Berlin Heidelberg New York, 2006.

The rectangular meshes (linear in the 1D case) used to analyse the problem are created by modifing the MATLAB script written by prof. Paola Gervasio[^2] and suitably converted into `csv` files form which the code extrapolates the necessary information to compute the approximated solution of the problem considered.

[^2]: [Paola Gervasio](gervasio@ing.unibs.it), Department of Mathematics, University of Brescia, 25133 Brescia (Italy). 

The solution is then written onto a `.vtp` file and can be handled with visualization softwares such as ParaView.

It's worth mentioning the premises of the Spectral Element method, i.e. the compuation of the solution on the mesh is carried out in a **matrix free** way: while considering the diffent elements and computing their specific contribution on the structure of the gloabl system matrix, we avoid explicitly building the local matrix but, instead, we exploit the local stiffness and mass matrix of the reference element $[-1,1]\times[-1,1]$ and then apply the transformation to the specific element of the mesh considered (such trasformation is performed as the algebraic product between matrixes).

<br/>

## Requirements
---
<br/>

Since the **code heavily relies on the** `Eigen` **library**, to deal with algrabraic structures, as  well **as the** `VTK` **library**, to export the solution to a file (for later visualization and processing), to properly run the program, these libraries must be in scope when compiling using the incorporated CMake file.
Finally, `Boost` **libraries are required** for option parsing to properly start the execution of the code .

<br/>

## Assumptions
---
<br/>

Some assuptions and generality restrictions  where taken along the way, to simplify the coding process:

1.  The project was carried out with the intention of preserving the generality in dimensions, degree and number of quadrature nodes per cell considered - this is why most of the classes and objects defined are strongly dependent on some global constant parameters.
Yet along the development of the code, some assumpions where taken:
    - Even if different classes are endowed with a member (typically an array) to store the degree of the finite element space used along each considered dimension indipendenlty, the <span style="color:darkorange">degree **r** </span> of the finite element space used is supposed to be the same for both $x$ an $y$ directions (for the 2D case): it is possible to change it in the `variables.hpp` file, along with the <span style="color:darkorange">dimension of the problem **DIM**</span>.

    - The same holds true for the number of quadrature points along the two directions, set to be equal to $r$ (A total of $r+1$ points along each side of each cell).

2. As for the imposition of boundary conditions, in all the test cases considered, Dirichelet conditions are inposed at the bounday of the considered domain. This is done explicitly inside the methods of the solver class, which means that imposing different conditions at the boundary would require to change the src code direcly - the changes would not be extensive, and only limited to the source files containing the impelementation of the solver, i.e. solver.cpp, psolver.cpp and the likes.

<br/>

## How to run the code
---

<br/>

The code comes with a CMake file to make the compilation process easier: to run the code, follow these steps form terminal inside the directory of the project

```
mkdir build
cd build
cmake ..
make
```

Finally launch the project by typing `./A_Finite_Element_project`

To execute the code it is necessary to specify some options. To visualize the list of options available, type when launching the project:

```
-- help
```

After launching the code with the desired options, the code will compute the solution of the problem considered on the mesh specified in the `mesh.csv` file inside the `build` directory.


<br/>

## Test case
---

<br/>


The test cases considered are the following:

### **1D Test case**

Let $\Omega = (0, 1)$, and let us consider the following Poisson problem:

$$ 
\begin{cases}
    - (\mu (x) u' (x))' = f && x\space\epsilon\space\Omega \\ 
     u(0) = u(1) = 0
\end{cases}
$$

where $\ mu(x) = 1$ and $f(x)$ is defined as:

$$ 
f(x) = 
\begin{cases}
    0  && \textrm{ if } x \leq \frac{1}{2} \\ 
    - \sqrt{x - 1/2} && \textrm{ if } x > \frac{1}{2}
\end{cases}
$$

The exact solution is

$$ 
u_ex(x) = 
\begin{cases}
    A && \textrm{ if } x \leq \frac{1}{2} \\ 
    Ax + \frac{4}{15}(x - \frac{1}{2})^\frac{5}{2} && \textrm{ if } x > \frac{1}{2} 
\end{cases}
$$

where $A = - \frac{4}{15}(\frac{1}{2})^\frac{5}{2}$


### **2D Test Case**

Let $\Omega = (0, 1) × (0, 1)$, and let us consider the following Poisson problem with homogeneous Dirichlet boundary conditions:

$$ 
\begin{cases}
    \nabla(\mu \nabla u) + \sigma u = f && x\space\epsilon\space\Omega \\ 
     u = 0 && \textrm{on}\space\delta\Omega
\end{cases}
$$

where $x = (x, y)^T$ , $\mu(x) = 1$, $\sigma = 1$ and
$$f(x)=(20π2 +1)sin(2πx)sin(4πy)$$
The exact solution to this problem is
$$u_{ex}(x, y) = sin(2\pi x) sin(4\pi y)$$
<br/>


## The parallel code
---
The parallelisation was carried out using OpenMP:
<span style="color:darkorange">_**three main approaches**_</span> where considered when parallelising the code, all of them centered around the assembly phase of global system matrix.

<span style="color:orange">**PLEASE NOTE**</span>: For the parallel version, the solving process was carried out in a parallel way too, i.e. the parallel version of the corresponding Eigen iterative solving method was used: no actual programming effort went into it.


1.  The <span style="color:darkorange">**first approach**</span> would require using some _**local support**_ varibles, specifically <span style="color:darkorange">_**unordered maps**_</span>, private to each thread: each of them creates an unordered map storing the non-zero elements. To do so, the keys of such maps are pairs of unsigned integers, representing the indexes of the global matrix. After the computation of the local matrixes in parallel is completed, each thread "loads its contribution" onto the global matrix serially form the private unordered maps.\
This last passage is necessarly carried out in serial, as to avoid **race conditions** on the elements of the global matrix.


2. The <span style="color:darkorange">**second approach**</span> is based on the use of <span style="color:darkorange">_**OpenMP locks**_</span>: when the single thread computes a subset of the local matrixes, it may happen that two or more of them are to write on the same entry of the matrix. To avoid data races, locks where used.
Idealy it would be better to have a single lock for each coefficinet of the matrix, but doing so would result in terrible performances due to the use of lock's overheads.\
To mitigate such phenomena, the idea is to define a **resticted number of locks**: specifically, we define a set of lock to partition the global matrix into sub-matrixes (partition of the rows of the global matrix).
Each time a thread gains access to an element of one of these blocks, it acquires the lock for the specifc sub-matrix and then returns it when it has finished updating the coefficient.

3. Finally the <span style="color:darkorange">**third approac**h</span> revolves around the <span style="color:darkorange">_**coloring algorithm**_</span> based on graph theory: the elements assigned to some color, in such a way that elements having some sharing points do not have the same color. Having divided the elements in such a way, the ones belonging to the same group can be processed in parallel, as their contribution to the global matrix do not map on the same coefficients, which means that **there can be no possibility of race conditions occurring**. \
It is important to notice that this approach implies the necessity to **group the elements of the mesh per color**. This is done in the setting process of the solver, and has a considerable impact on the performances of the code:
>_Even though this is appealing, the algorithm has its overhead. This includes the already discussed need to establish the colouring. However, only the inner loop can be parallelized. All threads should finish processing the element’s inner loop be-
fore processing elements for the next colour. This requires synchronization. Finally, there is also an overhead connected to the creation and termination of threads for each colour. In the colouring algorithm, the individual threads are created and terminated inside the outer loop over individual colours, but this overhead is typically very small, considering typical assembly times for real problems. Eventually, it can have some minor impact on overhead costs._[^3]

[^3]: _M. Bošanský, B. Patzák, ["Parallelization of assembly operations in finite element Method"](https://doi.org/10.14311/AP.2020.60.0025)_

It is worth noticing that almost all of the parallel code has been implemented with a specific intent to **reduce to the minimum the number of memory accesses**, in an attempt to make up for the overheads due to the managing of threads.
For such reason, it isn't uncommon to find the parallel version of the code filled with more control pa than the  serial counterpart

## References
---
<br/>

* _A. Quarteroni, [Numerical Models for Differential Problems, MS&A 16](https://doi.org/10.1007/978-3-319-49316-9)_

* _P. Gervasio, [Spectral Element Library - CHQZ lib](https://paola-gervasio.unibs.it/software/)_

* _M. Bošanský, B. Patzák, [Parallelization of assembly operations in finite element Method](https://doi.org/10.14311/AP.2020.60.0025)_


## Thanks
---

<br/>

Special thanks to PHD student [Matteo Caldana](matteo.caldana@polimi.it) and prof. [Luca Formaggia](luca.formaggia@polimi.it), MOX Politecnico di Milano, Department of Mathematics, Via Bonardi n9, 20133 Milano

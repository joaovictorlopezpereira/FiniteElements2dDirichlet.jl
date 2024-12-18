# Problem Formulation

This package solves the following 2D ordinary differential equation using the finite elements method:

## Strong Formulation

Let:
- ``\Omega = [0,1] \times [0,1] \subset \mathbb{R}^2``;
- ``\Gamma`` be the boundary of ``\Omega``;
- ``\hat{\Omega} = \Gamma \cup \Omega``;
- ``\Delta u(x,y) = u_{xx}(x,y) + u_{yy}(x,y)``.

Given ``\alpha > 0``, ``\beta \geq 0``, and a function ``f : \Omega \to \mathbb{R}``, we want to find ``u : \hat{\Omega} \to \mathbb{R}`` such that:

$\begin{cases}
  -\alpha \Delta u(x,y) + \beta u(x,y) = f(x,y) & \forall (x,y) \in \Omega \\
  \\
  u(x,y) = 0 & \forall (x,y) \in \Gamma.
\end{cases}$


## Weak Formulation

Let:
- ``V`` represent the space of test functions;
- ``H`` represent the space of solution functions that satisfy the given variational formulation;
- ``v \in V`` be a test function;
- ``u \in H`` be a solution function;
- ``\displaystyle\kappa(f,g) = \alpha \int_\Omega \Delta f(x) \Delta g(x) d\Omega + \beta \int_\Omega f(x)g(x) d\Omega``;
- ``\displaystyle(f,g) = \int_\Omega f(x) g(x) d\Omega``.

Given ``\alpha > 0``, ``\beta \geq 0``, and a function ``f : \Omega \to \mathbb{R}``, we want to find `` u \in H, u : \hat{\Omega} \to \mathbb{R}`` such that:

$\begin{cases}
  \kappa(u,v) = (f,v) & \forall (x,y) \in \Omega \\
  \\
  u(x,y) = 0 & \forall (x,y) \in \Gamma.
\end{cases}$


## Approximate Problem

Let:
- ``V_h`` represent the finite-dimensional space of approximate test functions;
- ``H_h`` represent the finite-dimensional space of approximate solution functions that satisfy the discrete variational formulation;
- ``v_h \in V_h`` be an approximate test function;
- ``u_h \in H_h`` be an approximate solution function.

Given ``\alpha > 0``, ``\beta \geq 0``, and a function ``f : \Omega \to \mathbb{R}``, we want to find `` u_h \in H_h, u_h : \hat{\Omega} \to \mathbb{R}`` such that:

$\begin{cases}
  \kappa(u_h,v_h) = (f,v_h) & \forall (x,y) \in \Omega \\
  \\
  u_h(x,y) = 0 & \forall (x,y) \in \Gamma.
\end{cases}$


## Matrix Formulation

Let
- ``\varphi_i`` be the basis functions of the finite-dimensional space ``H_h = V_h``.

Given a matrix ``K`` and a vector ``F``, we want to find ``C`` such that ``KC = F``

```math
\begin{align*} \begin{bmatrix}
\kappa\big(\varphi_1,\varphi_1\big)& \dots& \kappa\big(\varphi_m,\varphi_1\big) \\
\vdots& \ddots& \vdots \\
\kappa\big(\varphi_1,\varphi_m\big)& \dots& \kappa\big(\varphi_m,\varphi_m\big) \\
\end{bmatrix} \begin{bmatrix}
c_1\\ \vdots\\ c_m
\end{bmatrix} = \begin{bmatrix}
(f,\varphi_1)\\ \vdots\\ (f,\varphi_m)\\
\end{bmatrix}. \end{align*}
```

More details on the calculations can be found at [Two Dimensional Analysis and Calculations](https://github.com/joaovictorlopezpereira/Finite-Elements-Method/blob/main/Analysis%20and%20Calculations/two-dim.pdf).


---


# Function Signatures

Functions signatures from the FiniteElements2dDirichlet package.

## EQ, LG and m Initializers

```@docs
init_LG_matrix(Nx, Ny)
```
```@docs
init_EQ_vector_and_m(Nx, Ny)
```


## Fe and F Initializers

```@docs
init_Fe_vector!(f, Xs, Ys, Fe, P, W)
```

```@docs
init_F_vector(f, X_matrix, Y_matrix, m, EQ, LG)
```


## Ke and K Initializers

```@docs
init_Ke_matrix!(alpha, beta, Xs, Ys, Ke, P, W)
```
```@docs
init_K_matrix(alpha, beta, X_matrix, Y_matrix, m, EQ, LG)
```


## Mesh Initializer

```@docs
init_mesh(Nx, Ny; ns=false, plot=false)
```


## System Solver

```@docs
solve_system(alpha, beta, f, Nx, Ny; EQLG=false, XY_matrix=false, noise=false)
```


## Error Convergence

```@docs
error_convergence(lb, ub, alpha, beta, u, f; see_plot=false, ns=false)
```


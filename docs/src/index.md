# Functions Signatures from FiniteElements2dDirichlet


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


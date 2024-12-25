# FiniteElements2dDirichlet

## Information About the Package

This repository contains a solver for a 2D partial differential equation using the finite element method.

In this case, we are solving a two-dimensional stationary equation with homogeneous Dirichlet boundary conditions.

-----

## Equation Being Solved by the Package

Let:
  - $\Omega = [0,1] \times [0,1] \subset \mathbb{R}^2$;
  - $\Gamma$ be the border of $\Omega$;
  - $\hat{\Omega} = \Gamma \cup \Omega$;
  - $\Delta u(x,y) = u_{xx}(x,y) + u_{yy}(x,y)$.

Given $\alpha > 0$, $\beta \geq 0$ and a function $f : \Omega \to \mathbb{R}$, this system finds $u : \hat{\Omega} \to \mathbb{R}$ such that:

$$
\begin{cases}
  -\alpha \Delta u(x,y) + \beta u(x,y) = f(x,y) & \forall (x,y) \in \Omega \\
  \\ \\
  u(x,y) = 0 & \forall (x,y) \in \Gamma.
\end{cases}
$$

More information about the package and the mathematical formulation can be found at the ``index.html`` file in ``/docs/build``.

-----

## The Future of the Package

At this moment, the package solves only a single equation. In the future, it is planned to extend its capabilities to solve multiple equations, allowing the same functions to be used across different equations.

[![Build Status](https://github.com/joaovictorlopezpereira/FiniteElements2dDirichlet.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/joaovictorlopezpereira/FiniteElements2dDirichlet.jl/actions/workflows/CI.yml?query=branch%3Amaster)

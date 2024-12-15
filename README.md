# FiniteElements2dDirichlet

This repository contains a solver of a 2D ordinary differential equation by using the finite elements method.

Let $\Omega \subset \mathbb{R}^2$, $\Gamma$ be the border of $\Omega$, $\hat{\Omega} = \Gamma \union \Omega$, and $\Delta u(x,y) = u_{xx}(x) + u_{yy}(x)$.

Given $\alpha > 0$, $\beta \geq 0$ and a function $f : \Omega \to \mathbb{R}$, this system finds $u : \hat{\Omega} \to \mathbb{R}$ such that:

$$
\begin{cases}
  -\alpha \Delta u(x, y) + \beta u(x, y) = f(x, y) & \forall x \in \Omega \\
  \\
  u(x) = 0 & \forall x \in \Gamma.
\end{cases}
$$

Some of the calculations used in this implementation can be found at [Two Dimensional Analysis and Calculations](https://github.com/joaovictorlopezpereira/Finite-Elements-Method/blob/main/Analysis%20and%20Calculations/two-dim.pdf).


[![Build Status](https://github.com/joaovictorlopezpereira/FiniteElements2dDirichlet.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/joaovictorlopezpereira/FiniteElements2dDirichlet.jl/actions/workflows/CI.yml?query=branch%3Amaster)

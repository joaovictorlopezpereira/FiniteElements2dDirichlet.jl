
# X function-vector
x = [
    (xi1, xi2, Xs) -> sum(Xs[k] * phi[k](xi1, xi2) for k in 1:4),
    (xi1, xi2, Ys) -> sum(Ys[k] * phi[k](xi1, xi2) for k in 1:4)
]

# Derivative of the x function-vector
d_xi_x = [
    [(xi1, xi2, coords) -> sum(coords[k] * d_xi_phi[1][k](xi1, xi2) for k in 1:4),
    (xi1, xi2, coords) -> sum(coords[k] * d_xi_phi[2][k](xi1, xi2) for k in 1:4)],
    [(xi1, xi2, coords) -> sum(coords[k] * d_xi_phi[1][k](xi1, xi2) for k in 1:4),
    (xi1, xi2, coords) -> sum(coords[k] * d_xi_phi[2][k](xi1, xi2) for k in 1:4)]
]
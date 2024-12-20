
# Phi function-vector
phi = [
    (xi1, xi2) -> (1 - xi1) * (1 - xi2) * (1 / 4);
    (xi1, xi2) -> (1 + xi1) * (1 - xi2) * (1 / 4);
    (xi1, xi2) -> (1 + xi1) * (1 + xi2) * (1 / 4);
    (xi1, xi2) -> (1 - xi1) * (1 + xi2) * (1 / 4)
]

# Derivative of the phi function-vector
d_xi_phi = [
    [((xi1, xi2) -> (-1 / 4) * (1 - xi2)),
    ((xi1, xi2) -> ( 1 / 4) * (1 - xi2)),
    ((xi1, xi2) -> ( 1 / 4) * (1 + xi2)),
    ((xi1, xi2) -> (-1 / 4) * (1 + xi2))],
    [((xi1, xi2) -> (-1 / 4) * (1 - xi1)),
    ((xi1, xi2) -> (-1 / 4) * (1 + xi1)),
    ((xi1, xi2) -> ( 1 / 4) * (1 + xi1)),
    ((xi1, xi2) -> ( 1 / 4) * (1 - xi1))]
]
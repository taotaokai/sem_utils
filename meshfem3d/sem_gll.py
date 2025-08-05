#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import warnings

import numpy as np

# n-point Gauss–Lobatto quadrature on interval [-1, 1]
# nodes {x_i, i = 1,...,n}: roots of (1-x^2) * P'_{n-1}(x)
# weights {w_i, i = 1,...,n}: 2 / (n*(n-1)) / (P_{n-1}(x_i))^2


def gll_nodes_weights(n):
    # Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    x = np.cos(np.pi * np.arange(n) / (n-1))
    # n1 = n - 1
    # The Legendre Vandermonde Matrix
    P = np.zeros((n, n))
    # Compute P_(N) using the recursion relation
    # Compute its first and second derivatives and
    # update x using the Newton-Raphson method.
    xold = 2 * np.ones_like(x)
    while max(abs(x - xold)) > 1e-5:
        xold = x
        P[:, 0] = 1
        P[:, 1] = x
        for k in range(1, n-1):
            P[:, k + 1] = ((2 * k - 1) * x * P[:, k] - (k - 1) * P[:, k - 1]) / k

        x = xold - (x * P[:, -1] - P[:, -2]) / (n * P[:, -1])
        print(x)

    w = 2 / (n * (n -1) * P[:, -1] ** 2)

    return x, w


def dLagIdzJ(nodes):
    """
    derivative of Lagrange polynomials evaluated at n-point GLL nodes
    """
    z = nodes
    n = len(nodes)
    # a[i,j] = dLag_i/dz(z_j)
    dLag = np.zeros((n, n))
    ind = np.arange(n)
    for i in range(n):
        denumerator = np.prod((z[i] - z[ind != i]))
        for j in range(n):
            numerator = 0
            for k in range(n):
                if k == i:
                    continue
                mask = (ind != i) & (ind != k)
                numerator += np.prod(z[j] - z[mask])
            dLag[i, j] = numerator / denumerator
    return dLag

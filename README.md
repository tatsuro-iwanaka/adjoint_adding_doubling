# Adjoint Adding-Doubling Radiative Transfer Model
This repository contains a C++ implementation of the Adjoint Adding-Doubling Method for plane-parallel radiative transfer problems.
Designed for inverse problems in planetary atmospheric science, this code computes the exact analytical gradients of the cost function with respect to atmospheric parameters.
It utilizes the Eigen library for matrix operations.
Supports both scattering and thermal emission.

# Forward Adding-Doubling Method Implementation

## Overview
This repository implements the **Adding-Doubling method** for radiative transfer calculations in a plane-parallel atmosphere. The core algorithm computes the reflection and transmission matrices of a combined layer by interacting two sub-layers.

The implementation explicitly handles the **direct beam attenuation** (matrix $E$) and accounts for the **azimuthal dependence** of the phase function.

## Numerical Integration (Gauss-Radau Quadrature)

The code employs **Gauss-Radau quadrature** for discretizing the zenith angle $\mu$. Unlike the standard Gauss-Legendre quadrature, the Gauss-Radau method fixes one integration node at the boundary of the interval.

### Mathematical Formulation
The integral of a function $f(\mu)$ over $[0, 1]$ is approximated as:

$$
\int_0^1 f(\mu) d\mu \approx \sum_{i=1}^{n} w_i f(\mu_i)
$$

where the last node is fixed at $\mu_n = 1.0$. This is particularly advantageous for radiative transfer problems where the radiation field at the vertical boundary ($\mu=1$) is of specific interest and should be calculated explicitly without extrapolation.

### Algorithm
The nodes $x_i$ (in the standard interval $[-1, 1]$) and weights $w_i^{gr}$ are computed using the **Golub-Welsch algorithm** adapted for Radau quadrature:

1.  **Jacobi Matrix Construction**: A symmetric tridiagonal matrix $J$ is constructed using the recurrence relation coefficients of the Jacobi polynomials $P_{n-1}^{(1, 0)}(x)$ (fixing the endpoint at $+1$).
2.  **Eigenvalue Decomposition**: The eigenvalues of $J$ correspond to the internal nodes $x_i$.
3.  **Mapping**: The nodes and weights are linearly mapped from $[-1, 1]$ to the physical interval $[0, 1]$:　$\mu_i = \frac{x_i + 1}{2}, \quad w_i = \frac{w_i^{gr}}{2}$

This method integrates polynomials of degree up to $2n-2$ exactly.

## Matrix Formulation of Doubling Equations

The doubling process combines two identical layers (top and bottom) with optical thickness $\tau$ to form a layer of thickness $2\tau$.

### 1. Direct Attenuation Matrix ($E$)
The exponential attenuation of the direct beam for each quadrature angle $\mu_i$ is represented by a diagonal matrix $E$:

$$
E_{ij} = \delta_{ij} \exp\left(-\frac{\tau}{\mu_i}\right)
$$

### 2. Multiple Scattering Interaction
The interaction between the two layers involves infinite reflections. This is solved using a Neumann series (geometric series of matrices):

$$
\begin{aligned}
Q_1 &= c R_{\text{bot}} W R_{\text{top}} \\
Q_2 &= c Q_1 W \\
S &= Q_1 + Q_2 Q_1 + Q_2 Q_2 Q_1 + \dots = (I - Q_2)^{-1} Q_1
\end{aligned}
$$

* $W$: Diagonal matrix of quadrature weights.
* $c$: Azimuthal weighting factor ($c=2$ for $m=0$, $c=1$ for $m>0$).
* $S$: The source term representing the summation of multiple scattering between layers.

### 3. Internal Radiation Fields
The downward ($D$) and upward ($U$) diffuse intensities at the interface between the two layers are calculated as:

$$
\begin{aligned}
D &= c S W T_{\text{top}} + T_{\text{top}} + S E \\
U &= c R_{\text{top}} W D + R_{\text{top}} E
\end{aligned}
$$

Here, terms like $SE$ and $R_{\text{top}}E$ represent the contribution from the direct beam scattered into the diffuse field at the interface.

### 4. Resulting Layer Properties
Finally, the reflection ($R^{\text{res}}$) and transmission ($T^{\text{res}}$) matrices for the combined layer are:

$$
\begin{aligned}
R_{\text{top}}^{\text{res}} &= c T_{\text{bot}} W U + R_{\text{top}} + E U \\
T_{\text{top}}^{\text{res}} &= c T_{\text{top}} W D + T_{\text{top}} E + E D
\end{aligned}
$$

For the doubling method (where top and bottom layers are identical), $R_{\text{bot}}^{\text{res}} = R_{\text{top}}^{\text{res}}$ and $T_{\text{bot}}^{\text{res}} = T_{\text{top}}^{\text{res}}$.
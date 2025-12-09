# Adjoint Adding-Doubling Radiative Transfer Model (AAD-RTM)

This repository contains a C++ implementation of the **Adjoint Adding-Doubling Method** for plane-parallel radiative transfer problems.
Designed for inverse problems in planetary atmospheric science, this code computes the exact analytical gradients of the cost function with respect to atmospheric parameters.
It utilizes the Eigen library for matrix operations.
Supports both scattering and thermal emission.

# Forward Adding-Doubling Method Implementation

## Numerical Integration (Gauss-Radau Quadrature)

The code employs **Gauss-Radau quadrature** for discretizing the zenith angle $\mu$. Unlike the standard Gauss-Legendre quadrature, the Gauss-Radau method fixes one integration node at the boundary of the interval.

### Mathematical Formulation
The integral of a function $f(\mu)$ over $[0, 1]$ is approximated as:

$$
\int_0^1 f(\mu) \mu d\mu \approx \sum_{i=1}^{n} w_i f(\mu_i)
$$

where the last node is fixed at $\mu_n = 1.0$. This is particularly advantageous for radiative transfer problems where the radiation field at the vertical boundary ($\mu=1$) is of specific interest and should be calculated explicitly without extrapolation.

### Algorithm
The nodes $x_i$ (in the standard interval $[-1, 1]$) and weights $w_i^{gr}$ are computed using the **Golub-Welsch algorithm** adapted for Radau quadrature:

1.  **Jacobi Matrix Construction**: A symmetric tridiagonal matrix $J$ is constructed using the recurrence relation coefficients of the Jacobi polynomials $P_{n-1}^{(1, 0)}(x)$ (fixing the endpoint at $+1$).
2.  **Eigenvalue Decomposition**: The eigenvalues of $J$ correspond to the internal nodes $x_i$.
3.  **Mapping**: The nodes and weights are linearly mapped from $[-1, 1]$ to the physical interval $[0, 1]$.

$$
\mu_i = \frac{x_i + 1}{2}, \quad w_i = \frac{\mu_i w_i^{gr}}{2}
$$

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
Q_1 &= c R_{\mathrm{bot}} W R_{\mathrm{top}} \\
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
D &= c S W T_{\mathrm{top}} + T_{\mathrm{top}} + S E \\
U &= c R_{\mathrm{top}} W D + R_{\mathrm{top}} E
\end{aligned}
$$

Here, terms like $SE$ and $R_{\mathrm{top}}E$ represent the contribution from the direct beam scattered into the diffuse field at the interface.

### 4. Resulting Layer Properties
Finally, the reflection ($R^{\mathrm{res}}$) and transmission ($T^{\mathrm{res}}$) matrices for the combined layer are:

$$
\begin{aligned}
R_{\mathrm{top}}^{\mathrm{res}} &= c T_{\mathrm{bot}} W U + R_{\mathrm{top}} + E U \\
T_{\mathrm{top}}^{\mathrm{res}} &= c T_{\mathrm{top}} W D + T_{\mathrm{top}} E + E D
\end{aligned}
$$

For the doubling method (where top and bottom layers are identical), $R_{\mathrm{bot}}^{\mathrm{res}} = R_{\mathrm{top}}^{\mathrm{res}}$ and $T_{\mathrm{bot}}^{\mathrm{res}} = T_{\mathrm{top}}^{\mathrm{res}}$.

## Matrix Formulation of Adding Equations

This section describes the combination of two dissimilar layers: an upper layer (superscript $\mathrm{top}$) and a lower layer (superscript $\mathrm{bot}$).

**Notation Convention:**
* **Superscripts** ($\mathrm{top}, \mathrm{bot}$): Denote the physical layer (Upper layer, Lower layer).
* **Subscripts** ($\mathrm{top}, \mathrm{bot}$): Denote the direction of the incident light (incident from above, incident from below).

### 1. Direct Attenuation Matrix ($E$)
The exponential attenuation of the direct beam for the upper and lower layers is given by the diagonal matrices $E^\mathrm{top}$ and $E^\mathrm{bot}$:

$$
E^\mathrm{top}_{ij} = \delta_{ij} \exp\left(-\frac{\tau^\mathrm{top}}{\mu_i}\right),\ E^\mathrm{bot}_{ij} = \delta_{ij} \exp\left(-\frac{\tau^\mathrm{bot}}{\mu_i}\right)
$$

### 2. Multiple Scattering Interaction
The interaction between the two layers involves infinite reflections at the interface. This is solved using a Neumann series:

$$
\begin{aligned}
Q_1 &= c R^\mathrm{top}_{\mathrm{bot}} W R^\mathrm{bot}_{\mathrm{top}} \\
Q_2 &= c Q_1 W \\
S &= Q_1 + Q_2 Q_1 + Q_2 Q_2 Q_1 + \dots = (I - Q_2)^{-1} Q_1
\end{aligned}
$$

* $R^\mathrm{top}_{\mathrm{bot}}$: Reflection of the **upper layer** for light incident from **below**.
* $R^\mathrm{bot}_{\mathrm{top}}$: Reflection of the **lower layer** for light incident from **above**.
* $S$: The source term representing the summation of multiple scattering between layers.

### 3. Internal Radiation Fields
The downward ($D$) and upward ($U$) diffuse intensities at the interface between the two layers are calculated as:

$$
\begin{aligned}
D &= c S W T^\mathrm{top}_{\mathrm{top}} + T^\mathrm{top}_{\mathrm{top}} + S E^\mathrm{top} \\
U &= c R^\mathrm{bot}_{\mathrm{top}} W D + R^\mathrm{bot}_{\mathrm{top}} E^\mathrm{top}
\end{aligned}
$$

### 4. Resulting Layer Properties
Finally, the reflection ($R^{\mathrm{res}}_\mathrm{top}$) and transmission ($T^{\mathrm{res}}_\mathrm{top}$) matrices for the combined layer (incident from space) are:

$$
\begin{aligned}
R_{\mathrm{top}}^{\mathrm{res}} &= c T^\mathrm{top}_{\mathrm{bot}} W U + R^\mathrm{top}_{\mathrm{top}} + E^\mathrm{top} U \\
T_{\mathrm{top}}^{\mathrm{res}} &= c T^\mathrm{bot}_{\mathrm{top}} W D + T^\mathrm{bot}_{\mathrm{top}} E^\mathrm{top} + E^\mathrm{bot} D
\end{aligned}
$$

> **Note:**
> As seen in the equations above, the reflection and transmission matrices for light incident from below the lower layer ($R^\mathrm{bot}_{\mathrm{bot}}$ and $T^\mathrm{bot}_{\mathrm{bot}}$) are not required to compute $R_{\mathrm{top}}^{\mathrm{res}}$ and $T_{\mathrm{top}}^{\mathrm{res}}$. Therefore, for simulations considering only satellite observations, calculating the response for light incident from the bottom of the combined layer is unnecessary.
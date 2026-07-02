# Background

## Mechanical Homogenization
Mechanical homogenization seeks to find the macroscopic response of heterogeneous microstructures in order to replace them with equivalent homogeneous materials that reproduce this response.

Heterogeneous microstructures are defined on a unit cell ``\Omega \subset \mathbb{R}^d``, with ``d\in\lbrace 2, 3\rbrace``, typically representing a rectangular domain ``\Omega = \prod_{i=1}^d [0, L_i]``. Within the unit cell, the material behaviour is described by a spatially varying material tensor ``\mathbb{C}(\boldsymbol{x})``, ``\boldsymbol{x} \in \Omega``.

The macroscopic response of a microstructure is represented by the cell average defined as:
```math
\langle \mathcal{Q} \rangle = \frac{1}{\mid\Omega\mid}\int_{\Omega} \mathcal{Q}(\boldsymbol{x})\, d\boldsymbol{x},
```
where ``\mathcal{Q}(\boldsymbol{x})`` denotes a spatially varying quantity. Accordingly, the macroscopic stress and strain are given by ``\langle \boldsymbol{\sigma} \rangle`` and ``\langle \boldsymbol{\varepsilon} \rangle``, respectively, and represent the effective (macroscopic) response of the microstructure. Computing these quantities requires solving the micromechanical boundary value problem:
```math
\begin{align*}
    &\nabla \cdot \boldsymbol{\sigma}(\boldsymbol{x}) = \boldsymbol{0} \quad \forall \boldsymbol{x} \in \Omega\\
    &\boldsymbol{\sigma}(\boldsymbol{x}) = \mathbb{C}(\boldsymbol{x}) : \boldsymbol{\varepsilon}\big(\boldsymbol{u}(\boldsymbol{x})\big)\\
    &\boldsymbol{\varepsilon}\big(\boldsymbol{u}(\boldsymbol{x})\big) = \frac{1}{2}\Big(  \nabla\boldsymbol{u}(\boldsymbol{x}) + \big(\nabla\boldsymbol{u}(\boldsymbol{x})\big)^T  \Big)\\
    &\text{boundary conditions on } \partial\Omega.
\end{align*}
```
Boundary conditions can be prescribed displacements, prescribed tractions, or combinations thereof, including periodic boundary conditions.

`FFTMechHomo.jl` implements a solver for this class of micromechanical boundary value problems based on a Fast Fourier Transform (FFT) formulation. The implementation is restricted to small strains. An example demonstrating the computation of effective linear elastic properties is provided in the accompanying tutorials.

## Lippmann-Schwinger Formalism
The FFT-based homogenization approach implemented in `FFTMechHomo.jl` follows the seminal work of Moulinec and Suquet [moulinec_fast_1994,moulinec_numerical_1998](@cite). The formulation assumes a periodic unit cell and decomposes the displacement and strain fields into macroscopic and fluctuating parts,

```math
\boldsymbol{u}(\boldsymbol{x}) = \tilde{\boldsymbol{u}}(\boldsymbol{x}) + \boldsymbol{E}\boldsymbol{x}, \quad \boldsymbol{\varepsilon}(\boldsymbol{x}) = \boldsymbol{\varepsilon}\big(\tilde{\boldsymbol{u}}(\boldsymbol{x})\big)+ \boldsymbol{E},
```

where ``\boldsymbol{E}`` denotes a prescribed macroscopic strain and ``\tilde{\boldsymbol{u}}`` is periodic. Additionally, a homogeneous linear elastic reference material with stiffness tensor ``\mathbb{C}^0`` is introduced to rewrite the constitutive law as:

```math
\boldsymbol{\sigma}(\boldsymbol{x}) = \mathbb{C}^0 : \boldsymbol{\varepsilon}(\boldsymbol{x})+ \boldsymbol{\tau}(\boldsymbol{x}),
```

where ``\boldsymbol{\tau}(\boldsymbol{x}) = \big(\mathbb{C}(\boldsymbol{x}) - \mathbb{C}^0\big) : \boldsymbol{\varepsilon}(\boldsymbol{x})`` is a polarization field. Note that, although the reference material is linear elastic, the constitutive behaviour of the actual material phases is not restricted to be linear.

The boundary value problem can now be formulated as:

```math
\begin{align*}
    &\nabla \cdot \boldsymbol{\sigma}(\boldsymbol{x}) = \boldsymbol{0} \quad \forall \boldsymbol{x} \in \Omega\\
    &\boldsymbol{\sigma}(\boldsymbol{x}) = \mathbb{C}^0 : \boldsymbol{\varepsilon}(\boldsymbol{x}) + \boldsymbol{\tau}(\boldsymbol{x})\\
    &\boldsymbol{\varepsilon}(\boldsymbol{x}) = \frac{1}{2}\Big(  \nabla\boldsymbol{u}(\boldsymbol{x}) + \big(\nabla\boldsymbol{u}(\boldsymbol{x})\big)^T  \Big)\\
    &\tilde{\boldsymbol{u}}(\boldsymbol{x}) \text{ periodic}\\
    &\boldsymbol{\sigma}(\boldsymbol{x})\cdot\boldsymbol{n} \text{ anti-periodic}.
\end{align*}
```

This problem can be solved by introducing the periodic Green operator ``\boldsymbol{\Gamma}^0`` associated with the reference material. The Green operator maps the polarization field to the compatible strain field and yields the integral equation:

```math
\boldsymbol{\varepsilon}(\boldsymbol{x}) = -\boldsymbol{\Gamma}^0(\boldsymbol{x}) * \boldsymbol{\tau}(\boldsymbol{x}) + \boldsymbol{E},
```

which is commonly referred to as the Lippmann–Schwinger equation. Since the displacement fluctuations are periodic, all fields can be represented by Fourier series. In this representation, periodicity is satisfied by construction and the convolution appearing in the Lippmann–Schwinger equation reduces to a pointwise multiplication in Fourier space:

```math
\hat{\boldsymbol{\varepsilon}}(\boldsymbol{\xi}) = -\hat{\boldsymbol{\Gamma}}^0(\boldsymbol{\xi}) : \hat{\boldsymbol{\tau}}(\boldsymbol{\xi}) \quad \forall \boldsymbol{\xi} \neq \boldsymbol{0}, \; \hat{\boldsymbol{\varepsilon}}(\boldsymbol{0})=\boldsymbol{E}.
```

The Green operator in Fourier space ``\boldsymbol{\Gamma}^0`` can be derived analytically for linear elastic reference materials. 

In practice, the Fourier coefficients of the involved fields are computed efficiently using FFTs, and the Green operator is applied pointwise in Fourier space.

In real applications, microstructures are typically provided in a discretized form (e.g. voxel data). Accordingly, the underlying fields must be represented on a discrete grid. The original Moulinec–Suquet scheme employs a regular voxel-wise representation for all fields.

Since the polarization field depends on the unknown strain field, the Lippmann–Schwinger equation must be solved iteratively. 

Consequently, an FFT-based homogenization method is characterized by two fundamental ingredients: a discretization of the underlying fields and an iterative solution algorithm. Numerous improvements have been proposed for both aspects, including alternative discretizations and advanced iterative solvers; an extensive review is given in Schneider [schneider_review_2021](@cite).

The discretizations and solvers currently implemented in `FFTMechHomo.jl` are described in [Discretization](@ref) and [Solver](@ref), respectively.
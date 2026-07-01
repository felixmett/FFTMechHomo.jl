# Mechanical Homogenization
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
    &\boldsymbol{\sigma}(\boldsymbol{x}) = \mathbb{C}(\boldsymbol{x}) : \boldsymbol{\varepsilon}(\boldsymbol{x})\\
    &\boldsymbol{\varepsilon}(\boldsymbol{x}) = \frac{1}{2}\Big(  \nabla\boldsymbol{u}(\boldsymbol{x}) + \big(\nabla\boldsymbol{u}(\boldsymbol{x})\big)^T  \Big)\\
    &\text{boundary conditions on } \partial\Omega.
\end{align*}
```
Boundary conditions can be prescribed displacements, prescribed tractions, or combinations thereof, including periodic boundary conditions.

`FFTMechHomo.jl` implements a solver for this class of micromechanical boundary value problems based on a Fast Fourier Transform (FFT) formulation. The implementation is restricted to small strains and infinitesimal kinematics. An example demonstrating the computation of effective linear elastic properties is provided in the accompanying tutorials.

# Lippmann-Schwinger Formalism
- Assume a periodic microstructure
- Local strain field ``\boldsymbol{\varepsilon}\big(\boldsymbol{u}(\boldsymbol{x})\big)`` is split into its average ``\boldsymbol{E}`` and a periodic fluctuation term ``\boldsymbol{\varepsilon}\big(\tilde{\boldsymbol{u}}(\boldsymbol{x})\big)``
- Local strain ``\boldsymbol{\varepsilon}\big(\boldsymbol{u}(\boldsymbol{x})\big) = \boldsymbol{\varepsilon}\big(\tilde{\boldsymbol{u}}(\boldsymbol{x})\big) + \boldsymbol{E}``
- Local displacement ``\boldsymbol{u}(\boldsymbol{x}) = \tilde{\boldsymbol{u}}(\boldsymbol{x}) + \boldsymbol{E}\boldsymbol{x}``
- Refo
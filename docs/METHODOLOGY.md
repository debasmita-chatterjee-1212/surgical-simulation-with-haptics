# Computational Methodology

## 1. Overview

The surgical simulation framework investigates computational
methods for interactive anatomical tissue deformation,
instrument–tissue interaction and force estimation.

The framework consists of three simulation applications:

- Rigid temporal bone probing.
- Deformable soft-organ palpation.
- Soft-tissue surface incision.

The computational approach combines anatomical geometry,
linear-elastic tissue modelling, tensor-mass discretization,
explicit numerical integration and instrument contact mechanics.

The graphical applications are implemented in C++ using
CHAI3D and GLFW.

## 2. Constitutive Modelling

### 2.1. Linear Elasticity

Soft tissue is approximated as a homogeneous, isotropic,
linear-elastic material within each anatomical model.

Under the small-strain assumption, the strain tensor is:

```math
\boldsymbol{\varepsilon}
= \frac{1}{2}
\left(\nabla\mathbf{u}+(\nabla\mathbf{u})^T\right)
```

where the displacement field is represented by u.

The stress–strain relationship is expressed as:

```math
\boldsymbol{\sigma}
= \lambda\,\mathrm{tr}(\boldsymbol{\varepsilon})\mathbf{I}
+ 2\mu\boldsymbol{\varepsilon}
```

where:

- σ is the Cauchy stress tensor.
- ε is the infinitesimal strain tensor.
- λ and μ are the Lamé parameters.
- I is the identity tensor.

### 2.2. Material Parameters

The Lamé parameters are calculated from Young's modulus
and Poisson's ratio:

```math
\lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}
```

```math
\mu = \frac{E}{2(1+\nu)}
```

Material-specific parameters allow the framework to
represent different mechanical responses across
anatomical tissues.

| Anatomical tissue | Young's modulus |
|---|---:|
| Pancreas | 1.4 kPa |
| Liver | 15 kPa |
| Kidney | 19.6 kPa |
| Ear cartilage | 1.66 MPa |
| Temporal bone | 10 GPa |

The material values represent modelling parameters rather
than patient-specific measurements.

Further details regarding material properties, assumptions
and parameter sources are provided in
[Material Parameters](MATERIALS.md).

## 3. Tensor-Mass Formulation

### 3.1. Model Selection

The tensor-mass formulation combines the local connectivity
of mass-based discrete models with stiffness tensors derived
from continuum elasticity.

The formulation follows the approach described by
Cotin, Delingette and Ayache (2000).

Unlike conventional scalar spring-mass formulations,
the interaction between connected vertices is represented
using matrix-valued stiffness coefficients.

This allows the model to incorporate directional
coupling derived from the elastic constitutive law.

### 3.2. Tetrahedral Discretization

The anatomical tissue volume is discretized using
linear tetrahedral elements.

Each tetrahedron contains four vertices, with a
three-dimensional displacement vector associated
with each vertex.

The displacement of vertex i is defined as:

```math
\mathbf{u}_i = \mathbf{P}_i-\mathbf{P}_i^0
```

where:

- P is the current vertex position.
- P⁰ is the corresponding rest position.

The linear tetrahedral shape functions have constant
spatial gradients within each element.

### 3.3. Element Stiffness Tensors

For a tetrahedral element T, the stiffness block
associated with local vertices k and l is:

```math
\mathbf{K}_{kl}^{T}
=
V_T\left[
\lambda\mathbf{b}_k\mathbf{b}_l^T
+\mu\mathbf{b}_l\mathbf{b}_k^T
+\mu(\mathbf{b}_k\cdot\mathbf{b}_l)\mathbf{I}
\right]
```

where:

- V is the tetrahedral element volume.
- b is the gradient of the corresponding shape function.
- λ and μ are the Lamé parameters.

Each stiffness block is a 3 × 3 matrix.

The element contributions are assembled into
vertex-associated and edge-associated stiffness tensors.

### 3.4. Internal Force Computation

The restoring force acting on vertex i is calculated as:

```math
\mathbf{F}_i^{\mathrm{int}}
=
-\left(
\mathbf{K}_{ii}\mathbf{u}_i
+\sum_{j\in\mathcal{N}(i)}
\mathbf{K}_{ij}\mathbf{u}_j
\right)
```

where N(i) represents the neighbouring vertices
connected to vertex i.

The stiffness tensors are precomputed from the
undeformed tetrahedral geometry.

The symmetry relationship

```math
\mathbf{K}_{ji}=\mathbf{K}_{ij}^{T}
```

allows a single stiffness tensor to be stored
for each connected edge.

This avoids repeatedly assembling element stiffness
matrices during deformation updates.

### 3.5. Lumped Mass Representation

The mass of each tetrahedral element is distributed
equally among its four vertices.

The resulting nodal mass is:

```math
m_i=\sum_{T\ni i}\frac{\rho V_T}{4}
```

where ρ is the material density.

The lumped mass representation permits independent
nodal position updates during explicit integration.

### 3.6. Boundary Conditions

Selected vertices are constrained to their rest positions
to represent anatomical attachment and support regions.

External instrument forces are applied to unconstrained
vertices through the tissue interaction model.

## 4. Numerical Integration

### 4.1. Governing Equation

The motion of each unconstrained vertex is governed by
Newton's second law with viscous damping:

```math
m_i\ddot{\mathbf{P}}_i
+c_i\dot{\mathbf{P}}_i
=
\mathbf{F}_i^{\mathrm{int}}
+\mathbf{F}_i^{\mathrm{ext}}
```

where:

- m is the lumped nodal mass.
- c is the non-negative damping coefficient.
- F_int is the internal restoring force.
- F_ext is the externally applied force.

The damping coefficient is proportional to nodal mass:

```math
c_i=\alpha m_i
```

where α is the damping coefficient per unit mass.

### 4.2. Explicit Central-Difference Scheme

Velocity and acceleration are approximated using
central finite differences:

```math
\dot{\mathbf{P}}_i^n
\approx
\frac{\mathbf{P}_i^{n+1}-\mathbf{P}_i^{n-1}}
{2\Delta t}
```

```math
\ddot{\mathbf{P}}_i^n
\approx
\frac{\mathbf{P}_i^{n+1}-2\mathbf{P}_i^n
+\mathbf{P}_i^{n-1}}{\Delta t^2}
```

Substitution into the governing equation gives the
explicit position update:

```math
\mathbf{P}_i^{n+1}
=
\frac{
\mathbf{F}_i^n
+\frac{2m_i}{\Delta t^2}\mathbf{P}_i^n
-B_i\mathbf{P}_i^{n-1}
}{A_i}
```

where:

```math
A_i=\frac{m_i}{\Delta t^2}+\frac{c_i}{2\Delta t}
```

```math
B_i=\frac{m_i}{\Delta t^2}-\frac{c_i}{2\Delta t}
```

The explicit formulation avoids solving a global
linear system at every integration step.

### 4.3. Numerical Stability

Explicit central-difference integration is conditionally
stable.

For the undamped linear system, the critical timestep
is related to the highest natural frequency:

```math
\Delta t_{\mathrm{crit}}
=\frac{2}{\omega_{\max}}
```

The maximum natural frequency is related to the
largest eigenvalue of the mass-normalized stiffness
operator:

```math
\omega_{\max}^2
=\lambda_{\max}(\mathbf{M}^{-1}\mathbf{K})
```

The implementation estimates the dominant eigenvalue
using power iteration.

A safety factor is applied to the estimated critical
timestep, with additional consideration of contact
stiffness and configured timestep limits.

Increasing material stiffness or reducing element
size generally imposes more restrictive timestep
requirements.

This is particularly significant when simulating
rigid anatomical structures such as temporal bone.

### 4.4. Computational Considerations

The implementation incorporates:

- Precomputed vertex and edge stiffness tensors.
- Lumped nodal masses.
- Explicit nodal position updates.
- Multithreaded internal force computation.
- Numerical timestep estimation.
- Fixed-vertex boundary conditions.

These design choices reduce the computational
requirements of interactive tissue deformation.

No specific hardware feedback frequency or computational
performance benchmark is claimed here.

## 5. Numerical Verification

The tensor-mass formulation was investigated through
analytical comparisons involving:

- Uniaxial compression.
- Cantilever beam bending.
- Distributed loading of an elastic half-space.
- Integrated proxy contact.

The experiments examine displacement accuracy,
mesh sensitivity and the influence of Poisson's ratio.

Detailed experimental configurations, analytical
references, numerical results and limitations are
documented in the
[Numerical Verification Report](../validation/README.md).

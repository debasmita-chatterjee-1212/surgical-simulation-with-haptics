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

## 6. Anatomical Geometry and Spatial Discretization

### 6.1. Anatomical Models

The simulation framework incorporates anatomical geometries
representing both rigid and deformable biological tissues.

The anatomical geometry provides the spatial domain for
instrument interaction, deformation and material removal.

The computational representation is selected according to
the mechanical characteristics of the tissue and the
requirements of the surgical interaction.

### 6.2. Volumetric Mesh Representation

Deformable anatomical structures are represented using
volumetric tetrahedral meshes.

The computational mesh contains:

- Three-dimensional vertex coordinates.
- Tetrahedral element connectivity.
- Material-specific mechanical parameters.
- Fixed vertices representing support regions.

The volumetric representation allows internal mechanical
forces to be computed throughout the tissue rather than
only on its external surface.

Geometry preprocessing utilities are retained in the
[tools directory](../tools/).

### 6.3. Boundary Surface Extraction

The visible tissue surface is extracted from the active
tetrahedral mesh.

Each tetrahedron contributes four triangular faces.
Faces shared by two active tetrahedra are classified
as internal faces.

Faces belonging to only one active tetrahedron form
the external boundary surface.

The boundary triangles are oriented such that their
surface normals point away from the active tissue volume.

This representation supports:

- Anatomical visualization.
- Surface contact detection.
- Instrument interaction.
- Visualization of newly exposed tissue after cutting.

### 6.4. Surface Smoothing

The discretized boundary may exhibit staircase-like
geometrical artefacts.

Taubin smoothing is used to improve the representation
of the exposed tissue surface.

The smoothing offsets are maintained separately from
the underlying simulated vertex positions.

This allows surface geometry to be improved without
directly replacing the mechanical degrees of freedom
of the tensor-mass model.


## 7. Rigid Temporal Bone Simulation

### 7.1. Motivation for Rigid Modelling

Temporal bone has a substantially higher elastic modulus
than soft biological tissue.

Direct explicit integration of a highly stiff material
requires extremely small numerical timesteps.

Furthermore, the deformation of bone under typical
probing forces is considerably smaller than the
deformation observed in soft tissue.

A rigid approximation therefore provides a computationally
efficient representation of intact bone for interactive
instrument contact.

### 7.2. Rigid–Deformable Modelling Strategy

The framework distinguishes between:

- Rigid anatomical structures, for which elastic
  deformation is neglected.
- Deformable anatomical structures, for which nodal
  displacements are computed using the tensor-mass model.

In rigid mode, the uncut tissue is not advanced through
the explicit elastic integration scheme.

Instrument interaction is instead represented through
contact detection and reaction-force estimation.

### 7.3. Penalty-Based Contact

Penalty-based contact approximates the reaction force
generated when an instrument intersects the tissue
contact surface.

A general linear penalty formulation is:

```math
F_n=k_c\delta
```

where:

- F_n is the normal contact force.
- k_c is the contact stiffness.
- δ is the penetration depth.

The area-weighted contact formulation distributes
contact contributions across surface triangles.

A representative local force contribution is:

```math
\mathbf{F}_c=k_A A_s\delta\mathbf{n}
```

where:

- k_A is the contact stiffness per unit area.
- A_s is the associated surface area.
- n is the contact normal.

Force contributions are distributed to the associated
surface vertices.

The reaction force is limited to a configured maximum
to prevent excessively large force estimates.

The penalty parameters are computational modelling
parameters and are not experimentally calibrated
bone contact constants.

### 7.4. Bone Material Removal

The volumetric representation also supports progressive
material removal.

Tetrahedral elements within the cutting region may be
deactivated, exposing new boundary surfaces.

The instrument is constrained relative to the remaining
intact material.

The estimated reaction force represents resistance
associated with the contact model.

It does not constitute a validated drilling-force
prediction accounting for burr rotational speed,
feed rate or bone microstructure.

The rigid bone application is documented in:

[Temporal Bone Simulation](../examples/42-bone-probe/)


## 8. Deformable Soft-Tissue Interaction

### 8.1. Proxy-Based Contact

Soft-tissue palpation uses a proxy-based contact approach
inspired by the god-object method used in haptic rendering.

The interaction model distinguishes between:

- The desired instrument position, referred to as
  the hand position.
- The constrained instrument position, referred to
  as the proxy position.

The proxy is represented by a spherical contact geometry.

The hand position may move towards the tissue interior,
while the proxy is constrained relative to the
deforming tissue surface.

This separation reduces the penetration artefacts
associated with direct penalty-based interaction.

### 8.2. Proxy Position Update

During each physics update, the proxy moves towards
the desired hand position.

The displacement is limited according to a configured
maximum instrument velocity.

Contact constraints are evaluated using nearby
surface triangles.

The closest point on each candidate triangle is
computed to determine whether the spherical proxy
intersects the surface.

The proxy position is corrected to resolve detected
penetration.

Previously identified contact normals are also used
to restrict motion into the tissue surface while
allowing sliding and retreating motion.

### 8.3. Virtual Coupling

The hand and proxy positions are coupled through
a virtual spring.

The unconstrained spring force is expressed as:

```math
\mathbf{F}_s=k_s(\mathbf{x}_p-\mathbf{x}_h)
```

where:

- k_s is the virtual coupling stiffness.
- x_p is the proxy position.
- x_h is the desired hand position.

The implemented interaction uses the compressive
normal component of the spring force and applies
a configured force limit.

The resulting force represents the estimated
reaction experienced by the instrument.

An opposing load is applied to the tissue model.

### 8.4. Distributed Surface Loading

Applying the entire instrument force to a single
mesh vertex can produce unrealistic localized
deformation, particularly when the instrument
radius is smaller than the mesh resolution.

The contact force is therefore distributed across
a neighbourhood of exposed surface vertices.

The load distribution radius is determined from
both instrument size and mesh resolution:

```math
r_L=\max(aR,b\bar{l})
```

where:

- R is the probe radius.
- l̄ is the mean surface edge length.
- a and b are configurable scaling parameters.

The radial weighting function decreases quadratically
with distance from the contact region.

A mass-weighted representation is:

```math
w_i=m_i\left(1-\frac{d_i}{r_L}\right)^2
```

for vertices within the loading radius.

Here, d_i is the distance from the contact region
and m_i is the lumped nodal mass.

The normalized force contribution is:

```math
\mathbf{F}_i=
\frac{w_i}{\sum_j w_j}
\mathbf{F}_{\mathrm{tissue}}
```

This distributes the applied load while preserving
the total force across the selected vertices.

The approach is a numerical approximation of
distributed instrument loading rather than a
complete analytical contact mechanics solution.

### 8.5. Coupled Deformation

The external contact loads are incorporated into
the governing equations of the tensor-mass model.

The tissue deformation is advanced through explicit
numerical integration.

The updated tissue geometry subsequently influences
the proxy contact configuration.

The resulting interaction forms a coupled loop
between instrument position, contact force and
tissue deformation.

### 8.6. Contact Stability

Virtual coupling introduces an additional stiffness
into the interaction system.

The characteristic contact timestep is estimated
using a spring–mass approximation:

```math
\Delta t_{\mathrm{contact}}
\approx 2\sqrt{\frac{m_{\mathrm{eff}}}{k_s}}
```

where m_eff represents an estimated effective
mass of the loaded surface region.

The simulation timestep accounts for both mesh
stability and contact-related stability constraints.

This estimate is a practical numerical safeguard
rather than a formal guarantee of stability for
every nonlinear contact configuration.

### 8.7. Numerical Assessment

The isolated tissue deformation model was evaluated
using distributed loading of an elastic half-space.

The integrated proxy contact model was additionally
investigated through force–indentation experiments.

The latter experiments identified discrepancies
between proxy indentation and maximum tissue
displacement, together with unsteady contact
in selected configurations.

These findings are discussed in the
[Numerical Verification Report](../validation/README.md).

The deformable tissue application is documented in:

[Soft-Organ Palpation](../examples/43-soft-organ-probe/)


## 9. Surgical Cutting and Material Removal

### 9.1. Objective

Surgical cutting extends instrument interaction
beyond reversible tissue deformation.

Unlike palpation, cutting introduces permanent
changes to the computational representation
of the tissue.

The framework investigates these changes through
material removal and updates to the exposed
tissue geometry.

### 9.2. Progressive Material Removal

The shared simulation framework supports an
energy-based material removal approach.

A characteristic energy requirement for removing
a tetrahedral element is expressed as:

```math
W_T=q_cV_T
```

where:

- W_T is the characteristic removal energy.
- q_c is the configured cutting energy per unit volume.
- V_T is the tetrahedral element volume.

The available cutting power determines the energy
budget over a simulation interval:

```math
\Delta W=P_c\Delta t
```

where P_c is the configured cutting power.

The cutting energy density and power are prototype
parameters rather than experimentally validated
biological tissue properties.

### 9.3. Topological Updates

Material removal changes the connectivity of
the active computational domain.

Following tetrahedral element deactivation,
the simulation updates the affected mechanical
and geometrical representations.

These updates include:

- Active tetrahedral elements.
- Nodal masses.
- Internal stiffness contributions.
- Exposed boundary triangles.
- Instrument contact geometry.
- Numerical stability estimates.

Fragments disconnected from the constrained
support region can also be removed.

The resulting surface represents the remaining
active tissue volume.

### 9.4. Scalpel Incision

A separate simulation application investigates
soft-tissue interaction using a scalpel.

The application extends the surgical simulation
framework towards incision-oriented interaction.

The scalpel-specific formulation is documented
separately from the volumetric material removal
model to distinguish the two cutting approaches.

[Soft-Organ Scalpel Simulation](../examples/44-soft-organ-scalpel/)

### 9.5. Cutting Model Limitations

The computational cutting models are intended
to investigate interactive surgical manipulation.

They do not constitute a complete fracture
mechanics description of biological tissue.

Quantitative prediction of surgical cutting forces
would require additional experimental calibration
and consideration of tissue-specific failure
mechanisms.


## 10. Computational Architecture

The framework separates the underlying mechanical
simulation from the graphical application.

### 10.1. Core Computational Components

| Component | Responsibility |
|---|---|
| Tensor-mass model | Mesh connectivity and stiffness tensors |
| Dynamics solver | Nodal masses, forces and explicit integration |
| Surgical interaction model | Contact and material removal |
| Geometry representation | Boundary extraction and surface updates |
| Graphical application | Visualization and instrument interaction |

### 10.2. Computational Optimizations

Several implementation strategies are employed
to reduce computational overhead:

- Precomputation of stiffness tensors.
- Local vertex and edge connectivity.
- Multithreaded internal force computation.
- Lumped mass representation.
- Localized contact candidate selection.
- Warm-started timestep estimation after cutting.
- Rigid approximation for highly stiff tissues.

These strategies support interactive simulation
without requiring a global implicit system solve
during every physics update.

### 10.3. Force Estimation and Haptic Interface

The computational interaction model accepts
a desired instrument position and returns an
estimated reaction-force vector.

This separation provides a computational interface
for potential integration with a physical
haptic device.

The original research investigated surgical
simulation using a Phantom Omni device.

The original hardware interface implementation
is not included in the present repository.

Consequently, the available implementation
documents instrument interaction and force
estimation rather than a complete operational
hardware feedback system.


## 11. Modelling Assumptions and Limitations

The principal limitations of the framework are:

- Homogeneous and isotropic material assumptions.
- Linear elasticity and infinitesimal strain.
- Constant-strain tetrahedral discretization.
- Mesh-dependent bending and indentation accuracy.
- Volumetric locking at higher Poisson's ratios.
- Conditional stability of explicit integration.
- Approximate instrument contact mechanics.
- Simplified material removal and cutting criteria.
- Absence of patient-specific mechanical calibration.
- Absence of physical tissue validation.
- Absence of the original haptic hardware interface
  from the archived implementation.

The framework is a computational research prototype.

The numerical experiments establish selected
properties of the underlying deformation model,
but do not establish clinical validity or
surgical training effectiveness.


## 12. Theoretical Foundations

The computational approach draws upon established
research in deformable tissue simulation, continuum
mechanics and haptic interaction.

### Tensor-Mass Modelling

Cotin, S., Delingette, H., and Ayache, N. (2000).
A hybrid elastic model for real-time cutting,
deformations, and force feedback for surgery training
and simulation.
The Visual Computer, 16, 437–452.

### Proxy-Based Haptic Interaction

Zilles, C. B., and Salisbury, J. K. (1995).
A constraint-based god-object method for haptic display.
IEEE/RSJ International Conference on Intelligent
Robots and Systems.

Ruspini, D. C., Kolarov, K., and Khatib, O. (1997).
The haptic display of complex graphical environments.
Proceedings of SIGGRAPH.

### Additional Documentation

- [Material Parameters](MATERIALS.md)
- [Numerical Verification](../validation/README.md)
- [Anatomical Geometry Processing](../tools/)
- [Simulation Applications](../examples/)

# System Architecture

## Overview

The surgical simulation framework separates anatomical
geometry processing, mechanical modelling, instrument
interaction and graphical visualization.

This document describes the computational architecture
and information flow between the principal components.

For mathematical formulations, refer to the
[Computational Methodology](METHODOLOGY.md).

## 1. Overall System Architecture

```mermaid
flowchart TD
    A[Anatomical Geometry] --> B[Geometry Preprocessing]
    B --> C[Volumetric Mesh]
    M[Material Parameters] --> D[Mechanical Model]
    C --> D
    D --> E[Rigid Bone Model]
    D --> F[Tensor-Mass Model]
    I[Instrument Input] --> G[Contact and Cutting]
    E --> G
    F <--> G
    G --> H[Reaction Force Estimation]
    G --> J[Surface Updates]
    J --> K[CHAI3D Visualization]
    H --> K
```

The framework adopts a rigid approximation for intact
bone and a deformable tensor-mass representation for
soft tissues.

Instrument interaction is coupled to the mechanical
model through contact forces and geometry updates.

## 2. Anatomical Geometry Processing

```mermaid
flowchart TD
    A[Anatomical Geometry] --> B[Geometry Preprocessing]
    B --> C[Volumetric Discretization]
    C --> D[Tetrahedral Connectivity]
    D --> E[Vertex Coordinates]
    D --> F[Element Connectivity]
    D --> G[Boundary Constraints]
    E --> H[Computational Mesh]
    F --> H
    G --> H
    H --> I[Boundary Surface Extraction]
    I --> J[Visualization and Contact Geometry]
```

The computational representation contains vertex
coordinates, tetrahedral connectivity and constrained
boundary vertices.

The exposed surface provides the geometry required
for visualization and instrument contact.

Geometry preprocessing utilities are retained in
the [tools directory](../tools/).

## 3. Tensor-Mass Solver

```mermaid
flowchart TD
    A[Rest Mesh and Material Properties] --> B[Lamé Parameters]
    B --> C[Element Stiffness Tensors]
    C --> D[Vertex and Edge Assembly]
    A --> E[Lumped Nodal Masses]
    D --> F[Internal Force Computation]
    E --> G[Timestep Estimation]
    D --> G
    H[External Instrument Forces] --> I[Total Nodal Forces]
    F --> I
    I --> J[Central Difference Integration]
    G --> J
    J --> K[Updated Vertex Positions]
    K --> L[Boundary Surface Update]
    K --> F
```

The mechanical solver uses precomputed stiffness
tensors, lumped nodal masses and explicit numerical
integration.

Internal force computation and nodal position updates
form the principal deformation loop.

## 4. Instrument–Tissue Interaction

```mermaid
flowchart TD
    A[Desired Instrument Position] --> B[Proxy Position Update]
    B --> C[Surface Contact Detection]
    C --> D{Contact Detected?}
    D -->|No| E[Free Instrument Motion]
    D -->|Yes| F[Proxy Position Correction]
    F --> G[Virtual Spring Force]
    G --> H[Normal Force and Force Limit]
    H --> I[Estimated Reaction Force]
    H --> J[Distributed Tissue Loading]
    J --> K[Tensor-Mass Dynamics]
    K --> L[Updated Tissue Surface]
    L --> C
    E --> M[Graphical Visualization]
    I --> M
    L --> M
```

The deformable tissue interaction model distinguishes
between the desired instrument position and its
constrained proxy position.

The estimated contact force is distributed across
the tissue surface. The resulting deformation
subsequently influences the contact configuration.

The original Phantom Omni hardware interface is
not included in this repository.

## 5. Numerical Verification Framework

```mermaid
flowchart TD
    A[Numerical Verification] --> B[Solid Deformation Model]
    A --> C[Integrated Contact Model]
    B --> D[Uniaxial Compression]
    B --> E[Cantilever Bending]
    B --> F[Distributed Indentation]
    D --> G[Axial Elasticity Solution]
    E --> H[Timoshenko Beam Theory]
    F --> I[Elastic Half-Space Solution]
    G --> J[Displacement Error]
    H --> J
    I --> J
    J --> K[Mesh Sensitivity Analysis]
    C --> L[Force Indentation Assessment]
    L --> M[Proxy and Tissue Displacement]
    M --> N[Contact Steadiness]
    K --> O[Results and Limitations]
    N --> O
```

The solid deformation model is assessed through
comparisons with analytical reference solutions.

The integrated contact model is assessed separately
through force–indentation consistency and contact
steadiness.

The verification framework operates independently
of the graphical application.

Detailed results are available in the
[Numerical Verification Report](../validation/README.md).

## 6. Principal Implementation Components

| Component | Responsibility |
|---|---|
| `tensor_mass.h` | Tensor-mass formulation and stiffness assembly |
| `tensor_mass_dynamics.h` | Force computation and numerical integration |
| `surgical_sim.h` | Instrument interaction and geometry updates |
| `tools/` | Anatomical geometry preprocessing |
| `examples/` | Graphical surgical simulation applications |
| `validation/` | Numerical verification and experimental records |

## Related Documentation

- [Computational Methodology](METHODOLOGY.md)
- [Material Parameters](MATERIALS.md)
- [Numerical Verification](../validation/README.md)
- [Simulation Applications](../examples/)
- [Project Overview](../README.md)

# FVCOM-MP: Microplastic Transport Model

## Overview

This is a microplastic-adapted version of the **Finite Volume Community Ocean Model (FVCOM)**, originally developed for simulating unstructured grid, finite-volume 3D coastal ocean circulation. FVCOM provides robust capabilities for modeling complex coastal and estuarine systems with irregular boundaries.

## Project Information

**Development**: WaterPACT Project  
**Funding**: U.S. Department of Energy Water Power Technologies Office (WPTO) and Pacific Northwest National Laboratory (PNNL) Laboratory Directed Research and Development (LDRD)  
**Institution**: Pacific Northwest National Laboratory

## Features

This enhanced FVCOM version includes comprehensive microplastic transport modeling capabilities:

### Core Microplastic Module (`mod_plast.F`)

**Completed Features:**
- **Transport Processes**: Advection-diffusion with settling and resuspension
- **Size-Dependent Physics**: Multi-class plastic particles with varying properties
- **Sediment Interactions**: Bed deposition, erosion, and active layer dynamics
- **Toxic Vector Modeling**: Embedded and free toxin transport with leaching kinetics
- **Shape-Dependent Processes**: Particle orientation effects on transport
- **Infiltration Module**: Plastic transport into sediment layers with DLVO+strain filtration theory, including retention/blockage kinetics and burial processes

**Recent Developments:**
- **Microplastic-Sediment Flocculation**: Integration of MP-floc aggregation/disaggregation processes with size-dependent exclusion mechanisms, parameterized using experimental data from N. Wu et al. (2023)
- Advanced infiltration physics with mechanistic filtration theory
- Spatially-varying diffusion in sediment matrix
- Multi-layer sediment column transport

### Upcoming Development

**Planned Enhancements:**
- Further refinement of flocculation kinetics under varying turbulence conditions
- Integration with biofilm growth processes

## Model Capabilities

- Unstructured triangular grids for complex coastal geometries
- Multiple plastic size classes and shapes
- Coupled hydrodynamic-plastic transport
- Comprehensive source/sink terms
- NetCDF-based I/O for analysis and visualization

## Usage

This model is designed for researchers studying microplastic fate and transport in coastal and estuarine systems, with applications in environmental impact assessment and pollution source tracking.

---

*For technical details, see the comprehensive documentation in `mod_plast.F` header sections.*
# FVCOM-MP Microplastic Infiltration Implementation Plan

## Overview
This document provides a comprehensive analysis of the FVCOM-MP codebase and outlines the plan for implementing 1D infiltration of microplastics (MP) and associated toxics into sediment columns after bottom settling.

## Current Code Structure Analysis

### Main Module: mod_plast.F (37,333 tokens)
**Primary Location**: `mod_plast.F:119-end`
**Key Dependencies**: 
- `mod_par.F`, `mod_prec.F`, `mod_types.F`, `mod_wd.F`
- `Control` module for `TOXBDIS`, `TOXFDIS`, `PLASTDIS`
- `all_vars` for `CNSTNT`, `UNIFORM`

**Current Features**:
- Microplastic transport and fate modeling 
- Toxic-vector functionality with bonded/free toxic states
- Bottom boundary interaction (5 modes including concrete_bottom, active_bottom, etc.)
- Settling velocity calculations using Mehta (2021) method
- Current bottom treatment modes (line 33-46 in generic_plastic.inp):
  1. `levitation`: MP levitated with no settling/resuspension
  2. `concrete_bottom`: MP can settle/resuspend but no sediment bed interaction
  3. `active_bottom`: MP interaction with sediment in active layer
  4. `bedload_bottom`: Horizontal MP movement with sediment through bedload
  5. `burial_bottom`: **Target for infiltration implementation**

### Key Subroutines Identified (mod_plast.F)
1. `Setup_PLAST` (line 392): Parameter setup and initialization
2. `Advance_plast` (line 749): Main time-stepping routine
3. `Calc_Deposition` (line 1674): Current settling implementation
4. `Calc_Wset_Plastic` (line 1767): Settling velocity calculations

### Relevant Files for Infiltration Implementation

#### Core FVCOM Files
1. **mod_plast.F** - Main plastic module (primary modification target)
2. **mod_main.F** - Global variables and namelist structures
3. **mod_input.F** - Input module for namelist reading
4. **mod_ncdio.F** - NetCDF output for concentration/bottom/bed fields
5. **mod_startup.F** - Variable setup before computation
6. **fvcom.F** - Main program calling plastic setup
7. **internal_step.F** - Calls plastic advance routine
8. **init_plast.F** - Plastic model initialization

#### Sediment-Related Files (for reference and interaction)
1. **mod_sed.F** - Sediment transport module
2. **mod_bbl.F** - Bottom boundary layer module
3. **init_sed.F** - Sediment initialization
4. **mod_cstms_vars.F** - Community Sediment Transport Model variables
5. **mod_sed_cstms.F** - CSTMS implementation

#### Configuration Files
1. **FVCOM_MP_run/waterPACT_run.nml** - Main namelist file
   - Line 269-276: Current sediment/plastic model settings
   - Line 275: `PLASTIC_MODEL = T`
   - Line 276: `PLASTIC_MODEL_FILE = 'generic_plastic.inp'`

2. **FVCOM_MP_run/generic_plastic.inp** - Plastic model parameters
   - Line 60: `PLAST_BEDSET = concrete_bottom` (needs upgrade to burial_bottom)
   - Line 33-46: Bottom treatment options documentation

3. **FVCOM_MP_run/makefile** - Compilation directives
   - Line 76: `mod_plast.F` included in MODS
   - Line 93: `init_plast.F` included in MAIN

### C Preprocessor Directives Usage
**Primary Directive**: `#if defined (PLASTIC)` (line 121 in mod_plast.F)
**Pattern**: Mixed C/Fortran compilation using CPP preprocessor
- `.F` files processed through CPP before Fortran compilation
- Conditional compilation for different model configurations
- `WET_DRY` directive for dry node handling (line 1692)

## Current Settling Implementation and Infiltration Integration

### Calc_Deposition Subroutine (mod_plast.F:1674-1760)
**Current Process**:
1. Loop over plastic classes and grid nodes
2. Skip dry nodes if WET_DRY is defined
3. Calculate settling velocity using `Calc_Wset_Plastic`
4. Apply flux-limited settling equation with CFL condition
5. Store depositional flux (`dflx`) and floating flux (`fflx`)
6. Handle bonded toxic deposition if toxic-vector is enabled

**Key Variables and Integration Points**:
- `dflx`: depositional flux into bed [kg m^-2] → **FEED TO INFILTRATION**
- `iflx`: infiltration flux into sediment [kg m^-2] → **NEW VARIABLE ADDED** ✓
- `fflx`: floating flux
- `mass`: plastic mass density at bottom [kg m^-2]
- `wset`: settling velocity [m/s]
- `tiflx_b`, `tiflx_f`: toxic infiltration fluxes → **NEW VARIABLES ADDED** ✓

**Infiltration Integration Strategy**:
The settled mass from `dflx` becomes the source for infiltration processes, tracked through the new `iflx` arrays for each species.

## Proposed Infiltration Implementation Plan

### Phase 1: Data Structure Implementation (COMPLETED)

#### 1.1 Infiltration Variables Added to plast_type (mod_plast.F:171-267)

The user has successfully implemented a streamlined, species-focused approach for infiltration variables:

**Core Infiltration Variables (lines 238, 242, 256, 261, 266-267):**
```fortran
! Species-specific sediment concentrations for each MP type
real(sp), allocatable :: conc_sed(:,:)      ! [node, layer] MP concentration in sediment
real(sp), allocatable :: tconc_b_sed(:,:)   ! [node, layer] Embedded toxic concentration in sediment
real(sp), allocatable :: tconc_f_sed(:,:)   ! [node, layer] Free toxic concentration in sediment

! Species-specific infiltration fluxes
real(sp), allocatable :: iflx(:)            ! [node] MP infiltration flux (+into bed)
real(sp), allocatable :: tiflx_b(:)         ! [node] Embedded toxic infiltration flux 
real(sp), allocatable :: tiflx_f(:)         ! [node] Free toxic infiltration flux
real(sp), allocatable :: tiflx_f_deep(:)    ! [node] Deep free toxic infiltration flux
```

**DLVO Physics Parameters (lines 231-233):**
```fortran
! Added to each plastic species for individual parameterization
real(sp) :: A_Hamaker                       ! [-] Hamaker constant for DLVO theory
real(sp) :: ws_sed                           ! [m/s] MP settling velocity in sediment matrix
real(sp) :: D_plast                          ! [m^2/s] Species-specific diffusivity
```

**Key Design Features:**
- **Species-Specific Approach**: Each MP type (iplast=1,Nplast) has individual infiltration parameters
- **I/O Ready**: All sediment concentration arrays designed for NetCDF output 
- **Flux Tracking**: Separate flux variables for MP, embedded toxics, and free toxics
- **Physics Integration**: Hamaker constant and settling velocity per species for DLVO calculations
- **Memory Efficient**: Streamlined variable set focused on essential infiltration physics

#### 1.2 Global Infiltration Switch (IMPLEMENTED)

**Added to Public Variables (mod_plast.F:289):**
```fortran
! plast_infil    : switch of plastic infiltration
logical :: plast_infil = .false.
```

**Namelist Integration Required:**
The user has prepared for namelist integration by adding the `plast_infil` switch. The following parameters should be added to `generic_plastic.inp`:

```
!===============================================================================
! INFILTRATION MODULE PARAMETERS  
!===============================================================================

! Basic infiltration control
PLAST_INFIL = T                        ! Enable infiltration module

! Sediment layer configuration
INFILTRATION_LAYERS = 10               ! Number of vertical sediment layers
INFILTRATION_MAX_DEPTH = 1.0           ! Maximum infiltration depth [m]
BOUNDARY_CONDITION_MODE = 1            ! 1=natural, 2=laboratory

! Global sediment properties (per-species parameters in plast_type)
SEDIMENT_GRAIN_SIZE = 0.0005           ! d_g sediment grain size [m]
SEDIMENT_POROSITY = 0.4                ! phi sediment porosity [-]
SEDIMENT_TORTUOSITY = 2.0              ! tau tortuosity factor [-]
MOLECULAR_DIFFUSIVITY = 1.0e-9         ! D_m molecular diffusivity [m^2/s]
MECHANICAL_DISPERSIVITY = 0.01         ! c mechanical diffusivity coefficient [-]

! Seepage flow (uniform or spatially-varying)
SEEPAGE_FLOW_TYPE = 'uniform'          ! 'uniform' or 'spatially_varying'
SEEPAGE_VELOCITY = -1.0e-6             ! Uniform Darcy velocity [m/s] (negative=downward)

! Species-specific parameters set per plastic type:
! A_Hamaker    : Set individually for each MP species [J]
! ws_sed       : Set individually for each MP species [m/s]  
! D_plast      : Set individually for each MP species [m^2/s]
```

#### 1.3 Setup_PLAST Extensions (PENDING IMPLEMENTATION)
Extend the existing setup routine to handle the streamlined infiltration arrays:

```fortran
! Add to Setup_PLAST subroutine after existing allocations
if (plast_infil) then
    ! Allocate sediment concentration arrays (species-specific)
    allocate(plast(iplast)%conc_sed(m, n_infil_layers))
    if (tox_vector) then
        allocate(plast(iplast)%tconc_b_sed(m, n_infil_layers))
        allocate(plast(iplast)%tconc_f_sed(m, n_infil_layers))
    endif
    
    ! Allocate infiltration flux arrays (nodal)
    allocate(plast(iplast)%iflx(m))
    if (tox_vector) then
        allocate(plast(iplast)%tiflx_b(m))
        allocate(plast(iplast)%tiflx_f(m))
        allocate(plast(iplast)%tiflx_f_deep(m))
    endif
    
    ! Initialize all arrays to zero
    plast(iplast)%conc_sed = 0.0
    plast(iplast)%iflx = 0.0
    if (tox_vector) then
        plast(iplast)%tconc_b_sed = 0.0
        plast(iplast)%tconc_f_sed = 0.0
        plast(iplast)%tiflx_b = 0.0
        plast(iplast)%tiflx_f = 0.0
        plast(iplast)%tiflx_f_deep = 0.0
    endif
    
    ! Initialize species-specific parameters (read from namelist)
    ! plast(iplast)%A_Hamaker already allocated as scalar
    ! plast(iplast)%ws_sed already allocated as scalar
    ! plast(iplast)%D_plast already allocated as scalar
endif
```

#### 1.4 Read_PLAST_Params Extensions (PENDING IMPLEMENTATION)
Add streamlined infiltration parameter reading to existing parameter input routine:

```fortran
! Add to Read_PLAST_Params subroutine
if (plast_infil) then
    ! Global infiltration parameters
    call FREAD(INFILTRATION_LAYERS, FILE_NAME,'INFILTRATION_LAYERS')
    call FREAD(INFILTRATION_MAX_DEPTH, FILE_NAME,'INFILTRATION_MAX_DEPTH') 
    call FREAD(BOUNDARY_CONDITION_MODE, FILE_NAME,'BOUNDARY_CONDITION_MODE')
    call FREAD(SEEPAGE_FLOW_TYPE, FILE_NAME,'SEEPAGE_FLOW_TYPE')
    call FREAD(SEEPAGE_VELOCITY, FILE_NAME,'SEEPAGE_VELOCITY')
    call FREAD(SEDIMENT_GRAIN_SIZE, FILE_NAME,'SEDIMENT_GRAIN_SIZE')
    call FREAD(SEDIMENT_POROSITY, FILE_NAME,'SEDIMENT_POROSITY')
    
    ! Species-specific parameters (read per plastic type)
    ! A_Hamaker, ws_sed, D_plast already in plast_type structure
    ! These will be read in the existing plastic species loop
endif
```

### Phase 2: Core Implementation

#### 2.1 New Subroutine: Calc_Infiltration
**Location**: Add after `Calc_Deposition` in mod_plast.F
**Purpose**: Calculate 1D infiltration into sediment column
**Key Features**:
- Diffusion-based or advection-based infiltration
- Consideration of sediment porosity
- Interaction with existing sediment layers
- Mass balance conservation

#### 2.2 Integration with Advance_plast
**Modification Point**: Line 815-818 in mod_plast.F
**Current Flow**:
```
!calculated depositional flux from settling (settle flux function is called here)
call calc_deposition
!calculate erosive flux (resuspension function is called here)  
call calc_erosion
```

**New Flow**:
```
!calculated depositional flux from settling
call calc_deposition
!calculate infiltration into sediment column
if (PLAST_BEDSET == 'burial_bottom') then
    call calc_infiltration
endif
!calculate erosive flux (resuspension function is called here)
call calc_erosion
```

#### 2.3 Modify Bottom Boundary Conditions
**Target**: Enhance bottom treatment in settling calculations
**Integration**: Update `Calc_Deposition` to pass settled mass to infiltration module

### Phase 3: Input/Output Enhancements

#### 3.1 NetCDF Output Extension
**File**: mod_ncdio.F
**Add Variables**:
- `mp_infil_depth`: Current infiltration depth
- `mp_infil_mass_profile`: Mass profile in sediment layers
- `mp_infil_rate`: Infiltration rate field

#### 3.2 Initialization Support
**File**: init_plast.F
**Enhancements**:
- Initialize infiltration arrays
- Set up sediment layer structure
- Configure infiltration parameters

### Phase 4: Advanced Features

#### 4.1 Toxic Infiltration
**Extension**: Apply similar infiltration to bonded toxics
**Consideration**: Different infiltration rates for toxics vs. MPs
**Integration**: Extend toxic-vector functionality

#### 4.2 Sediment Interaction
**Coupling**: Interface with mod_sed.F for sediment dynamics
**Features**: 
- Sediment layer evolution affects infiltration
- MP infiltration affects sediment properties

## Implementation Strategy (UPDATED)

### Phase 1: Data Structure Setup (COMPLETED ✓)
- ✓ Extended plast_type with species-specific infiltration variables
- ✓ Added infiltration flux tracking (iflx, tiflx_b, tiflx_f, tiflx_f_deep)
- ✓ Added sediment concentration arrays (conc_sed, tconc_b_sed, tconc_f_sed)
- ✓ Added DLVO physics parameters per species (A_Hamaker, ws_sed, D_plast)
- ✓ Added global infiltration switch (plast_infil)

### Phase 2: Namelist and Initialization (NEXT PRIORITY)
- Add infiltration parameters to generic_plastic.inp
- Extend Setup_PLAST routine for array allocation
- Update Read_PLAST_Params for parameter reading
- Extend init_plast.F for sediment layer initialization

### Phase 3: Core Physics Implementation
- Implement Calc_Infiltration subroutine based on governing equations
- Add 1D transport solver for movable/blocked phases
- Implement DLVO retention mechanisms using species-specific parameters
- Ensure mass conservation across water column-sediment interface

### Phase 4: Integration and Testing
- Modify Advance_plast to call infiltration when plast_infil=T
- Test with burial_bottom boundary treatment
- Validate mass balance and flux conservation
- Test species-specific parameterization

### Phase 5: I/O and Advanced Features
- Extend NetCDF output for sediment concentration profiles
- Add infiltration flux visualization
- Implement toxic leaching mechanisms
- Performance optimization and parallel compatibility

## Detailed Infiltration Physics Model

### Seepage Flow Input
Two implementation options for vertical seepage velocity (Darcy velocity):
1. **Uniform Option**: Single `w_D` value for entire domain via namelist parameter
2. **Spatial Varying**: Node-specific `w_D` values read from input file

### Governing Equations for Microplastic Infiltration

#### Primary MP Transport Equation
$$\frac{\partial p_{MP}}{\partial t} + (w_D + w_{ss}) \frac{\partial p_{MP}}{\partial z} = \frac{\partial}{\partial z}\left[D_L \frac{\partial p_{MP}}{\partial z}\right] - \lambda|w_D + w_{ss}|p_{MP} - S_{leach}^p$$

#### Blocked MP Equation  
$$\frac{\partial s_{MP}}{\partial t} = \lambda|w_D + w_{ss}|p_{MP} - S_{leach}^s$$

#### Total Concentration
$$c_{MP} = p_{MP} + s_{MP}$$

Where:
- $p_{MP}$: Movable MP concentration in sediment
- $s_{MP}$: Blocked MP concentration (immobilized by solid matrix)
- $w_D$: Darcy velocity (seepage flow, positive/negative)
- $w_{ss}$: MP settling velocity within solid matrix
- $z$: Vertical coordinate (sigma coordinate system)
- $\lambda$: Effective percolation blockage rate
- $D_L$: Mechanical seepage diffusivity

### Mechanical Diffusivity Model
$$D_L = c \cdot d_g |w_D + w_{ss}| + \frac{D_m}{\tau}$$

Parameters:
- $c$: Solid-matrix mechanical diffusivity coefficient
- $d_g$: Solid-matrix grain size
- $D_m$: Molecular diffusivity  
- $\tau$: Tortuosity in porous media

### Percolation Blockage Rate
$$\lambda = \lambda_{DLVO} + \lambda_{strain}$$

#### DLVO Theory Component
$$\lambda_{DLVO} = \frac{3(1-\phi)}{2d_g}[\eta_D + \eta_I + \eta_G]\alpha_{primary}$$

Where $\alpha_{primary} = 1$ when Hamaker constant $A \geq 0$, otherwise $\alpha_{primary} = 0$.

**Retention Efficiency Components:**

1. **Diffusion-induced retention** (van der Waals + Brownian motion):
   $$\eta_D = 2.4 A_S^{1/3} N_R^{-0.081} N_{Pe}^{-0.715} N_{vdW}^{0.052}$$

2. **Mechanical interception retention**:
   $$\eta_I = 0.55 A_S N_R^{1.675} N_A^{0.125}$$

3. **Gravitational retention** (buoyancy-induced trajectory changes):
   $$\eta_G = 0.22 N_R^{-0.24} N_G^{1.11} N_{vdW}^{0.053}$$

### Dimensionless Parameters

#### Porosity-dependent parameter:
$$A_S = \frac{2(1-\gamma^5)}{2-3\gamma+3\gamma^5-2\gamma^6}$$
$$\gamma = (1-\phi)^{1/3}$$

#### Size ratio:
$$N_R = \frac{d_{MP}}{d_g}$$

#### Peclet number (convection vs diffusion):
$$N_{Pe} = \frac{|w_D|d_g}{D_L}$$

#### Van der Waals number:
$$N_{vdW} = \frac{A}{kT}$$

#### Attraction number:
$$N_A = \frac{A}{12\pi\mu(d_{MP}/2)^2|w_D|}$$

#### Gravity number:
$$N_G = \frac{|w_{ss}|}{|w_D|}$$

Where:
- $\phi$: Porosity
- $d_{MP}$: Microplastic particle diameter
- $A$: Hamaker constant
- $k$: Boltzmann constant
- $T$: Temperature
- $\mu$: Dynamic viscosity of water

#### Strain-induced Blockage Component
$$\lambda_{strain} = \frac{k_{str}\psi_{str}}{|w_D + w_{ss}|}$$

$$k_{str} = 269.7\left(\frac{d_{MP}}{d_g}\right)^{1.42}$$

$$\psi_{str} = \left(\frac{d_{MP}+z}{d_{MP}}\right)^{-0.43}$$

Where $z$ is measured from sediment-water interface to the layer considered.

## Embedded Toxic Transport Equations

### Governing Equations for Embedded Toxics
Embedded toxics are carried by microplastics and follow identical transport dynamics with the same parameters and leaching kinetics.

#### Primary Embedded Toxic Transport Equation
$$\frac{\partial p_T^e}{\partial t} + (w_D + w_{ss}) \frac{\partial p_T^e}{\partial z} = \frac{\partial}{\partial z}\left[D_L \frac{\partial p_T^e}{\partial z}\right] - \lambda|w_D + w_{ss}|p_T^e - S_{leach}^p$$

#### Blocked Embedded Toxic Equation  
$$\frac{\partial s_T^e}{\partial t} = \lambda|w_D + w_{ss}|p_T^e - S_{leach}^s$$

#### Total Embedded Toxic Concentration
$$c_T^e = p_T^e + s_T^e$$

**Variable Definitions:**
- $p_T^e$: Movable embedded toxic concentration in sediment
- $s_T^e$: Blocked embedded toxic concentration (immobilized by solid matrix)
- $c_T^e$: Total embedded toxic concentration in sediment

**Key Features:**
- **All transport parameters** ($w_D$, $w_{ss}$, $D_L$, $\lambda$) **identical to MP** since toxics are carried by MPs
- **Leaching kinetics** ($S_{leach}^p$, $S_{leach}^s$) same as MP system
- **Same retention mechanisms** apply (DLVO + strain-induced)
- **Coupled with MP transport** through identical physical processes

## Free Toxic Transport Equations

### Governing Equation for Free Toxics
Free toxics are not subject to blocking mechanisms and follow simplified transport dynamics.

#### Free Toxic Transport Equation
$$\frac{\partial c_T^f}{\partial t} + w_D \frac{\partial c_T^f}{\partial z} = \frac{\partial}{\partial z}\left[D_L \frac{\partial c_T^f}{\partial z}\right] + S_{leach}$$

**Variable Definitions:**
- $c_T^f$: Mass concentration of free toxics in sediment
- $w_D$: Darcy velocity (seepage flow, identical to MP/embedded toxic equations)
- $D_L$: Mechanical seepage diffusivity (identical to MP/embedded toxic equations)
- $S_{leach}$: Combined leaching source from blocked and movable phases

#### Leaching Source Term
The leaching source combines toxin release from both blocked and movable MP phases:
$$S_{leach} = S_{leach}^p + S_{leach}^s$$

Where:
- $S_{leach}^p$: Leaching from movable MP phase
- $S_{leach}^s$: Leaching from blocked MP phase

**Key Features:**
- **No settling velocity** ($w_{ss} = 0$) - free toxics don't settle independently
- **No blocking mechanisms** - free toxics are not retained by DLVO or strain effects
- **Pure advection-diffusion** transport with seepage flow
- **Source-driven dynamics** from MP leaching processes
- **Same physical diffusivity** ($D_L$) as other phases

## Boundary Conditions

The model supports two distinct simulation modes with different boundary condition treatments.

### Mode 1: Natural Penetration Mode
Simulates natural infiltration of MP and toxics from surface water into sediment column.

#### Upper Boundary Conditions (z = 0, sediment-water interface)

**Movable MP and Embedded Toxics** - Robin boundary condition:
$$w_D c_{MP}|_{bot} = \left[-D_L \frac{\partial p_{MP}}{\partial z} + w_D p_{MP}\right]|_{z=0}$$

$$w_D c_T^e|_{bot} = \left[-D_L \frac{\partial p_T^e}{\partial z} + w_D p_T^e\right]|_{z=0}$$

**Free Toxics** - Dirichlet boundary condition:
$$c_T^f|_{z=0} = c_T^f|_{bot}$$

**Blocked MP and Embedded Toxics** - Neumann boundary condition:
$$\frac{\partial s_{MP}}{\partial z}\bigg|_{z=0} = 0$$

$$\frac{\partial s_T^e}{\partial z}\bigg|_{z=0} = 0$$

#### Lower Boundary Conditions (z = D_s, bottom outlet)

**All tracers** - Neumann boundary condition:
$$\frac{\partial \phi}{\partial z}\bigg|_{z=D_s} = 0 \quad \text{for all tracers } \phi$$

### Mode 2: Laboratory Leaching Experiment Mode
Simulates controlled lab experiments with water inflow and toxic outflow, with MPs pre-distributed in column.

#### Upper Boundary Conditions (z = 0, inlet)

**Free Toxics** - Dirichlet boundary condition:
$$c_T^f|_{z=0} = 0$$

**All other tracers** - Neumann boundary condition:
$$\frac{\partial \phi}{\partial z}\bigg|_{z=0} = 0 \quad \text{for } p_{MP}, s_{MP}, p_T^e, s_T^e$$

#### Lower Boundary Conditions (z = D_s, outlet)

**Movable MP and Embedded Toxics** - Robin boundary condition (no mass flux):
$$\left[-D_L \frac{\partial p_{MP}}{\partial z} + w_D p_{MP}\right]|_{z=D_s} = 0$$

$$\left[-D_L \frac{\partial p_T^e}{\partial z} + w_D p_T^e\right]|_{z=D_s} = 0$$

**Free Toxics** - Robin boundary condition (free outflow):
$$w_D c_T^f|_{z=D_s} = \left[-D_L \frac{\partial c_T^f}{\partial z} + w_D c_T^f\right]|_{z=D_s}$$

**Blocked MP and Embedded Toxics** - Neumann boundary condition:
$$\frac{\partial s_{MP}}{\partial z}\bigg|_{z=D_s} = 0$$

$$\frac{\partial s_T^e}{\partial z}\bigg|_{z=D_s} = 0$$

### Variable Definitions for Boundary Conditions
- $c_{MP}|_{bot}$: MP concentration in bottom layer of surface water
- $c_T^e|_{bot}$: Embedded toxic concentration in bottom layer of surface water  
- $c_T^f|_{bot}$: Free toxic concentration in bottom layer of surface water
- $z = 0$: Sediment-water interface (uppermost sediment layer)
- $D_s$: Total thickness of sediment matrix layer
- Robin conditions: Relate flux to concentration gradients
- Dirichlet conditions: Specify concentration values
- Neumann conditions: Specify zero gradient (no flux)

## Initial Conditions Framework

### Initialization Options
The model supports multiple initialization approaches for sediment column concentrations:

#### 1. Uniform Concentration Initialization
Read from `generic_plastic.inp` parameters:
```
! Sediment column initial concentrations (uniform)
INIT_MP_CONC_SEDIMENT = 0.0001        ! kg/m^3
INIT_ETOX_CONC_SEDIMENT = 0.00001     ! kg/m^3  
INIT_FTOX_CONC_SEDIMENT = 0.0         ! kg/m^3
INIT_SMASS_DISTRIBUTION = 0.5         ! fraction in blocked phase
```

#### 2. NetCDF File Initialization
Read spatially-varying initial conditions from NetCDF files via namelist:
```
&NML_ADDITIONAL_MODELS
INIT_INFILTRATION_TYPE = 'netcdf'
INIT_INFILTRATION_FILE = 'initial_sediment_profiles.nc'
```

**NetCDF File Structure:**
- `mp_sediment_conc(node, layer, nplast)`: Movable MP concentrations
- `mp_sediment_blocked(node, layer, nplast)`: Blocked MP concentrations  
- `etox_sediment_conc(node, layer, nplast)`: Embedded toxic concentrations
- `etox_sediment_blocked(node, layer, nplast)`: Blocked embedded toxics
- `ftox_sediment_conc(node, layer)`: Free toxic concentrations

#### 3. Implementation in init_plast.F
**Current Structure** (lines 69-97):
- `LOAD_INIT_PLAST_MASS()`: Bottom mass initialization
- `LOAD_INIT_TOXIN_MASS()`: Toxic mass initialization

**Extensions Needed**:
- `LOAD_INIT_SEDIMENT_PROFILES()`: 3D concentration profiles
- Layer-by-layer initialization for each tracer type
- Boundary condition coupling with water column

## Numerical Treatment Framework

### Existing Numerical Methods (Reference: mod_plast.F)
The infiltration module will leverage existing proven numerical schemes:

#### 1. Flux-Limited Settling Algorithm
**Reference**: `Settle_Flux_new()` subroutine (mod_plast.F:1829-1982)
- **CFL-limited time stepping**: `DTmax = Settle_CFL*minval(dx/wset)`  
- **Multi-cycle approach**: `mcyc = int(DTplast/DTmax + 1)`
- **Flux limiters**: Support for upwind, minmod, Van Leer limiters
- **Bottom boundary handling**: Different treatments for boundary conditions

**Adaptation for Infiltration**:
- Replace settling velocity `wset` with combined seepage-settling: `w_D + w_ss`
- Modify boundary conditions for sediment interface
- Add blocking term integration: `-λ|w_D + w_ss|p_MP`

#### 2. Vertical Diffusion Treatment  
**Reference**: `vdif_scal()` calls (mod_plast.F:863, 907)
- **Implicit scheme**: Handles stiff diffusion problems
- **Tridiagonal solver**: Efficient 1D solution
- **Multiple tracers**: Supports different diffusivity values

**Adaptation for Infiltration**:
- Use mechanical diffusivity `D_L` instead of turbulent diffusivity
- Apply to all infiltration tracers: `p_MP`, `s_MP`, `p_T^e`, `s_T^e`, `c_T^f`

#### 3. Time Integration Strategy
**Current Approach** (mod_plast.F:749-950):
- **Operator splitting**: Advection → Diffusion → Reactions
- **Multi-step updates**: `cnew` → `conc` progression  
- **Source term integration**: Explicit treatment

**Infiltration Integration**:
```fortran
! Infiltration time-stepping (pseudo-code)
call calc_infiltration_advection(p_MP, w_D + w_ss, DT_infil)
call calc_infiltration_diffusion(p_MP, D_L, DT_infil)  
call calc_blocking_reactions(p_MP, s_MP, lambda, DT_infil)
call calc_toxic_leaching(p_T_e, s_T_e, c_T_f, S_leach, DT_infil)
```

### Proposed Infiltration Subroutines

#### 1. Main Integration: `Calc_Infiltration()`
**Location**: Add after `Calc_Deposition` in mod_plast.F
**Function**: Coordinate all infiltration processes
**Structure**:
- Loop over nodes and plastic classes
- Apply boundary conditions (Mode 1 vs Mode 2)
- Call individual process subroutines
- Update concentration arrays

#### 2. Transport Solver: `Infiltrate_1D_Transport()`
**Function**: Solve 1D advection-diffusion-reaction system
**Method**: Based on existing `Settle_Flux_new()` framework
**Features**:
- CFL-limited time stepping
- Flux limiters for stability
- Robin/Dirichlet/Neumann boundary conditions

#### 3. Blocking Kinetics: `Calc_Retention_Kinetics()`
**Function**: Handle DLVO + strain blocking mechanisms  
**Method**: Explicit integration of blocking rates
**Components**:
- Calculate all dimensionless parameters
- Compute retention efficiencies
- Apply blocking/release rates

#### 4. Boundary Coupling: `Apply_Infiltration_BC()`
**Function**: Interface with water column concentrations
**Features**:
- Mode selection (natural vs laboratory)
- Robin condition implementation
- Mass balance enforcement

## Key Considerations

### Scientific Aspects
1. **Multi-phase Transport**: Coupled movable/blocked MP concentrations
2. **Physical Blocking Mechanisms**: DLVO theory + strain-induced retention  
3. **Seepage-settling Coupling**: Combined Darcy flow and gravitational settling
4. **Mass Conservation**: Total MP conservation across phase transitions

### Technical Aspects
1. **Memory Management**: New arrays for sediment layer profiles
2. **Performance**: Minimize computational overhead
3. **Parallelization**: Ensure new code works with MPI parallelization
4. **Backward Compatibility**: Don't break existing functionality

### Testing Strategy
1. **Unit Tests**: Test infiltration algorithms independently
2. **Integration Tests**: Test with existing FVCOM-MP functionality
3. **Validation**: Compare against analytical solutions where possible
4. **Performance Tests**: Ensure acceptable computational cost

## File Modification Summary

### Primary Files to Modify:
1. **mod_plast.F**: Core infiltration implementation
2. **generic_plastic.inp**: New parameters
3. **mod_ncdio.F**: Output variables
4. **init_plast.F**: Initialization support

### Reference Files:
1. **mod_sed.F**: For sediment layer structure reference
2. **mod_main.F**: For global variable patterns
3. **waterPACT_run.nml**: For testing configuration

## Summary of User's Streamlined Implementation

### Key Advantages of the Implemented Approach

1. **Species-Specific Focus**: Each microplastic species (iplast=1,Nplast) has individual DLVO parameters (A_Hamaker, ws_sed, D_plast), enabling realistic multi-species infiltration modeling.

2. **Memory Efficient Design**: Streamlined variable set with only essential arrays (conc_sed, iflx, tiflx_*) rather than comprehensive physics parameters, reducing computational overhead.

3. **I/O Ready Structure**: All sediment concentration arrays designed for NetCDF output capability, enabling direct visualization and analysis of infiltration profiles.

4. **Modular Control**: Simple global switch (plast_infil) for easy enable/disable of infiltration functionality without affecting existing FVCOM-MP operations.

5. **Flux Conservation**: Clear separation of depositional flux (dflx) and infiltration flux (iflx) ensures proper mass balance tracking across water column-sediment interface.

6. **Toxic Integration**: Complete toxic-vector support with separate arrays for embedded (tiflx_b) and free (tiflx_f, tiflx_f_deep) toxic infiltration.

### Next Development Priority

Based on the completed data structure implementation, the immediate next steps are:
1. Namelist parameter integration in generic_plastic.inp
2. Array allocation in Setup_PLAST routine
3. Core Calc_Infiltration subroutine development
4. Integration with Advance_plast main loop

This streamlined approach provides a solid foundation for implementing the complete infiltration physics while maintaining code efficiency and scientific rigor.
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

## Implementation Strategy (UPDATED)

### Phase 1: Data Structure Setup
- ✓ Extended plast_type with species-specific infiltration variables
- ✓ Added infiltration flux tracking (iflx, tiflx_b, tiflx_f, tiflx_f_deep)
- ✓ Added sediment concentration arrays (conc_sed, tconc_b_sed, tconc_f_sed)
- ✓ Added DLVO physics parameters per species (A_Hamaker, ws_sed, D_plast)
- ✓ Added global infiltration switch (plast_infil)

### Phase 2: Namelist and Parameter Implementation (NEXT PRIORITY)
- ✓ Add infiltration parameters to generic_plastic.inp following existing patterns
- ✓ Extend Setup_PLAST routine for array allocation within species loop
- ✓ Update Read_PLAST_Params using Get_Val calls consistent with existing code
- ✓ Add global variable declarations to mod_plast.F public section

### Phase 3: Core Physics Implementation
- ✓ Implement Calc_Infiltration subroutine based on governing equations
- ✓ Implement a new vertical diffusion scalar solver
- ✓ Adapt 1D transport solver for movable/blocked phases
- ✓ Implement DLVO+strain retention mechanisms using species-specific parameters
- ✓ Modify Advance_plast to call infiltration when plast_infil=T
- ✓ Implement toxic leaching mechanisms for deep sediment layers
- ✓ Remove burial_bottom boundary treatment
- ✓ Ensure mass conservation across water column-sediment interface at bottom update

### Phase 4: I/O and Advanced Features
- ✓ Extend initialization file input to consider new bed properties
- ✓ Extend initialization file input to initialize plastic/toxin inside bed layers
- Extend NetCDF output for concentration profiles (conc_sed, tconc_b_sed, tconc_f_sed)
# FVCOM-MP Infiltration Development Notes

## Project Overview
Developing 1D infiltration functionality for FVCOM-MP (microplastic transport and fate module) to model the infiltration of bottom-settling microplastics and toxics into sediment columns.

## Code Structure Analysis Completed (2025-08-20)

### Primary Module: mod_plast.F 
- **Size**: 37,333 tokens (large file requiring section-by-section reading)
- **Key Features**: Microplastic transport, toxic-vector functionality, settling dynamics
- **Current Bottom Treatment**: 5 modes including `burial_bottom` (target for infiltration)
- **Main Subroutines**:
  - `Setup_PLAST`: Parameter initialization  
  - `Advance_plast`: Main time-stepping (line 749)
  - `Calc_Deposition`: Current settling implementation (line 1674)
  - `Calc_Wset_Plastic`: Settling velocity calculations (line 1767)

### Key Dependencies
- **Core**: `mod_par.F`, `mod_prec.F`, `mod_types.F`, `mod_wd.F`
- **Sediment Interaction**: `mod_sed.F`, `mod_bbl.F`, `mod_cstms_vars.F`
- **I/O**: `mod_ncdio.F`, `mod_input.F`
- **Main Program**: `fvcom.F`, `internal_step.F`

### Configuration Files
- **Main Namelist**: `FVCOM_MP_run/waterPACT_run.nml`
  - Plastic model enabled: `PLASTIC_MODEL = T`
  - Config file: `PLASTIC_MODEL_FILE = 'generic_plastic.inp'`
- **Plastic Parameters**: `FVCOM_MP_run/generic_plastic.inp`
  - Current setting: `PLAST_BEDSET = concrete_bottom`
  - Target upgrade: `PLAST_BEDSET = burial_bottom`

### C Preprocessor Integration
- **Primary Directive**: `#if defined (PLASTIC)`
- **Compilation**: `.F` files processed through CPP before Fortran compilation
- **Conditional Features**: WET_DRY handling, MULTIPROCESSOR support

### Current Settling Implementation
- **Location**: `Calc_Deposition` subroutine (mod_plast.F:1674-1760)
- **Process**: Flux-limited settling with CFL condition
- **Key Variables**:
  - `dflx`: depositional flux [kg m^-2]
  - `mass`: bottom plastic mass density [kg m^-2] 
  - `wset`: settling velocity [m/s]

## Implementation Plan
Comprehensive plan documented in `infiltration_plan.md` covering:

### Phase 1: Data Structure Extension
- Add infiltration variables to `plast_type`
- Extend namelist parameters
- Initialize sediment layer arrays

### Phase 2: Core Algorithm
- Implement `Calc_Infiltration` subroutine
- Add 1D infiltration physics (diffusion/advection-based)
- Integrate with existing `Advance_plast` routine

### Phase 3: I/O Enhancement  
- Extend NetCDF output for infiltration profiles
- Add visualization support
- Update initialization routines

### Phase 4: Advanced Features
- Toxic infiltration coupling
- Sediment-MP interaction
- Performance optimization

## Next Steps
1. **Implement Data Structures**: Extend `plast_type` with infiltration arrays
2. **Develop Core Algorithm**: Create `Calc_Infiltration` subroutine with mass-conserving 1D infiltration
3. **Integration Testing**: Modify `Advance_plast` to call infiltration functionality
4. **Validation**: Test against simplified analytical cases

## File Locations Reference
- **Main Module**: `mod_plast.F` (37,333 tokens)
- **Configuration**: `FVCOM_MP_run/generic_plastic.inp`
- **Namelist**: `FVCOM_MP_run/waterPACT_run.nml`
- **Plan Document**: `infiltration_plan.md`
- **Makefile**: `FVCOM_MP_run/makefile`
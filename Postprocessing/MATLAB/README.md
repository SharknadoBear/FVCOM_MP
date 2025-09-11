# FVCOM Plastic Initialization Tools

This directory contains MATLAB tools for creating NetCDF initialization files for the FVCOM plastic transport model.

## 📁 Files Overview

### Main Functions
- **`create_init_plast_conc.m`** - Creates water column plastic concentration initialization files
- **`create_init_plast_sed_conc.m`** - Creates sediment plastic concentration initialization files

### Example Scripts
- **`example_create_plastic_water_init.m`** - Examples for water column initialization
- **`example_create_plastic_sed_init.m`** - Examples for sediment initialization (if available)

### Testing
- **`test_plastic_init_tools.m`** - Comprehensive test suite for both tools

## 🚀 Quick Start

### 1. Water Column Initialization

```matlab
% Basic usage
create_init_plast_conc('water_init.nc');

% With custom parameters
create_init_plast_conc('water_init.nc', ...
    'num_nodes', 1500, ...              % Match your mesh
    'num_sigma_layers', 20, ...         % Match KB in FVCOM
    'num_plastic_types', 3, ...         % Match NPLAST
    'init_conc', 0.05);                 % kg/m³
```

### 2. Sediment Initialization

```matlab
% Basic usage
create_init_plast_sed_conc('sediment_init.nc');

% With custom parameters
create_init_plast_sed_conc('sediment_init.nc', ...
    'num_nodes', 1500, ...              % Match your mesh
    'num_sed_layers', 5, ...            % Match N_sed_layer
    'num_plastic_types', 3, ...         % Match NPLAST
    'conc_total', 0.02);                % kg/m³
```

## 📊 File Structures

### Water Column NetCDF Structure
```
Dimensions:
  node     = N_nodes      (mesh nodes)
  ksl      = N_sigma      (sigma levels, KB)
  nplast   = N_types      (plastic classes)

Variables:
  init_conc(node, ksl, nplast) - Initial plastic concentration [kg/m³]
```

### Sediment NetCDF Structure
```
Dimensions:
  node         = N_nodes      (mesh nodes)
  N_sed_layer  = N_layers     (sediment layers)
  nplast       = N_types      (plastic classes)

Variables:
  init_conc_sed(node, N_sed_layer, nplast)   - Total plastic in sediment [kg/m³]
  init_conc_sed_s(node, N_sed_layer, nplast) - Stopped plastic [kg/m³]
  init_conc_sed_p(node, N_sed_layer, nplast) - Passing plastic [kg/m³]
```

## 🔧 Parameters

### Common Parameters
- **`num_nodes`** - Number of mesh nodes (must match your FVCOM mesh)
- **`num_plastic_types`** - Number of plastic classes (must match NPLAST in FVCOM)
- **`plastic_names`** - Cell array of plastic type names (optional)
- **`units`** - Concentration units (default: 'kg/m^3')
- **`verbose`** - Display progress messages (default: true)

### Water Column Specific
- **`num_sigma_layers`** - Number of sigma layers (must match KB in FVCOM)
- **`init_conc`** - Initial concentration (scalar or array)
- **`depth_profile`** - Vertical distribution: 'uniform', 'surface', 'exponential'
- **`surface_factor`** - Surface enhancement factor
- **`decay_scale`** - Decay scale for exponential profiles

### Sediment Specific
- **`num_sed_layers`** - Number of sediment layers (must match N_sed_layer)
- **`conc_total`** - Total plastic concentration (scalar or array)
- **`conc_stopped`** - Stopped plastic concentration (scalar or array)
- **`conc_passing`** - Passing plastic concentration (scalar or array)

## 📈 Vertical Profiles (Water Column)

### Uniform Distribution
```matlab
'depth_profile', 'uniform'
```
Equal concentration at all depths.

### Surface Enhanced
```matlab
'depth_profile', 'surface', 'surface_factor', 3.0
```
Linear decrease from surface with specified enhancement factor.

### Exponential Decay
```matlab
'depth_profile', 'exponential', 'decay_scale', 0.2
```
Exponential decay from surface (typical for buoyant plastics).

### Custom Arrays
```matlab
% Create your own 3D concentration array
conc_array = your_concentration_data;  % Size: [nodes × layers × types]
create_init_plast_conc('file.nc', 'init_conc', conc_array);
```

## 💡 Usage Examples

### Example 1: Basic Regional Model
```matlab
% Set up for a regional model with 2500 nodes
create_init_plast_conc('regional_water.nc', ...
    'num_nodes', 2500, ...
    'num_sigma_layers', 15, ...
    'num_plastic_types', 4, ...
    'init_conc', 0.03, ...
    'plastic_names', {'Microplastics', 'Bottles', 'Bags', 'Fragments'});

create_init_plast_sed_conc('regional_sediment.nc', ...
    'num_nodes', 2500, ...
    'num_sed_layers', 3, ...
    'num_plastic_types', 4, ...
    'conc_total', 0.01, ...
    'plastic_names', {'Microplastics', 'Bottles', 'Bags', 'Fragments'});
```

### Example 2: Pollution Hotspot
```matlab
% Surface-enhanced distribution near pollution sources
create_init_plast_conc('hotspot_water.nc', ...
    'num_nodes', 1000, ...
    'init_conc', 0.1, ...
    'depth_profile', 'surface', ...
    'surface_factor', 5.0);
```

### Example 3: River Discharge Scenario
```matlab
% Custom spatial distribution
nodes = 500; layers = 12; types = 2;
conc_array = create_river_plume_distribution(nodes, layers, types);
create_init_plast_conc('river_discharge.nc', ...
    'num_nodes', nodes, 'init_conc', conc_array);
```

## 🔗 FVCOM Integration

### 1. Set Input Files
In your FVCOM namelist, specify:
```fortran
&NML_PLASTIC
  INIT_PLAST_CONC_FILE     = 'water_init.nc'
  INIT_PLAST_SED_CONC_FILE = 'sediment_init.nc'
  ...
/
```

### 2. Ensure Dimension Consistency
- **Nodes**: Must match your mesh file
- **Sigma layers**: Must match KB parameter
- **Sediment layers**: Must match N_sed_layer parameter  
- **Plastic types**: Must match NPLAST parameter

### 3. Verify Files
Use `ncdump -h filename.nc` to inspect file structure before running FVCOM.

## 🧪 Testing

Run the test suite to verify everything works:
```matlab
test_plastic_init_tools
```

## 📋 Requirements

- MATLAB with NetCDF support
- Statistics and Signal Processing Toolbox (for some advanced examples)

## 🐛 Troubleshooting

### Common Issues

1. **Dimension mismatch errors**
   - Ensure `num_nodes` matches your FVCOM mesh
   - Ensure `num_sigma_layers` matches KB in your model
   - Ensure `num_plastic_types` matches NPLAST

2. **NetCDF creation fails**
   - Check write permissions in target directory
   - Ensure sufficient disk space
   - Try with smaller test arrays first

3. **FVCOM cannot read files**
   - Verify file structure with `ncdump -h filename.nc`
   - Check that variable names match expected: `init_conc`, `init_conc_sed`, etc.
   - Ensure dimensions are named correctly: `node`, `ksl`/`N_sed_layer`, `nplast`

4. **Mass balance warnings (sediment)**
   - This is normal when total ≠ stopped + passing
   - The tools automatically enforce mass balance
   - Check the verbose output to see what adjustments were made

### Getting Help

1. Run the example scripts to see working configurations
2. Check the test suite for validation examples  
3. Use `help function_name` in MATLAB for detailed parameter descriptions

## 📝 Notes

- All concentration values are in kg/m³ by default
- Sigma coordinates: k=1 (surface) to k=KB (bottom)
- Sediment layers: k=1 (top) to k=N_sed_layer (bottom)
- Mass balance is automatically enforced for sediment concentrations
- Files are CF-1.6 compliant NetCDF format

## 🔄 Version History

- v1.0 - Initial release with water column and sediment initialization tools
- Compatible with FVCOM plastic transport model subroutines:
  - `LOAD_INIT_PLAST_CONC` (water column)
  - `LOAD_INIT_PLAST_SED_CONC` (sediment)

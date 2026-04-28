# waterPACT Experiment b NetCDF Structure

Representative source on Kestrel:
`/kfs3/scratch/yhuang168/waterPACT_MP_floc/OUTPUT_b/waterPACT_b_0001.nc`

This note records the NetCDF structure reported for experiment `b`, the
floc-enabled reference case. It is intended as a practical reference for
writing post-processing and mass-budget diagnostic scripts.

## File Metadata

| Item | Value |
| --- | --- |
| Format | 64-bit NetCDF |
| title | `waterPACT` |
| institution | `School for Marine Science and Technology` |
| source | `FVCOM_4.3.1` |
| history | `model started at: 27/04/2026   14:18` |
| Conventions | `CF-1.0` |
| CoordinateSystem | `GeoReferenced` |
| CoordinateProjection | `none` |
| Tidal forcing | `run2019` |
| Rivers | 9 rivers: `DR_1` to `DR_6`, `SR_1` to `SR_3` |
| Groundwater forcing | off |
| Surface heat forcing | off |
| Surface wind forcing | `waterPACT_DRE_wind_2019.nc`, WRF structured grid, starts `2019-01-01 00:00:00` |
| Surface precipitation/evaporation forcing | off |

## Dimensions

| Dimension | Size | Meaning |
| --- | ---: | --- |
| `nele` | 126722 | triangular elements / faces |
| `node` | 68416 | mesh nodes |
| `siglay` | 30 | sigma layer centers |
| `siglev` | 31 | sigma layer interfaces |
| `three` | 3 | triangle vertex count |
| `time` | 160 | unlimited output time records in this stack |
| `maxnode` | 11 | max surrounding nodes per node |
| `maxelem` | 9 | max surrounding elements per node |
| `four` | 4 | auxiliary grid metric dimension |

## Coordinate and Mesh Variables

| Variable | Dimensions | Type | Notes |
| --- | --- | --- | --- |
| `nprocs` | scalar | int32 | number of processors |
| `partition` | `nele` | int32 | element partition |
| `x`, `y` | `node` | single | nodal Cartesian coordinates, meters |
| `lon`, `lat` | `node` | single | nodal longitude/latitude, degrees |
| `xc`, `yc` | `nele` | single | element-center Cartesian coordinates, meters |
| `lonc`, `latc` | `nele` | single | element-center longitude/latitude, degrees |
| `siglay` | `node,siglay` | single | node sigma layers, positive up |
| `siglev` | `node,siglev` | single | node sigma levels, positive up |
| `siglay_center` | `nele,siglay` | single | element-center sigma layers |
| `siglev_center` | `nele,siglev` | single | element-center sigma levels |
| `h` | `node` | single | nodal bathymetry, m, positive down |
| `h_center` | `nele` | single | element-center bathymetry, m, positive down |
| `nv` | `nele,three` | int32 | nodes surrounding each element |
| `nbe` | `nele,three` | int32 | elements surrounding each element |
| `ntsn` | `node` | int32 | number of nodes surrounding each node |
| `nbsn` | `node,maxnode` | int32 | nodes surrounding each node |
| `ntve` | `node` | int32 | number of elements surrounding each node |
| `nbve` | `node,maxelem` | int32 | elements surrounding each node |
| `a1u`, `a2u` | `nele,four` | single | element metric coefficients |
| `aw0`, `awx`, `awy` | `nele,three` | single | element interpolation/gradient weights |
| `art1` | `node` | single | node-based control volume area |
| `art2` | `node` | single | area of elements around each node |

## Time Variables

| Variable | Dimensions | Type | Units / Notes |
| --- | --- | --- | --- |
| `iint` | `time` | int32 | internal mode iteration number |
| `time` | `time` | single | days since `0.0`, time zone `none` |
| `Itime` | `time` | int32 | days since `0.0`, time zone `none` |
| `Itime2` | `time` | int32 | msec since `00:00:00`, time zone `none` |

## Hydrodynamic and Turbulence Variables

| Variable | Dimensions | Type | Units | Location |
| --- | --- | --- | --- | --- |
| `zeta` | `node,time` | single | meters | node |
| `u` | `nele,siglay,time` | single | meters s-1 | face |
| `v` | `nele,siglay,time` | single | meters s-1 | face |
| `tauc` | `nele,time` | single | m2 s-2 | face |
| `temp` | `node,siglay,time` | single | degrees C | node |
| `salinity` | `node,siglay,time` | single | 1e-3 | node |
| `viscofm` | `nele,siglay,time` | single | m2 s-1 | face |
| `viscofh` | `node,siglay,time` | single | m2 s-1 | node |
| `km` | `node,siglev,time` | single | m2 s-1 | node/interface |
| `kh` | `node,siglev,time` | single | m2 s-1 | node/interface |
| `kq` | `node,siglev,time` | single | m2 s-1 | node/interface |
| `q2` | `node,siglev,time` | single | m2 s-2 | node/interface |
| `q2l` | `node,siglev,time` | single | m3 s-2 | node/interface |
| `l` | `node,siglev,time` | single | m3 s-2 | node/interface |

## Wetting and Drying Variables

| Variable | Dimensions | Type | Location |
| --- | --- | --- | --- |
| `wet_nodes` | `node,time` | int32 | node |
| `wet_cells` | `nele,time` | int32 | face |
| `wet_nodes_prev_int` | `node,time` | int32 | node |
| `wet_cells_prev_int` | `nele,time` | int32 | face |
| `wet_cells_prev_ext` | `nele,time` | int32 | face |
| `inundation_cells` | `nele,time` | int32 | face |

## Sediment Variables

| Variable | Dimensions | Type | Units | Notes |
| --- | --- | --- | --- | --- |
| `coarse_sand_1` | `node,siglay,time` | single | g/L | suspended sediment concentration |
| `coarse_sand_1_bedfrac` | `node,time` | single | `-` | surface-layer bed fraction |

## Plastic and Floc Variables

These fields are the primary targets for later plastic mass-budget diagnostics.
For case `b`, `PLAST_FLOC = T`, so the split floc variables are present.

| Variable | Dimensions | Type | Units | Interpretation |
| --- | --- | --- | --- | --- |
| `taub_total` | `node,time` | single | m2 s-2 | total bed stress magnitude used to force sediment/plastic exchange |
| `mp1` | `node,siglay,time` | single | kg/m3 | total suspended plastic concentration |
| `mp1_agg` | `node,siglay,time` | single | kg/m3 | aggregated plastic concentration, `conc_a` |
| `mp1_dis` | `node,siglay,time` | single | kg/m3 | disaggregated plastic concentration, `conc_d` |
| `settle_vel_mp1` | `node,siglay,time` | single | m/s | free/disaggregated plastic settling velocity |
| `settle_vel_floc_mp1` | `node,siglay,time` | single | m/s | floc/aggregated plastic settling velocity |
| `lambda_a_mp1` | `node,siglay,time` | single | 1/s | aggregation rate |
| `lambda_d_mp1` | `node,siglay,time` | single | 1/s | disaggregation rate |
| `dep_flux_mp1` | `node,time` | single | kg/m2 | deposition mass density/flux output |
| `ero_flux_mp1` | `node,time` | single | kg/m2 | erosion mass density/flux output |
| `bot_mass_mp1` | `node,time` | single | kg/m2 | total bottom surface plastic mass density |
| `bot_mass_agg_mp1` | `node,time` | single | kg/m2 | aggregated bottom plastic mass density |
| `bot_mass_dis_mp1` | `node,time` | single | kg/m2 | disaggregated bottom plastic mass density |
| `tau_ce_mp1` | `node,time` | single | N m-2 | critical shear stress for plastic erosion |
| `frac_actv_mp1` | `node,time` | single | dimensionless | active-layer mass fraction |
| `actv_thck_mp1` | `node,time` | single | m | active-layer thickness |
| `ambient_rho` | `node,siglay,time` | single | kg/m3 | ambient water density |
| `bottom_stress` | `node,time` | single | kg m-1 s-1 | bottom current shear stress |

## Script-Writing Notes

- Use node-centered arrays for plastic and sediment mass budgets:
  `mp1`, `mp1_agg`, `mp1_dis`, `coarse_sand_1`, `bot_mass_mp1`,
  `bot_mass_agg_mp1`, and `bot_mass_dis_mp1`.
- Use `art1(node)` as the natural node-control-volume area for nodal
  inventory integrals.
- Use `h(node)`, `zeta(node,time)`, and `siglay/siglev` or layer thicknesses
  derived from sigma spacing to convert water-column concentrations to mass.
- For a layer-integrated plastic inventory, the likely form is
  `sum_node sum_layer C(node,layer,time) * volume(node,layer,time)`, with
  `volume` based on `art1`, total water depth, and sigma-layer thickness.
- Bed inventory should be area-integrated with `sum_node bot_mass(node,time) *
  art1(node)`.
- In floc-enabled outputs, `mp1` should be checked against
  `mp1_agg + mp1_dis`. This is the first low-cost consistency check before a
  full mass budget.
- For experiment `a`, floc-only variables such as `mp1_agg`, `mp1_dis`,
  `lambda_a_mp1`, and `bot_mass_agg_mp1` are absent because `PLAST_FLOC = F`.
  See `waterPACT_netcdf_structure_a.md` for the no-floc baseline structure.

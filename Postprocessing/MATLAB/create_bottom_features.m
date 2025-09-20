function create_bottom_features(output_filename, varargin)
% CREATE_BOTTOM_FEATURES - Creates a NetCDF file for initializing bottom
% sediment features for FVCOM plastic transport model
%
% SYNTAX:
%   create_bottom_features(output_filename)
%   create_bottom_features(output_filename, 'Parameter', Value, ...)
%
% INPUTS:
%   output_filename - String, name of the NetCDF file to create
%
% OPTIONAL PARAMETERS:
%   'num_nodes'    - Number of nodes in the mesh (default: 100)
%   'sd50'         - Sediment median grain size [m] - array or scalar (default: 1e-4)
%   'dens'         - Sediment bulk density [kg/m^3] - array or scalar (default: 2650)
%   'mudp'         - Mud percentage [0-1] - array or scalar (default: 0.1)
%   'poro'         - Sediment porosity [0-1] - array or scalar (default: 0.4)
%   'Tsed'         - Sediment temperature [°C] - array or scalar (default: 15.0)
%   'wDcy'         - Vertical Darcy velocity [m/s] - array or scalar (default: 1e-6)
%   'Lsed'         - Sediment layer thickness [m] - array or scalar (default: 0.1)
%   'spatial_variation' - Apply spatial variation ('uniform', 'random', 'gradient')
%   'variation_factor'  - Factor for spatial variation (default: 0.1)
%   'verbose'      - Display progress messages (default: true)
%
% EXAMPLES:
%   % Basic usage with default values
%   create_bottom_features('bottom_features.nc');
%
%   % Specify custom dimensions and uniform values
%   create_bottom_features('bottom_features.nc', ...
%       'num_nodes', 500, 'sd50', 2e-4, 'poro', 0.35, 'dens', 2700);
%
%   % Add spatial variation
%   create_bottom_features('bottom_features.nc', ...
%       'spatial_variation', 'random', 'variation_factor', 0.2, ...
%       'sd50', 1.5e-4, 'poro', 0.4);
%
%   % Use custom arrays
%   nodes = 100;
%   sd50_array = rand(nodes, 1) * 2e-4 + 1e-4;
%   poro_array = rand(nodes, 1) * 0.2 + 0.3;
%   create_bottom_features('bottom_features.nc', ...
%       'num_nodes', nodes, 'sd50', sd50_array, 'poro', poro_array);
%
% NOTES:
%   - Compatible with FVCOM's LOAD_BOTTOM_FEATURES subroutine
%   - All arrays should have dimension: (num_nodes,)
%   - Default values represent typical sediment properties
%   - Spatial variation options:
%     * 'uniform': No variation (default)
%     * 'random': Random variation around base values
%     * 'gradient': Linear gradient across domain
%
% VARIABLE DESCRIPTIONS:
%   sd50 : Sediment median grain size [m] (typical: 1e-5 to 1e-3)
%   dens : Sediment bulk density [kg/m^3] (typical: 2000-3000)
%   mudp : Mud percentage [0-1] (0=pure sand, 1=pure mud)
%   poro : Sediment porosity [0-1] (typical: 0.3-0.6)
%   Tsed : Sediment temperature [°C] (affects toxin kinetics)
%   wDcy : Vertical Darcy velocity [m/s] (pore water flow)
%   Lsed : Sediment layer thickness [m]
%
% Author: Generated for FVCOM plastic transport model
% Date: September 2025

% Parse input arguments
p = inputParser;
addRequired(p, 'output_filename', @ischar);
addParameter(p, 'num_nodes', 100, @(x) isnumeric(x) && x > 0);
addParameter(p, 'sd50', 1e-4, @isnumeric);
addParameter(p, 'dens', 2650, @isnumeric);
addParameter(p, 'mudp', 0.1, @isnumeric);
addParameter(p, 'poro', 0.4, @isnumeric);
addParameter(p, 'Tsed', 15.0, @isnumeric);
addParameter(p, 'wDcy', 1e-6, @isnumeric);
addParameter(p, 'Lsed', 0.1, @isnumeric);
addParameter(p, 'spatial_variation', 'uniform', @(x) ismember(x, {'uniform', 'random', 'gradient'}));
addParameter(p, 'variation_factor', 0.1, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'verbose', true, @islogical);

parse(p, output_filename, varargin{:});

% Extract parsed values
num_nodes = p.Results.num_nodes;
sd50_input = p.Results.sd50;
dens_input = p.Results.dens;
mudp_input = p.Results.mudp;
poro_input = p.Results.poro;
Tsed_input = p.Results.Tsed;
wDcy_input = p.Results.wDcy;
Lsed_input = p.Results.Lsed;
spatial_variation = p.Results.spatial_variation;
variation_factor = p.Results.variation_factor;
verbose = p.Results.verbose;

if verbose
    fprintf('Creating NetCDF file for FVCOM bottom features...\n');
    fprintf('  Output file: %s\n', output_filename);
    fprintf('  Number of nodes: %d\n', num_nodes);
    fprintf('  Spatial variation: %s\n', spatial_variation);
end

% Process each variable
target_dims = [num_nodes, 1];

sd50 = process_bottom_variable(sd50_input, target_dims, spatial_variation, ...
    variation_factor, num_nodes, 'sd50', verbose);
dens = process_bottom_variable(dens_input, target_dims, spatial_variation, ...
    variation_factor, num_nodes, 'dens', verbose);
mudp = process_bottom_variable(mudp_input, target_dims, spatial_variation, ...
    variation_factor, num_nodes, 'mudp', verbose);
poro = process_bottom_variable(poro_input, target_dims, spatial_variation, ...
    variation_factor, num_nodes, 'poro', verbose);
Tsed = process_bottom_variable(Tsed_input, target_dims, spatial_variation, ...
    variation_factor, num_nodes, 'Tsed', verbose);
wDcy = process_bottom_variable(wDcy_input, target_dims, spatial_variation, ...
    variation_factor, num_nodes, 'wDcy', verbose);
Lsed = process_bottom_variable(Lsed_input, target_dims, spatial_variation, ...
    variation_factor, num_nodes, 'Lsed', verbose);

% Apply physical constraints
poro = max(0.1, min(0.8, poro));  % Porosity bounds
mudp = max(0.0, min(1.0, mudp));  % Mud percentage bounds
sd50 = max(1e-6, min(1e-2, sd50)); % Grain size bounds
dens = max(1500, min(4000, dens)); % Density bounds

if verbose
    fprintf('  Variable ranges after processing:\n');
    fprintf('    sd50: %.2e to %.2e m\n', min(sd50), max(sd50));
    fprintf('    dens: %.1f to %.1f kg/m^3\n', min(dens), max(dens));
    fprintf('    mudp: %.3f to %.3f (fraction)\n', min(mudp), max(mudp));
    fprintf('    poro: %.3f to %.3f (fraction)\n', min(poro), max(poro));
    fprintf('    Tsed: %.1f to %.1f °C\n', min(Tsed), max(Tsed));
    fprintf('    wDcy: %.2e to %.2e m/s\n', min(wDcy), max(wDcy));
    fprintf('    Lsed: %.3f to %.3f m\n', min(Lsed), max(Lsed));
end

% Create NetCDF file
if verbose
    fprintf('  Writing NetCDF file...\n');
end

% Delete existing file if it exists
if exist(output_filename, 'file')
    delete(output_filename);
end

% Create file
ncid = netcdf.create(output_filename, 'CLOBBER');

try
    % Define dimensions
    node_dimid = netcdf.defDim(ncid, 'node', num_nodes);

    % Define variables
    sd50_varid = netcdf.defVar(ncid, 'sd50', 'double', node_dimid);
    dens_varid = netcdf.defVar(ncid, 'dens', 'double', node_dimid);
    mudp_varid = netcdf.defVar(ncid, 'mudp', 'double', node_dimid);
    poro_varid = netcdf.defVar(ncid, 'poro', 'double', node_dimid);
    Tsed_varid = netcdf.defVar(ncid, 'Tsed', 'double', node_dimid);
    wDcy_varid = netcdf.defVar(ncid, 'wDcy', 'double', node_dimid);
    Lsed_varid = netcdf.defVar(ncid, 'Lsed', 'double', node_dimid);

    % Add attributes for sd50
    netcdf.putAtt(ncid, sd50_varid, 'long_name', 'Sediment median grain size');
    netcdf.putAtt(ncid, sd50_varid, 'units', 'm');
    netcdf.putAtt(ncid, sd50_varid, 'description', ...
        'Median grain size (D50) of bottom sediment');
    netcdf.putAtt(ncid, sd50_varid, 'valid_range', [1e-6, 1e-2]);

    % Add attributes for dens
    netcdf.putAtt(ncid, dens_varid, 'long_name', 'Sediment bulk density');
    netcdf.putAtt(ncid, dens_varid, 'units', 'kg/m^3');
    netcdf.putAtt(ncid, dens_varid, 'description', ...
        'Bulk density of bottom sediment (dry weight per unit volume)');
    netcdf.putAtt(ncid, dens_varid, 'valid_range', [1500, 4000]);

    % Add attributes for mudp
    netcdf.putAtt(ncid, mudp_varid, 'long_name', 'Mud percentage');
    netcdf.putAtt(ncid, mudp_varid, 'units', 'dimensionless');
    netcdf.putAtt(ncid, mudp_varid, 'description', ...
        'Fraction of sediment that is mud (0=pure sand, 1=pure mud)');
    netcdf.putAtt(ncid, mudp_varid, 'valid_range', [0.0, 1.0]);

    % Add attributes for poro
    netcdf.putAtt(ncid, poro_varid, 'long_name', 'Sediment porosity');
    netcdf.putAtt(ncid, poro_varid, 'units', 'dimensionless');
    netcdf.putAtt(ncid, poro_varid, 'description', ...
        'Porosity of bottom sediment (void space fraction)');
    netcdf.putAtt(ncid, poro_varid, 'valid_range', [0.1, 0.8]);

    % Add attributes for Tsed
    netcdf.putAtt(ncid, Tsed_varid, 'long_name', 'Sediment temperature');
    netcdf.putAtt(ncid, Tsed_varid, 'units', 'degrees_C');
    netcdf.putAtt(ncid, Tsed_varid, 'description', ...
        'Temperature of bottom sediment (affects chemical kinetics)');

    % Add attributes for wDcy
    netcdf.putAtt(ncid, wDcy_varid, 'long_name', 'Vertical Darcy velocity');
    netcdf.putAtt(ncid, wDcy_varid, 'units', 'm/s');
    netcdf.putAtt(ncid, wDcy_varid, 'description', ...
        'Vertical pore water velocity in sediment (positive upward)');

    % Add attributes for Lsed
    netcdf.putAtt(ncid, Lsed_varid, 'long_name', 'Sediment layer thickness');
    netcdf.putAtt(ncid, Lsed_varid, 'units', 'm');
    netcdf.putAtt(ncid, Lsed_varid, 'description', ...
        'Thickness of active sediment layer');

    % Add global attributes
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'title', ...
        'FVCOM Bottom Sediment Features');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'source', ...
        'Generated by MATLAB script for FVCOM plastic transport model');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'creation_date', datestr(now));
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'conventions', 'CF-1.6');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'dimensions_info', ...
        sprintf('nodes=%d', num_nodes));
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'spatial_variation', spatial_variation);

    if ~strcmp(spatial_variation, 'uniform')
        netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'variation_factor', variation_factor);
    end

    % Add variable information
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'variables_info', ...
        'sd50, dens, mudp, poro, Tsed, wDcy, Lsed');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'compatible_subroutine', ...
        'LOAD_BOTTOM_FEATURES');

    % End definition mode
    netcdf.endDef(ncid);

    % Write data (ensure column vectors)
    netcdf.putVar(ncid, sd50_varid, sd50(:));
    netcdf.putVar(ncid, dens_varid, dens(:));
    netcdf.putVar(ncid, mudp_varid, mudp(:));
    netcdf.putVar(ncid, poro_varid, poro(:));
    netcdf.putVar(ncid, Tsed_varid, Tsed(:));
    netcdf.putVar(ncid, wDcy_varid, wDcy(:));
    netcdf.putVar(ncid, Lsed_varid, Lsed(:));

    % Close file
    netcdf.close(ncid);

    if verbose
        fprintf('  Successfully created NetCDF file: %s\n', output_filename);

        % Display file information
        fprintf('\nFile Information:\n');
        fprintf('  Variables created:\n');
        fprintf('    sd50 : Sediment median grain size [m]\n');
        fprintf('    dens : Sediment bulk density [kg/m^3]\n');
        fprintf('    mudp : Mud percentage [0-1]\n');
        fprintf('    poro : Sediment porosity [0-1]\n');
        fprintf('    Tsed : Sediment temperature [°C]\n');
        fprintf('    wDcy : Vertical Darcy velocity [m/s]\n');
        fprintf('    Lsed : Sediment layer thickness [m]\n');

        % Display representative values
        fprintf('\n  Representative values:\n');
        fprintf('    Node 1: sd50=%.2e, dens=%.1f, poro=%.3f, Tsed=%.1f\n', ...
            sd50(1), dens(1), poro(1), Tsed(1));
        if num_nodes > 1
            mid_node = round(num_nodes/2);
            fprintf('    Node %d: sd50=%.2e, dens=%.1f, poro=%.3f, Tsed=%.1f\n', ...
                mid_node, sd50(mid_node), dens(mid_node), poro(mid_node), Tsed(mid_node));
            fprintf('    Node %d: sd50=%.2e, dens=%.1f, poro=%.3f, Tsed=%.1f\n', ...
                num_nodes, sd50(end), dens(end), poro(end), Tsed(end));
        end
    end

catch ME
    % Close file if error occurs
    if exist('ncid', 'var')
        netcdf.close(ncid);
    end
    rethrow(ME);
end

end

function var_array = process_bottom_variable(var_input, target_dims, ...
    spatial_variation, variation_factor, num_nodes, var_name, verbose)
% Helper function to process bottom variable input and apply spatial variation

if isscalar(var_input)
    % Create base variable array
    base_value = var_input;
    var_array = ones(target_dims) * base_value;

    % Apply spatial variation if not uniform
    if ~strcmp(spatial_variation, 'uniform')
        if verbose
            fprintf('  Applying %s spatial variation to %s...\n', spatial_variation, var_name);
        end

        switch spatial_variation
            case 'random'
                % Random variation around base value
                variation = (rand(num_nodes, 1) - 0.5) * 2 * variation_factor;
                var_array = base_value * (1 + variation);

            case 'gradient'
                % Linear gradient across domain
                gradient = linspace(-variation_factor, variation_factor, num_nodes)';
                var_array = base_value * (1 + gradient);
        end
    end

else
    % Validate array dimensions
    if ~isequal(size(var_input), target_dims)
        % Try to handle different orientations
        if isvector(var_input) && length(var_input) == num_nodes
            var_array = var_input(:);  % Force column vector
        else
            error('Array %s has dimensions [%s], but expected [%s]', ...
                var_name, sprintf('%d ', size(var_input)), ...
                sprintf('%d ', target_dims));
        end
    else
        var_array = var_input;
    end
end

end
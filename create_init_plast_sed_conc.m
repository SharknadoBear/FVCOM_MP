function create_init_plast_sed_conc(output_filename, varargin)
% CREATE_INIT_PLAST_SED_CONC - Creates a NetCDF file for initializing plastic 
% concentration in sediment layers for FVCOM plastic transport model
%
% SYNTAX:
%   create_init_plast_sed_conc(output_filename)
%   create_init_plast_sed_conc(output_filename, 'Parameter', Value, ...)
%
% INPUTS:
%   output_filename - String, name of the NetCDF file to create
%
% OPTIONAL PARAMETERS:
%   'num_nodes'         - Number of nodes in the mesh (default: 100)
%   'num_sed_layers'    - Number of sediment layers (default: 3)
%   'num_plastic_types' - Number of plastic classes (default: 5)
%   'conc_total'        - Total plastic concentration array or scalar (default: 0.0)
%   'conc_stopped'      - Stopped plastic concentration array or scalar (default: 0.0)
%   'conc_passing'      - Passing plastic concentration array or scalar (default: 0.0)
%   'plastic_names'     - Cell array of plastic type names (optional)
%   'units'             - Units for concentration (default: 'kg/m^3')
%   'verbose'           - Display progress messages (default: true)
%
% EXAMPLES:
%   % Basic usage with default values
%   create_init_plast_sed_conc('init_plastic_sed.nc');
%
%   % Specify custom dimensions
%   create_init_plast_sed_conc('init_plastic_sed.nc', ...
%       'num_nodes', 500, 'num_sed_layers', 5, 'num_plastic_types', 3);
%
%   % Set uniform concentration values
%   create_init_plast_sed_conc('init_plastic_sed.nc', ...
%       'conc_total', 0.1, 'conc_stopped', 0.05, 'conc_passing', 0.05);
%
%   % Use custom concentration arrays
%   nodes = 100; layers = 3; types = 5;
%   conc_array = rand(nodes, layers, types) * 0.1;
%   create_init_plast_sed_conc('init_plastic_sed.nc', ...
%       'num_nodes', nodes, 'conc_total', conc_array);
%
% NOTES:
%   - The NetCDF file created is compatible with FVCOM's LOAD_INIT_PLAST_SED_CONC subroutine
%   - Concentration arrays should have dimensions: (num_nodes, num_sed_layers, num_plastic_types)
%   - If scalar values are provided for concentrations, they will be expanded to full arrays
%   - The sum constraint: conc_total = conc_stopped + conc_passing is enforced
%
% Author: Generated for FVCOM plastic transport model
% Date: September 2025

% Parse input arguments
p = inputParser;
addRequired(p, 'output_filename', @ischar);
addParameter(p, 'num_nodes', 100, @(x) isnumeric(x) && x > 0);
addParameter(p, 'num_sed_layers', 3, @(x) isnumeric(x) && x > 0);
addParameter(p, 'num_plastic_types', 5, @(x) isnumeric(x) && x > 0);
addParameter(p, 'conc_total', 0.0, @isnumeric);
addParameter(p, 'conc_stopped', 0.0, @isnumeric);
addParameter(p, 'conc_passing', 0.0, @isnumeric);
addParameter(p, 'plastic_names', {}, @iscell);
addParameter(p, 'units', 'kg/m^3', @ischar);
addParameter(p, 'verbose', true, @islogical);

parse(p, output_filename, varargin{:});

% Extract parsed values
num_nodes = p.Results.num_nodes;
num_sed_layers = p.Results.num_sed_layers;
num_plastic_types = p.Results.num_plastic_types;
conc_total_input = p.Results.conc_total;
conc_stopped_input = p.Results.conc_stopped;
conc_passing_input = p.Results.conc_passing;
plastic_names = p.Results.plastic_names;
units = p.Results.units;
verbose = p.Results.verbose;

if verbose
    fprintf('Creating NetCDF file for FVCOM plastic sediment initialization...\n');
    fprintf('  Output file: %s\n', output_filename);
    fprintf('  Dimensions: %d nodes, %d sediment layers, %d plastic types\n', ...
        num_nodes, num_sed_layers, num_plastic_types);
end



% Set default plastic names if not provided
if isempty(plastic_names)
    plastic_names = cell(num_plastic_types, 1);
    for i = 1:num_plastic_types
        plastic_names{i} = sprintf('Plastic_Type_%d', i);
    end
end

% Validate plastic names
if length(plastic_names) ~= num_plastic_types
    error('Length of plastic_names (%d) must match num_plastic_types (%d)', ...
        length(plastic_names), num_plastic_types);
end

% Process concentration inputs
target_dims = [num_nodes, num_sed_layers, num_plastic_types];

% Function to expand scalar or validate array
expand_or_validate = @(input_val, var_name) process_concentration_input(input_val, target_dims, var_name);

conc_total = expand_or_validate(conc_total_input, 'conc_total');
conc_stopped = expand_or_validate(conc_stopped_input, 'conc_stopped');
conc_passing = expand_or_validate(conc_passing_input, 'conc_passing');

% Enforce mass balance constraint: total = stopped + passing
% If total is zero and others are non-zero, calculate total
if all(conc_total(:) == 0) && (any(conc_stopped(:) ~= 0) || any(conc_passing(:) ~= 0))
    conc_total = conc_stopped + conc_passing;
    if verbose
        fprintf('  Calculated total concentration from stopped + passing components\n');
    end
% If total is non-zero but components are zero, distribute evenly
elseif any(conc_total(:) ~= 0) && all(conc_stopped(:) == 0) && all(conc_passing(:) == 0)
    conc_stopped = conc_total * 0.5;
    conc_passing = conc_total * 0.5;
    if verbose
        fprintf('  Distributed total concentration evenly between stopped and passing\n');
    end
% If all are specified, check consistency (with tolerance)
elseif any(conc_total(:) ~= 0) && (any(conc_stopped(:) ~= 0) || any(conc_passing(:) ~= 0))
    calculated_total = conc_stopped + conc_passing;
    tolerance = 1e-10;
    if max(abs(conc_total(:) - calculated_total(:))) > tolerance
        warning('Total concentration does not equal stopped + passing. Adjusting to maintain mass balance.');
        conc_total = calculated_total;
    end
end

if verbose
    fprintf('  Concentration ranges:\n');
    fprintf('    Total: %.6f to %.6f %s\n', min(conc_total(:)), max(conc_total(:)), units);
    fprintf('    Stopped: %.6f to %.6f %s\n', min(conc_stopped(:)), max(conc_stopped(:)), units);
    fprintf('    Passing: %.6f to %.6f %s\n', min(conc_passing(:)), max(conc_passing(:)), units);
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
    sed_layer_dimid = netcdf.defDim(ncid, 'N_sed_layer', num_sed_layers);
    nplast_dimid = netcdf.defDim(ncid, 'nplast', num_plastic_types);
    
    % Define variables for concentrations
    init_conc_sed_varid = netcdf.defVar(ncid, 'init_conc_sed', 'double', ...
        [node_dimid, sed_layer_dimid, nplast_dimid]);
    init_conc_sed_s_varid = netcdf.defVar(ncid, 'init_conc_sed_s', 'double', ...
        [node_dimid, sed_layer_dimid, nplast_dimid]);
    init_conc_sed_p_varid = netcdf.defVar(ncid, 'init_conc_sed_p', 'double', ...
        [node_dimid, sed_layer_dimid, nplast_dimid]);
    

    
    % Add attributes for main variables
    netcdf.putAtt(ncid, init_conc_sed_varid, 'long_name', 'Total plastic concentration in sediment');
    netcdf.putAtt(ncid, init_conc_sed_varid, 'units', units);
    netcdf.putAtt(ncid, init_conc_sed_varid, 'description', 'Total concentration of plastic in each sediment layer');
    
    netcdf.putAtt(ncid, init_conc_sed_s_varid, 'long_name', 'Stopped plastic concentration in sediment');
    netcdf.putAtt(ncid, init_conc_sed_s_varid, 'units', units);
    netcdf.putAtt(ncid, init_conc_sed_s_varid, 'description', 'Concentration of stopped/settled plastic in each sediment layer');
    
    netcdf.putAtt(ncid, init_conc_sed_p_varid, 'long_name', 'Passing plastic concentration in sediment');
    netcdf.putAtt(ncid, init_conc_sed_p_varid, 'units', units);
    netcdf.putAtt(ncid, init_conc_sed_p_varid, 'description', 'Concentration of passing/mobile plastic in each sediment layer');
    
    % Add global attributes
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'title', ...
        'FVCOM Initial Plastic Concentration in Sediment');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'source', ...
        'Generated by MATLAB script for FVCOM plastic transport model');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'creation_date', datestr(now));
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'conventions', 'CF-1.6');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'dimensions_info', ...
        sprintf('nodes=%d, sediment_layers=%d, plastic_types=%d', ...
        num_nodes, num_sed_layers, num_plastic_types));
    
    % Add plastic type names as attributes
    for i = 1:num_plastic_types
        attr_name = sprintf('plastic_type_%d_name', i);
        netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), attr_name, plastic_names{i});
    end
    
    % End definition mode
    netcdf.endDef(ncid);
    
    % Write data
    netcdf.putVar(ncid, init_conc_sed_varid, conc_total);
    netcdf.putVar(ncid, init_conc_sed_s_varid, conc_stopped);
    netcdf.putVar(ncid, init_conc_sed_p_varid, conc_passing);
    
    % Close file
    netcdf.close(ncid);
    
    if verbose
        fprintf('  Successfully created NetCDF file: %s\n', output_filename);
        
        % Display file information
        fprintf('\nFile Information:\n');
        fprintf('  Variables created:\n');
        fprintf('    init_conc_sed    : Total plastic concentration [%s]\n', units);
        fprintf('    init_conc_sed_s  : Stopped plastic concentration [%s]\n', units);
        fprintf('    init_conc_sed_p  : Passing plastic concentration [%s]\n', units);
        fprintf('\n  Plastic types: ');
        fprintf('%s ', plastic_names{:});
        fprintf('\n');
    end
    
catch ME
    % Close file if error occurs
    if exist('ncid', 'var')
        netcdf.close(ncid);
    end
    rethrow(ME);
end

end

% Helper function to process concentration inputs
function output_array = process_concentration_input(input_val, target_dims, var_name)
    if isscalar(input_val)
        % Expand scalar to full array
        output_array = ones(target_dims) * input_val;
    else
        % Validate array dimensions
        if ~isequal(size(input_val), target_dims)
            error('Array %s has dimensions [%s], but expected [%s]', ...
                var_name, ...
                sprintf('%d ', size(input_val)), ...
                sprintf('%d ', target_dims));
        end
        output_array = input_val;
    end
end

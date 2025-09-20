function create_init_tox_sed_conc(output_filename, varargin)
% CREATE_INIT_TOX_SED_CONC - Creates a NetCDF file for initializing toxin
% concentration in sediment layers for FVCOM toxin transport model
%
% SYNTAX:
%   create_init_tox_sed_conc(output_filename)
%   create_init_tox_sed_conc(output_filename, 'Parameter', Value, ...)
%
% INPUTS:
%   output_filename - String, name of the NetCDF file to create
%
% OPTIONAL PARAMETERS:
%   'num_nodes'         - Number of nodes in the mesh (default: 100)
%   'num_sed_layers'    - Number of sediment layers (default: 3)
%   'num_plastic_types' - Number of plastic classes (default: 5)
%   'init_tconc_b_sed'  - Initial bound toxin concentration in sediment array or scalar (default: 0.0)
%   'init_tconc_b_sed_s' - Initial bound toxin stopped concentration array or scalar (default: 0.0)
%   'init_tconc_b_sed_p' - Initial bound toxin passing concentration array or scalar (default: 0.0)
%   'init_tconc_f_sed'  - Initial free toxin concentration in sediment array or scalar (default: 0.0)
%   'plastic_names'     - Cell array of plastic type names (optional)
%   'units'             - Units for concentration (default: 'kg/m^3')
%   'depth_profile'     - Vertical distribution profile ('uniform', 'exponential', 'custom')
%   'decay_factor'      - Exponential decay factor with depth (default: 0.5)
%   'verbose'           - Display progress messages (default: true)
%
% EXAMPLES:
%   % Basic usage with default values
%   create_init_tox_sed_conc('init_toxin_sed.nc');
%
%   % Specify custom dimensions and uniform concentrations
%   create_init_tox_sed_conc('init_toxin_sed.nc', ...
%       'num_nodes', 500, 'num_sed_layers', 5, 'num_plastic_types', 3, ...
%       'init_tconc_b_sed', 0.02, 'init_tconc_f_sed', 0.01);
%
%   % Exponential decay profile with depth
%   create_init_tox_sed_conc('init_toxin_sed.nc', ...
%       'depth_profile', 'exponential', 'decay_factor', 0.3, ...
%       'init_tconc_b_sed', 0.05, 'init_tconc_f_sed', 0.02);
%
%   % Use custom concentration arrays
%   nodes = 100; layers = 3; types = 5;
%   tconc_b_array = rand(nodes, layers, types) * 0.02;
%   tconc_f_array = rand(nodes, layers, types) * 0.01;
%   create_init_tox_sed_conc('init_toxin_sed.nc', ...
%       'num_nodes', nodes, 'init_tconc_b_sed', tconc_b_array, ...
%       'init_tconc_f_sed', tconc_f_array);
%
% NOTES:
%   - Compatible with FVCOM's LOAD_INIT_TOXIN_SED_CONC subroutine
%   - Concentration arrays have dimensions: (num_nodes, num_sed_layers, num_plastic_types)
%   - The 'N_sed_layer' dimension in NetCDF corresponds to sediment layers
%   - Sediment layers are numbered from 1 (top) to N_sed_layer (bottom)
%   - Creates four variables: 'init_tconc_b_sed', 'init_tconc_b_sed_s',
%     'init_tconc_b_sed_p', and 'init_tconc_f_sed'
%   - Constraint: init_tconc_b_sed = init_tconc_b_sed_s + init_tconc_b_sed_p
%
% Author: Generated for FVCOM toxin transport model
% Date: September 2025

% Parse input arguments
p = inputParser;
addRequired(p, 'output_filename', @ischar);
addParameter(p, 'num_nodes', 100, @(x) isnumeric(x) && x > 0);
addParameter(p, 'num_sed_layers', 3, @(x) isnumeric(x) && x > 0);
addParameter(p, 'num_plastic_types', 5, @(x) isnumeric(x) && x > 0);
addParameter(p, 'init_tconc_b_sed', 0.0, @isnumeric);
addParameter(p, 'init_tconc_b_sed_s', 0.0, @isnumeric);
addParameter(p, 'init_tconc_b_sed_p', 0.0, @isnumeric);
addParameter(p, 'init_tconc_f_sed', 0.0, @isnumeric);
addParameter(p, 'plastic_names', {}, @iscell);
addParameter(p, 'units', 'kg/m^3', @ischar);
addParameter(p, 'depth_profile', 'uniform', @(x) ismember(x, {'uniform', 'exponential', 'custom'}));
addParameter(p, 'decay_factor', 0.5, @(x) isnumeric(x) && x > 0);
addParameter(p, 'verbose', true, @islogical);

parse(p, output_filename, varargin{:});

% Extract parsed values
num_nodes = p.Results.num_nodes;
num_sed_layers = p.Results.num_sed_layers;
num_plastic_types = p.Results.num_plastic_types;
init_tconc_b_sed_input = p.Results.init_tconc_b_sed;
init_tconc_b_sed_s_input = p.Results.init_tconc_b_sed_s;
init_tconc_b_sed_p_input = p.Results.init_tconc_b_sed_p;
init_tconc_f_sed_input = p.Results.init_tconc_f_sed;
plastic_names = p.Results.plastic_names;
units = p.Results.units;
depth_profile = p.Results.depth_profile;
decay_factor = p.Results.decay_factor;
verbose = p.Results.verbose;

if verbose
    fprintf('Creating NetCDF file for FVCOM toxin sediment initialization...\n');
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

% Process bound toxin total concentration
init_tconc_b_sed = process_concentration_input(init_tconc_b_sed_input, target_dims, ...
    depth_profile, decay_factor, num_nodes, num_sed_layers, num_plastic_types, ...
    'bound toxin total', verbose);

% Process bound toxin stopped concentration
init_tconc_b_sed_s = process_concentration_input(init_tconc_b_sed_s_input, target_dims, ...
    depth_profile, decay_factor, num_nodes, num_sed_layers, num_plastic_types, ...
    'bound toxin stopped', verbose);

% Process bound toxin passing concentration
init_tconc_b_sed_p = process_concentration_input(init_tconc_b_sed_p_input, target_dims, ...
    depth_profile, decay_factor, num_nodes, num_sed_layers, num_plastic_types, ...
    'bound toxin passing', verbose);

% Process free toxin concentration
init_tconc_f_sed = process_concentration_input(init_tconc_f_sed_input, target_dims, ...
    depth_profile, decay_factor, num_nodes, num_sed_layers, num_plastic_types, ...
    'free toxin', verbose);

% Ensure constraint: init_tconc_b_sed = init_tconc_b_sed_s + init_tconc_b_sed_p
if isscalar(init_tconc_b_sed_input) && isscalar(init_tconc_b_sed_s_input) && isscalar(init_tconc_b_sed_p_input)
    % Automatically adjust if using scalar inputs
    total_input = init_tconc_b_sed_input;
    stopped_input = init_tconc_b_sed_s_input;
    passing_input = init_tconc_b_sed_p_input;

    if total_input == 0 && (stopped_input > 0 || passing_input > 0)
        % Set total to sum of parts
        init_tconc_b_sed = init_tconc_b_sed_s + init_tconc_b_sed_p;
        if verbose
            fprintf('  Auto-adjusting total bound toxin to match sum of stopped + passing\n');
        end
    elseif total_input > 0 && stopped_input == 0 && passing_input == 0
        % Distribute total equally between stopped and passing
        init_tconc_b_sed_s = init_tconc_b_sed * 0.5;
        init_tconc_b_sed_p = init_tconc_b_sed * 0.5;
        if verbose
            fprintf('  Auto-distributing total bound toxin equally between stopped and passing\n');
        end
    end
end

if verbose
    fprintf('  Bound toxin total range: %.6f to %.6f %s\n', ...
        min(init_tconc_b_sed(:)), max(init_tconc_b_sed(:)), units);
    fprintf('  Bound toxin stopped range: %.6f to %.6f %s\n', ...
        min(init_tconc_b_sed_s(:)), max(init_tconc_b_sed_s(:)), units);
    fprintf('  Bound toxin passing range: %.6f to %.6f %s\n', ...
        min(init_tconc_b_sed_p(:)), max(init_tconc_b_sed_p(:)), units);
    fprintf('  Free toxin range: %.6f to %.6f %s\n', ...
        min(init_tconc_f_sed(:)), max(init_tconc_f_sed(:)), units);

    if ~strcmp(depth_profile, 'uniform')
        fprintf('  Surface/bottom ratio (bound total): %.2f\n', ...
            mean(init_tconc_b_sed(:, 1, :), 'all') / mean(init_tconc_b_sed(:, end, :), 'all'));
        fprintf('  Surface/bottom ratio (free): %.2f\n', ...
            mean(init_tconc_f_sed(:, 1, :), 'all') / mean(init_tconc_f_sed(:, end, :), 'all'));
    end
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
    N_sed_layer_dimid = netcdf.defDim(ncid, 'N_sed_layer', num_sed_layers);
    nplast_dimid = netcdf.defDim(ncid, 'nplast', num_plastic_types);

    % Define variables for toxin concentrations
    init_tconc_b_sed_varid = netcdf.defVar(ncid, 'init_tconc_b_sed', 'double', ...
        [node_dimid, N_sed_layer_dimid, nplast_dimid]);
    init_tconc_b_sed_s_varid = netcdf.defVar(ncid, 'init_tconc_b_sed_s', 'double', ...
        [node_dimid, N_sed_layer_dimid, nplast_dimid]);
    init_tconc_b_sed_p_varid = netcdf.defVar(ncid, 'init_tconc_b_sed_p', 'double', ...
        [node_dimid, N_sed_layer_dimid, nplast_dimid]);
    init_tconc_f_sed_varid = netcdf.defVar(ncid, 'init_tconc_f_sed', 'double', ...
        [node_dimid, N_sed_layer_dimid, nplast_dimid]);

    % Add attributes for bound toxin total variable
    netcdf.putAtt(ncid, init_tconc_b_sed_varid, 'long_name', 'Initial bound toxin concentration in sediment');
    netcdf.putAtt(ncid, init_tconc_b_sed_varid, 'units', units);
    netcdf.putAtt(ncid, init_tconc_b_sed_varid, 'description', ...
        'Total initial concentration of bound toxin in sediment layers');
    netcdf.putAtt(ncid, init_tconc_b_sed_varid, 'coordinates', 'node N_sed_layer nplast');

    % Add attributes for bound toxin stopped variable
    netcdf.putAtt(ncid, init_tconc_b_sed_s_varid, 'long_name', 'Initial stopped bound toxin concentration in sediment');
    netcdf.putAtt(ncid, init_tconc_b_sed_s_varid, 'units', units);
    netcdf.putAtt(ncid, init_tconc_b_sed_s_varid, 'description', ...
        'Initial concentration of stopped bound toxin in sediment layers (blocked fraction)');
    netcdf.putAtt(ncid, init_tconc_b_sed_s_varid, 'coordinates', 'node N_sed_layer nplast');

    % Add attributes for bound toxin passing variable
    netcdf.putAtt(ncid, init_tconc_b_sed_p_varid, 'long_name', 'Initial passing bound toxin concentration in sediment');
    netcdf.putAtt(ncid, init_tconc_b_sed_p_varid, 'units', units);
    netcdf.putAtt(ncid, init_tconc_b_sed_p_varid, 'description', ...
        'Initial concentration of passing bound toxin in sediment layers (mobile fraction)');
    netcdf.putAtt(ncid, init_tconc_b_sed_p_varid, 'coordinates', 'node N_sed_layer nplast');

    % Add attributes for free toxin variable
    netcdf.putAtt(ncid, init_tconc_f_sed_varid, 'long_name', 'Initial free toxin concentration in sediment');
    netcdf.putAtt(ncid, init_tconc_f_sed_varid, 'units', units);
    netcdf.putAtt(ncid, init_tconc_f_sed_varid, 'description', ...
        'Initial concentration of free toxin in sediment layers');
    netcdf.putAtt(ncid, init_tconc_f_sed_varid, 'coordinates', 'node N_sed_layer nplast');

    % Add global attributes
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'title', ...
        'FVCOM Initial Toxin Concentration in Sediment');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'source', ...
        'Generated by MATLAB script for FVCOM toxin transport model');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'creation_date', datestr(now));
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'conventions', 'CF-1.6');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'dimensions_info', ...
        sprintf('nodes=%d, sediment_layers=%d, plastic_types=%d', ...
        num_nodes, num_sed_layers, num_plastic_types));
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'depth_profile', depth_profile);

    if ~strcmp(depth_profile, 'uniform')
        netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'decay_factor', decay_factor);
    end

    % Add plastic type names as attributes
    for i = 1:num_plastic_types
        attr_name = sprintf('plastic_type_%d_name', i);
        netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), attr_name, plastic_names{i});
    end

    % Add layer information
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'sediment_layer_info', ...
        'Sediment layers: 1 (top/surface) to N_sed_layer (bottom)');

    % End definition mode
    netcdf.endDef(ncid);

    % Write data
    netcdf.putVar(ncid, init_tconc_b_sed_varid, init_tconc_b_sed);
    netcdf.putVar(ncid, init_tconc_b_sed_s_varid, init_tconc_b_sed_s);
    netcdf.putVar(ncid, init_tconc_b_sed_p_varid, init_tconc_b_sed_p);
    netcdf.putVar(ncid, init_tconc_f_sed_varid, init_tconc_f_sed);

    % Close file
    netcdf.close(ncid);

    if verbose
        fprintf('  Successfully created NetCDF file: %s\n', output_filename);

        % Display file information
        fprintf('\nFile Information:\n');
        fprintf('  Variables created:\n');
        fprintf('    init_tconc_b_sed   : Initial bound toxin concentration in sediment [%s]\n', units);
        fprintf('    init_tconc_b_sed_s : Initial stopped bound toxin in sediment [%s]\n', units);
        fprintf('    init_tconc_b_sed_p : Initial passing bound toxin in sediment [%s]\n', units);
        fprintf('    init_tconc_f_sed   : Initial free toxin concentration in sediment [%s]\n', units);
        fprintf('  Sediment layers: 1 (top/surface) to %d (bottom)\n', num_sed_layers);
        fprintf('  Plastic types: ');
        fprintf('%s ', plastic_names{:});
        fprintf('\n');

        % Display layer-averaged concentrations
        if verbose && num_sed_layers <= 10  % Only show for reasonable number of layers
            fprintf('\n  Layer-averaged concentrations (bound toxin total):\n');
            for k = 1:num_sed_layers
                layer_avg = mean(mean(init_tconc_b_sed(:, k, :), 1), 3);
                if k == 1
                    fprintf('    k=%2d (top)   : %.6f %s\n', k, layer_avg, units);
                elseif k == num_sed_layers
                    fprintf('    k=%2d (bottom): %.6f %s\n', k, layer_avg, units);
                else
                    fprintf('    k=%2d         : %.6f %s\n', k, layer_avg, units);
                end
            end

            fprintf('\n  Layer-averaged concentrations (free toxin):\n');
            for k = 1:num_sed_layers
                layer_avg = mean(mean(init_tconc_f_sed(:, k, :), 1), 3);
                if k == 1
                    fprintf('    k=%2d (top)   : %.6f %s\n', k, layer_avg, units);
                elseif k == num_sed_layers
                    fprintf('    k=%2d (bottom): %.6f %s\n', k, layer_avg, units);
                else
                    fprintf('    k=%2d         : %.6f %s\n', k, layer_avg, units);
                end
            end
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

function conc_array = process_concentration_input(conc_input, target_dims, ...
    depth_profile, decay_factor, num_nodes, num_sed_layers, num_plastic_types, ...
    conc_type, verbose)
% Helper function to process concentration input and apply depth profiles

if isscalar(conc_input)
    % Create base concentration array
    base_conc = conc_input;
    conc_array = ones(target_dims) * base_conc;

    % Apply depth profile if not uniform
    if ~strcmp(depth_profile, 'uniform')
        if verbose
            fprintf('  Applying %s depth profile to %s...\n', depth_profile, conc_type);
        end

        % Create normalized depth coordinate (0 = top, 1 = bottom)
        depth_coord = linspace(0, 1, num_sed_layers)';

        switch depth_profile
            case 'exponential'
                % Exponential decay with depth
                profile_factor = exp(-depth_coord / decay_factor);

            case 'custom'
                % User can modify this section for custom profiles
                profile_factor = ones(num_sed_layers, 1);
                warning('Custom profile selected but not implemented. Using uniform distribution.');
        end

        % Apply profile to all nodes and plastic types
        for inode = 1:num_nodes
            for iplastic = 1:num_plastic_types
                conc_array(inode, :, iplastic) = base_conc * profile_factor';
            end
        end
    end

else
    % Validate array dimensions
    if ~isequal(size(conc_input), target_dims)
        error('Array %s has dimensions [%s], but expected [%s]', ...
            conc_type, sprintf('%d ', size(conc_input)), ...
            sprintf('%d ', target_dims));
    end
    conc_array = conc_input;
end

end
function create_init_tox_conc(output_filename, varargin)
% CREATE_INIT_TOX_CONC - Creates a NetCDF file for initializing toxin
% concentration in the water column and sediment for FVCOM toxin transport model
%
% SYNTAX:
%   create_init_tox_conc(output_filename)
%   create_init_tox_conc(output_filename, 'Parameter', Value, ...)
%
% INPUTS:
%   output_filename - String, name of the NetCDF file to create
%
% OPTIONAL PARAMETERS:
%   'num_nodes'         - Number of nodes in the mesh (default: 100)
%   'num_sigma_layers'  - Number of sigma layers (kb) (default: 10)
%   'num_plastic_types' - Number of plastic classes (default: 5)
%   'init_tconc_b'      - Initial toxin concentration in bound form array or scalar (default: 0.0)
%   'init_tconc_f'      - Initial toxin concentration in free form array or scalar (default: 0.0)
%   'plastic_names'     - Cell array of plastic type names (optional)
%   'units'             - Units for concentration (default: 'kg/m^3')
%   'depth_profile'     - Vertical distribution profile ('uniform', 'surface', 'exponential', 'custom')
%   'surface_factor'    - Factor for surface concentration enhancement (default: 2.0)
%   'decay_scale'       - Decay scale for exponential profile (default: 0.3)
%   'verbose'           - Display progress messages (default: true)
%
% EXAMPLES:
%   % Basic usage with default values
%   create_init_tox_conc('init_toxin_water.nc');
%
%   % Specify custom dimensions and uniform concentration
%   create_init_tox_conc('init_toxin_water.nc', ...
%       'num_nodes', 500, 'num_sigma_layers', 15, 'num_plastic_types', 3, ...
%       'init_tconc_b', 0.05, 'init_tconc_f', 0.02);
%
%   % Surface-enhanced distribution
%   create_init_tox_conc('init_toxin_water.nc', ...
%       'depth_profile', 'surface', 'surface_factor', 5.0, ...
%       'init_tconc_b', 0.02, 'init_tconc_f', 0.01);
%
%   % Custom concentration arrays
%   nodes = 100; layers = 10; types = 5;
%   tconc_b_array = rand(nodes, layers, types) * 0.1;
%   tconc_f_array = rand(nodes, layers, types) * 0.05;
%   create_init_tox_conc('init_toxin_water.nc', ...
%       'num_nodes', nodes, 'init_tconc_b', tconc_b_array, 'init_tconc_f', tconc_f_array);
%
% NOTES:
%   - Compatible with FVCOM's LOAD_INIT_TOXIN_CONC subroutine
%   - Concentration arrays have dimensions: (num_nodes, num_sigma_layers, num_plastic_types)
%   - The 'ksl' dimension in NetCDF corresponds to sigma layers (kb in FVCOM)
%   - Sigma coordinates: surface (k=1) to bottom (k=kb)
%   - Creates two variables: 'init_tconc_b' (bound toxin) and 'init_tconc_f' (free toxin)
%
% Author: Generated for FVCOM toxin transport model
% Date: September 2025

% Parse input arguments
p = inputParser;
addRequired(p, 'output_filename', @ischar);
addParameter(p, 'num_nodes', 100, @(x) isnumeric(x) && x > 0);
addParameter(p, 'num_sigma_layers', 10, @(x) isnumeric(x) && x > 0);
addParameter(p, 'num_plastic_types', 5, @(x) isnumeric(x) && x > 0);
addParameter(p, 'init_tconc_b', 0.0, @isnumeric);
addParameter(p, 'init_tconc_f', 0.0, @isnumeric);
addParameter(p, 'plastic_names', {}, @iscell);
addParameter(p, 'units', 'kg/m^3', @ischar);
addParameter(p, 'depth_profile', 'uniform', @(x) ismember(x, {'uniform', 'surface', 'exponential', 'custom'}));
addParameter(p, 'surface_factor', 2.0, @(x) isnumeric(x) && x > 0);
addParameter(p, 'decay_scale', 0.3, @(x) isnumeric(x) && x > 0);
addParameter(p, 'verbose', true, @islogical);

parse(p, output_filename, varargin{:});

% Extract parsed values
num_nodes = p.Results.num_nodes;
num_sigma_layers = p.Results.num_sigma_layers;
num_plastic_types = p.Results.num_plastic_types;
init_tconc_b_input = p.Results.init_tconc_b;
init_tconc_f_input = p.Results.init_tconc_f;
plastic_names = p.Results.plastic_names;
units = p.Results.units;
depth_profile = p.Results.depth_profile;
surface_factor = p.Results.surface_factor;
decay_scale = p.Results.decay_scale;
verbose = p.Results.verbose;

if verbose
    fprintf('Creating NetCDF file for FVCOM toxin water column initialization...\n');
    fprintf('  Output file: %s\n', output_filename);
    fprintf('  Dimensions: %d nodes, %d sigma layers, %d plastic types\n', ...
        num_nodes, num_sigma_layers, num_plastic_types);
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
target_dims = [num_nodes, num_sigma_layers, num_plastic_types];

% Process bound toxin concentration
init_tconc_b = process_concentration_input(init_tconc_b_input, target_dims, ...
    depth_profile, surface_factor, decay_scale, num_nodes, num_sigma_layers, ...
    num_plastic_types, 'bound toxin', verbose);

% Process free toxin concentration
init_tconc_f = process_concentration_input(init_tconc_f_input, target_dims, ...
    depth_profile, surface_factor, decay_scale, num_nodes, num_sigma_layers, ...
    num_plastic_types, 'free toxin', verbose);

if verbose
    fprintf('  Bound toxin concentration range: %.6f to %.6f %s\n', ...
        min(init_tconc_b(:)), max(init_tconc_b(:)), units);
    fprintf('  Free toxin concentration range: %.6f to %.6f %s\n', ...
        min(init_tconc_f(:)), max(init_tconc_f(:)), units);
    if ~strcmp(depth_profile, 'uniform')
        fprintf('  Surface/bottom ratio (bound): %.2f\n', ...
            mean(init_tconc_b(:, 1, :), 'all') / mean(init_tconc_b(:, end, :), 'all'));
        fprintf('  Surface/bottom ratio (free): %.2f\n', ...
            mean(init_tconc_f(:, 1, :), 'all') / mean(init_tconc_f(:, end, :), 'all'));
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
    ksl_dimid = netcdf.defDim(ncid, 'ksl', num_sigma_layers);  % ksl = k sigma levels
    nplast_dimid = netcdf.defDim(ncid, 'nplast', num_plastic_types);

    % Define variables for initial toxin concentrations
    init_tconc_b_varid = netcdf.defVar(ncid, 'init_tconc_b', 'double', ...
        [node_dimid, ksl_dimid, nplast_dimid]);
    init_tconc_f_varid = netcdf.defVar(ncid, 'init_tconc_f', 'double', ...
        [node_dimid, ksl_dimid, nplast_dimid]);

    % Add attributes for bound toxin variable
    netcdf.putAtt(ncid, init_tconc_b_varid, 'long_name', 'Initial bound toxin concentration in water column');
    netcdf.putAtt(ncid, init_tconc_b_varid, 'units', units);
    netcdf.putAtt(ncid, init_tconc_b_varid, 'description', ...
        'Initial concentration of bound toxin in water column at each sigma layer');
    netcdf.putAtt(ncid, init_tconc_b_varid, 'coordinates', 'node ksl nplast');

    % Add attributes for free toxin variable
    netcdf.putAtt(ncid, init_tconc_f_varid, 'long_name', 'Initial free toxin concentration in water column');
    netcdf.putAtt(ncid, init_tconc_f_varid, 'units', units);
    netcdf.putAtt(ncid, init_tconc_f_varid, 'description', ...
        'Initial concentration of free toxin in water column at each sigma layer');
    netcdf.putAtt(ncid, init_tconc_f_varid, 'coordinates', 'node ksl nplast');

    % Add global attributes
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'title', ...
        'FVCOM Initial Toxin Concentration in Water Column');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'source', ...
        'Generated by MATLAB script for FVCOM toxin transport model');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'creation_date', datestr(now));
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'conventions', 'CF-1.6');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'dimensions_info', ...
        sprintf('nodes=%d, sigma_layers=%d, plastic_types=%d', ...
        num_nodes, num_sigma_layers, num_plastic_types));
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'vertical_profile', depth_profile);

    if ~strcmp(depth_profile, 'uniform')
        netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'surface_factor', surface_factor);
        if strcmp(depth_profile, 'exponential')
            netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'decay_scale', decay_scale);
        end
    end

    % Add plastic type names as attributes
    for i = 1:num_plastic_types
        attr_name = sprintf('plastic_type_%d_name', i);
        netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), attr_name, plastic_names{i});
    end

    % Add coordinate information
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'sigma_coordinate_info', ...
        'Sigma layers: k=1 (surface) to k=kb (bottom)');

    % End definition mode
    netcdf.endDef(ncid);

    % Write data
    netcdf.putVar(ncid, init_tconc_b_varid, init_tconc_b);
    netcdf.putVar(ncid, init_tconc_f_varid, init_tconc_f);

    % Close file
    netcdf.close(ncid);

    if verbose
        fprintf('  Successfully created NetCDF file: %s\n', output_filename);

        % Display file information
        fprintf('\nFile Information:\n');
        fprintf('  Variables created:\n');
        fprintf('    init_tconc_b : Initial bound toxin concentration in water column [%s]\n', units);
        fprintf('    init_tconc_f : Initial free toxin concentration in water column [%s]\n', units);
        fprintf('  Sigma layers: 1 (surface) to %d (bottom)\n', num_sigma_layers);
        fprintf('  Plastic types: ');
        fprintf('%s ', plastic_names{:});
        fprintf('\n');

        % Display layer-averaged concentrations
        if verbose && num_sigma_layers <= 20  % Only show for reasonable number of layers
            fprintf('\n  Layer-averaged concentrations (bound toxin):\n');
            for k = 1:num_sigma_layers
                layer_avg = mean(mean(init_tconc_b(:, k, :), 1), 3);
                if k == 1
                    fprintf('    k=%2d (surface): %.6f %s\n', k, layer_avg, units);
                elseif k == num_sigma_layers
                    fprintf('    k=%2d (bottom) : %.6f %s\n', k, layer_avg, units);
                elseif k <= 5 || k > num_sigma_layers - 3
                    fprintf('    k=%2d          : %.6f %s\n', k, layer_avg, units);
                elseif k == 6 && num_sigma_layers > 8
                    fprintf('    ...           \n');
                end
            end

            fprintf('\n  Layer-averaged concentrations (free toxin):\n');
            for k = 1:num_sigma_layers
                layer_avg = mean(mean(init_tconc_f(:, k, :), 1), 3);
                if k == 1
                    fprintf('    k=%2d (surface): %.6f %s\n', k, layer_avg, units);
                elseif k == num_sigma_layers
                    fprintf('    k=%2d (bottom) : %.6f %s\n', k, layer_avg, units);
                elseif k <= 5 || k > num_sigma_layers - 3
                    fprintf('    k=%2d          : %.6f %s\n', k, layer_avg, units);
                elseif k == 6 && num_sigma_layers > 8
                    fprintf('    ...           \n');
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
    depth_profile, surface_factor, decay_scale, num_nodes, num_sigma_layers, ...
    num_plastic_types, conc_type, verbose)
% Helper function to process concentration input and apply vertical profiles

if isscalar(conc_input)
    % Create base concentration array
    base_conc = conc_input;
    conc_array = ones(target_dims) * base_conc;

    % Apply vertical profile if not uniform
    if ~strcmp(depth_profile, 'uniform')
        if verbose
            fprintf('  Applying %s vertical profile to %s...\n', depth_profile, conc_type);
        end

        % Create sigma coordinate (normalized depth: 0 = surface, 1 = bottom)
        sigma = linspace(0, 1, num_sigma_layers)';

        switch depth_profile
            case 'surface'
                % Enhanced surface concentration, linear decrease
                profile_factor = surface_factor * (1 - sigma) + 1;

            case 'exponential'
                % Exponential decay from surface
                profile_factor = exp(-sigma / decay_scale);

            case 'custom'
                % User can modify this section for custom profiles
                profile_factor = ones(num_sigma_layers, 1);
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
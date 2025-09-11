function create_init_plast_conc(output_filename, varargin)
% CREATE_INIT_PLAST_CONC - Creates a NetCDF file for initializing plastic 
% concentration in the water column for FVCOM plastic transport model
%
% SYNTAX:
%   create_init_plast_conc(output_filename)
%   create_init_plast_conc(output_filename, 'Parameter', Value, ...)
%
% INPUTS:
%   output_filename - String, name of the NetCDF file to create
%
% OPTIONAL PARAMETERS:
%   'num_nodes'         - Number of nodes in the mesh (default: 100)
%   'num_sigma_layers'  - Number of sigma layers (kb) (default: 10)
%   'num_plastic_types' - Number of plastic classes (default: 5)
%   'init_conc'         - Initial plastic concentration array or scalar (default: 0.0)
%   'plastic_names'     - Cell array of plastic type names (optional)
%   'units'             - Units for concentration (default: 'kg/m^3')
%   'depth_profile'     - Vertical distribution profile ('uniform', 'surface', 'exponential', 'custom')
%   'surface_factor'    - Factor for surface concentration enhancement (default: 2.0)
%   'decay_scale'       - Decay scale for exponential profile (default: 0.3)
%   'verbose'           - Display progress messages (default: true)
%
% EXAMPLES:
%   % Basic usage with default values
%   create_init_plast_conc('init_plastic_water.nc');
%
%   % Specify custom dimensions and uniform concentration
%   create_init_plast_conc('init_plastic_water.nc', ...
%       'num_nodes', 500, 'num_sigma_layers', 15, 'num_plastic_types', 3, ...
%       'init_conc', 0.05);
%
%   % Surface-enhanced distribution
%   create_init_plast_conc('init_plastic_water.nc', ...
%       'depth_profile', 'surface', 'surface_factor', 5.0, 'init_conc', 0.02);
%
%   % Custom concentration array
%   nodes = 100; layers = 10; types = 5;
%   conc_array = rand(nodes, layers, types) * 0.1;
%   create_init_plast_conc('init_plastic_water.nc', ...
%       'num_nodes', nodes, 'init_conc', conc_array);
%
% NOTES:
%   - Compatible with FVCOM's LOAD_INIT_PLAST_CONC subroutine
%   - Concentration arrays have dimensions: (num_nodes, num_sigma_layers, num_plastic_types)
%   - The 'ksl' dimension in NetCDF corresponds to sigma layers (kb in FVCOM)
%   - Sigma coordinates: surface (k=1) to bottom (k=kb)
%
% Author: Generated for FVCOM plastic transport model
% Date: September 2025

% Parse input arguments
p = inputParser;
addRequired(p, 'output_filename', @ischar);
addParameter(p, 'num_nodes', 100, @(x) isnumeric(x) && x > 0);
addParameter(p, 'num_sigma_layers', 10, @(x) isnumeric(x) && x > 0);
addParameter(p, 'num_plastic_types', 5, @(x) isnumeric(x) && x > 0);
addParameter(p, 'init_conc', 0.0, @isnumeric);
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
init_conc_input = p.Results.init_conc;
plastic_names = p.Results.plastic_names;
units = p.Results.units;
depth_profile = p.Results.depth_profile;
surface_factor = p.Results.surface_factor;
decay_scale = p.Results.decay_scale;
verbose = p.Results.verbose;

if verbose
    fprintf('Creating NetCDF file for FVCOM plastic water column initialization...\n');
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

% Process concentration input
target_dims = [num_nodes, num_sigma_layers, num_plastic_types];

if isscalar(init_conc_input)
    % Create base concentration array
    base_conc = init_conc_input;
    init_conc = ones(target_dims) * base_conc;
    
    % Apply vertical profile if not uniform
    if ~strcmp(depth_profile, 'uniform')
        if verbose
            fprintf('  Applying %s vertical profile...\n', depth_profile);
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
                init_conc(inode, :, iplastic) = base_conc * profile_factor';
            end
        end
    end
    
else
    % Validate array dimensions
    if ~isequal(size(init_conc_input), target_dims)
        error('Array init_conc has dimensions [%s], but expected [%s]', ...
            sprintf('%d ', size(init_conc_input)), ...
            sprintf('%d ', target_dims));
    end
    init_conc = init_conc_input;
end

if verbose
    fprintf('  Concentration range: %.6f to %.6f %s\n', ...
        min(init_conc(:)), max(init_conc(:)), units);
    if ~strcmp(depth_profile, 'uniform')
        fprintf('  Surface/bottom ratio: %.2f\n', ...
            mean(init_conc(:, 1, :), 'all') / mean(init_conc(:, end, :), 'all'));
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
    
    % Define main variable for initial concentration
    init_conc_varid = netcdf.defVar(ncid, 'init_conc', 'double', ...
        [node_dimid, ksl_dimid, nplast_dimid]);
    
    % Add attributes for main variable
    netcdf.putAtt(ncid, init_conc_varid, 'long_name', 'Initial plastic concentration in water column');
    netcdf.putAtt(ncid, init_conc_varid, 'units', units);
    netcdf.putAtt(ncid, init_conc_varid, 'description', ...
        'Initial concentration of plastic in water column at each sigma layer');
    netcdf.putAtt(ncid, init_conc_varid, 'coordinates', 'node ksl nplast');
    
    % Add global attributes
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'title', ...
        'FVCOM Initial Plastic Concentration in Water Column');
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'source', ...
        'Generated by MATLAB script for FVCOM plastic transport model');
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
    netcdf.putVar(ncid, init_conc_varid, init_conc);
    
    % Close file
    netcdf.close(ncid);
    
    if verbose
        fprintf('  Successfully created NetCDF file: %s\n', output_filename);
        
        % Display file information
        fprintf('\nFile Information:\n');
        fprintf('  Variable created:\n');
        fprintf('    init_conc : Initial plastic concentration in water column [%s]\n', units);
        fprintf('  Sigma layers: 1 (surface) to %d (bottom)\n', num_sigma_layers);
        fprintf('  Plastic types: ');
        fprintf('%s ', plastic_names{:});
        fprintf('\n');
        
        % Display layer-averaged concentrations
        if verbose && num_sigma_layers <= 20  % Only show for reasonable number of layers
            fprintf('\n  Layer-averaged concentrations:\n');
            for k = 1:num_sigma_layers
                layer_avg = mean(mean(init_conc(:, k, :), 1), 3);
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

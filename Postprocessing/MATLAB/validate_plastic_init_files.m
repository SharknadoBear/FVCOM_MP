function validate_plastic_init_files(varargin)
% VALIDATE_PLASTIC_INIT_FILES - Validates FVCOM plastic initialization files
%
% SYNTAX:
%   validate_plastic_init_files()
%   validate_plastic_init_files('water_file', filename)
%   validate_plastic_init_files('sediment_file', filename)
%   validate_plastic_init_files('water_file', water_file, 'sediment_file', sed_file)
%   validate_plastic_init_files(..., 'mesh_file', mesh_file)
%
% INPUTS:
%   'water_file'    - Path to water column initialization NetCDF file
%   'sediment_file' - Path to sediment initialization NetCDF file  
%   'mesh_file'     - Path to FVCOM mesh file (optional, for node count validation)
%   'expected_kb'   - Expected number of sigma layers (optional)
%   'expected_nplast' - Expected number of plastic types (optional)
%
% EXAMPLES:
%   % Validate both files
%   validate_plastic_init_files('water_file', 'water_init.nc', ...
%                              'sediment_file', 'sed_init.nc');
%
%   % Check against mesh file
%   validate_plastic_init_files('water_file', 'water_init.nc', ...
%                              'mesh_file', 'my_mesh.nc');
%
%   % Quick validation of single file
%   validate_plastic_init_files('water_file', 'water_init.nc');
%
% Author: Generated for FVCOM plastic transport model
% Date: September 2025

% Parse inputs
p = inputParser;
addParameter(p, 'water_file', '', @ischar);
addParameter(p, 'sediment_file', '', @ischar);
addParameter(p, 'mesh_file', '', @ischar);
addParameter(p, 'expected_kb', [], @(x) isnumeric(x) && x > 0);
addParameter(p, 'expected_nplast', [], @(x) isnumeric(x) && x > 0);

parse(p, varargin{:});

water_file = p.Results.water_file;
sediment_file = p.Results.sediment_file;
mesh_file = p.Results.mesh_file;
expected_kb = p.Results.expected_kb;
expected_nplast = p.Results.expected_nplast;

fprintf('=== FVCOM Plastic Initialization File Validator ===\n\n');

% Check if any files provided
if isempty(water_file) && isempty(sediment_file)
    error('Must provide at least one file to validate');
end

% Validation results
validation_passed = true;
warnings_count = 0;

%% Validate mesh file if provided
mesh_info = [];
if ~isempty(mesh_file)
    fprintf('Validating mesh file: %s\n', mesh_file);
    try
        mesh_info = ncinfo(mesh_file);
        
        % Try to find node dimension in mesh
        mesh_dims = {mesh_info.Dimensions.Name};
        if ismember('node', mesh_dims)
            mesh_nodes = mesh_info.Dimensions(strcmp(mesh_dims, 'node')).Length;
            fprintf('  Mesh nodes: %d\n', mesh_nodes);
        elseif ismember('nnode', mesh_dims)
            mesh_nodes = mesh_info.Dimensions(strcmp(mesh_dims, 'nnode')).Length;
            fprintf('  Mesh nodes: %d\n', mesh_nodes);
        else
            mesh_nodes = [];
            fprintf('  WARNING: Could not determine number of nodes in mesh file\n');
            warnings_count = warnings_count + 1;
        end
        
    catch ME
        fprintf('  ERROR: Cannot read mesh file: %s\n', ME.message);
        validation_passed = false;
        mesh_nodes = [];
    end
    fprintf('\n');
end

%% Validate water column file
water_info = [];
if ~isempty(water_file)
    fprintf('Validating water column file: %s\n', water_file);
    
    try
        % Check file exists
        if ~exist(water_file, 'file')
            fprintf('  ERROR: File does not exist\n');
            validation_passed = false;
        else
            water_info = ncinfo(water_file);
            
            % Check dimensions
            water_dims = {water_info.Dimensions.Name};
            required_dims = {'node', 'ksl', 'nplast'};
            
            for i = 1:length(required_dims)
                if ismember(required_dims{i}, water_dims)
                    dim_size = water_info.Dimensions(strcmp(water_dims, required_dims{i})).Length;
                    fprintf('  %s: %d\n', required_dims{i}, dim_size);
                    
                    % Store important dimensions
                    switch required_dims{i}
                        case 'node'
                            water_nodes = dim_size;
                        case 'ksl'
                            water_kb = dim_size;
                        case 'nplast'
                            water_nplast = dim_size;
                    end
                else
                    fprintf('  ERROR: Missing dimension ''%s''\n', required_dims{i});
                    validation_passed = false;
                end
            end
            
            % Check variables
            water_vars = {water_info.Variables.Name};
            if ismember('init_conc', water_vars)
                var_info = water_info.Variables(strcmp(water_vars, 'init_conc'));
                fprintf('  init_conc variable: dimensions [%s]\n', ...
                    strjoin({var_info.Dimensions.Name}, ', '));
                
                % Check data range
                try
                    conc_data = ncread(water_file, 'init_conc');
                    fprintf('  Concentration range: %.6f to %.6f\n', ...
                        min(conc_data(:)), max(conc_data(:)));
                    
                    if any(conc_data(:) < 0)
                        fprintf('  WARNING: Negative concentrations found\n');
                        warnings_count = warnings_count + 1;
                    end
                    
                catch ME
                    fprintf('  WARNING: Cannot read concentration data: %s\n', ME.message);
                    warnings_count = warnings_count + 1;
                end
            else
                fprintf('  ERROR: Missing variable ''init_conc''\n');
                validation_passed = false;
            end
            
            % Validate against mesh
            if exist('mesh_nodes', 'var') && ~isempty(mesh_nodes)
                if water_nodes ~= mesh_nodes
                    fprintf('  ERROR: Node count (%d) does not match mesh (%d)\n', ...
                        water_nodes, mesh_nodes);
                    validation_passed = false;
                end
            end
            
            % Validate against expected values
            if ~isempty(expected_kb) && exist('water_kb', 'var')
                if water_kb ~= expected_kb
                    fprintf('  WARNING: KB (%d) does not match expected (%d)\n', ...
                        water_kb, expected_kb);
                    warnings_count = warnings_count + 1;
                end
            end
        end
        
    catch ME
        fprintf('  ERROR: Cannot validate water file: %s\n', ME.message);
        validation_passed = false;
    end
    fprintf('\n');
end

%% Validate sediment file
sediment_info = [];
if ~isempty(sediment_file)
    fprintf('Validating sediment file: %s\n', sediment_file);
    
    try
        % Check file exists
        if ~exist(sediment_file, 'file')
            fprintf('  ERROR: File does not exist\n');
            validation_passed = false;
        else
            sediment_info = ncinfo(sediment_file);
            
            % Check dimensions
            sed_dims = {sediment_info.Dimensions.Name};
            required_dims = {'node', 'N_sed_layer', 'nplast'};
            
            for i = 1:length(required_dims)
                if ismember(required_dims{i}, sed_dims)
                    dim_size = sediment_info.Dimensions(strcmp(sed_dims, required_dims{i})).Length;
                    fprintf('  %s: %d\n', required_dims{i}, dim_size);
                    
                    % Store important dimensions
                    switch required_dims{i}
                        case 'node'
                            sed_nodes = dim_size;
                        case 'N_sed_layer'
                            sed_layers = dim_size;
                        case 'nplast'
                            sed_nplast = dim_size;
                    end
                else
                    fprintf('  ERROR: Missing dimension ''%s''\n', required_dims{i});
                    validation_passed = false;
                end
            end
            
            % Check variables
            sed_vars = {sediment_info.Variables.Name};
            required_vars = {'init_conc_sed', 'init_conc_sed_s', 'init_conc_sed_p'};
            
            for i = 1:length(required_vars)
                if ismember(required_vars{i}, sed_vars)
                    fprintf('  %s variable: OK\n', required_vars{i});
                else
                    fprintf('  ERROR: Missing variable ''%s''\n', required_vars{i});
                    validation_passed = false;
                end
            end
            
            % Check mass balance
            try
                conc_total = ncread(sediment_file, 'init_conc_sed');
                conc_stopped = ncread(sediment_file, 'init_conc_sed_s');
                conc_passing = ncread(sediment_file, 'init_conc_sed_p');
                
                fprintf('  Total concentration range: %.6f to %.6f\n', ...
                    min(conc_total(:)), max(conc_total(:)));
                
                % Check mass balance
                calculated_total = conc_stopped + conc_passing;
                tolerance = 1e-10;
                max_error = max(abs(conc_total(:) - calculated_total(:)));
                
                if max_error < tolerance
                    fprintf('  Mass balance: OK (max error: %.2e)\n', max_error);
                else
                    fprintf('  WARNING: Mass balance error: %.2e\n', max_error);
                    warnings_count = warnings_count + 1;
                end
                
                % Check for negative values
                if any(conc_total(:) < 0) || any(conc_stopped(:) < 0) || any(conc_passing(:) < 0)
                    fprintf('  WARNING: Negative concentrations found\n');
                    warnings_count = warnings_count + 1;
                end
                
            catch ME
                fprintf('  WARNING: Cannot validate concentration data: %s\n', ME.message);
                warnings_count = warnings_count + 1;
            end
            
            % Validate against mesh
            if exist('mesh_nodes', 'var') && ~isempty(mesh_nodes)
                if sed_nodes ~= mesh_nodes
                    fprintf('  ERROR: Node count (%d) does not match mesh (%d)\n', ...
                        sed_nodes, mesh_nodes);
                    validation_passed = false;
                end
            end
        end
        
    catch ME
        fprintf('  ERROR: Cannot validate sediment file: %s\n', ME.message);
        validation_passed = false;
    end
    fprintf('\n');
end

%% Cross-file validation
if ~isempty(water_file) && ~isempty(sediment_file) && ...
   exist('water_nodes', 'var') && exist('sed_nodes', 'var')
    
    fprintf('Cross-file validation:\n');
    
    % Check node consistency
    if water_nodes == sed_nodes
        fprintf('  Node count consistency: OK (%d nodes)\n', water_nodes);
    else
        fprintf('  ERROR: Node count mismatch (water: %d, sediment: %d)\n', ...
            water_nodes, sed_nodes);
        validation_passed = false;
    end
    
    % Check plastic type consistency
    if exist('water_nplast', 'var') && exist('sed_nplast', 'var')
        if water_nplast == sed_nplast
            fprintf('  Plastic type consistency: OK (%d types)\n', water_nplast);
        else
            fprintf('  ERROR: Plastic type count mismatch (water: %d, sediment: %d)\n', ...
                water_nplast, sed_nplast);
            validation_passed = false;
        end
    end
    
    % Validate against expected values
    if ~isempty(expected_nplast) && exist('water_nplast', 'var')
        if water_nplast ~= expected_nplast
            fprintf('  WARNING: NPLAST (%d) does not match expected (%d)\n', ...
                water_nplast, expected_nplast);
            warnings_count = warnings_count + 1;
        end
    end
    
    fprintf('\n');
end

%% Final summary
fprintf('=== Validation Summary ===\n');
if validation_passed && warnings_count == 0
    fprintf('✅ VALIDATION PASSED: All files are valid and ready for FVCOM\n');
elseif validation_passed && warnings_count > 0
    fprintf('⚠️  VALIDATION PASSED WITH WARNINGS: %d warnings found\n', warnings_count);
    fprintf('   Files should work with FVCOM but check warnings above\n');
else
    fprintf('❌ VALIDATION FAILED: Critical errors found\n');
    fprintf('   Files need to be fixed before using with FVCOM\n');
end

%% Recommendations
fprintf('\n=== Recommendations ===\n');
if exist('water_info', 'var') && ~isempty(water_info)
    fprintf('Water column file is ready for FVCOM use.\n');
    fprintf('Set: INIT_PLAST_CONC_FILE = ''%s''\n', water_file);
end

if exist('sediment_info', 'var') && ~isempty(sediment_info)
    fprintf('Sediment file is ready for FVCOM use.\n');
    fprintf('Set: INIT_PLAST_SED_CONC_FILE = ''%s''\n', sediment_file);
end

if exist('water_kb', 'var')
    fprintf('Ensure KB = %d in your FVCOM configuration\n', water_kb);
end

if exist('water_nplast', 'var')
    fprintf('Ensure NPLAST = %d in your FVCOM plastic configuration\n', water_nplast);
end

if exist('sed_layers', 'var')
    fprintf('Ensure N_sed_layer = %d in your FVCOM plastic configuration\n', sed_layers);
end

end

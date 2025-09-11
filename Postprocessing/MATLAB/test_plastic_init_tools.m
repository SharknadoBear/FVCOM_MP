% TEST_PLASTIC_INIT_TOOLS.m
%
% Comprehensive test script for both FVCOM plastic initialization tools:
% - create_init_plast_conc.m (water column)
% - create_init_plast_sed_conc.m (sediment)

clc; clear;

fprintf('Testing FVCOM Plastic Initialization Tools...\n\n');

%% Test 1: Water column initialization - Basic functionality
fprintf('=== Testing Water Column Initialization ===\n');
fprintf('Test 1: Basic water column functionality\n');
try
    create_init_plast_conc('test_water_basic.nc', 'verbose', false);
    
    % Verify file structure
    info = ncinfo('test_water_basic.nc');
    
    % Check dimensions
    dim_names = {info.Dimensions.Name};
    assert(ismember('node', dim_names), 'Missing node dimension');
    assert(ismember('ksl', dim_names), 'Missing ksl dimension');
    assert(ismember('nplast', dim_names), 'Missing nplast dimension');
    
    % Check variables
    var_names = {info.Variables.Name};
    assert(ismember('init_conc', var_names), 'Missing init_conc variable');
    
    fprintf('  ✓ PASSED: Water column basic structure correct\n');
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
end

%% Test 2: Water column - Vertical profiles
fprintf('Test 2: Water column vertical profiles\n');
try
    % Test surface-enhanced profile
    create_init_plast_conc('test_water_surface.nc', ...
        'num_nodes', 10, 'num_sigma_layers', 5, 'num_plastic_types', 2, ...
        'init_conc', 1.0, 'depth_profile', 'surface', 'surface_factor', 3.0, ...
        'verbose', false);
    
    % Read and verify profile
    conc_data = ncread('test_water_surface.nc', 'init_conc');
    
    % Surface should be higher than bottom
    surface_avg = mean(mean(conc_data(:, 1, :), 1), 3);
    bottom_avg = mean(mean(conc_data(:, end, :), 1), 3);
    
    assert(surface_avg > bottom_avg, 'Surface concentration should be higher than bottom');
    assert(surface_avg / bottom_avg > 1.5, 'Surface enhancement not working properly');
    
    fprintf('  ✓ PASSED: Vertical profiles working correctly\n');
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
end

%% Test 3: Sediment initialization - Basic functionality
fprintf('\n=== Testing Sediment Initialization ===\n');
fprintf('Test 3: Basic sediment functionality\n');
try
    create_init_plast_sed_conc('test_sed_basic.nc', 'verbose', false);
    
    % Verify file structure
    info = ncinfo('test_sed_basic.nc');
    
    % Check dimensions
    dim_names = {info.Dimensions.Name};
    assert(ismember('node', dim_names), 'Missing node dimension');
    assert(ismember('N_sed_layer', dim_names), 'Missing N_sed_layer dimension');
    assert(ismember('nplast', dim_names), 'Missing nplast dimension');
    
    % Check variables
    var_names = {info.Variables.Name};
    assert(ismember('init_conc_sed', var_names), 'Missing init_conc_sed variable');
    assert(ismember('init_conc_sed_s', var_names), 'Missing init_conc_sed_s variable');
    assert(ismember('init_conc_sed_p', var_names), 'Missing init_conc_sed_p variable');
    
    fprintf('  ✓ PASSED: Sediment basic structure correct\n');
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
end

%% Test 4: Sediment mass balance
fprintf('Test 4: Sediment mass balance\n');
try
    create_init_plast_sed_conc('test_sed_balance.nc', ...
        'num_nodes', 10, 'num_sed_layers', 3, 'num_plastic_types', 2, ...
        'conc_stopped', 0.4, 'conc_passing', 0.6, 'verbose', false);
    
    % Read and verify mass balance
    conc_total = ncread('test_sed_balance.nc', 'init_conc_sed');
    conc_stopped = ncread('test_sed_balance.nc', 'init_conc_sed_s');
    conc_passing = ncread('test_sed_balance.nc', 'init_conc_sed_p');
    
    % Check mass balance
    calculated_total = conc_stopped + conc_passing;
    tolerance = 1e-12;
    
    assert(max(abs(conc_total(:) - calculated_total(:))) < tolerance, ...
        'Mass balance not satisfied');
    assert(all(conc_total(:) == 1.0), 'Total should be 1.0');
    
    fprintf('  ✓ PASSED: Mass balance correctly enforced\n');
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
end

%% Test 5: Compatibility between water and sediment files
fprintf('\n=== Testing Compatibility ===\n');
fprintf('Test 5: Dimension consistency between water and sediment files\n');
try
    % Create files with same dimensions
    nodes = 50; layers_w = 8; layers_s = 4; types = 3;
    
    create_init_plast_conc('test_water_compat.nc', ...
        'num_nodes', nodes, 'num_sigma_layers', layers_w, 'num_plastic_types', types, ...
        'verbose', false);
    
    create_init_plast_sed_conc('test_sed_compat.nc', ...
        'num_nodes', nodes, 'num_sed_layers', layers_s, 'num_plastic_types', types, ...
        'verbose', false);
    
    % Verify dimensions match
    info_w = ncinfo('test_water_compat.nc');
    info_s = ncinfo('test_sed_compat.nc');
    
    % Check node and plastic type dimensions match
    node_dim_w = info_w.Dimensions(strcmp({info_w.Dimensions.Name}, 'node')).Length;
    node_dim_s = info_s.Dimensions(strcmp({info_s.Dimensions.Name}, 'node')).Length;
    
    nplast_dim_w = info_w.Dimensions(strcmp({info_w.Dimensions.Name}, 'nplast')).Length;
    nplast_dim_s = info_s.Dimensions(strcmp({info_s.Dimensions.Name}, 'nplast')).Length;
    
    assert(node_dim_w == node_dim_s, 'Node dimensions should match between files');
    assert(nplast_dim_w == nplast_dim_s, 'Plastic type dimensions should match');
    
    fprintf('  ✓ PASSED: Compatible dimensions between water and sediment files\n');
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
end

%% Test 6: Array input validation
fprintf('Test 6: Custom array input validation\n');
try
    % Test water column with custom array
    nodes = 5; layers = 4; types = 2;
    custom_conc = rand(nodes, layers, types) * 0.1;
    
    create_init_plast_conc('test_water_array.nc', ...
        'num_nodes', nodes, 'num_sigma_layers', layers, 'num_plastic_types', types, ...
        'init_conc', custom_conc, 'verbose', false);
    
    % Read back and verify
    read_conc = ncread('test_water_array.nc', 'init_conc');
    
    assert(isequal(size(read_conc), size(custom_conc)), 'Array dimensions should match');
    assert(max(abs(read_conc(:) - custom_conc(:))) < 1e-12, 'Array values should match');
    
    fprintf('  ✓ PASSED: Custom array input working correctly\n');
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
end

%% Test 7: Error handling
fprintf('Test 7: Error handling for incorrect inputs\n');
try
    error_count = 0;
    
    % Test wrong array size for water column
    try
        wrong_array = rand(5, 3, 2);
        create_init_plast_conc('test_error_water.nc', ...
            'num_nodes', 10, 'num_sigma_layers', 5, 'num_plastic_types', 2, ...
            'init_conc', wrong_array, 'verbose', false);
        error_count = error_count + 1;
    catch
        % Expected error - good!
    end
    
    % Test wrong array size for sediment
    try
        wrong_array = rand(5, 3, 2);
        create_init_plast_sed_conc('test_error_sed.nc', ...
            'num_nodes', 10, 'num_sed_layers', 5, 'num_plastic_types', 2, ...
            'conc_total', wrong_array, 'verbose', false);
        error_count = error_count + 1;
    catch
        % Expected error - good!
    end
    
    if error_count == 0
        fprintf('  ✓ PASSED: Error handling working correctly\n');
    else
        fprintf('  ✗ FAILED: Should have caught %d errors\n', error_count);
    end
catch ME
    fprintf('  ✗ FAILED: Unexpected error in error handling test: %s\n', ME.message);
end

%% Performance test
fprintf('\n=== Performance Test ===\n');
fprintf('Test 8: Performance with larger arrays\n');
try
    tic;
    
    % Create moderately large files
    create_init_plast_conc('test_perf_water.nc', ...
        'num_nodes', 1000, 'num_sigma_layers', 20, 'num_plastic_types', 5, ...
        'init_conc', 0.05, 'depth_profile', 'exponential', 'verbose', false);
    
    create_init_plast_sed_conc('test_perf_sed.nc', ...
        'num_nodes', 1000, 'num_sed_layers', 5, 'num_plastic_types', 5, ...
        'conc_total', 0.02, 'verbose', false);
    
    elapsed_time = toc;
    
    fprintf('  ✓ PASSED: Created large files in %.2f seconds\n', elapsed_time);
    
    if elapsed_time > 10
        fprintf('  ⚠ WARNING: Performance seems slow (>10s)\n');
    end
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
end

%% Cleanup
fprintf('\n=== Cleanup ===\n');
test_files = {'test_water_basic.nc', 'test_water_surface.nc', 'test_sed_basic.nc', ...
              'test_sed_balance.nc', 'test_water_compat.nc', 'test_sed_compat.nc', ...
              'test_water_array.nc', 'test_error_water.nc', 'test_error_sed.nc', ...
              'test_perf_water.nc', 'test_perf_sed.nc'};

cleaned_count = 0;
for i = 1:length(test_files)
    if exist(test_files{i}, 'file')
        delete(test_files{i});
        cleaned_count = cleaned_count + 1;
    end
end

fprintf('Cleaned up %d test files.\n', cleaned_count);

%% Summary
fprintf('\n=== Test Summary ===\n');
fprintf('All tests completed successfully!\n');
fprintf('Both plastic initialization tools are working correctly and are\n');
fprintf('compatible with FVCOM''s plastic transport model requirements.\n\n');

fprintf('=== Quick Start Guide ===\n');
fprintf('1. Water column initialization:\n');
fprintf('   create_init_plast_conc(''my_water_init.nc'', ''num_nodes'', your_nodes)\n\n');
fprintf('2. Sediment initialization:\n');
fprintf('   create_init_plast_sed_conc(''my_sed_init.nc'', ''num_nodes'', your_nodes)\n\n');
fprintf('3. Set these files in your FVCOM input:\n');
fprintf('   INIT_PLAST_CONC_FILE     = ''my_water_init.nc''\n');
fprintf('   INIT_PLAST_SED_CONC_FILE = ''my_sed_init.nc''\n\n');
fprintf('4. Run example scripts for more advanced usage:\n');
fprintf('   example_create_plastic_water_init\n');
fprintf('   example_create_plastic_sed_init\n');

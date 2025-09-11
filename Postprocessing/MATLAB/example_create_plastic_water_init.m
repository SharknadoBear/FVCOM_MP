% EXAMPLE_CREATE_PLASTIC_WATER_INIT.m
% 
% Example script demonstrating how to use create_init_plast_conc.m
% to create NetCDF initialization files for FVCOM plastic transport model
% (water column concentrations)

clc; clear; close all;

fprintf('=== FVCOM Plastic Water Column Initialization Examples ===\n\n');

%% Example 1: Basic usage with default values
fprintf('Example 1: Creating file with default parameters\n');
create_init_plast_conc('example1_water_default.nc');

%% Example 2: Custom dimensions with uniform distribution
fprintf('\nExample 2: Custom dimensions and uniform concentration\n');
create_init_plast_conc('example2_water_uniform.nc', ...
    'num_nodes', 200, ...
    'num_sigma_layers', 15, ...
    'num_plastic_types', 3, ...
    'init_conc', 0.08, ...          % kg/m³ uniform concentration
    'plastic_names', {'Microplastics', 'Bottles', 'Bags'});

%% Example 3: Surface-enhanced distribution
fprintf('\nExample 3: Surface-enhanced plastic distribution\n');
create_init_plast_conc('example3_water_surface.nc', ...
    'num_nodes', 150, ...
    'num_sigma_layers', 12, ...
    'num_plastic_types', 4, ...
    'init_conc', 0.05, ...          % Base concentration
    'depth_profile', 'surface', ... % Surface enhancement
    'surface_factor', 3.0, ...      % 3x higher at surface
    'plastic_names', {'PET_fragments', 'PE_pellets', 'PS_foam', 'Mixed_debris'});

%% Example 4: Exponential decay from surface (typical for buoyant plastics)
fprintf('\nExample 4: Exponential decay profile (buoyant plastics)\n');
create_init_plast_conc('example4_water_exponential.nc', ...
    'num_nodes', 180, ...
    'num_sigma_layers', 20, ...
    'num_plastic_types', 2, ...
    'init_conc', 0.1, ...           % Surface concentration
    'depth_profile', 'exponential', ...
    'decay_scale', 0.2, ...         % Rapid decay (mostly near surface)
    'plastic_names', {'Floating_plastics', 'Neutrally_buoyant'});

%% Example 5: Spatially and vertically varying concentrations
fprintf('\nExample 5: Custom spatially and vertically varying concentrations\n');

% Define dimensions
num_nodes = 100;
num_sigma_layers = 10;
num_plastic_types = 3;

% Create realistic concentration distribution
init_conc = zeros(num_nodes, num_sigma_layers, num_plastic_types);

% Simulate different plastic types with different vertical distributions
for inode = 1:num_nodes
    % Create some spatial variation (could be based on distance from source, bathymetry, etc.)
    spatial_factor = 0.5 + 0.5 * sin(2*pi*inode/num_nodes); % Sinusoidal variation
    
    for k = 1:num_sigma_layers
        % Sigma coordinate (0 = surface, 1 = bottom)
        sigma = (k - 1) / (num_sigma_layers - 1);
        
        for iplastic = 1:num_plastic_types
            switch iplastic
                case 1  % Floating microplastics - surface concentration
                    profile = exp(-sigma / 0.15) * 0.08;
                case 2  % Neutrally buoyant fragments - more uniform
                    profile = 0.04 * (1 - 0.5 * sigma);
                case 3  % Sinking particles - bottom enhanced
                    profile = 0.02 * (0.5 + sigma);
            end
            
            init_conc(inode, k, iplastic) = profile * spatial_factor;
        end
    end
end

create_init_plast_conc('example5_water_custom.nc', ...
    'num_nodes', num_nodes, ...
    'num_sigma_layers', num_sigma_layers, ...
    'num_plastic_types', num_plastic_types, ...
    'init_conc', init_conc, ...
    'plastic_names', {'Floating_microplastics', 'Neutral_fragments', 'Sinking_particles'});

%% Example 6: Pollution plume scenario
fprintf('\nExample 6: River discharge plume with plastic pollution\n');

num_nodes = 250;
num_sigma_layers = 8;
num_plastic_types = 2;

% Simulate river discharge at node 25 with surface plume
source_node = 25;
plume_extent = 50; % nodes downstream

init_conc = zeros(num_nodes, num_sigma_layers, num_plastic_types);

for inode = 1:num_nodes
    % Distance from river source
    distance = abs(inode - source_node);
    
    % River plume mixing (exponential decay downstream)
    if inode >= source_node && distance <= plume_extent
        mixing_factor = exp(-distance / 20); % Decay length = 20 nodes
    else
        mixing_factor = 0.1; % Background level
    end
    
    for k = 1:num_sigma_layers
        sigma = (k - 1) / (num_sigma_layers - 1);
        
        % River water is fresher and tends to stay near surface
        vertical_mixing = exp(-sigma / 0.25); % Surface stratified
        
        for iplastic = 1:num_plastic_types
            base_conc = [0.15, 0.08]; % Different source concentrations
            init_conc(inode, k, iplastic) = base_conc(iplastic) * mixing_factor * vertical_mixing;
        end
    end
end

create_init_plast_conc('example6_water_plume.nc', ...
    'num_nodes', num_nodes, ...
    'num_sigma_layers', num_sigma_layers, ...
    'num_plastic_types', num_plastic_types, ...
    'init_conc', init_conc, ...
    'plastic_names', {'River_plastics', 'Urban_runoff'});

%% Display summary
fprintf('\n=== Summary of Created Files ===\n');
files = {'example1_water_default.nc', 'example2_water_uniform.nc', 'example3_water_surface.nc', ...
         'example4_water_exponential.nc', 'example5_water_custom.nc', 'example6_water_plume.nc'};

for i = 1:length(files)
    if exist(files{i}, 'file')
        info = ncinfo(files{i});
        fprintf('%s:\n', files{i});
        fprintf('  Dimensions: %d nodes, %d sigma layers, %d plastic types\n', ...
            info.Dimensions(1).Length, info.Dimensions(2).Length, info.Dimensions(3).Length);
        
        % Read and display concentration ranges
        init_conc = ncread(files{i}, 'init_conc');
        fprintf('  Concentration range: %.6f to %.6f kg/m³\n', ...
            min(init_conc(:)), max(init_conc(:)));
        
        % Surface vs bottom comparison
        surface_avg = mean(init_conc(:, 1, :), 'all');
        bottom_avg = mean(init_conc(:, end, :), 'all');
        if bottom_avg > 0
            fprintf('  Surface/bottom ratio: %.2f\n', surface_avg / bottom_avg);
        end
    end
end

fprintf('\nAll water column initialization files created successfully!\n');

%% Visualization example
fprintf('\n=== Creating visualizations ===\n');

try
    % Load data from surface-enhanced example
    filename = 'example3_water_surface.nc';
    if exist(filename, 'file')
        conc_data = ncread(filename, 'init_conc');
        
        figure('Name', 'Plastic Water Column Concentration Profiles', 'Position', [100 100 1200 800]);
        
        % Plot 1: Vertical profiles for different plastic types
        subplot(2,2,1);
        sigma_levels = 1:size(conc_data, 2);
        for itype = 1:min(3, size(conc_data, 3))
            profile = squeeze(mean(conc_data(:, :, itype), 1));
            plot(profile, sigma_levels, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
            hold on;
        end
        set(gca, 'YDir', 'reverse');
        xlabel('Average Concentration (kg/m³)');
        ylabel('Sigma Level (1=surface, kb=bottom)');
        title('Vertical Concentration Profiles');
        legend('Type 1', 'Type 2', 'Type 3', 'Location', 'best');
        grid on;
        
        % Plot 2: Surface distribution
        subplot(2,2,2);
        node_ids = 1:size(conc_data, 1);
        surface_conc = squeeze(conc_data(:, 1, 1)); % First plastic type at surface
        plot(node_ids, surface_conc, 'b-', 'LineWidth', 2);
        xlabel('Node ID');
        ylabel('Surface Concentration (kg/m³)');
        title('Surface Layer Distribution (Type 1)');
        grid on;
        
        % Plot 3: Depth-integrated concentration
        subplot(2,2,3);
        integrated_conc = squeeze(sum(conc_data, 2)); % Sum over all layers
        plot(node_ids, integrated_conc(:, 1), 'r-', 'LineWidth', 2);
        xlabel('Node ID');
        ylabel('Integrated Concentration (kg/m²)');
        title('Depth-Integrated Concentration');
        grid on;
        
        % Plot 4: 3D visualization of first nodes
        subplot(2,2,4);
        [X, Y] = meshgrid(1:min(20, size(conc_data, 1)), sigma_levels);
        Z = squeeze(conc_data(1:size(X,1), :, 1))';
        surf(X, Y, Z);
        set(gca, 'YDir', 'reverse');
        xlabel('Node ID');
        ylabel('Sigma Level');
        zlabel('Concentration (kg/m³)');
        title('3D Concentration Distribution');
        colorbar;
        
        fprintf('Visualization created successfully!\n');
    end
catch ME
    fprintf('Visualization failed: %s\n', ME.message);
    fprintf('(This is optional - NetCDF files were created successfully)\n');
end

fprintf('\n=== Usage Instructions ===\n');
fprintf('To use these files with FVCOM:\n');
fprintf('1. Set INIT_PLAST_CONC_FILE in your FVCOM input file\n');
fprintf('2. Ensure the number of nodes matches your mesh\n');
fprintf('3. Ensure num_sigma_layers matches KB in your model\n');
fprintf('4. Ensure num_plastic_types matches NPLAST in your model\n');

% EXTRACT_SPATIAL_SNAPSHOTS_B_C_DAY90
% Extract spatial snapshot and station time-series fields from the case-b
% and case-c FVCOM-MP NetCDF outputs centred on simulation day 90 (Memo 05).
%
% For each node and each record within ±12 h of day 90 this script saves:
%   mp1_surf        -- surface-layer (k=1) mp1 concentration       [kg/m3]
%   mp1_bot         -- bottom-layer (k=siglay) mp1 concentration    [kg/m3]
%   mp1_davg        -- depth-averaged (sigma-weighted) mp1          [kg/m3]
%   mp1_agg_surf    -- aggregated mp1, surface layer                [kg/m3]
%   mp1_agg_bot     -- aggregated mp1, bottom layer                 [kg/m3]
%   mp1_agg_davg    -- aggregated mp1, depth-averaged               [kg/m3]
%   mp1_dis_surf    -- disaggregated mp1, surface layer             [kg/m3]
%   mp1_dis_bot     -- disaggregated mp1, bottom layer              [kg/m3]
%   mp1_dis_davg    -- disaggregated mp1, depth-averaged            [kg/m3]
%   bot_mass        -- total bed plastic mass density               [kg/m2]
%   bot_mass_agg    -- aggregated bed plastic mass density          [kg/m2]
%   bot_mass_dis    -- disaggregated bed plastic mass density       [kg/m2]
%   sed_surf        -- surface-layer coarse_sand_1 concentration    [g/L]
%   sed_bot         -- bottom-layer coarse_sand_1 concentration     [g/L]
%   sed_davg        -- depth-averaged coarse_sand_1                 [g/L]
%   sed_bedfrac     -- coarse_sand_1 surface bed fraction           [-]
%
% In addition, full time-series profiles at 8 monitoring stations are
% extracted by nearest-node lookup.
%
% Default NetCDF inputs:
%   <test_root>/OUTPUT_b/backup_data/waterPACT_b_0009.nc
%   <test_root>/OUTPUT_c/backup_data/waterPACT_c_0009.nc
%
% Default MAT outputs:
%   <test_root>/OUTPUT_b/spatial_snapshots_b_day90.mat  (v7.3 / HDF5)
%   <test_root>/OUTPUT_c/spatial_snapshots_c_day90.mat  (v7.3 / HDF5)
%
% HPC usage (run from any directory on Kestrel data-analysis node):
%   cd /scratch/yhuang168/waterPACT_MP_floc
%   matlab -nodisplay -batch "run('MATLAB/extract_spatial_snapshots_b_c_day90.m')"
%
% Optional environment overrides:
%   WATERPACT_CASE_IDS      Comma-separated list of case IDs.  Default: b,c
%   WATERPACT_OUTPUT_ROOT   Root dir containing OUTPUT_b, OUTPUT_c.
%   WATERPACT_OUTPUT_DIRS   Comma-separated explicit output directories.
%   WATERPACT_NC_FILES      Comma-separated exact NC file paths.
%   WATERPACT_TARGET_DAY    Centre of the extraction window in simulation days.
%                           Default: 90
%   WATERPACT_LAST_NHOURS   Width of the extraction window in hours (centred on
%                           WATERPACT_TARGET_DAY).  Default: 24

%% -------------------------------------------------------------------------
%  User Settings
%  Absolute Kestrel paths are set as defaults.  Override via environment
%  variables for portability to other systems.
% --------------------------------------------------------------------------
case_ids     = {'b', 'c'};
output_root  = '/kfs3/scratch/yhuang168/waterPACT_MP_floc';
output_dirs  = { ...
    '/kfs3/scratch/yhuang168/waterPACT_MP_floc/OUTPUT_b'; ...
    '/kfs3/scratch/yhuang168/waterPACT_MP_floc/OUTPUT_c' };
nc_files     = { ...
    '/kfs3/scratch/yhuang168/waterPACT_MP_floc/OUTPUT_b/backup_data/waterPACT_b_0009.nc'; ...
    '/kfs3/scratch/yhuang168/waterPACT_MP_floc/OUTPUT_c/backup_data/waterPACT_c_0009.nc' };
target_day   = 90;
last_n_hours = 24;

%% Environment overrides
case_ids_env = getenv('WATERPACT_CASE_IDS');
if ~isempty(case_ids_env)
    case_ids = strtrim(strsplit(case_ids_env, ','));
end
output_root  = getenv_default('WATERPACT_OUTPUT_ROOT', output_root);
output_dirs_env = getenv('WATERPACT_OUTPUT_DIRS');
if ~isempty(output_dirs_env)
    output_dirs = strtrim(strsplit(output_dirs_env, ','));
end
nc_files_env = getenv('WATERPACT_NC_FILES');
if ~isempty(nc_files_env)
    nc_files = strtrim(strsplit(nc_files_env, ','));
end
target_day   = getenv_number('WATERPACT_TARGET_DAY',  target_day);
last_n_hours = getenv_number('WATERPACT_LAST_NHOURS', last_n_hours);

%% Resolve script root
% WATERPACT_TEST_ROOT overrides the auto-detected root.  Set this when
% MATLAB is launched via run() with a relative path and mfilename may
% return an empty string (e.g. on Kestrel data-analysis nodes).
%   export WATERPACT_TEST_ROOT=/kfs3/scratch/yhuang168/waterPACT_MP_floc
test_root_env = getenv_default('WATERPACT_TEST_ROOT', '');
if ~isempty(test_root_env)
    test_root = test_root_env;
else
    script_dir = fileparts(mfilename('fullpath'));
    if isempty(script_dir)
        script_dir = pwd;
    end
    test_root = fileparts(script_dir);
end

if isempty(output_root)
    output_root = test_root;
end

if ~isempty(output_dirs) && numel(output_dirs) ~= numel(case_ids)
    error('WATERPACT_OUTPUT_DIRS must match number of cases.');
end
if ~isempty(nc_files) && numel(nc_files) ~= numel(case_ids)
    error('WATERPACT_NC_FILES must match number of cases.');
end

%% -------------------------------------------------------------------------
%  Station Table  (8 monitoring points, Delaware River / Bay)
% --------------------------------------------------------------------------
station_names = { ...
    'Chesapeake Canal'; ...
    'Main stream 1'; ...
    'Main stream 2'; ...
    'Main stream 3'; ...
    'Main stream 4'; ...
    'Main stream 5'; ...
    'Upstream Delaware'; ...
    'Upstream Schuylkill' };

station_lat = [39.5350; 39.5644; 39.7826; 39.8788; 40.0147; 40.1244; 40.1956; 39.9256];
station_lon = [-75.7867; -75.5417; -75.4562; -75.1951; -75.0386; -74.8056; -74.7594; -75.2121];
n_stations  = numel(station_names);

%% -------------------------------------------------------------------------
%  Loop over cases
% --------------------------------------------------------------------------
fprintf('Starting day-90 spatial snapshot extraction for %d case(s).\n\n', numel(case_ids));

for icase = 1:numel(case_ids)
    case_id = char(case_ids{icase});

    % Resolve output directory
    if isempty(output_dirs)
        case_output_dir = fullfile(output_root, ['OUTPUT_' case_id]);
    else
        case_output_dir = char(output_dirs{icase});
    end

    % Resolve NetCDF file
    if isempty(nc_files)
        nc_file = fullfile(case_output_dir, 'backup_data', ...
            sprintf('waterPACT_%s_0009.nc', case_id));
    else
        nc_file_raw = char(nc_files{icase});
        if is_absolute_path(nc_file_raw)
            nc_file = nc_file_raw;
        else
            nc_file = fullfile(case_output_dir, 'backup_data', nc_file_raw);
        end
    end

    out_mat = fullfile(case_output_dir, ['spatial_snapshots_' case_id '_day90.mat']);

    fprintf('=== Case %s  (day-90 snapshot) ===\n', case_id);
    fprintf('  NetCDF : %s\n', nc_file);
    fprintf('  MAT out: %s\n', out_mat);

    if ~isfile(nc_file)
        error('NetCDF file not found: %s', nc_file);
    end

    extract_and_save(case_id, nc_file, out_mat, target_day, last_n_hours, ...
        station_names, station_lat, station_lon, n_stations);

    fprintf('\n');
end

fprintf('All cases complete.\n');

%% =========================================================================
%  Core extraction function
% ==========================================================================

function extract_and_save(case_id, nc_file, out_mat, target_day, last_n_hours, ...
        station_names, station_lat, station_lon, n_stations)

%% Read grid once
fprintf('  Reading grid ...\n');
lon    = double(ncread(nc_file, 'lon'));      lon    = lon(:);
lat    = double(ncread(nc_file, 'lat'));      lat    = lat(:);
lonc   = double(ncread(nc_file, 'lonc'));     lonc   = lonc(:);
latc   = double(ncread(nc_file, 'latc'));     latc   = latc(:);
h      = double(ncread(nc_file, 'h'));        h      = h(:);
art1   = double(ncread(nc_file, 'art1'));     art1   = art1(:);
nv     = double(ncread(nc_file, 'nv'));       % (nele, 3)
siglev = double(ncread(nc_file, 'siglev'));   % (node, siglev)

n_nodes  = numel(lon);
dzfrac   = abs(diff(siglev, 1, 2));           % (node, siglay)
n_siglay = size(dzfrac, 2);

%% Check which optional variables are present
vars_avail  = inspect_variables(nc_file);
has_agg     = has_var(vars_avail, 'mp1_agg');
has_dis     = has_var(vars_avail, 'mp1_dis');
has_bot_agg = has_var(vars_avail, 'bot_mass_agg_mp1');
has_bot_dis = has_var(vars_avail, 'bot_mass_dis_mp1');
has_sed     = has_var(vars_avail, 'coarse_sand_1');
has_bedfrac = has_var(vars_avail, 'coarse_sand_1_bedfrac');

fprintf('  mp1_agg/dis: %d/%d  bot_agg/dis: %d/%d  sed: %d  bedfrac: %d\n', ...
    has_agg, has_dis, has_bot_agg, has_bot_dis, has_sed, has_bedfrac);

%% Station nearest-node lookup
fprintf('  Locating nearest nodes for %d stations ...\n', n_stations);
station_node_idx     = zeros(n_stations, 1, 'int32');
station_node_lat     = zeros(n_stations, 1);
station_node_lon     = zeros(n_stations, 1);
station_node_dist_km = zeros(n_stations, 1);

for is = 1:n_stations
    d_km = haversine_km(lat, lon, station_lat(is), station_lon(is));
    [min_d, idx] = min(d_km);
    station_node_idx(is)     = int32(idx);
    station_node_lat(is)     = lat(idx);
    station_node_lon(is)     = lon(idx);
    station_node_dist_km(is) = min_d;
    fprintf('    %-22s  (%.4f, %.4f) -> node %d  (%.4f, %.4f)  %.3f km\n', ...
        station_names{is}, station_lat(is), station_lon(is), ...
        idx, lat(idx), lon(idx), min_d);
end

%% Select records within ±(last_n_hours/2) h of target_day
all_time = double(ncread(nc_file, 'time'));   % absolute days (MJD)
n_total  = numel(all_time);

% Convert to simulation-relative days: file _000N starts at (N-1)*10 sim-days.
tok = regexp(nc_file, '_(\d{4})\.nc', 'tokens', 'once');
if ~isempty(tok)
    file_num = str2double(tok{1});
else
    file_num = 1;
    warning('Could not parse file number from %s; assuming file 1.', nc_file);
end
rel_time = all_time - all_time(1) + (file_num - 1) * 10;

half_win = last_n_hours / 2 / 24.0;
snap_idx = find(rel_time >= (target_day - half_win) & ...
                rel_time <= (target_day + half_win));

if isempty(snap_idx)
    error('No records within %.1f h of day %.1f (available range: %.4f .. %.4f sim-days).', ...
        last_n_hours/2, target_day, rel_time(1), rel_time(end));
end
n_snaps   = numel(snap_idx);
snap_time = rel_time(snap_idx);   % simulation-relative days

fprintf('  Snapshot window: day %.1f ± %.1f h -> %d records (%.4f .. %.4f sim-days)\n', ...
    target_day, last_n_hours/2, n_snaps, snap_time(1), snap_time(end));

%% Allocate snapshot arrays  (node x n_snaps)
mk = @(r,c) zeros(r, c, 'single');
mp1_surf     = mk(n_nodes, n_snaps);
mp1_bot      = mk(n_nodes, n_snaps);
mp1_davg     = mk(n_nodes, n_snaps);
bot_mass_out = mk(n_nodes, n_snaps);

if has_agg
    mp1_agg_surf = mk(n_nodes, n_snaps);
    mp1_agg_bot  = mk(n_nodes, n_snaps);
    mp1_agg_davg = mk(n_nodes, n_snaps);
end
if has_dis
    mp1_dis_surf = mk(n_nodes, n_snaps);
    mp1_dis_bot  = mk(n_nodes, n_snaps);
    mp1_dis_davg = mk(n_nodes, n_snaps);
end
if has_bot_agg
    bot_mass_agg_out = mk(n_nodes, n_snaps);
end
if has_bot_dis
    bot_mass_dis_out = mk(n_nodes, n_snaps);
end
if has_sed
    sed_surf = mk(n_nodes, n_snaps);
    sed_bot  = mk(n_nodes, n_snaps);
    sed_davg = mk(n_nodes, n_snaps);
end
if has_bedfrac
    sed_bedfrac_out = mk(n_nodes, n_snaps);
end

%% Station time series  (n_stations x n_siglay x n_total)
station_time    = rel_time;           % simulation-relative days
station_mp1     = zeros(n_stations, n_siglay, n_total, 'single');
if has_agg
    station_mp1_agg = zeros(n_stations, n_siglay, n_total, 'single');
end
if has_dis
    station_mp1_dis = zeros(n_stations, n_siglay, n_total, 'single');
end
if has_sed
    station_sed     = zeros(n_stations, n_siglay, n_total, 'single');
end

%% Read loop -- station profiles (all records)
fprintf('  Extracting station profiles for all %d records ...\n', n_total);
for trec = 1:n_total
    mp1_3d = single(ncread(nc_file, 'mp1', [1 1 trec], [Inf Inf 1]));
    if has_agg
        agg_3d = single(ncread(nc_file, 'mp1_agg', [1 1 trec], [Inf Inf 1]));
    end
    if has_dis
        dis_3d = single(ncread(nc_file, 'mp1_dis', [1 1 trec], [Inf Inf 1]));
    end
    if has_sed
        sed_3d = single(ncread(nc_file, 'coarse_sand_1', [1 1 trec], [Inf Inf 1]));
    end
    for is = 1:n_stations
        nd = station_node_idx(is);
        station_mp1(is, :, trec) = mp1_3d(nd, :);
        if has_agg; station_mp1_agg(is, :, trec) = agg_3d(nd, :); end
        if has_dis; station_mp1_dis(is, :, trec) = dis_3d(nd, :); end
        if has_sed; station_sed(is,  :, trec)    = sed_3d(nd, :); end
    end
    if mod(trec, 50) == 0 || trec == n_total
        fprintf('    station profiles: %d / %d\n', trec, n_total);
    end
end

%% Read loop -- spatial snapshots (window around day 90)
fprintf('  Extracting %d spatial snapshots ...\n', n_snaps);
for isnap = 1:n_snaps
    trec = snap_idx(isnap);

    zeta  = double(ncread(nc_file, 'zeta', [1 trec], [Inf 1]));
    depth = max(h + zeta(:), 0.0);
    base_vol = art1 .* depth;
    vol_lay  = bsxfun(@times, base_vol, dzfrac);   % (node, siglay)

    % mp1
    mp1_3d = double(ncread(nc_file, 'mp1', [1 1 trec], [Inf Inf 1]));
    mp1_surf(:, isnap) = single(mp1_3d(:, 1));
    mp1_bot(:,  isnap) = single(mp1_3d(:, n_siglay));
    mp1_davg(:, isnap) = single(depth_avg(mp1_3d, vol_lay));

    % bot_mass
    bm = double(ncread(nc_file, 'bot_mass_mp1', [1 trec], [Inf 1]));
    bot_mass_out(:, isnap) = single(bm(:));

    % mp1_agg
    if has_agg
        agg_3d = double(ncread(nc_file, 'mp1_agg', [1 1 trec], [Inf Inf 1]));
        mp1_agg_surf(:, isnap) = single(agg_3d(:, 1));
        mp1_agg_bot(:,  isnap) = single(agg_3d(:, n_siglay));
        mp1_agg_davg(:, isnap) = single(depth_avg(agg_3d, vol_lay));
    end

    % mp1_dis
    if has_dis
        dis_3d = double(ncread(nc_file, 'mp1_dis', [1 1 trec], [Inf Inf 1]));
        mp1_dis_surf(:, isnap) = single(dis_3d(:, 1));
        mp1_dis_bot(:,  isnap) = single(dis_3d(:, n_siglay));
        mp1_dis_davg(:, isnap) = single(depth_avg(dis_3d, vol_lay));
    end

    % bot_mass_agg / dis
    if has_bot_agg
        ba = double(ncread(nc_file, 'bot_mass_agg_mp1', [1 trec], [Inf 1]));
        bot_mass_agg_out(:, isnap) = single(ba(:));
    end
    if has_bot_dis
        bd = double(ncread(nc_file, 'bot_mass_dis_mp1', [1 trec], [Inf 1]));
        bot_mass_dis_out(:, isnap) = single(bd(:));
    end

    % coarse_sand_1
    if has_sed
        sed_3d = double(ncread(nc_file, 'coarse_sand_1', [1 1 trec], [Inf Inf 1]));
        sed_surf(:, isnap) = single(sed_3d(:, 1));
        sed_bot(:,  isnap) = single(sed_3d(:, n_siglay));
        sed_davg(:, isnap) = single(depth_avg(sed_3d, vol_lay));
    end
    if has_bedfrac
        bf = double(ncread(nc_file, 'coarse_sand_1_bedfrac', [1 trec], [Inf 1]));
        sed_bedfrac_out(:, isnap) = single(bf(:));
    end

    if mod(isnap, 6) == 0 || isnap == n_snaps
        fprintf('    snap %d / %d  (day %.4f)\n', isnap, n_snaps, snap_time(isnap));
    end
end

%% Pack output struct
out = struct();
out.config.case_id      = case_id;
out.config.nc_file      = nc_file;
out.config.out_mat      = out_mat;
out.config.target_day   = target_day;
out.config.last_n_hours = last_n_hours;
out.config.created_by   = mfilename;
out.config.created_on   = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

out.grid.lon    = single(lon);
out.grid.lat    = single(lat);
out.grid.lonc   = single(lonc);
out.grid.latc   = single(latc);
out.grid.h      = single(h);
out.grid.art1   = single(art1);
out.grid.nv     = int32(nv);
out.grid.dzfrac = single(dzfrac);
out.grid.n_nodes  = n_nodes;
out.grid.n_siglay = n_siglay;

out.snaps.snap_idx  = int32(snap_idx(:));
out.snaps.snap_time = snap_time(:);
out.snaps.n_snaps   = n_snaps;

out.snaps.mp1_surf    = mp1_surf;
out.snaps.mp1_bot     = mp1_bot;
out.snaps.mp1_davg    = mp1_davg;
out.snaps.bot_mass    = bot_mass_out;

if has_agg
    out.snaps.mp1_agg_surf = mp1_agg_surf;
    out.snaps.mp1_agg_bot  = mp1_agg_bot;
    out.snaps.mp1_agg_davg = mp1_agg_davg;
end
if has_dis
    out.snaps.mp1_dis_surf = mp1_dis_surf;
    out.snaps.mp1_dis_bot  = mp1_dis_bot;
    out.snaps.mp1_dis_davg = mp1_dis_davg;
end
if has_bot_agg; out.snaps.bot_mass_agg = bot_mass_agg_out; end
if has_bot_dis; out.snaps.bot_mass_dis = bot_mass_dis_out; end
if has_sed
    out.snaps.sed_surf = sed_surf;
    out.snaps.sed_bot  = sed_bot;
    out.snaps.sed_davg = sed_davg;
end
if has_bedfrac; out.snaps.sed_bedfrac = sed_bedfrac_out; end

out.snaps.vars_available = vars_avail;
out.snaps.units.mp1_surf      = 'kg/m3';
out.snaps.units.mp1_bot       = 'kg/m3';
out.snaps.units.mp1_davg      = 'kg/m3 (sigma-weighted depth-average)';
out.snaps.units.mp1_agg_surf  = 'kg/m3';
out.snaps.units.mp1_agg_bot   = 'kg/m3';
out.snaps.units.mp1_agg_davg  = 'kg/m3 (sigma-weighted depth-average)';
out.snaps.units.mp1_dis_surf  = 'kg/m3';
out.snaps.units.mp1_dis_bot   = 'kg/m3';
out.snaps.units.mp1_dis_davg  = 'kg/m3 (sigma-weighted depth-average)';
out.snaps.units.bot_mass      = 'kg/m2';
out.snaps.units.bot_mass_agg  = 'kg/m2';
out.snaps.units.bot_mass_dis  = 'kg/m2';
out.snaps.units.sed_surf      = 'g/L';
out.snaps.units.sed_bot       = 'g/L';
out.snaps.units.sed_davg      = 'g/L (sigma-weighted depth-average)';
out.snaps.units.sed_bedfrac   = 'dimensionless (surface bed fraction)';

out.stations.names        = station_names;
out.stations.target_lat   = station_lat;
out.stations.target_lon   = station_lon;
out.stations.node_idx     = station_node_idx;
out.stations.node_lat     = station_node_lat;
out.stations.node_lon     = station_node_lon;
out.stations.node_dist_km = station_node_dist_km;

out.stations.time    = station_time(:);
out.stations.mp1     = station_mp1;
out.stations.units.mp1 = 'kg/m3 (all sigma layers)';
if has_agg
    out.stations.mp1_agg = station_mp1_agg;
    out.stations.units.mp1_agg = 'kg/m3 (all sigma layers)';
end
if has_dis
    out.stations.mp1_dis = station_mp1_dis;
    out.stations.units.mp1_dis = 'kg/m3 (all sigma layers)';
end
if has_sed
    out.stations.sed = station_sed;
    out.stations.units.sed = 'g/L (all sigma layers)';
end

out.notes = { ...
    'mp1_surf/bot/davg: combined plastic conc at surface/bottom/depth-avg.'; ...
    'mp1_agg/dis: aggregated and disaggregated components (PLAST_FLOC=T).'; ...
    'bot_mass_agg/dis: aggregated and disaggregated bed mass density.'; ...
    'sed_surf/bot/davg: coarse_sand_1 concentration (g/L).'; ...
    'sed_bedfrac: coarse_sand_1 surface bed fraction (dimensionless).'; ...
    'Station profiles span the full model time (stations.time).'; ...
    sprintf('Snapshot arrays span a ±%.1f-h window centred on day %.0f.', last_n_hours/2, target_day) ...
    };

out_dir_path = fileparts(out_mat);
if ~isempty(out_dir_path) && ~exist(out_dir_path, 'dir')
    mkdir(out_dir_path);
end
save(out_mat, 'out', '-v7.3');
fprintf('  Saved: %s\n', out_mat);
fprintf('  Snapshot grid: %d nodes x %d records\n', n_nodes, n_snaps);
fprintf('  Station profiles: %d stations x %d layers x %d records\n', ...
    n_stations, n_siglay, n_total);
end

%% =========================================================================
%  Shared helper functions
% ==========================================================================

function davg = depth_avg(conc_3d, vol_lay)
% Sigma-weighted depth average.  conc_3d: (node, siglay); vol_lay: (node, siglay).
wc_mass = nansum(conc_3d .* vol_lay, 2);
wc_vol  = nansum(vol_lay, 2);
davg = wc_mass ./ max(wc_vol, eps);
davg(wc_vol < eps) = 0.0;
end

function value = getenv_default(name, default_value)
value = getenv(name);
if isempty(value)
    value = default_value;
end
end

function value = getenv_number(name, default_value)
txt = getenv(name);
if isempty(txt)
    value = default_value;
    return;
end
value = str2double(txt);
if isnan(value)
    error('Environment variable %s must be numeric. Got: %s', name, txt);
end
end

function tf = is_absolute_path(p)
p = char(p);
tf = startsWith(p, '/') || startsWith(p, '\') || ...
    ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once'));
end

function vars_avail = inspect_variables(ncfile)
info = ncinfo(ncfile);
names = {info.Variables.Name};
vars_avail = struct();
for i = 1:numel(names)
    vars_avail.(matlab.lang.makeValidName(names{i})) = true;
end
end

function tf = has_var(vars_avail, varname)
tf = isfield(vars_avail, matlab.lang.makeValidName(varname)) && ...
    vars_avail.(matlab.lang.makeValidName(varname));
end

function d_km = haversine_km(lat_vec, lon_vec, lat0, lon0)
R = 6371.0;
dlat = deg2rad(lat_vec - lat0);
dlon = deg2rad(lon_vec - lon0);
a = sin(dlat/2).^2 + cos(deg2rad(lat0)) .* cos(deg2rad(lat_vec)) .* sin(dlon/2).^2;
d_km = 2 * R * asin(min(1, sqrt(a)));
end

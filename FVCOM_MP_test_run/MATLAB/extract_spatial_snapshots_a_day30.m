% EXTRACT_SPATIAL_SNAPSHOTS_A_DAY30
% Extract spatial snapshot and station time-series fields from the case-a
% FVCOM-MP NetCDF output centred on simulation day 30 (Memo 05).
%
% For each node and each record within ±12 h of day 30 this script saves:
%   mp1_surf   -- surface-layer (k=1) mp1 concentration  [kg/m3]
%   mp1_bot    -- bottom-layer (k=siglay) mp1 concentration  [kg/m3]
%   mp1_davg   -- depth-averaged (sigma-weighted) mp1 concentration  [kg/m3]
%   bot_mass   -- bed plastic mass density  [kg/m2]
%
% In addition, for each of 8 monitoring stations the full time series of
% every mp1 vertical layer is extracted by nearest-node lookup.
%
% Default NetCDF input:
%   <test_root>/OUTPUT_a/backup_data/waterPACT_a_0003.nc
%
% Default MAT output:
%   <test_root>/OUTPUT_a/spatial_snapshots_a_day30.mat  (v7.3 / HDF5)
%
% HPC usage (run from any directory on Kestrel data-analysis node):
%   cd /scratch/yhuang168/waterPACT_MP_floc
%   matlab -nodisplay -batch "run('MATLAB/extract_spatial_snapshots_a_day30.m')"
%
% Optional environment overrides:
%   WATERPACT_CASE_ID       Case label.  Default: a
%   WATERPACT_OUTPUT_DIR    Directory that holds backup_data/ and the output MAT.
%   WATERPACT_NC_FILE       Exact path to NetCDF file.
%   WATERPACT_OUT_MAT       MAT file to write.
%   WATERPACT_TARGET_DAY    Centre of the extraction window in simulation days.
%                           Default: 30
%   WATERPACT_LAST_NHOURS   Width of the extraction window in hours (centred on
%                           WATERPACT_TARGET_DAY).  Default: 24

%% -------------------------------------------------------------------------
%  User Settings
%  Absolute Kestrel paths are set as defaults.  Override via environment
%  variables for portability to other systems.
% --------------------------------------------------------------------------
case_id       = 'a';
output_dir    = '/kfs3/scratch/yhuang168/waterPACT_MP_floc/OUTPUT_a';
nc_file       = '/kfs3/scratch/yhuang168/waterPACT_MP_floc/OUTPUT_a/backup_data/waterPACT_a_0003.nc';
out_mat       = '/kfs3/scratch/yhuang168/waterPACT_MP_floc/OUTPUT_a/spatial_snapshots_a_day30.mat';
target_day    = 30;
last_n_hours  = 24;

%% Environment overrides
case_id      = getenv_default('WATERPACT_CASE_ID',    case_id);
output_dir   = getenv_default('WATERPACT_OUTPUT_DIR', output_dir);
nc_file      = getenv_default('WATERPACT_NC_FILE',    nc_file);
out_mat      = getenv_default('WATERPACT_OUT_MAT',    out_mat);
target_day   = getenv_number ('WATERPACT_TARGET_DAY',  target_day);
last_n_hours = getenv_number ('WATERPACT_LAST_NHOURS', last_n_hours);

%% Resolve paths
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

if isempty(output_dir)
    output_dir = fullfile(test_root, ['OUTPUT_' case_id]);
end
if isempty(nc_file)
    nc_file = fullfile(output_dir, 'backup_data', ...
        sprintf('waterPACT_%s_0003.nc', case_id));
end
if isempty(out_mat)
    out_mat = fullfile(output_dir, ['spatial_snapshots_' case_id '_day30.mat']);
end

if ~isfile(nc_file)
    error('NetCDF file not found: %s', nc_file);
end

fprintf('Case %s  (day-30 snapshot)\n', case_id);
fprintf('  NetCDF : %s\n', nc_file);
fprintf('  MAT out: %s\n', out_mat);

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
%  Read grid once
% --------------------------------------------------------------------------
fprintf('Reading grid ...\n');
lon    = double(ncread(nc_file, 'lon'));      lon    = lon(:);
lat    = double(ncread(nc_file, 'lat'));      lat    = lat(:);
lonc   = double(ncread(nc_file, 'lonc'));     lonc   = lonc(:);
latc   = double(ncread(nc_file, 'latc'));     latc   = latc(:);
h      = double(ncread(nc_file, 'h'));        h      = h(:);
art1   = double(ncread(nc_file, 'art1'));     art1   = art1(:);
nv     = double(ncread(nc_file, 'nv'));       % (nele, 3)
siglev = double(ncread(nc_file, 'siglev'));   % (node, siglev)

n_nodes  = numel(lon);
dzfrac   = abs(diff(siglev, 1, 2));          % (node, siglay)
n_siglay = size(dzfrac, 2);

%% -------------------------------------------------------------------------
%  Station nearest-node lookup  (haversine)
% --------------------------------------------------------------------------
fprintf('Locating nearest nodes for %d stations ...\n', n_stations);
station_node_idx = zeros(n_stations, 1, 'int32');
station_node_lat = zeros(n_stations, 1);
station_node_lon = zeros(n_stations, 1);
station_node_dist_km = zeros(n_stations, 1);

for is = 1:n_stations
    d_km = haversine_km(lat, lon, station_lat(is), station_lon(is));
    [min_d, idx] = min(d_km);
    station_node_idx(is)      = int32(idx);
    station_node_lat(is)      = lat(idx);
    station_node_lon(is)      = lon(idx);
    station_node_dist_km(is)  = min_d;
    fprintf('  %-22s  target (%.4f, %.4f) -> node %d  (%.4f, %.4f)  %.3f km\n', ...
        station_names{is}, station_lat(is), station_lon(is), ...
        idx, lat(idx), lon(idx), min_d);
end

%% -------------------------------------------------------------------------
%  Select records within ±(last_n_hours/2) h of target_day
% --------------------------------------------------------------------------
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

half_win    = last_n_hours / 2 / 24.0;        % half-window in days
snap_idx    = find(rel_time >= (target_day - half_win) & ...
                   rel_time <= (target_day + half_win));

if isempty(snap_idx)
    error('No records found within %.1f h of day %.1f (available range: %.4f .. %.4f sim-days).', ...
        last_n_hours/2, target_day, rel_time(1), rel_time(end));
end
n_snaps = numel(snap_idx);
fprintf('Snapshot window: day %.1f ± %.1f h -> %d records (%.4f .. %.4f sim-days)\n', ...
    target_day, last_n_hours/2, n_snaps, rel_time(snap_idx(1)), rel_time(snap_idx(end)));

%% -------------------------------------------------------------------------
%  Allocate snapshot arrays
% --------------------------------------------------------------------------
snap_time    = rel_time(snap_idx);            % (n_snaps,) simulation-relative days
mp1_surf     = zeros(n_nodes, n_snaps, 'single');
mp1_bot      = zeros(n_nodes, n_snaps, 'single');
mp1_davg     = zeros(n_nodes, n_snaps, 'single');
bot_mass_out = zeros(n_nodes, n_snaps, 'single');

%% -------------------------------------------------------------------------
%  Allocate station arrays  (full time series)
% --------------------------------------------------------------------------
station_time = rel_time;                      % (n_total,) simulation-relative days
% station_mp1_profile: n_stations x n_siglay x n_total
station_mp1  = zeros(n_stations, n_siglay, n_total, 'single');

%% -------------------------------------------------------------------------
%  Read loop -- station time series (all records, all sigma layers)
% --------------------------------------------------------------------------
fprintf('Extracting station profiles for all %d records ...\n', n_total);
for trec = 1:n_total
    mp1_3d = single(ncread(nc_file, 'mp1', [1 1 trec], [Inf Inf 1]));
    % mp1_3d is (node, siglay)
    for is = 1:n_stations
        nd = station_node_idx(is);
        station_mp1(is, :, trec) = mp1_3d(nd, :);
    end
    if mod(trec, 50) == 0 || trec == n_total
        fprintf('  station profiles: %d / %d\n', trec, n_total);
    end
end

%% -------------------------------------------------------------------------
%  Read loop -- spatial snapshots (window around day 30)
% --------------------------------------------------------------------------
fprintf('Extracting spatial snapshots for %d records ...\n', n_snaps);
for isnap = 1:n_snaps
    trec = snap_idx(isnap);

    % 3-D mp1  (node, siglay)
    mp1_3d = double(ncread(nc_file, 'mp1', [1 1 trec], [Inf Inf 1]));

    % Depth  (needed for sigma-weighted depth-average)
    zeta  = double(ncread(nc_file, 'zeta', [1 trec], [Inf 1]));
    zeta  = zeta(:);
    depth = max(h + zeta, 0.0);         % (node,)

    % Volume per layer  (node, siglay)
    base_vol = art1 .* depth;
    vol_lay  = bsxfun(@times, base_vol, dzfrac);   % (node, siglay)

    % Surface layer (k=1) and bottom layer (k=n_siglay)
    mp1_surf(:, isnap) = single(mp1_3d(:, 1));
    mp1_bot(:,  isnap) = single(mp1_3d(:, n_siglay));

    % Depth-average: sum(C * dz_frac * depth) / sum(dz_frac * depth)
    %   = sum(C * vol) / water_volume    (wet nodes only)
    wc_mass  = nansum(mp1_3d .* vol_lay, 2);    % (node,)
    wc_vol   = nansum(vol_lay, 2);               % (node,)
    davg_raw = wc_mass ./ max(wc_vol, eps);
    davg_raw(wc_vol < eps) = 0.0;
    mp1_davg(:, isnap) = single(davg_raw);

    % Bed mass
    bot = double(ncread(nc_file, 'bot_mass_mp1', [1 trec], [Inf 1]));
    bot_mass_out(:, isnap) = single(bot(:));

    if mod(isnap, 6) == 0 || isnap == n_snaps
        fprintf('  snapshots: %d / %d  (day %.4f)\n', isnap, n_snaps, snap_time(isnap));
    end
end

%% -------------------------------------------------------------------------
%  Pack output struct and save
% --------------------------------------------------------------------------
out = struct();

% Config / provenance
out.config.case_id       = case_id;
out.config.nc_file       = nc_file;
out.config.out_mat       = out_mat;
out.config.target_day    = target_day;
out.config.last_n_hours  = last_n_hours;
out.config.created_by    = mfilename;
out.config.created_on    = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

% Grid
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

% Snapshot info
out.snaps.snap_idx    = int32(snap_idx(:));
out.snaps.snap_time   = snap_time(:);
out.snaps.n_snaps     = n_snaps;

% Snapshot fields  (node x n_snaps)
out.snaps.mp1_surf    = mp1_surf;
out.snaps.mp1_bot     = mp1_bot;
out.snaps.mp1_davg    = mp1_davg;
out.snaps.bot_mass    = bot_mass_out;

% Snapshot units
out.snaps.units.mp1_surf  = 'kg/m3';
out.snaps.units.mp1_bot   = 'kg/m3';
out.snaps.units.mp1_davg  = 'kg/m3 (sigma-weighted depth-average)';
out.snaps.units.bot_mass  = 'kg/m2';

% Station lookup
out.stations.names         = station_names;
out.stations.target_lat    = station_lat;
out.stations.target_lon    = station_lon;
out.stations.node_idx      = station_node_idx;
out.stations.node_lat      = station_node_lat;
out.stations.node_lon      = station_node_lon;
out.stations.node_dist_km  = station_node_dist_km;

% Station time series  (n_stations x n_siglay x n_total)
out.stations.time     = station_time(:);
out.stations.mp1      = station_mp1;
out.stations.units.mp1 = 'kg/m3 (all sigma layers)';

% Notes
out.notes = { ...
    'mp1_surf: surface sigma-layer (k=1) concentration at each node.'; ...
    'mp1_bot: bottom sigma-layer (k=n_siglay) concentration at each node.'; ...
    'mp1_davg: sigma-volume-weighted depth-averaged concentration.'; ...
    'bot_mass: bed plastic mass density at each node.'; ...
    'Station profiles span the full model time (station_time).'; ...
    sprintf('Snapshot arrays span a ±%.1f-h window centred on day %.0f.', last_n_hours/2, target_day); ...
    'Nearest-node lookup used for station extraction (haversine distance).' ...
    };

out_dir = fileparts(out_mat);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

save(out_mat, 'out', '-v7.3');
fprintf('\nSaved: %s\n', out_mat);
fprintf('  Snapshot grid: %d nodes x %d records\n', n_nodes, n_snaps);
fprintf('  Station profiles: %d stations x %d sigma layers x %d records\n', ...
    n_stations, n_siglay, n_total);

%% =========================================================================
%  Local helper functions
% ==========================================================================

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

function d_km = haversine_km(lat_vec, lon_vec, lat0, lon0)
% Haversine distance in km from a single point (lat0,lon0) to all (lat_vec,lon_vec).
R = 6371.0;  % Earth radius, km
dlat = deg2rad(lat_vec - lat0);
dlon = deg2rad(lon_vec - lon0);
a = sin(dlat/2).^2 + cos(deg2rad(lat0)) .* cos(deg2rad(lat_vec)) .* sin(dlon/2).^2;
d_km = 2 * R * asin(min(1, sqrt(a)));
end

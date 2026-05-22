function extract_spatial_distribution_v2_core(target_day, nc_stack)
%EXTRACT_SPATIAL_DISTRIBUTION_V2_CORE
% Shared day-average spatial extractor for FVCOM-MP case c/d comparisons.
%
% This function is called by the day-specific scripts:
%   extract_spatial_distribution_v2_day10.m
%   extract_spatial_distribution_v2_day30.m
%   extract_spatial_distribution_v2_day90.m
%
% Run the day scripts from inside OUTPUT_c or OUTPUT_d. The case ID is
% detected from the current directory name unless WATERPACT_CASE_ID is set.

if nargin ~= 2
    error('Usage: extract_spatial_distribution_v2_core(target_day, nc_stack)');
end

rho_floc_eff = 1300.0;   % kg/m3, PLAST_FLOC_RHO_EFF default for case d
grav = 9.81;             % m/s2

output_dir = pwd;
case_id = detect_case_id(output_dir);

nc_override = getenv('WATERPACT_NC_FILE');
if isempty(nc_override)
    nc_file = fullfile(output_dir, sprintf('waterPACT_%s_%04d.nc', case_id, nc_stack));
    if ~exist(nc_file, 'file')
        nc_file = fullfile(output_dir, 'backup_data', ...
            sprintf('waterPACT_%s_%04d.nc', case_id, nc_stack));
    end
else
    nc_file = nc_override;
end
if ~exist(nc_file, 'file')
    error('NetCDF file not found: %s', nc_file);
end

out_override = getenv('WATERPACT_OUT_MAT');
if isempty(out_override)
    out_mat = fullfile(output_dir, ...
        sprintf('spatial_distribution_v2_%s_day%d.mat', case_id, target_day));
else
    out_mat = out_override;
end

fprintf('Spatial distribution v2 extraction\n');
fprintf('  case       : %s\n', case_id);
fprintf('  target day : %d\n', target_day);
fprintf('  nc stack   : %04d\n', nc_stack);
fprintf('  nc file    : %s\n', nc_file);
fprintf('  out mat    : %s\n', out_mat);

vars_avail = inspect_variables(nc_file);
required = {'mp1', 'coarse_sand_1', 'settle_vel_floc_mp1', 'ambient_rho', 'zeta'};
assert_required_vars(vars_avail, required, nc_file);

%% Grid
lon    = double(ncread(nc_file, 'lon'));      lon = lon(:);
lat    = double(ncread(nc_file, 'lat'));      lat = lat(:);
lonc   = double(ncread(nc_file, 'lonc'));     lonc = lonc(:);
latc   = double(ncread(nc_file, 'latc'));     latc = latc(:);
h      = double(ncread(nc_file, 'h'));        h = h(:);
art1   = double(ncread(nc_file, 'art1'));     art1 = art1(:);
nv     = double(ncread(nc_file, 'nv'));
siglev = double(ncread(nc_file, 'siglev'));

n_nodes = numel(lon);
dzfrac = abs(diff(siglev, 1, 2));
if size(dzfrac, 1) ~= n_nodes && size(dzfrac, 2) == n_nodes
    dzfrac = dzfrac';
end
n_siglay = size(dzfrac, 2);

%% Time window: final 24 h ending at target_day
all_time = double(ncread(nc_file, 'time'));
rel_time = all_time(:) - all_time(1) + (nc_stack - 1) * 10.0;
window_start = target_day - 1.0;
window_end = target_day;
rec_idx = find(rel_time >= window_start & rel_time <= window_end);
if isempty(rec_idx)
    error('No records found for day-%d diurnal window [%.3f, %.3f]. Available %.3f..%.3f.', ...
        target_day, window_start, window_end, rel_time(1), rel_time(end));
end
n_records = numel(rec_idx);
fprintf('  records    : %d records over simulation days %.4f..%.4f\n', ...
    n_records, rel_time(rec_idx(1)), rel_time(rec_idx(end)));

%% Fields
field3d = { ...
    'mp1',                  'mp1'; ...
    'mp1_agg',              'mp1_agg'; ...
    'mp1_dis',              'mp1_dis'; ...
    'coarse_sand_1',        'sed'; ...
    'settle_vel_mp1',       'wset_mp1'; ...
    'settle_vel_floc_mp1',  'wset_floc'; ...
    'lambda_a_mp1',         'lambda_a'; ...
    'lambda_d_mp1',         'lambda_d'; ...
    'ambient_rho',          'ambient_rho'; ...
    'temp',                 'temp' ...
    };

field2d = { ...
    'bot_mass_mp1',       'bot_mass'; ...
    'bot_mass_agg_mp1',   'bot_mass_agg'; ...
    'bot_mass_dis_mp1',   'bot_mass_dis' ...
    };

acc = struct();
present = struct();
for ii = 1:size(field3d, 1)
    src = field3d{ii, 1};
    tag = field3d{ii, 2};
    present.(tag) = has_var(vars_avail, src);
    if present.(tag)
        acc.(tag).surf = zeros(n_nodes, 1);
        acc.(tag).bot  = zeros(n_nodes, 1);
        acc.(tag).davg = zeros(n_nodes, 1);
    end
end

for ii = 1:size(field2d, 1)
    src = field2d{ii, 1};
    tag = field2d{ii, 2};
    present.(tag) = has_var(vars_avail, src);
    if present.(tag)
        acc.(tag).mean = zeros(n_nodes, 1);
    end
end

acc.d_floc_eff_m.surf = zeros(n_nodes, 1);
acc.d_floc_eff_m.bot  = zeros(n_nodes, 1);
acc.d_floc_eff_m.davg = zeros(n_nodes, 1);

temp_fallback_used = ~present.temp;
if temp_fallback_used
    fprintf('  temp       : missing; using 20 degC fallback for Stokes inversion.\n');
else
    fprintf('  temp       : present; using model temperature for Stokes inversion.\n');
end

%% Read and average
for ir = 1:n_records
    trec = rec_idx(ir);

    zeta = read_node_var(nc_file, 'zeta', trec, n_nodes);
    depth = max(h + zeta, 0.0);
    base_vol = art1 .* depth;
    vol_lay = bsxfun(@times, base_vol, dzfrac);

    for ii = 1:size(field3d, 1)
        src = field3d{ii, 1};
        tag = field3d{ii, 2};
        if ~present.(tag)
            continue;
        end
        val = read_3d_var(nc_file, src, trec, n_nodes, n_siglay);
        acc.(tag).surf = acc.(tag).surf + val(:, 1);
        acc.(tag).bot  = acc.(tag).bot  + val(:, n_siglay);
        acc.(tag).davg = acc.(tag).davg + depth_avg(val, vol_lay);
    end

    for ii = 1:size(field2d, 1)
        src = field2d{ii, 1};
        tag = field2d{ii, 2};
        if ~present.(tag)
            continue;
        end
        val = read_node_var(nc_file, src, trec, n_nodes);
        acc.(tag).mean = acc.(tag).mean + val;
    end

    wset_floc = read_3d_var(nc_file, 'settle_vel_floc_mp1', trec, n_nodes, n_siglay);
    ambient_rho = read_3d_var(nc_file, 'ambient_rho', trec, n_nodes, n_siglay);
    if present.temp
        temp_c = read_3d_var(nc_file, 'temp', trec, n_nodes, n_siglay);
    else
        temp_c = 20.0 * ones(n_nodes, n_siglay);
    end
    d_floc_eff = calc_stokes_equiv_size(wset_floc, ambient_rho, temp_c, ...
        rho_floc_eff, grav);
    acc.d_floc_eff_m.surf = acc.d_floc_eff_m.surf + d_floc_eff(:, 1);
    acc.d_floc_eff_m.bot  = acc.d_floc_eff_m.bot  + d_floc_eff(:, n_siglay);
    acc.d_floc_eff_m.davg = acc.d_floc_eff_m.davg + depth_avg(d_floc_eff, vol_lay);

    if mod(ir, 12) == 0 || ir == n_records
        fprintf('    record %d / %d (sim day %.4f)\n', ir, n_records, rel_time(trec));
    end
end

%% Pack output
out = struct();
out.config.case_id = case_id;
out.config.target_day = target_day;
out.config.nc_stack = nc_stack;
out.config.nc_file = nc_file;
out.config.out_mat = out_mat;
out.config.window_start_day = window_start;
out.config.window_end_day = window_end;
out.config.rho_floc_eff = rho_floc_eff;
out.config.grav = grav;
out.config.temp_fallback_used = temp_fallback_used;
out.config.created_by = mfilename;
out.config.created_on = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

out.grid.lon = single(lon);
out.grid.lat = single(lat);
out.grid.lonc = single(lonc);
out.grid.latc = single(latc);
out.grid.h = single(h);
out.grid.art1 = single(art1);
out.grid.nv = int32(nv);
out.grid.dzfrac = single(dzfrac);
out.grid.n_nodes = n_nodes;
out.grid.n_siglay = n_siglay;

out.snaps.target_day = target_day;
out.snaps.record_idx = int32(rec_idx(:));
out.snaps.record_time = single(rel_time(rec_idx));
out.snaps.n_records = n_records;
out.snaps.vars_available = vars_avail;

for ii = 1:size(field3d, 1)
    tag = field3d{ii, 2};
    if ~present.(tag)
        continue;
    end
    out.snaps.([tag '_surf']) = single(acc.(tag).surf / n_records);
    out.snaps.([tag '_bot'])  = single(acc.(tag).bot  / n_records);
    out.snaps.([tag '_davg']) = single(acc.(tag).davg / n_records);
end

for ii = 1:size(field2d, 1)
    tag = field2d{ii, 2};
    if ~present.(tag)
        continue;
    end
    out.snaps.(tag) = single(acc.(tag).mean / n_records);
end

d_m_surf = acc.d_floc_eff_m.surf / n_records;
d_m_bot  = acc.d_floc_eff_m.bot  / n_records;
d_m_davg = acc.d_floc_eff_m.davg / n_records;
out.snaps.d_floc_eff_m_surf = single(d_m_surf);
out.snaps.d_floc_eff_m_bot  = single(d_m_bot);
out.snaps.d_floc_eff_m_davg = single(d_m_davg);
out.snaps.d_floc_eff_um_surf = single(d_m_surf * 1.0e6);
out.snaps.d_floc_eff_um_bot  = single(d_m_bot  * 1.0e6);
out.snaps.d_floc_eff_um_davg = single(d_m_davg * 1.0e6);

out.snaps.units.mp1_surf = 'kg/m3';
out.snaps.units.mp1_bot = 'kg/m3';
out.snaps.units.mp1_davg = 'kg/m3, depth averaged';
out.snaps.units.mp1_agg_surf = 'kg/m3';
out.snaps.units.mp1_agg_bot = 'kg/m3';
out.snaps.units.mp1_agg_davg = 'kg/m3, depth averaged';
out.snaps.units.mp1_dis_surf = 'kg/m3';
out.snaps.units.mp1_dis_bot = 'kg/m3';
out.snaps.units.mp1_dis_davg = 'kg/m3, depth averaged';
out.snaps.units.sed_surf = 'g/L or kg/m3, as written by FVCOM sediment output';
out.snaps.units.sed_bot = 'g/L or kg/m3, as written by FVCOM sediment output';
out.snaps.units.sed_davg = 'g/L or kg/m3, depth averaged';
out.snaps.units.wset_mp1_surf = 'm/s';
out.snaps.units.wset_mp1_bot = 'm/s';
out.snaps.units.wset_mp1_davg = 'm/s, depth averaged';
out.snaps.units.wset_floc_surf = 'm/s';
out.snaps.units.wset_floc_bot = 'm/s';
out.snaps.units.wset_floc_davg = 'm/s, depth averaged';
out.snaps.units.lambda_a_surf = '1/s';
out.snaps.units.lambda_a_bot = '1/s';
out.snaps.units.lambda_a_davg = '1/s, depth averaged';
out.snaps.units.lambda_d_surf = '1/s';
out.snaps.units.lambda_d_bot = '1/s';
out.snaps.units.lambda_d_davg = '1/s, depth averaged';
out.snaps.units.ambient_rho_surf = 'kg/m3';
out.snaps.units.ambient_rho_bot = 'kg/m3';
out.snaps.units.ambient_rho_davg = 'kg/m3, depth averaged';
out.snaps.units.temp_surf = 'degC';
out.snaps.units.temp_bot = 'degC';
out.snaps.units.temp_davg = 'degC, depth averaged';
out.snaps.units.bot_mass = 'kg/m2';
out.snaps.units.bot_mass_agg = 'kg/m2';
out.snaps.units.bot_mass_dis = 'kg/m2';
out.snaps.units.d_floc_eff_m_surf = 'm';
out.snaps.units.d_floc_eff_m_bot = 'm';
out.snaps.units.d_floc_eff_m_davg = 'm, depth averaged';
out.snaps.units.d_floc_eff_um_surf = 'micron';
out.snaps.units.d_floc_eff_um_bot = 'micron';
out.snaps.units.d_floc_eff_um_davg = 'micron, depth averaged';

out.notes = { ...
    'Each field is a simple arithmetic mean over records in the final 24 h ending at target_day.'; ...
    'Depth averages are sigma-layer volume weighted using art1*(h+zeta)*dzfrac.'; ...
    'd_floc_eff is computed layer-by-layer from settle_vel_floc_mp1, ambient_rho, and temp.'; ...
    'No station time series are extracted in v2.' ...
    };

save(out_mat, 'out', '-v7.3');
fprintf('  Saved %s\n', out_mat);

end

function case_id = detect_case_id(output_dir)
env_case = getenv('WATERPACT_CASE_ID');
if ~isempty(env_case)
    case_id = lower(strtrim(env_case));
    return;
end
[~, folder] = fileparts(output_dir);
tok = regexp(folder, '^OUTPUT_([A-Za-z0-9]+)$', 'tokens', 'once');
if isempty(tok)
    error('Cannot detect case ID from current directory: %s. Set WATERPACT_CASE_ID.', output_dir);
end
case_id = lower(tok{1});
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
fname = matlab.lang.makeValidName(varname);
tf = isfield(vars_avail, fname) && vars_avail.(fname);
end

function assert_required_vars(vars_avail, required, ncfile)
missing = {};
for i = 1:numel(required)
    if ~has_var(vars_avail, required{i})
        missing{end+1} = required{i}; %#ok<AGROW>
    end
end
if ~isempty(missing)
    error('Missing required variable(s) in %s: %s', ncfile, strjoin(missing, ', '));
end
end

function val = read_node_var(ncfile, varname, trec, n_nodes)
raw = squeeze(double(ncread(ncfile, varname, [1 trec], [Inf 1])));
if numel(raw) ~= n_nodes
    error('Variable %s record %d has %d values, expected %d nodes.', ...
        varname, trec, numel(raw), n_nodes);
end
val = raw(:);
end

function val = read_3d_var(ncfile, varname, trec, n_nodes, n_siglay)
raw = squeeze(double(ncread(ncfile, varname, [1 1 trec], [Inf Inf 1])));
if isequal(size(raw), [n_nodes n_siglay])
    val = raw;
elseif isequal(size(raw), [n_siglay n_nodes])
    val = raw';
else
    error('Variable %s record %d has size [%s], expected [%d %d].', ...
        varname, trec, num2str(size(raw)), n_nodes, n_siglay);
end
end

function davg = depth_avg(field3d, vol_lay)
valid = isfinite(field3d) & isfinite(vol_lay) & vol_lay > 0.0;
weighted = zeros(size(field3d));
vol_valid = zeros(size(field3d));
weighted(valid) = field3d(valid) .* vol_lay(valid);
vol_valid(valid) = vol_lay(valid);
den = sum(vol_valid, 2);
davg = sum(weighted, 2) ./ max(den, eps);
davg(den <= eps) = NaN;
end

function d_floc = calc_stokes_equiv_size(wset_floc, ambient_rho, temp_c, rho_floc_eff, grav)
nu = 1.0e-6 .* exp(-0.025 .* (temp_c - 20.0));
nu = max(0.5e-6, min(2.0e-6, nu));
mu = ambient_rho .* nu;
rho_diff = rho_floc_eff - ambient_rho;

d_floc = NaN(size(wset_floc));
valid = isfinite(wset_floc) & isfinite(mu) & isfinite(rho_diff) & ...
    wset_floc >= 0.0 & rho_diff > 0.0;
d_floc(valid) = sqrt(18.0 .* mu(valid) .* wset_floc(valid) ./ ...
    (grav .* rho_diff(valid)));
end

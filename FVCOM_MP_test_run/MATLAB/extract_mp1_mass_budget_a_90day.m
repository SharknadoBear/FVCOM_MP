% extract_mp1_mass_budget_a_90day
% Integrate FVCOM-MP mp1 water and bed mass across three consecutive output
% files (_0001 through _0009) to cover a ~90-day period.
%
% This diagnostic concatenates waterPACT_a_0001.nc, waterPACT_a_0002.nc, and
% waterPACT_a_0003.nc from the OUTPUT_a directory, computes domain-integrated
% water-column and bed mass time series, and saves a compact MAT file for
% local analysis. Raw 3-D concentration fields are not saved.
%
% This is a regular script. Edit the User Settings block below, or override
% paths from the shell with environment variables before launching MATLAB.
%
% Example HPC usage after editing paths below:
%   cd /kfs3/scratch/yhuang168/waterPACT_MP_floc/FVCOM_MP_test_run/MATLAB
%   matlab -batch "run('extract_mp1_mass_budget_a_90day.m')"
%
% Optional environment overrides:
%   WATERPACT_CASE_ID       Case label. Default: a
%   WATERPACT_OUTPUT_DIR    Directory containing the target NetCDF file
%   WATERPACT_NC_FILE       Exact NetCDF file name or absolute path
%   WATERPACT_OUT_MAT       MAT file to write
%   WATERPACT_TIME_STRIDE   Read every Nth output record
%   WATERPACT_READ_SPLIT    true/false for mp1_agg/mp1_dis if present
%   WATERPACT_READ_FLUX     true/false for dep_flux_mp1/ero_flux_mp1
%   WATERPACT_SAVE_GRID     true/false for saving art1/h/sigma fractions

%% User Settings
% Edit these values for a normal script run.
%  Absolute Kestrel paths are set as defaults.  Override via environment
%  variables for portability to other systems.
case_id = 'a';
output_dir = '/kfs3/scratch/yhuang168/waterPACT_MP_floc/OUTPUT_a';
out_mat = '';
nc_file = '';
time_stride = 1;
read_split = true;
read_flux_fields = true;
save_grid_metrics = true;
cap_threshold_kgm3 = 100.0;
cap_tolerance_kgm3 = 1.0e-6;
% File suffix labels to search.  Files that do not exist are skipped with a
% warning, so a partial local set (e.g. only _0003.nc) works fine.
%   e.g. set WATERPACT_FILE_LABELS=0003 to use only the third file.
file_labels = {'0001', '0002', '0003', '0004', '0005', '0006', '0007', '0008', '0009'};

% Optional batch-job overrides. If an environment variable is empty, the
% editable value above is kept.
case_id = getenv_default('WATERPACT_CASE_ID', case_id);
output_dir = getenv_default('WATERPACT_OUTPUT_DIR', output_dir);
out_mat = getenv_default('WATERPACT_OUT_MAT', out_mat);
nc_file = getenv_default('WATERPACT_NC_FILE', nc_file);
time_stride = getenv_number('WATERPACT_TIME_STRIDE', time_stride);
read_split = getenv_logical('WATERPACT_READ_SPLIT', read_split);
read_flux_fields = getenv_logical('WATERPACT_READ_FLUX', read_flux_fields);
save_grid_metrics = getenv_logical('WATERPACT_SAVE_GRID', save_grid_metrics);
cap_threshold_kgm3 = getenv_number('WATERPACT_CAP_THRESHOLD_KGM3', cap_threshold_kgm3);
cap_tolerance_kgm3 = getenv_number('WATERPACT_CAP_TOLERANCE_KGM3', cap_tolerance_kgm3);
file_labels = getenv_list('WATERPACT_FILE_LABELS', file_labels);

%% Resolve root path
% WATERPACT_TEST_ROOT overrides auto-detection so the script works when
% MATLAB is launched via run() in batch mode and mfilename returns empty.
%   export WATERPACT_TEST_ROOT=/kfs3/scratch/yhuang168/waterPACT_MP_floc
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
test_root_env = getenv_default('WATERPACT_TEST_ROOT', '');
if ~isempty(test_root_env)
    test_root = test_root_env;
else
    test_root = fileparts(script_dir);
end

if isempty(output_dir)
    output_dir = fullfile(test_root, ['OUTPUT_' case_id]);
else
    case_id = infer_case_id(output_dir);
end

% Build list of output files from file_labels; skip any that don't exist.
files = {};
for ilab = 1:numel(file_labels)
    fname = sprintf('waterPACT_%s_%s.nc', case_id, char(file_labels{ilab}));
    fpath = fullfile(output_dir, 'backup_data', fname);
    if isfile(fpath)
        files{end+1} = fpath; %#ok<AGROW>
    else
        fprintf('  [skip] Not found: %s\n', fpath);
    end
end
if isempty(files)
    error('No readable NetCDF files found in %s/backup_data/ for case %s.', output_dir, case_id);
end

if isempty(out_mat)
    out_mat = fullfile(output_dir, ['mp1_mass_budget_' case_id '_90day.mat']);
end

config = struct();
config.output_dir = char(output_dir);
config.nc_file = char(files{1});   % first file stored for reference
config.out_mat = char(out_mat);
config.time_stride = time_stride;
config.read_split = logical(read_split);
config.read_flux_fields = logical(read_flux_fields);
config.save_grid_metrics = logical(save_grid_metrics);
config.cap_threshold_kgm3 = cap_threshold_kgm3;
config.cap_tolerance_kgm3 = cap_tolerance_kgm3;
config.case_id = case_id;
config.created_by = mfilename;
config.created_on = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

config.time_stride = max(1, round(config.time_stride));

fprintf('Reading %d NetCDF files (~90-day concatenation):\n', numel(files));
for ilab = 1:numel(files)
    fprintf('  %s\n', files{ilab});
end

first_file = files{1};
vars_available = inspect_variables(first_file);
required_vars = {'art1', 'h', 'siglev', 'zeta', 'mp1', 'bot_mass_mp1', 'time'};
assert_required_vars(vars_available, required_vars, first_file);

art1 = double(ncread(first_file, 'art1'));
h = double(ncread(first_file, 'h'));
siglev = double(ncread(first_file, 'siglev'));
art1 = art1(:);
h = h(:);

if size(siglev, 1) ~= numel(h)
    error('siglev first dimension (%d) does not match node count (%d).', size(siglev, 1), numel(h));
end

dzfrac = abs(diff(siglev, 1, 2));
sigma_sum = sum(dzfrac, 2);
grid_summary = struct();
grid_summary.node_count = numel(h);
grid_summary.siglay_count = size(dzfrac, 2);
grid_summary.domain_area_m2 = nansum_local(art1);
grid_summary.h_min_m = finite_min(h);
grid_summary.h_max_m = finite_max(h);
grid_summary.h_mean_m = nansum_local(h .* art1) / max(nansum_local(art1), eps);
grid_summary.sigma_fraction_sum_min = finite_min(sigma_sum);
grid_summary.sigma_fraction_sum_max = finite_max(sigma_sum);
grid_summary.sigma_fraction_sum_max_abs_error = finite_max(abs(sigma_sum - 1.0));

has_wet_nodes = has_var(vars_available, 'wet_nodes');
has_split_water = config.read_split && has_var(vars_available, 'mp1_agg') && has_var(vars_available, 'mp1_dis');
has_split_bed = config.read_split && has_var(vars_available, 'bot_mass_agg_mp1') ...
    && has_var(vars_available, 'bot_mass_dis_mp1');
has_dep_flux = config.read_flux_fields && has_var(vars_available, 'dep_flux_mp1');
has_ero_flux = config.read_flux_fields && has_var(vars_available, 'ero_flux_mp1');

record_plan = build_record_plan(files, config.time_stride);
nrec = numel(record_plan.file_index);
if nrec == 0
    error('No time records selected. Check TimeStride and input files.');
end

budget = initialize_budget(nrec, has_split_water, has_split_bed, has_dep_flux, has_ero_flux);

fprintf('Selected %d output record(s), reading every %d record(s).\n', nrec, config.time_stride);

out_idx = 0;
for ifile = 1:numel(files)
    ncfile = files{ifile};
    local_records = record_plan.records_by_file{ifile};
    if isempty(local_records)
        continue;
    end

    fprintf('  %s: %d selected record(s)\n', get_filename(ncfile), numel(local_records));

    file_time = double(ncread(ncfile, 'time'));   % absolute days (MJD)
    % Convert to simulation-relative days: file _000N starts at (N-1)*10 sim-days.
    tok_t = regexp(ncfile, '_(\d{4})\.nc', 'tokens', 'once');
    if ~isempty(tok_t)
        file_num_t = str2double(tok_t{1});
    else
        file_num_t = 1;
        warning('Could not parse file number from %s; assuming file 1.', ncfile);
    end
    rel_time = file_time - file_time(1) + (file_num_t - 1) * 10;
    file_iint = safe_read_vector(ncfile, 'iint');
    file_itime = safe_read_vector(ncfile, 'Itime');
    file_itime2 = safe_read_vector(ncfile, 'Itime2');

    for ilocal = 1:numel(local_records)
        tindex = local_records(ilocal);
        out_idx = out_idx + 1;

        zeta = double(ncread(ncfile, 'zeta', [1 tindex], [Inf 1]));
        zeta = zeta(:);
        depth = max(h + zeta, 0.0);

        if has_wet_nodes
            wet_nodes = double(ncread(ncfile, 'wet_nodes', [1 tindex], [Inf 1])) > 0;
            wet_nodes = wet_nodes(:);
        else
            wet_nodes = depth > 0.0;
        end

        base_volume = art1 .* depth .* double(wet_nodes);
        volume = bsxfun(@times, base_volume, dzfrac);

        mp1 = double(ncread(ncfile, 'mp1', [1 1 tindex], [Inf Inf 1]));
        mp1 = squeeze(mp1);
        bot_mass = double(ncread(ncfile, 'bot_mass_mp1', [1 tindex], [Inf 1]));
        bot_mass = bot_mass(:);

        budget.file_index(out_idx) = ifile;
        budget.record_index(out_idx) = tindex;
        budget.time_days(out_idx) = rel_time(tindex);  % simulation-relative days
        budget.iint(out_idx) = vector_value_or_nan(file_iint, tindex);
        budget.Itime(out_idx) = vector_value_or_nan(file_itime, tindex);
        budget.Itime2(out_idx) = vector_value_or_nan(file_itime2, tindex);

        budget.water_volume_m3(out_idx) = nansum_local(volume(:));
        budget.wet_area_m2(out_idx) = nansum_local(art1(wet_nodes));
        budget.zeta_min_m(out_idx) = finite_min(zeta);
        budget.zeta_max_m(out_idx) = finite_max(zeta);
        budget.zeta_area_mean_m(out_idx) = nansum_local(zeta .* art1) / max(nansum_local(art1), eps);

        mp1_mass_field = mp1 .* volume;
        budget.water_mass_kg(out_idx) = nansum_local(mp1_mass_field(:));
        budget.water_negative_mass_kg(out_idx) = -nansum_local(min(mp1, 0.0) .* volume);
        budget.water_valid_volume_m3(out_idx) = nansum_local(volume(isfinite(mp1)));
        budget.mp1_min_kgm3(out_idx) = finite_min(mp1(:));
        budget.mp1_max_kgm3(out_idx) = finite_max(mp1(:));
        budget.max_mp1_conc_kgm3(out_idx) = budget.mp1_max_kgm3(out_idx);
        cap_mask = isfinite(mp1) & mp1 >= config.cap_threshold_kgm3;
        near_cap_mask = isfinite(mp1) & ...
            mp1 >= config.cap_threshold_kgm3 - config.cap_tolerance_kgm3;
        budget.mp1_count_ge_cap(out_idx) = nnz(cap_mask);
        budget.mp1_count_near_or_ge_cap(out_idx) = nnz(near_cap_mask);
        budget.mp1_volume_near_or_ge_cap_m3(out_idx) = nansum_local(volume(near_cap_mask));
        budget.mp1_nan_count(out_idx) = nnz(~isfinite(mp1(:)));
        budget.mp1_negative_count(out_idx) = nnz(mp1(:) < 0.0);

        bot_area = art1(isfinite(bot_mass));
        budget.bed_mass_all_nodes_kg(out_idx) = nansum_local(bot_mass .* art1);
        budget.bed_mass_wet_nodes_kg(out_idx) = nansum_local(bot_mass(wet_nodes) .* art1(wet_nodes));
        budget.bed_valid_area_m2(out_idx) = nansum_local(bot_area);
        budget.bot_mass_min_kgm2(out_idx) = finite_min(bot_mass);
        budget.bot_mass_max_kgm2(out_idx) = finite_max(bot_mass);
        budget.bot_mass_nan_count(out_idx) = nnz(~isfinite(bot_mass));
        budget.bot_mass_negative_count(out_idx) = nnz(bot_mass < 0.0);
        budget.bot_mass_negative_kg(out_idx) = -nansum_local(min(bot_mass, 0.0) .* art1);

        if has_split_water
            mp1_agg = double(ncread(ncfile, 'mp1_agg', [1 1 tindex], [Inf Inf 1]));
            mp1_dis = double(ncread(ncfile, 'mp1_dis', [1 1 tindex], [Inf Inf 1]));
            mp1_agg = squeeze(mp1_agg);
            mp1_dis = squeeze(mp1_dis);
            budget.water_agg_mass_kg(out_idx) = nansum_local(mp1_agg(:) .* volume(:));
            budget.water_dis_mass_kg(out_idx) = nansum_local(mp1_dis(:) .* volume(:));
            budget.water_split_minus_total_kg(out_idx) = ...
                budget.water_agg_mass_kg(out_idx) + budget.water_dis_mass_kg(out_idx) - budget.water_mass_kg(out_idx);
            budget.mp1_split_minus_total_min_kgm3(out_idx) = finite_min(mp1_agg(:) + mp1_dis(:) - mp1(:));
            budget.mp1_split_minus_total_max_kgm3(out_idx) = finite_max(mp1_agg(:) + mp1_dis(:) - mp1(:));
        end

        if has_split_bed
            bot_agg = double(ncread(ncfile, 'bot_mass_agg_mp1', [1 tindex], [Inf 1]));
            bot_dis = double(ncread(ncfile, 'bot_mass_dis_mp1', [1 tindex], [Inf 1]));
            bot_agg = bot_agg(:);
            bot_dis = bot_dis(:);
            budget.bed_agg_mass_all_nodes_kg(out_idx) = nansum_local(bot_agg .* art1);
            budget.bed_dis_mass_all_nodes_kg(out_idx) = nansum_local(bot_dis .* art1);
            budget.bed_split_minus_total_kg(out_idx) = ...
                budget.bed_agg_mass_all_nodes_kg(out_idx) + budget.bed_dis_mass_all_nodes_kg(out_idx) ...
                - budget.bed_mass_all_nodes_kg(out_idx);
            budget.bot_split_minus_total_min_kgm2(out_idx) = finite_min(bot_agg + bot_dis - bot_mass);
            budget.bot_split_minus_total_max_kgm2(out_idx) = finite_max(bot_agg + bot_dis - bot_mass);
        end

        if has_dep_flux
            dep_flux = double(ncread(ncfile, 'dep_flux_mp1', [1 tindex], [Inf 1]));
            dep_flux = dep_flux(:);
            budget.dep_flux_area_integral_kg(out_idx) = nansum_local(dep_flux .* art1);
        end

        if has_ero_flux
            ero_flux = double(ncread(ncfile, 'ero_flux_mp1', [1 tindex], [Inf 1]));
            ero_flux = ero_flux(:);
            budget.ero_flux_area_integral_kg(out_idx) = nansum_local(ero_flux .* art1);
        end
    end
end

budget.water_mean_conc_kgm3 = safe_divide(budget.water_mass_kg, budget.water_valid_volume_m3);
budget.bed_mean_mass_density_kgm2 = safe_divide(budget.bed_mass_all_nodes_kg, budget.bed_valid_area_m2);
budget.total_mass_all_bed_kg = budget.water_mass_kg + budget.bed_mass_all_nodes_kg;
budget.total_mass_wet_bed_kg = budget.water_mass_kg + budget.bed_mass_wet_nodes_kg;
budget.total_equiv_conc_all_bed_kgm3 = safe_divide(budget.total_mass_all_bed_kg, budget.water_volume_m3);

budget.water_mass_change_kg = change_from_first(budget.water_mass_kg);
budget.bed_mass_change_kg = change_from_first(budget.bed_mass_all_nodes_kg);
budget.total_mass_change_kg = change_from_first(budget.total_mass_all_bed_kg);
budget.total_mass_change_percent = percent_change_from_first(budget.total_mass_all_bed_kg);
budget.total_mass_rate_kg_per_day = rate_from_time(budget.total_mass_all_bed_kg, budget.time_days);

summary = make_summary(budget, config, files, grid_summary, has_split_water, has_split_bed);

diagnostic = struct();
diagnostic.config = config;
diagnostic.files = files(:);
diagnostic.vars_available = vars_available;
diagnostic.grid_summary = grid_summary;
diagnostic.budget = budget;
diagnostic.summary = summary;
diagnostic.notes = { ...
    'water_mass_kg = integral mp1 * art1 * max(h+zeta,0) * sigma_layer_fraction over wet nodes'; ...
    'bed_mass_all_nodes_kg = integral bot_mass_mp1 * art1 over all nodes'; ...
    'bed_mass_wet_nodes_kg is also saved, but all-node bed mass is usually the cleaner conservation companion'; ...
    'dep_flux_mp1 and ero_flux_mp1 are area-integrated if present; FVCOM output units are kg/m^2'; ...
    'max_mp1_conc_kgm3 and mp1_count_near_or_ge_cap diagnose whether the 100 kg/m^3 concentration cap appears in output' ...
    };

if config.save_grid_metrics
    grid = struct();
    grid.art1_m2 = art1;
    grid.h_m = h;
    grid.sigma_layer_fraction = dzfrac;
    diagnostic.grid = grid;
end

out_dir = fileparts(config.out_mat);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

save(config.out_mat, 'diagnostic', '-v7.3');

fprintf('\nSaved MAT diagnostic:\n  %s\n', config.out_mat);
fprintf('Initial total mass: %.12g kg\n', summary.initial_total_mass_kg);
fprintf('Final total mass:   %.12g kg\n', summary.final_total_mass_kg);
fprintf('Final change:       %.12g kg (%.6g %%)\n', ...
    summary.final_total_mass_change_kg, summary.final_total_mass_change_percent);

function value = getenv_default(name, default_value)
value = getenv(name);
if isempty(value)
    value = default_value;
end
end

function value = getenv_number(name, default_value)
text_value = getenv(name);
if isempty(text_value)
    value = default_value;
    return;
end

value = str2double(text_value);
if isnan(value)
    error('Environment variable %s must be numeric. Got: %s', name, text_value);
end
end

function value = getenv_logical(name, default_value)
text_value = getenv(name);
if isempty(text_value)
    value = default_value;
    return;
end

switch lower(strtrim(text_value))
    case {'true', 't', '1', 'yes', 'y', 'on'}
        value = true;
    case {'false', 'f', '0', 'no', 'n', 'off'}
        value = false;
    otherwise
        error('Environment variable %s must be true/false. Got: %s', name, text_value);
end
end

function values = getenv_list(name, default_values)
text_value = getenv(name);
if isempty(text_value)
    values = default_values;
    return;
end
parts = regexp(text_value, ',', 'split');
values = {};
for i = 1:numel(parts)
    item = strtrim(parts{i});
    if ~isempty(item)
        values{end + 1} = item; %#ok<AGROW>
    end
end
end

function case_id = infer_case_id(output_dir)
[~, name] = fileparts(char(output_dir));
case_id = regexprep(name, '^OUTPUT_', '');
if isempty(case_id) || strcmp(case_id, name)
    case_id = getenv_default('WATERPACT_CASE_ID', 'custom');
end
end

function nc_file = resolve_nc_file(output_dir, nc_file, case_id)
if isempty(nc_file)
    nc_file = sprintf('waterPACT_%s_0001.nc', case_id);
end

if ~is_absolute_path(nc_file)
    nc_file = fullfile(output_dir, nc_file);
end
end

function tf = is_absolute_path(path_value)
path_value = char(path_value);
tf = startsWith(path_value, '/') || startsWith(path_value, '\') ...
    || ~isempty(regexp(path_value, '^[A-Za-z]:[\\/]', 'once'));
end

function vars_available = inspect_variables(ncfile)
info = ncinfo(ncfile);
names = {info.Variables.Name};
vars_available = struct();
for i = 1:numel(names)
    safe_name = matlab.lang.makeValidName(names{i});
    vars_available.(safe_name) = true;
end
end

function assert_required_vars(vars_available, required_vars, ncfile)
missing = {};
for i = 1:numel(required_vars)
    safe_name = matlab.lang.makeValidName(required_vars{i});
    if ~isfield(vars_available, safe_name) || ~vars_available.(safe_name)
        missing{end + 1} = required_vars{i}; %#ok<AGROW>
    end
end

if ~isempty(missing)
    error('Missing required variable(s) in %s: %s', ncfile, strjoin(missing, ', '));
end
end

function tf = has_var(vars_available, varname)
safe_name = matlab.lang.makeValidName(varname);
tf = isfield(vars_available, safe_name) && vars_available.(safe_name);
end

function record_plan = build_record_plan(files, time_stride)
nfiles = numel(files);
records_by_file = cell(nfiles, 1);
file_index = [];
record_index = [];

for ifile = 1:nfiles
    t = ncread(files{ifile}, 'time');
    records = 1:time_stride:numel(t);
    records_by_file{ifile} = records;
    file_index = [file_index; repmat(ifile, numel(records), 1)]; %#ok<AGROW>
    record_index = [record_index; records(:)]; %#ok<AGROW>
end

record_plan = struct();
record_plan.records_by_file = records_by_file;
record_plan.file_index = file_index;
record_plan.record_index = record_index;
end

function budget = initialize_budget(nrec, has_split_water, has_split_bed, has_dep_flux, has_ero_flux)
nan_col = nan(nrec, 1);
zero_col = zeros(nrec, 1);

budget = struct();
budget.file_index = zero_col;
budget.record_index = zero_col;
budget.time_days = nan_col;
budget.iint = nan_col;
budget.Itime = nan_col;
budget.Itime2 = nan_col;

budget.water_volume_m3 = nan_col;
budget.water_valid_volume_m3 = nan_col;
budget.wet_area_m2 = nan_col;
budget.zeta_min_m = nan_col;
budget.zeta_max_m = nan_col;
budget.zeta_area_mean_m = nan_col;

budget.water_mass_kg = nan_col;
budget.water_negative_mass_kg = nan_col;
budget.water_mean_conc_kgm3 = nan_col;
budget.mp1_min_kgm3 = nan_col;
budget.mp1_max_kgm3 = nan_col;
budget.max_mp1_conc_kgm3 = nan_col;
budget.mp1_count_ge_cap = zero_col;
budget.mp1_count_near_or_ge_cap = zero_col;
budget.mp1_volume_near_or_ge_cap_m3 = nan_col;
budget.mp1_nan_count = zero_col;
budget.mp1_negative_count = zero_col;

budget.bed_mass_all_nodes_kg = nan_col;
budget.bed_mass_wet_nodes_kg = nan_col;
budget.bed_valid_area_m2 = nan_col;
budget.bed_mean_mass_density_kgm2 = nan_col;
budget.bot_mass_min_kgm2 = nan_col;
budget.bot_mass_max_kgm2 = nan_col;
budget.bot_mass_nan_count = zero_col;
budget.bot_mass_negative_count = zero_col;
budget.bot_mass_negative_kg = nan_col;

budget.total_mass_all_bed_kg = nan_col;
budget.total_mass_wet_bed_kg = nan_col;
budget.total_equiv_conc_all_bed_kgm3 = nan_col;
budget.water_mass_change_kg = nan_col;
budget.bed_mass_change_kg = nan_col;
budget.total_mass_change_kg = nan_col;
budget.total_mass_change_percent = nan_col;
budget.total_mass_rate_kg_per_day = nan_col;

if has_split_water
    budget.water_agg_mass_kg = nan_col;
    budget.water_dis_mass_kg = nan_col;
    budget.water_split_minus_total_kg = nan_col;
    budget.mp1_split_minus_total_min_kgm3 = nan_col;
    budget.mp1_split_minus_total_max_kgm3 = nan_col;
end

if has_split_bed
    budget.bed_agg_mass_all_nodes_kg = nan_col;
    budget.bed_dis_mass_all_nodes_kg = nan_col;
    budget.bed_split_minus_total_kg = nan_col;
    budget.bot_split_minus_total_min_kgm2 = nan_col;
    budget.bot_split_minus_total_max_kgm2 = nan_col;
end

if has_dep_flux
    budget.dep_flux_area_integral_kg = nan_col;
end

if has_ero_flux
    budget.ero_flux_area_integral_kg = nan_col;
end
end

function values = safe_read_vector(ncfile, varname)
info = ncinfo(ncfile);
names = {info.Variables.Name};
if any(strcmp(names, varname))
    values = double(ncread(ncfile, varname));
    values = values(:);
else
    values = [];
end
end

function value = vector_value_or_nan(values, index)
if isempty(values) || index > numel(values)
    value = nan;
else
    value = values(index);
end
end

function value = nansum_local(values)
values = values(:);
values = values(isfinite(values));
if isempty(values)
    value = 0.0;
else
    value = sum(values);
end
end

function value = finite_min(values)
values = values(:);
values = values(isfinite(values));
if isempty(values)
    value = nan;
else
    value = min(values);
end
end

function value = finite_max(values)
values = values(:);
values = values(isfinite(values));
if isempty(values)
    value = nan;
else
    value = max(values);
end
end

function out = safe_divide(num, den)
out = nan(size(num));
mask = isfinite(num) & isfinite(den) & abs(den) > 0.0;
out(mask) = num(mask) ./ den(mask);
end

function out = change_from_first(values)
out = nan(size(values));
idx = find(isfinite(values), 1, 'first');
if ~isempty(idx)
    out = values - values(idx);
end
end

function out = percent_change_from_first(values)
out = nan(size(values));
idx = find(isfinite(values), 1, 'first');
if ~isempty(idx) && abs(values(idx)) > 0.0
    out = 100.0 * (values - values(idx)) ./ abs(values(idx));
end
end

function out = rate_from_time(values, time_days)
out = nan(size(values));
if numel(values) < 2
    return;
end
dt = diff(time_days);
dv = diff(values);
mask = isfinite(dt) & isfinite(dv) & abs(dt) > 0.0;
rate = nan(numel(dt), 1);
rate(mask) = dv(mask) ./ dt(mask);
out(2:end) = rate;
end

function summary = make_summary(budget, config, files, grid_summary, has_split_water, has_split_bed)
summary = struct();
summary.case_id = config.case_id;
summary.output_dir = config.output_dir;
summary.file_count = numel(files);
summary.record_count = numel(budget.time_days);
summary.time_start_days = first_finite(budget.time_days);
summary.time_end_days = last_finite(budget.time_days);
summary.initial_water_mass_kg = first_finite(budget.water_mass_kg);
summary.final_water_mass_kg = last_finite(budget.water_mass_kg);
summary.initial_bed_mass_kg = first_finite(budget.bed_mass_all_nodes_kg);
summary.final_bed_mass_kg = last_finite(budget.bed_mass_all_nodes_kg);
summary.initial_total_mass_kg = first_finite(budget.total_mass_all_bed_kg);
summary.final_total_mass_kg = last_finite(budget.total_mass_all_bed_kg);
summary.final_total_mass_change_kg = last_finite(budget.total_mass_change_kg);
summary.final_total_mass_change_percent = last_finite(budget.total_mass_change_percent);
summary.total_mass_min_kg = finite_min(budget.total_mass_all_bed_kg);
summary.total_mass_max_kg = finite_max(budget.total_mass_all_bed_kg);
summary.max_abs_total_mass_change_kg = finite_max(abs(budget.total_mass_change_kg));
summary.max_abs_total_mass_change_percent = finite_max(abs(budget.total_mass_change_percent));
summary.max_mp1_negative_mass_kg = finite_max(budget.water_negative_mass_kg);
summary.max_bot_negative_mass_kg = finite_max(budget.bot_mass_negative_kg);
summary.cap_threshold_kgm3 = config.cap_threshold_kgm3;
summary.cap_tolerance_kgm3 = config.cap_tolerance_kgm3;
summary.max_mp1_conc_kgm3 = finite_max(budget.max_mp1_conc_kgm3);
summary.max_mp1_count_ge_cap = finite_max(budget.mp1_count_ge_cap);
summary.max_mp1_count_near_or_ge_cap = finite_max(budget.mp1_count_near_or_ge_cap);
summary.max_mp1_volume_near_or_ge_cap_m3 = finite_max(budget.mp1_volume_near_or_ge_cap_m3);
summary.grid_summary = grid_summary;
summary.has_split_water = has_split_water;
summary.has_split_bed = has_split_bed;

if has_split_water
    summary.max_abs_water_split_minus_total_kg = finite_max(abs(budget.water_split_minus_total_kg));
end

if has_split_bed
    summary.max_abs_bed_split_minus_total_kg = finite_max(abs(budget.bed_split_minus_total_kg));
end
end

function value = first_finite(values)
idx = find(isfinite(values), 1, 'first');
if isempty(idx)
    value = nan;
else
    value = values(idx);
end
end

function value = last_finite(values)
idx = find(isfinite(values), 1, 'last');
if isempty(idx)
    value = nan;
else
    value = values(idx);
end
end

function name = get_filename(path_value)
[~, base, ext] = fileparts(path_value);
name = [base ext];
end

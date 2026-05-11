% EXTRACT_MP1_MASS_BUDGET_B_C_30DAY
% Integrate FVCOM-MP floc-case mp1 mass diagnostics for experiments b and c
% across three consecutive output files (_0001, _0002, _0003) to cover a
% ~30-day period.
%
% This regular script reads one designated FVCOM output NetCDF file per case,
% computes compact domain-integrated mass-budget time series, and saves them
% to MAT files.
% It is designed for the floc-enabled cases where these variables exist:
%
%   mp1                  combined water-column plastic concentration
%   mp1_agg              aggregated water-column plastic concentration
%   mp1_dis              disaggregated water-column plastic concentration
%   bot_mass_mp1         combined bed plastic mass density
%   bot_mass_agg_mp1     aggregated bed plastic mass density
%   bot_mass_dis_mp1     disaggregated bed plastic mass density
%
% The key consistency checks are:
%
%   water_split_minus_mp1_kg =
%       integral((mp1_agg + mp1_dis - mp1) * volume)
%
%   bed_split_minus_mp1_kg =
%       integral((bot_mass_agg_mp1 + bot_mass_dis_mp1 - bot_mass_mp1) * area)
%
%   total_split_minus_combined_kg =
%       water_split_minus_mp1_kg + bed_split_minus_mp1_kg
%
% Mass accumulation procedure:
%   1. Node/layer water volume is art1 * max(h + zeta, 0) * sigma_layer_frac.
%   2. Water mass is sum(concentration * node/layer volume) over wet nodes.
%   3. Bed mass is sum(bottom_mass_density * art1) over all nodes.
%   4. Combined total mass is water(mp1) + bed(bot_mass_mp1).
%   5. Split total mass is water(mp1_agg + mp1_dis)
%      + bed(bot_mass_agg_mp1 + bot_mass_dis_mp1).
%   6. Optional dep/ero flux fields are area-integrated per output record.
%      Their cumulative sums are output-sampled diagnostics only unless every
%      model/plastic step is written to NetCDF.
%
% Example HPC usage:
%   cd /kfs3/scratch/yhuang168/waterPACT_MP_floc/FVCOM_MP_test_run/MATLAB
%   matlab -batch "run('extract_mp1_mass_budget_b_c_30day.m')"
%
% Optional environment overrides:
%   WATERPACT_CASE_IDS       Comma-separated case IDs. Default: b,c
%   WATERPACT_OUTPUT_ROOT    Root containing OUTPUT_b and OUTPUT_c
%   WATERPACT_OUTPUT_DIRS    Comma-separated explicit output directories
%   WATERPACT_OUT_MAT        Combined MAT file to write
%   WATERPACT_NC_FILES       Comma-separated exact file names or absolute paths
%   WATERPACT_TIME_STRIDE    Read every Nth output record. Default: 1
%   WATERPACT_SAVE_GRID      true/false for saving art1/h/sigma fractions
%   WATERPACT_SAVE_CASE_MATS true/false for saving one MAT per case

%% User Settings
%  Absolute Kestrel paths are set as defaults.  Override via environment
%  variables for portability to other systems.
case_ids = {'b', 'c'};
output_root = '/kfs3/scratch/yhuang168/waterPACT_MP_floc';
output_dirs = {};
nc_files = {};
out_mat = '';
time_stride = 1;
save_grid_metrics = true;
save_case_mats = true;
cap_threshold_kgm3 = 100.0;
cap_tolerance_kgm3 = 1.0e-6;
% File suffix labels to search per case.  Files that do not exist are
% skipped with a warning, so you can safely list all labels even if only
% a subset is available locally.
%   e.g. set WATERPACT_FILE_LABELS=0003 to use only the third file.
file_labels = {'0001', '0002', '0003'};

%% Environment Overrides
case_ids = getenv_list('WATERPACT_CASE_IDS', case_ids);
output_root = getenv_default('WATERPACT_OUTPUT_ROOT', output_root);
output_dirs = getenv_list('WATERPACT_OUTPUT_DIRS', output_dirs);
nc_files = getenv_list('WATERPACT_NC_FILES', nc_files);
out_mat = getenv_default('WATERPACT_OUT_MAT', out_mat);
time_stride = getenv_number('WATERPACT_TIME_STRIDE', time_stride);
save_grid_metrics = getenv_logical('WATERPACT_SAVE_GRID', save_grid_metrics);
save_case_mats = getenv_logical('WATERPACT_SAVE_CASE_MATS', save_case_mats);
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

if isempty(output_root)
    output_root = test_root;
end

case_ids = case_ids(:);
if isempty(case_ids)
    error('No cases requested. Set case_ids or WATERPACT_CASE_IDS.');
end

if ~isempty(output_dirs) && numel(output_dirs) ~= numel(case_ids)
    error('WATERPACT_OUTPUT_DIRS must have the same number of entries as WATERPACT_CASE_IDS.');
end

if isempty(nc_files)
    nc_files = default_nc_files_for_cases(case_ids);
elseif numel(nc_files) ~= numel(case_ids)
    error('WATERPACT_NC_FILES must have the same number of entries as WATERPACT_CASE_IDS.');
end

default_output_dir = fullfile(script_dir, 'output');
if ~exist(default_output_dir, 'dir')
    mkdir(default_output_dir);
end

if isempty(out_mat)
    out_mat = fullfile(default_output_dir, ['mp1_mass_budget_' strjoin(case_ids(:).', '_') '_30day.mat']);
end

config = struct();
config.case_ids = case_ids;
config.output_root = char(output_root);
config.output_dirs = output_dirs;
config.nc_files = nc_files;
config.out_mat = char(out_mat);
config.time_stride = max(1, round(time_stride));
config.save_grid_metrics = logical(save_grid_metrics);
config.save_case_mats = logical(save_case_mats);
config.cap_threshold_kgm3 = cap_threshold_kgm3;
config.cap_tolerance_kgm3 = cap_tolerance_kgm3;
config.created_by = mfilename;
config.created_on = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

diagnostics = cell(numel(case_ids), 1);
case_mat_files = cell(numel(case_ids), 1);

fprintf('Starting floc-case mp1 mass extraction for %d case(s).\n', numel(case_ids));
fprintf('Output MAT: %s\n\n', config.out_mat);

for icase = 1:numel(case_ids)
    case_id = char(case_ids{icase});
    if isempty(output_dirs)
        case_output_dir = fullfile(config.output_root, ['OUTPUT_' case_id]);
    else
        case_output_dir = char(output_dirs{icase});
    end
    % Build the candidate file list from file_labels, then keep only those
    % that actually exist on disk.  Missing files are skipped with a warning
    % so that a partial set (e.g. only _0003.nc available locally) works.
    case_nc_files = {};
    for ilab = 1:numel(file_labels)
        fname = sprintf('waterPACT_%s_%s.nc', case_id, char(file_labels{ilab}));
        fpath = fullfile(case_output_dir, fname);
        if isfile(fpath)
            case_nc_files{end+1} = fpath; %#ok<AGROW>
        else
            fprintf('  [skip] Not found: %s\n', fpath);
        end
    end

    fprintf('Case %s\n', case_id);
    fprintf('  NetCDF directory: %s\n', case_output_dir);
    if isempty(case_nc_files)
        warning('No NetCDF files found for case %s in %s. Skipping.', case_id, case_output_dir);
        continue;
    end
    for ilab = 1:numel(case_nc_files)
        fprintf('  File %d: %s\n', ilab, case_nc_files{ilab});
    end

    diagnostics{icase} = extract_one_case(case_id, case_output_dir, case_nc_files, ...
        config.time_stride, config.save_grid_metrics, ...
        config.cap_threshold_kgm3, config.cap_tolerance_kgm3);

    if config.save_case_mats
        case_diagnostic = diagnostics{icase};
        case_mat_files{icase} = fullfile(case_output_dir, ...
            ['mp1_mass_budget_' case_id '_30day.mat']);
        save(case_mat_files{icase}, 'case_diagnostic', '-v7.3');
        fprintf('  Saved per-case MAT: %s\n', case_mat_files{icase});
    end

    print_case_summary(diagnostics{icase}.summary);
    fprintf('\n');
end

summary_table = build_summary_table(diagnostics);

combined = struct();
combined.config = config;
combined.diagnostics = diagnostics;
combined.case_mat_files = case_mat_files;
combined.summary_table = summary_table;
combined.notes = { ...
    'Water volume = art1 * max(h + zeta, 0) * abs(diff(siglev)) on wet nodes.'; ...
    'water_mp1_mass_kg = sum(mp1 * water_volume).'; ...
    'water_split_sum_mass_kg = sum((mp1_agg + mp1_dis) * water_volume).'; ...
    'bed_mp1_mass_kg = sum(bot_mass_mp1 * art1) over all nodes.'; ...
    'bed_split_sum_mass_kg = sum((bot_mass_agg_mp1 + bot_mass_dis_mp1) * art1).'; ...
    'total_combined_mass_kg = water_mp1_mass_kg + bed_mp1_mass_kg.'; ...
    'total_split_mass_kg = water_split_sum_mass_kg + bed_split_sum_mass_kg.'; ...
    'Flux cumulative fields are output-sampled cumsum values, not exact continuous integrals unless every plastic step is written.'; ...
    'max_*_conc_kgm3 and *_count_near_or_ge_cap diagnose whether the 100 kg/m^3 split concentration cap appears in output.' ...
    };

out_dir = fileparts(config.out_mat);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
save(config.out_mat, 'combined', '-v7.3');

fprintf('Saved combined MAT diagnostic:\n  %s\n', config.out_mat);

%% Local Functions
function diagnostic = extract_one_case(case_id, output_dir, nc_files_list, time_stride, save_grid_metrics, cap_threshold_kgm3, cap_tolerance_kgm3)
if ischar(nc_files_list)
    nc_files_list = {nc_files_list};
end
% Filter to existing files (missing files were already warned about in the
% main loop; re-filter here as a safety net).
files = {};
for ilab = 1:numel(nc_files_list)
    if isfile(nc_files_list{ilab})
        files{end+1} = nc_files_list{ilab}; %#ok<AGROW>
    else
        warning('extract_one_case: skipping missing file: %s', nc_files_list{ilab});
    end
end
if isempty(files)
    error('No readable NetCDF files for case %s.', case_id);
end

fprintf('  Reading %d NetCDF file(s) (~30-day concatenation).\n', numel(files));

first_file = files{1};
vars_available = inspect_variables(first_file);
required_vars = { ...
    'art1', 'h', 'siglev', 'zeta', 'time', ...
    'mp1', 'mp1_agg', 'mp1_dis', ...
    'bot_mass_mp1', 'bot_mass_agg_mp1', 'bot_mass_dis_mp1' ...
    };
assert_required_vars(vars_available, required_vars, first_file);

art1 = double(ncread(first_file, 'art1'));
h = double(ncread(first_file, 'h'));
siglev = double(ncread(first_file, 'siglev'));
art1 = art1(:);
h = h(:);

if size(siglev, 1) ~= numel(h)
    error('siglev first dimension (%d) does not match node count (%d).', ...
        size(siglev, 1), numel(h));
end

dzfrac = abs(diff(siglev, 1, 2));
sigma_sum = sum(dzfrac, 2);

grid_summary = struct();
grid_summary.node_count = numel(h);
grid_summary.siglay_count = size(dzfrac, 2);
grid_summary.domain_area_m2 = nansum_local(art1);
grid_summary.h_min_m = finite_min(h);
grid_summary.h_max_m = finite_max(h);
grid_summary.h_area_mean_m = nansum_local(h .* art1) / max(nansum_local(art1), eps);
grid_summary.sigma_fraction_sum_min = finite_min(sigma_sum);
grid_summary.sigma_fraction_sum_max = finite_max(sigma_sum);
grid_summary.sigma_fraction_sum_max_abs_error = finite_max(abs(sigma_sum - 1.0));

has_wet_nodes = has_var(vars_available, 'wet_nodes');
has_dep_flux = has_var(vars_available, 'dep_flux_mp1');
has_ero_flux = has_var(vars_available, 'ero_flux_mp1');

record_plan = build_record_plan(files, time_stride);
nrec = numel(record_plan.file_index);
if nrec == 0
    error('No time records selected. Check WATERPACT_TIME_STRIDE.');
end

budget = initialize_budget(nrec, has_dep_flux, has_ero_flux);
fprintf('  Selected %d output record(s), stride = %d.\n', nrec, time_stride);

out_idx = 0;
for ifile = 1:numel(files)
    ncfile = files{ifile};
    local_records = record_plan.records_by_file{ifile};
    if isempty(local_records)
        continue;
    end

    fprintf('    %s: %d selected record(s)\n', get_filename(ncfile), numel(local_records));

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

        mp1 = read_2d_layer_field(ncfile, 'mp1', tindex);
        mp1_agg = read_2d_layer_field(ncfile, 'mp1_agg', tindex);
        mp1_dis = read_2d_layer_field(ncfile, 'mp1_dis', tindex);

        bot_mass = read_1d_time_field(ncfile, 'bot_mass_mp1', tindex);
        bot_agg = read_1d_time_field(ncfile, 'bot_mass_agg_mp1', tindex);
        bot_dis = read_1d_time_field(ncfile, 'bot_mass_dis_mp1', tindex);

        water_split = mp1_agg + mp1_dis;
        water_split_diff = water_split - mp1;
        bed_split = bot_agg + bot_dis;
        bed_split_diff = bed_split - bot_mass;

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

        budget.water_mp1_mass_kg(out_idx) = nansum_local(mp1(:) .* volume(:));
        budget.water_agg_mass_kg(out_idx) = nansum_local(mp1_agg(:) .* volume(:));
        budget.water_dis_mass_kg(out_idx) = nansum_local(mp1_dis(:) .* volume(:));
        budget.water_split_sum_mass_kg(out_idx) = ...
            budget.water_agg_mass_kg(out_idx) + budget.water_dis_mass_kg(out_idx);
        budget.water_split_minus_mp1_kg(out_idx) = ...
            nansum_local(water_split_diff(:) .* volume(:));
        budget.water_split_minus_mp1_rel(out_idx) = safe_scalar_ratio( ...
            budget.water_split_minus_mp1_kg(out_idx), budget.water_mp1_mass_kg(out_idx));

        budget.mp1_min_kgm3(out_idx) = finite_min(mp1(:));
        budget.mp1_max_kgm3(out_idx) = finite_max(mp1(:));
        budget.max_mp1_conc_kgm3(out_idx) = budget.mp1_max_kgm3(out_idx);
        budget.mp1_agg_min_kgm3(out_idx) = finite_min(mp1_agg(:));
        budget.mp1_agg_max_kgm3(out_idx) = finite_max(mp1_agg(:));
        budget.max_mp1_agg_conc_kgm3(out_idx) = budget.mp1_agg_max_kgm3(out_idx);
        budget.mp1_dis_min_kgm3(out_idx) = finite_min(mp1_dis(:));
        budget.mp1_dis_max_kgm3(out_idx) = finite_max(mp1_dis(:));
        budget.max_mp1_dis_conc_kgm3(out_idx) = budget.mp1_dis_max_kgm3(out_idx);
        budget.water_split_diff_min_kgm3(out_idx) = finite_min(water_split_diff(:));
        budget.water_split_diff_max_kgm3(out_idx) = finite_max(water_split_diff(:));
        budget.water_split_diff_max_abs_kgm3(out_idx) = finite_max(abs(water_split_diff(:)));
        budget.water_split_diff_rms_kgm3(out_idx) = finite_rms(water_split_diff(:));
        mp1_cap_mask = isfinite(mp1) & mp1 >= cap_threshold_kgm3;
        mp1_near_cap_mask = isfinite(mp1) & mp1 >= cap_threshold_kgm3 - cap_tolerance_kgm3;
        agg_cap_mask = isfinite(mp1_agg) & mp1_agg >= cap_threshold_kgm3;
        agg_near_cap_mask = isfinite(mp1_agg) & mp1_agg >= cap_threshold_kgm3 - cap_tolerance_kgm3;
        dis_cap_mask = isfinite(mp1_dis) & mp1_dis >= cap_threshold_kgm3;
        dis_near_cap_mask = isfinite(mp1_dis) & mp1_dis >= cap_threshold_kgm3 - cap_tolerance_kgm3;
        budget.mp1_count_ge_cap(out_idx) = nnz(mp1_cap_mask);
        budget.mp1_count_near_or_ge_cap(out_idx) = nnz(mp1_near_cap_mask);
        budget.mp1_volume_near_or_ge_cap_m3(out_idx) = nansum_local(volume(mp1_near_cap_mask));
        budget.mp1_agg_count_ge_cap(out_idx) = nnz(agg_cap_mask);
        budget.mp1_agg_count_near_or_ge_cap(out_idx) = nnz(agg_near_cap_mask);
        budget.mp1_agg_volume_near_or_ge_cap_m3(out_idx) = nansum_local(volume(agg_near_cap_mask));
        budget.mp1_dis_count_ge_cap(out_idx) = nnz(dis_cap_mask);
        budget.mp1_dis_count_near_or_ge_cap(out_idx) = nnz(dis_near_cap_mask);
        budget.mp1_dis_volume_near_or_ge_cap_m3(out_idx) = nansum_local(volume(dis_near_cap_mask));
        budget.mp1_nan_count(out_idx) = nnz(~isfinite(mp1(:)));
        budget.mp1_agg_nan_count(out_idx) = nnz(~isfinite(mp1_agg(:)));
        budget.mp1_dis_nan_count(out_idx) = nnz(~isfinite(mp1_dis(:)));
        budget.mp1_negative_count(out_idx) = nnz(mp1(:) < 0.0);
        budget.mp1_agg_negative_count(out_idx) = nnz(mp1_agg(:) < 0.0);
        budget.mp1_dis_negative_count(out_idx) = nnz(mp1_dis(:) < 0.0);
        budget.water_mp1_negative_mass_kg(out_idx) = -nansum_local(min(mp1, 0.0) .* volume);
        budget.water_agg_negative_mass_kg(out_idx) = -nansum_local(min(mp1_agg, 0.0) .* volume);
        budget.water_dis_negative_mass_kg(out_idx) = -nansum_local(min(mp1_dis, 0.0) .* volume);

        budget.bed_mp1_mass_kg(out_idx) = nansum_local(bot_mass .* art1);
        budget.bed_agg_mass_kg(out_idx) = nansum_local(bot_agg .* art1);
        budget.bed_dis_mass_kg(out_idx) = nansum_local(bot_dis .* art1);
        budget.bed_split_sum_mass_kg(out_idx) = ...
            budget.bed_agg_mass_kg(out_idx) + budget.bed_dis_mass_kg(out_idx);
        budget.bed_split_minus_mp1_kg(out_idx) = nansum_local(bed_split_diff .* art1);
        budget.bed_split_minus_mp1_rel(out_idx) = safe_scalar_ratio( ...
            budget.bed_split_minus_mp1_kg(out_idx), budget.bed_mp1_mass_kg(out_idx));
        budget.bot_mass_min_kgm2(out_idx) = finite_min(bot_mass);
        budget.bot_mass_max_kgm2(out_idx) = finite_max(bot_mass);
        budget.bot_agg_min_kgm2(out_idx) = finite_min(bot_agg);
        budget.bot_agg_max_kgm2(out_idx) = finite_max(bot_agg);
        budget.bot_dis_min_kgm2(out_idx) = finite_min(bot_dis);
        budget.bot_dis_max_kgm2(out_idx) = finite_max(bot_dis);
        budget.bed_split_diff_min_kgm2(out_idx) = finite_min(bed_split_diff);
        budget.bed_split_diff_max_kgm2(out_idx) = finite_max(bed_split_diff);
        budget.bed_split_diff_max_abs_kgm2(out_idx) = finite_max(abs(bed_split_diff));
        budget.bed_split_diff_rms_kgm2(out_idx) = finite_rms(bed_split_diff);
        budget.bot_mass_nan_count(out_idx) = nnz(~isfinite(bot_mass));
        budget.bot_agg_nan_count(out_idx) = nnz(~isfinite(bot_agg));
        budget.bot_dis_nan_count(out_idx) = nnz(~isfinite(bot_dis));
        budget.bot_mass_negative_count(out_idx) = nnz(bot_mass < 0.0);
        budget.bot_agg_negative_count(out_idx) = nnz(bot_agg < 0.0);
        budget.bot_dis_negative_count(out_idx) = nnz(bot_dis < 0.0);
        budget.bed_mp1_negative_mass_kg(out_idx) = -nansum_local(min(bot_mass, 0.0) .* art1);
        budget.bed_agg_negative_mass_kg(out_idx) = -nansum_local(min(bot_agg, 0.0) .* art1);
        budget.bed_dis_negative_mass_kg(out_idx) = -nansum_local(min(bot_dis, 0.0) .* art1);

        if has_dep_flux
            dep_flux = read_1d_time_field(ncfile, 'dep_flux_mp1', tindex);
            budget.dep_flux_area_integral_kg(out_idx) = nansum_local(dep_flux .* art1);
        end

        if has_ero_flux
            ero_flux = read_1d_time_field(ncfile, 'ero_flux_mp1', tindex);
            budget.ero_flux_area_integral_kg(out_idx) = nansum_local(ero_flux .* art1);
        end
    end
end

budget.water_mean_mp1_kgm3 = safe_divide(budget.water_mp1_mass_kg, budget.water_volume_m3);
budget.water_mean_split_kgm3 = safe_divide(budget.water_split_sum_mass_kg, budget.water_volume_m3);
budget.bed_mean_mp1_kgm2 = safe_divide(budget.bed_mp1_mass_kg, grid_summary.domain_area_m2);
budget.bed_mean_split_kgm2 = safe_divide(budget.bed_split_sum_mass_kg, grid_summary.domain_area_m2);

budget.total_combined_mass_kg = budget.water_mp1_mass_kg + budget.bed_mp1_mass_kg;
budget.total_split_mass_kg = budget.water_split_sum_mass_kg + budget.bed_split_sum_mass_kg;
budget.total_split_minus_combined_kg = budget.total_split_mass_kg - budget.total_combined_mass_kg;
budget.total_split_minus_combined_rel = safe_divide( ...
    budget.total_split_minus_combined_kg, budget.total_combined_mass_kg);

budget.water_mp1_change_kg = change_from_first(budget.water_mp1_mass_kg);
budget.water_split_change_kg = change_from_first(budget.water_split_sum_mass_kg);
budget.bed_mp1_change_kg = change_from_first(budget.bed_mp1_mass_kg);
budget.bed_split_change_kg = change_from_first(budget.bed_split_sum_mass_kg);
budget.total_combined_change_kg = change_from_first(budget.total_combined_mass_kg);
budget.total_split_change_kg = change_from_first(budget.total_split_mass_kg);
budget.total_combined_change_percent = percent_change_from_first(budget.total_combined_mass_kg);
budget.total_split_change_percent = percent_change_from_first(budget.total_split_mass_kg);
budget.total_combined_rate_kg_per_day = rate_from_time(budget.total_combined_mass_kg, budget.time_days);
budget.total_split_rate_kg_per_day = rate_from_time(budget.total_split_mass_kg, budget.time_days);

if has_dep_flux
    budget.sampled_cumulative_dep_flux_kg = cumsum_nan_as_zero(budget.dep_flux_area_integral_kg);
end
if has_ero_flux
    budget.sampled_cumulative_ero_flux_kg = cumsum_nan_as_zero(budget.ero_flux_area_integral_kg);
end
if has_dep_flux && has_ero_flux
    budget.sampled_cumulative_bed_exchange_gain_kg = ...
        cumsum_nan_as_zero(budget.dep_flux_area_integral_kg - budget.ero_flux_area_integral_kg);
    budget.sampled_cumulative_water_exchange_gain_kg = ...
        cumsum_nan_as_zero(budget.ero_flux_area_integral_kg - budget.dep_flux_area_integral_kg);
end

summary = make_summary(case_id, output_dir, files, budget, grid_summary, ...
    has_dep_flux, has_ero_flux, cap_threshold_kgm3, cap_tolerance_kgm3);

diagnostic = struct();
diagnostic.case_id = case_id;
diagnostic.output_dir = char(output_dir);
diagnostic.nc_file = files(:);
diagnostic.files = files(:);
diagnostic.vars_available = vars_available;
diagnostic.grid_summary = grid_summary;
diagnostic.budget = budget;
diagnostic.summary = summary;
diagnostic.notes = { ...
    'mp1 is the combined concentration field written by FVCOM-MP.'; ...
    'mp1_agg and mp1_dis are conc_a and conc_d. Their sum should reproduce mp1 after final recombination.'; ...
    'bot_mass_agg_mp1 and bot_mass_dis_mp1 are the split bed pools. Their sum should reproduce bot_mass_mp1 if bed split bookkeeping is closed.'; ...
    'All bed mass diagnostics use all nodes, not only wet nodes, because the bed reservoir is not a water-column volume.'; ...
    'sampled_cumulative_* flux diagnostics are cumsum over output records and should not be interpreted as exact cumulative exchange unless output frequency equals model update frequency.'; ...
    'max_mp1_agg_conc_kgm3, max_mp1_dis_conc_kgm3, and count_near_or_ge_cap fields diagnose whether the 100 kg/m^3 split concentration cap appears in output.' ...
    };

if save_grid_metrics
    grid = struct();
    grid.art1_m2 = art1;
    grid.h_m = h;
    grid.sigma_layer_fraction = dzfrac;
    diagnostic.grid = grid;
end
end

function budget = initialize_budget(nrec, has_dep_flux, has_ero_flux)
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
budget.wet_area_m2 = nan_col;
budget.zeta_min_m = nan_col;
budget.zeta_max_m = nan_col;
budget.zeta_area_mean_m = nan_col;

budget.water_mp1_mass_kg = nan_col;
budget.water_agg_mass_kg = nan_col;
budget.water_dis_mass_kg = nan_col;
budget.water_split_sum_mass_kg = nan_col;
budget.water_split_minus_mp1_kg = nan_col;
budget.water_split_minus_mp1_rel = nan_col;
budget.water_mean_mp1_kgm3 = nan_col;
budget.water_mean_split_kgm3 = nan_col;
budget.mp1_min_kgm3 = nan_col;
budget.mp1_max_kgm3 = nan_col;
budget.max_mp1_conc_kgm3 = nan_col;
budget.mp1_agg_min_kgm3 = nan_col;
budget.mp1_agg_max_kgm3 = nan_col;
budget.max_mp1_agg_conc_kgm3 = nan_col;
budget.mp1_dis_min_kgm3 = nan_col;
budget.mp1_dis_max_kgm3 = nan_col;
budget.max_mp1_dis_conc_kgm3 = nan_col;
budget.water_split_diff_min_kgm3 = nan_col;
budget.water_split_diff_max_kgm3 = nan_col;
budget.water_split_diff_max_abs_kgm3 = nan_col;
budget.water_split_diff_rms_kgm3 = nan_col;
budget.mp1_count_ge_cap = zero_col;
budget.mp1_count_near_or_ge_cap = zero_col;
budget.mp1_volume_near_or_ge_cap_m3 = nan_col;
budget.mp1_agg_count_ge_cap = zero_col;
budget.mp1_agg_count_near_or_ge_cap = zero_col;
budget.mp1_agg_volume_near_or_ge_cap_m3 = nan_col;
budget.mp1_dis_count_ge_cap = zero_col;
budget.mp1_dis_count_near_or_ge_cap = zero_col;
budget.mp1_dis_volume_near_or_ge_cap_m3 = nan_col;
budget.mp1_nan_count = zero_col;
budget.mp1_agg_nan_count = zero_col;
budget.mp1_dis_nan_count = zero_col;
budget.mp1_negative_count = zero_col;
budget.mp1_agg_negative_count = zero_col;
budget.mp1_dis_negative_count = zero_col;
budget.water_mp1_negative_mass_kg = nan_col;
budget.water_agg_negative_mass_kg = nan_col;
budget.water_dis_negative_mass_kg = nan_col;

budget.bed_mp1_mass_kg = nan_col;
budget.bed_agg_mass_kg = nan_col;
budget.bed_dis_mass_kg = nan_col;
budget.bed_split_sum_mass_kg = nan_col;
budget.bed_split_minus_mp1_kg = nan_col;
budget.bed_split_minus_mp1_rel = nan_col;
budget.bed_mean_mp1_kgm2 = nan_col;
budget.bed_mean_split_kgm2 = nan_col;
budget.bot_mass_min_kgm2 = nan_col;
budget.bot_mass_max_kgm2 = nan_col;
budget.bot_agg_min_kgm2 = nan_col;
budget.bot_agg_max_kgm2 = nan_col;
budget.bot_dis_min_kgm2 = nan_col;
budget.bot_dis_max_kgm2 = nan_col;
budget.bed_split_diff_min_kgm2 = nan_col;
budget.bed_split_diff_max_kgm2 = nan_col;
budget.bed_split_diff_max_abs_kgm2 = nan_col;
budget.bed_split_diff_rms_kgm2 = nan_col;
budget.bot_mass_nan_count = zero_col;
budget.bot_agg_nan_count = zero_col;
budget.bot_dis_nan_count = zero_col;
budget.bot_mass_negative_count = zero_col;
budget.bot_agg_negative_count = zero_col;
budget.bot_dis_negative_count = zero_col;
budget.bed_mp1_negative_mass_kg = nan_col;
budget.bed_agg_negative_mass_kg = nan_col;
budget.bed_dis_negative_mass_kg = nan_col;

budget.total_combined_mass_kg = nan_col;
budget.total_split_mass_kg = nan_col;
budget.total_split_minus_combined_kg = nan_col;
budget.total_split_minus_combined_rel = nan_col;
budget.water_mp1_change_kg = nan_col;
budget.water_split_change_kg = nan_col;
budget.bed_mp1_change_kg = nan_col;
budget.bed_split_change_kg = nan_col;
budget.total_combined_change_kg = nan_col;
budget.total_split_change_kg = nan_col;
budget.total_combined_change_percent = nan_col;
budget.total_split_change_percent = nan_col;
budget.total_combined_rate_kg_per_day = nan_col;
budget.total_split_rate_kg_per_day = nan_col;

if has_dep_flux
    budget.dep_flux_area_integral_kg = nan_col;
    budget.sampled_cumulative_dep_flux_kg = nan_col;
end
if has_ero_flux
    budget.ero_flux_area_integral_kg = nan_col;
    budget.sampled_cumulative_ero_flux_kg = nan_col;
end
if has_dep_flux && has_ero_flux
    budget.sampled_cumulative_bed_exchange_gain_kg = nan_col;
    budget.sampled_cumulative_water_exchange_gain_kg = nan_col;
end
end

function summary = make_summary(case_id, output_dir, files, budget, grid_summary, ...
    has_dep_flux, has_ero_flux, cap_threshold_kgm3, cap_tolerance_kgm3)
summary = struct();
summary.case_id = case_id;
summary.output_dir = char(output_dir);
summary.file_count = numel(files);
summary.record_count = numel(budget.time_days);
summary.time_start_days = first_finite(budget.time_days);
summary.time_end_days = last_finite(budget.time_days);
summary.initial_total_combined_mass_kg = first_finite(budget.total_combined_mass_kg);
summary.final_total_combined_mass_kg = last_finite(budget.total_combined_mass_kg);
summary.final_total_combined_change_kg = last_finite(budget.total_combined_change_kg);
summary.final_total_combined_change_percent = last_finite(budget.total_combined_change_percent);
summary.initial_total_split_mass_kg = first_finite(budget.total_split_mass_kg);
summary.final_total_split_mass_kg = last_finite(budget.total_split_mass_kg);
summary.final_total_split_change_kg = last_finite(budget.total_split_change_kg);
summary.final_total_split_change_percent = last_finite(budget.total_split_change_percent);
summary.max_abs_total_split_minus_combined_kg = finite_max(abs(budget.total_split_minus_combined_kg));
summary.max_abs_water_split_minus_mp1_kg = finite_max(abs(budget.water_split_minus_mp1_kg));
summary.max_abs_bed_split_minus_mp1_kg = finite_max(abs(budget.bed_split_minus_mp1_kg));
summary.max_abs_water_split_diff_kgm3 = finite_max(budget.water_split_diff_max_abs_kgm3);
summary.max_abs_bed_split_diff_kgm2 = finite_max(budget.bed_split_diff_max_abs_kgm2);
summary.max_water_negative_mass_kg = finite_max([ ...
    budget.water_mp1_negative_mass_kg; ...
    budget.water_agg_negative_mass_kg; ...
    budget.water_dis_negative_mass_kg]);
summary.max_bed_negative_mass_kg = finite_max([ ...
    budget.bed_mp1_negative_mass_kg; ...
    budget.bed_agg_negative_mass_kg; ...
    budget.bed_dis_negative_mass_kg]);
summary.has_dep_flux = has_dep_flux;
summary.has_ero_flux = has_ero_flux;
summary.cap_threshold_kgm3 = cap_threshold_kgm3;
summary.cap_tolerance_kgm3 = cap_tolerance_kgm3;
summary.max_mp1_conc_kgm3 = finite_max(budget.max_mp1_conc_kgm3);
summary.max_mp1_agg_conc_kgm3 = finite_max(budget.max_mp1_agg_conc_kgm3);
summary.max_mp1_dis_conc_kgm3 = finite_max(budget.max_mp1_dis_conc_kgm3);
summary.max_mp1_count_ge_cap = finite_max(budget.mp1_count_ge_cap);
summary.max_mp1_count_near_or_ge_cap = finite_max(budget.mp1_count_near_or_ge_cap);
summary.max_mp1_agg_count_ge_cap = finite_max(budget.mp1_agg_count_ge_cap);
summary.max_mp1_agg_count_near_or_ge_cap = finite_max(budget.mp1_agg_count_near_or_ge_cap);
summary.max_mp1_dis_count_ge_cap = finite_max(budget.mp1_dis_count_ge_cap);
summary.max_mp1_dis_count_near_or_ge_cap = finite_max(budget.mp1_dis_count_near_or_ge_cap);
summary.grid_summary = grid_summary;
end

function print_case_summary(summary)
fprintf('  Time window: %.6g to %.6g days (%d records)\n', ...
    summary.time_start_days, summary.time_end_days, summary.record_count);
fprintf('  Combined total initial/final: %.12g -> %.12g kg\n', ...
    summary.initial_total_combined_mass_kg, summary.final_total_combined_mass_kg);
fprintf('  Combined total change: %.12g kg (%.6g %%)\n', ...
    summary.final_total_combined_change_kg, summary.final_total_combined_change_percent);
fprintf('  Split total initial/final: %.12g -> %.12g kg\n', ...
    summary.initial_total_split_mass_kg, summary.final_total_split_mass_kg);
fprintf('  Split total change: %.12g kg (%.6g %%)\n', ...
    summary.final_total_split_change_kg, summary.final_total_split_change_percent);
fprintf('  Max |split - combined| total: %.12g kg\n', ...
    summary.max_abs_total_split_minus_combined_kg);
fprintf('  Max |mp1_agg + mp1_dis - mp1| water: %.12g kg\n', ...
    summary.max_abs_water_split_minus_mp1_kg);
fprintf('  Max |bot split - bot total| bed: %.12g kg\n', ...
    summary.max_abs_bed_split_minus_mp1_kg);
end

function summary_table = build_summary_table(diagnostics)
n = numel(diagnostics);
case_id = cell(n, 1);
record_count = nan(n, 1);
time_start_days = nan(n, 1);
time_end_days = nan(n, 1);
final_total_combined_change_kg = nan(n, 1);
final_total_combined_change_percent = nan(n, 1);
max_abs_total_split_minus_combined_kg = nan(n, 1);
max_abs_water_split_minus_mp1_kg = nan(n, 1);
max_abs_bed_split_minus_mp1_kg = nan(n, 1);

for i = 1:n
    s = diagnostics{i}.summary;
    case_id{i} = s.case_id;
    record_count(i) = s.record_count;
    time_start_days(i) = s.time_start_days;
    time_end_days(i) = s.time_end_days;
    final_total_combined_change_kg(i) = s.final_total_combined_change_kg;
    final_total_combined_change_percent(i) = s.final_total_combined_change_percent;
    max_abs_total_split_minus_combined_kg(i) = s.max_abs_total_split_minus_combined_kg;
    max_abs_water_split_minus_mp1_kg(i) = s.max_abs_water_split_minus_mp1_kg;
    max_abs_bed_split_minus_mp1_kg(i) = s.max_abs_bed_split_minus_mp1_kg;
end

summary_table = table(case_id, record_count, time_start_days, time_end_days, ...
    final_total_combined_change_kg, final_total_combined_change_percent, ...
    max_abs_total_split_minus_combined_kg, ...
    max_abs_water_split_minus_mp1_kg, max_abs_bed_split_minus_mp1_kg);
end

function values = read_2d_layer_field(ncfile, varname, tindex)
values = double(ncread(ncfile, varname, [1 1 tindex], [Inf Inf 1]));
values = squeeze(values);
end

function values = read_1d_time_field(ncfile, varname, tindex)
values = double(ncread(ncfile, varname, [1 tindex], [Inf 1]));
values = values(:);
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

function value = finite_rms(values)
values = values(:);
values = values(isfinite(values));
if isempty(values)
    value = nan;
else
    value = sqrt(mean(values .^ 2));
end
end

function out = safe_divide(num, den)
out = nan(size(num));
if isscalar(den)
    mask = isfinite(num) & isfinite(den) & abs(den) > 0.0;
    out(mask) = num(mask) ./ den;
else
    mask = isfinite(num) & isfinite(den) & abs(den) > 0.0;
    out(mask) = num(mask) ./ den(mask);
end
end

function value = safe_scalar_ratio(num, den)
if isfinite(num) && isfinite(den) && abs(den) > 0.0
    value = num / den;
else
    value = nan;
end
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

function out = cumsum_nan_as_zero(values)
tmp = values;
tmp(~isfinite(tmp)) = 0.0;
out = cumsum(tmp);
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

function value = getenv_default(name, default_value)
value = getenv(name);
if isempty(value)
    value = default_value;
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

function name = get_filename(path_value)
[~, base, ext] = fileparts(path_value);
name = [base ext];
end

function nc_files = default_nc_files_for_cases(case_ids)
nc_files = cell(size(case_ids));
for i = 1:numel(case_ids)
    nc_files{i} = sprintf('waterPACT_%s_0001.nc', char(case_ids{i}));
end
end

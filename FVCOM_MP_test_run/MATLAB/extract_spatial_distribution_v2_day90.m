% Extract case c/d day-90 diurnal-average spatial fields for FVCOM-MP v2.
script_dir = fileparts(mfilename('fullpath'));
if ~isempty(script_dir)
    addpath(script_dir);
end
extract_spatial_distribution_v2_core(90, 9);

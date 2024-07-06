%% [Step 1] Call function & parameters - Step 1
clearvars; close all;
%==========================================================================
matlab_package_path = '/home/tfchengac/';   %/ Set the path of the folder where you store the matlab packages
%==========================================================================
addpath(genpath(fullfile(matlab_package_path,'MyMatlabPkg')));
addpath(genpath(fullfile(matlab_package_path,'MyMatlabFn')));
addpath(genpath(fullfile(matlab_package_path,'m_map1.4o')));
addpath(genpath(fullfile(matlab_package_path,'MCS_Codes')));
addpath(genpath(fullfile(matlab_package_path,'ITCC Codes')));
addpath(genpath(fullfile(matlab_package_path,'MJO_Codes')));
addpath(genpath(fullfile(matlab_package_path,'Flexpart_Codes')));
%/ Read nctoolbox for reading grib file, see also the website
%/ for solving the bug: https://www.mathworks.com/matlabcentral/answers/1635900-failed-to-setup-the-java-classpath
% setup_nctoolbox;

%/ Remove the path of CCToolbox to avoid confliction of the kmeans function!
rmpath(genpath('/home/tfchengac/MyMatlabPkg/CCToolbox/')); 

param = 'param_universal.m';
run(param); 

dataset = []; dataset.placeholder = [];

project_name    = 'CMIP6_MJO';
masterfolder    = '/disk/r059/tfchengac/CMIP6_MJO/';
data_folder     = strcat(masterfolder, 'prcssd_data_4plotting/');
plotting_folder = strcat(masterfolder, 'plots/');
year_list       = 1979:2014;   %/ 36 years but 35 cold seasons (Wang et al 2019)
str_years       = sprintf('_%d-%d', year_list(1), year_list(end));

%/ Load cmip6 models (good/poor)
model_list_all = {'ACCESS-CM2',    'ACCESS-ESM1-5', 'CESM2',         'CMCC-CM2-SR5',...
                  'CMCC-ESM2',     'CNRM-CM6-1',    'CNRM-CM6-1-HR', 'CNRM-ESM2-1',      'CanESM5',...
                  'EC-Earth3',     'EC-Earth3-CC',  'EC-Earth3-Veg', 'EC-Earth3-Veg-LR', 'FGOALS-g3',...
                  'GFDL-CM4',      'GISS-E2-1-G',   'IITM-ESM',      'INM-CM4-8',        'INM-CM5-0',...
                  'IPSL-CM6A-LR',  'MIROC-ES2H',    'MIROC-ES2L',    'MIROC6',           'MPI-ESM1-2-HR',...
                  'MPI-ESM1-2-LR', 'MRI-ESM2-0',    'NESM3',         'TaiESM1'};   %/ 28 models with full daily rlut data for historical, ssp245, ssp585

%==========================================================================
eval_scheme      = 7;
MJO_loc          = 'EIO';
select_field_MJO = 'daily_MJO_bw';  %/ The timescale for MJO events, but not necessarily for the variables to be plotted
%==========================================================================
str_eval_scheme  = sprintf('_EvalSch%d',eval_scheme);
interp_res       = 2;
hori_new         = 50:interp_res:160;
timelag_new      = -25:25;
obs_data         = {'NOAA_Interp'}; %/ Benchmark dataset
nmodel           = 28;

model_eval_filename = fullfile(data_folder, sprintf('model_eval%s_lon%d-%d_lag%d-%d_%s_against_%dmodels%s%s.mat',...
                               str_years, hori_new(1), hori_new(end), timelag_new(1), timelag_new(end), obs_data{:}, nmodel, str_years, str_eval_scheme));

if isfile(model_eval_filename)
    fprintf('*** Loading the model_eval data: %s ...***\n', model_eval_filename);
    model_eval = par_load(model_eval_filename, 'model_eval');
    disp('Poor models:')
    disp(model_eval.poor_models); 
    disp('Good models:')
    disp(model_eval.good_models);

    ens_list = read_CMIP_ens('project_name', project_name, 'model_list', model_list_all, 'var', 'rlut', 'exp', 'historical');
    if isequal(model_eval.ens_list, ens_list)   %/ [IMPORTANT] To verify if ens_list of the models stay the same
        fprintf('Verifying model_eval by ens_list: Most updated.')
    else
        warning('Verifying model_eval by ens_list: Not updated!')
        model_eval = [];  %/ Recalculate model_eval is needed [go to Step 5]
    end
else
    warning('!!! model_eval file %s not found !!!\n', model_eval_filename);
end

server_name      = getenv('HOSTNAME');          %/ get the server name
server_name_cell = strsplit(server_name, '.');
server_name      = server_name_cell{1};

%------------------------ Clim Study (1979-2020) -------------------------%
date_yyyymmdd_AllYr = date_array_gen('year_list', year_list, 'st_month', 1, 'st_day', 1, 'ed_month', 12, 'ed_day', 31, 'output_date_format', 'yyyymmdd');
fprintf('*** Year Period: %s ***\n', strrep(str_years, '_', ''));

date_mmdd_OneYr = date_array_gen('year_list', 1979, 'st_month', 1, 'st_day', 1, 'ed_month', 12, 'ed_day', 31, 'output_date_format', 'yyyymmdd');
date_mmdd_OneYr = mod(date_mmdd_OneYr, 1e4);  %/ 365 days

%/ Create string cells for pentads and the corresponding months.
str_mth = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
str_ptd_dates = cell(length(date_mmdd_OneYr)/5, 2);             %/ Skipped the leap day, as LinHo and Wang 2002 did.
for t = 1:length(str_ptd_dates)
    ptd_dates = date_mmdd_OneYr((1:5)+5*(t-1));
    str_ptd_dates{t, 1} = sprintf('P%02d_%d-%d', t, ptd_dates(1), ptd_dates(end));
    str_ptd_dates{t, 2} = str_mth{floor(ptd_dates(3)/1e2)};     %/ Record the month of the middle day in a given pentad. (fairest)
end

repmat_year = repmat(year_list, length(str_ptd_dates), 1);
repmat_year = reshape(repmat_year, [], 1);   
repmat_ptd  = repmat([1:length(str_ptd_dates)], 1, length(year_list))';
str_ptd_dates_AllYr = repmat(str_ptd_dates, length(year_list), 1);              %/ Convert into string for strcat().
str_ptd_dates_AllYr(:,1) = cellstr(strcat(string(repmat_year), '_', str_ptd_dates_AllYr(:,1)));     %/ cellstr(): convert string to cell.
str_ptd_dates_AllYr(:,3) = num2cell(repmat_year*1e2+repmat_ptd);           %/ For later convenience.
license('inuse');

%% [Step 2.2] read_CMIP 
fprintf('*** REMINDER: You may run read_CMIP_program.m in the background. ***\n')
model_list = []; exp_list = []; slct_ens = []; ori_field = []; MME_list = []; MME_name = []; compute_MME = 0; slct_level_range = [];
% CMIPdata = []; 

savemat             = 1;
recompute           = 0; %<- mind this
recompute_MME       = 0; %<- mind this
load_or_not         = 0; %<- mind this
delete_prob_files   = 0; %<- auto-delete all problmatic files
CMIP_folder         = {'/disk/v183.a/share/CMIP6_data'};        %/ CMIP folder
output_folder       = data_folder; 
output_field        = 'daily';   %/ 'daily', 'monthly'
units_level         = 'hPa';
MME_lon_res         = 2;   MME_lat_res   = 2;
interp_method       = 'linear';
NumWorkers          = 20;
noleap              = 1;         %/ toggle this on to ensure the lengths of daily data across all models are consistent!

% model_list = model_list_all;
model_list = model_eval.good_models;
% model_list = model_eval.poor_models;
% model_list = 'CESM2';
% model_list = 'GFDL-CM4';

%/ [WARNING] The signal will be much weaker for MME, do not do this!
% model_list = model_eval.good_models; compute_MME = 1; MME_name = sprintf('GoodMME%d', length(model_list)); 
% model_list = model_eval.poor_models; compute_MME = 1; MME_name = sprintf('PoorMME%d', length(model_list)); 

% exp_list = {'historical'}; slct_year_list = 1979:2014; 
% exp_list = {'historical', 'ssp585'}; slct_year_list = {1979:2014, 2064:2099}; %/ Some do not have output in 2100.
% exp_list = {'ssp245','ssp585'}; slct_year_list = {2064:2099, 2064:2099}; %/ Some do not have output in 2100.
% exp_list = {'historical', 'ssp585', 'ssp585'}; slct_year_list = {1979:2014, 2064:2099, 2015:2099}; %/ Some do not have output in 2100.
exp_list = {'historical', 'ssp245', 'ssp245', 'ssp585', 'ssp585'}; slct_year_list = {1979:2014, 2028:2063, 2064:2099, 2028:2063, 2064:2099}; %/ Some do not have output in 2100.

% var_list = {'rlut'};
% var_list = {'rlut'}; ori_field = 'daily'; output_field = 'monthly';
% var_list = {'Lcpr'};
% var_list = {'MC1000to850'}; output_field = 'daily'; %/ 1000-850 moisture convergence
% var_list = {'MC1000to850'}; output_field = 'monthly'; %/ 1000-850 moisture convergence

% pt_BL            = 850;
% pb_BL            = 1000;
% % var_list         = {'ua', 'va', 'hus'};  output_field = 'daily';
% var_list         = {'ua', 'va', 'hus'};  output_field = 'monthly';
% slct_level_range = {[pt_BL, pb_BL],[pt_BL, pb_BL],[pt_BL, pb_BL]};

% var_list = {'MC925'};  %/ 925-hPa moisture convergence
% var_list = {'ua', 'va', 'hus'}; slct_level_range = num2cell([850, 850, 850]);
% var_list = {'ua', 'va'}; slct_level_range = num2cell([850, 850]);
% var_list = {'hus', 'ta', 'zg'}; slct_level_range = repmat({[100, 1000]}, 1, length(var_list)); output_field = 'monthly';
% var_list = {'wap', 'hus','ua','va','ua','va'}; slct_level_range = num2cell([500, 850, 850, 850, 250, 250]); output_field = 'monthly';
% var_list = {'pr', 'rlut'};
% var_list = {'S500'}; output_field = 'monthly';
% var_list = {'hus', 'hus', 'ta', 'ta', 'ta'}; slct_level_range = num2cell([1000, 850, 250, 500, 700]); output_field = 'monthly'; %/ static stability at 500hPa
% var_list = {'S2_Wang1988', 'C0_Wang1988', 'alpha_Wang1988', 'q3_bar_Wang1988', 'q1_bar_Wang1988', 'dq_bar_Wang1988', 'I_Wang1988', 'C1_Wang1988'};  output_field = 'monthly';
% var_list = {'hus'}; slct_level_range = {[100, 1000]}; output_field = 'monthly';
% var_list = {'ta'}; slct_level_range = {[300]}; output_field = 'monthly';

% var_list = {'tos', 'pr', 'evspsbl', 'tas'};  output_field = 'monthly';
% var_list = {'ua','va'}; slct_level_range = num2cell([850, 850])
var_list = {'tos'}; output_field = 'monthly';
% var_list = {'tos'}; 
% var_list = {'hus500wap300'};
% var_list = {'wap'}; slct_level_range = num2cell([250]) %/ no hist-GHG, hist-aer, hist-nat for GFDL-ESM4!!!
% var_list = {'od550aer'};
% var_list = {'tas'};
% var_list = {'SS'};  %/ static stability (ta850 - ta500)
% var_list = {'SS2'};  %/ static stability (ta500 - ta300)
% var_list = {'div'}; slct_level_range = num2cell([500]) 

% clc;
CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model_list, 'var_list', var_list,...
                     'exp_list', exp_list, 'slct_year_list', slct_year_list, 'noleap', noleap, 'slct_level_range', slct_level_range, 'units_level', units_level,...
                     'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                     'compute_MME', compute_MME, 'MME_name', MME_name, 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                     'savemat', savemat, 'recompute', recompute, 'recompute_MME', recompute_MME, 'load_or_not', load_or_not, 'delete_prob_files', delete_prob_files);

% ncdisp(fullfile(CMIP_folder{2}, 'tos_Oday_EC-Earth3-CC_historical_r1i1p1f1_gn_18560101-18561231.nc'))
% unique(CMIPdata.ACCESS_ESM1_5.tos_historical)

% a = CMIPdata.ACCESS_ESM1_5.tos_historical(:,:,100);
% close all
% figure
% imagesc(flip(a', 1))

%/ Check CMIP6 model resolution
% for m = 1:length(model_list)
%     model = model_list{m};
%     model_strrep = strrep(model, '-', '_');
%     lon = CMIPdata.(model_strrep).lon;
%     lat = CMIPdata.(model_strrep).lat;
%     lon_res = abs(diff(lon(1:2)));
%     lat_res = abs(diff(lat(1:2)));
%     fprintf('%s (lonxlat): %.3f%sx%.3f%s (%d x %d) \n', model, lon_res, char(176), lat_res, char(176), length(lon), length(lat));
% end

%%
exit;
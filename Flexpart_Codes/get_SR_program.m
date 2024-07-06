%% [Step 1] Read data and FLEXPART-WaterSip experiments parameters (see also moisture_tracking.m)
clearvars;

%==========================================================================
matlab_package_path = '/home/tfchengac/';   %/ Set the path of the folder where you store the matlab packages
%==========================================================================

addpath(genpath(fullfile(matlab_package_path,'MyMatlabPkg')));
addpath(genpath(fullfile(matlab_package_path,'MyMatlabFn')));
addpath(genpath(fullfile(matlab_package_path,'m_map1.4o')));
addpath(genpath(fullfile(matlab_package_path,'MCS_Codes')));
addpath(genpath(fullfile(matlab_package_path,'Flexpart_Codes')));
WSV = []; dataset = [];  %/ clear up all to avoid using old data (e.g., sometimes mixing P in 2010 with P_LA in 1971-2010).

param = 'param_universal';
% param = 'param_flexpart';
run(param);
project_name     = 'TP';
masterfolder     = '/disk/r059/tfchengac/FLEXPART/';
server_name      = getenv('HOSTNAME');          %/ get the server name
server_name_cell = strsplit(server_name, '.');
server_name      = server_name_cell{1};
if ismember(server_name, {'hqlx145', 'hqlx146', 'hqlx147', 'hqlx148', 'hqlx149'})
    NumWorkers = 40;   %/ optimal for 128-core (better than 80 or 128).
elseif ismember(server_name, {'hqlx128'})
    NumWorkers = 10;
else
    NumWorkers = 20;   %/ optimal for 48-core     
end

%================ global or basin-wise / RHc_dqc_scheme ==================%
%/ 'from_basin':   0]: from global land
%/               1-3]: from global source hotspots
%/                 4]: from TP basins 
%/                 5]: from TP low-veg covers
%/                 6]: from TP grids (0.5x0.5)
%/                 7]: from Pakistan box region
%/     'expmnt': The name of your FLEXPART experiment (i.e., the folder name)
%/ 'output_res': The output WSV spatial resolution (in deg)
%/         'dt': FLEXPART time interval (in hr)

%/ See also 'read_RHc_dqc_scheme.m' and [Step 5.1]
expmnt = 'domfill_CERA_MPI'; ldirect = 'bwd'; from_basin = 6; output_res = 1; dt = 3; year_list = 1971:2010; optimal_rr = 0.99; RHc_dqc_scheme = 23; %/ *LATEST*
% expmnt = 'domfill_CERA_MPI'; ldirect = 'bwd'; from_basin = 5; output_res = 1; dt = 3; year_list = 1971:2010; optimal_rr = 0.99; RHc_dqc_scheme = 23; 
% expmnt = 'domfill_CERA_MPI'; ldirect = 'bwd'; from_basin = 4;  output_res = 1; dt = 3; year_list = 1971:2010; optimal_rr = 0.99; RHc_dqc_scheme = 23; 
% expmnt = 'domfill_CERA_MPI'; ldirect = 'bwd'; from_basin = 0; output_res = 1; dt = 3; year_list = 2010;      optimal_rr = 0.9;  RHc_dqc_scheme = 10; %/ The Base Scheme    

[str_RHc_dqc, str_remark, str_src_z, str_years, WSV_lon, WSV_lat,...
 maxtraj_day, str_BLH_factor, str_optimal, str_traj_rm_jump, dt_slct, data_folder,...
 plotting_folder, WSV_dir, str_expmntinfo, basin_catalog, str_domain, str_domain_trajfile, str_sharpcut, forcing] = ...
    load_tracking_param('expmnt', expmnt, 'ldirect', ldirect, 'output_res', output_res, 'dt', dt, 'year_list', year_list,...
                        'optimal_rr', optimal_rr, 'RHc_dqc_scheme', RHc_dqc_scheme, 'from_basin', from_basin, 'masterfolder', masterfolder);

% [str_RHc_dqc, ldirect, str_remark, optimal_tracking, traj_rm_jump, BLH_factor, str_src_z,...
%           str_years, WSV_lon, WSV_lat, maxtraj_day, str_BLH_factor, str_optimal, str_traj_rm_jump, dt_slct, data_folder,...
%           plotting_folder, WSV_dir, str_expmntinfo, basin_catalog, str_domain, ...
%           ldirect_fwd, str_remark_fwd, maxtraj_day_fwd, WSV_dir_fwd, str_expmntinfo_fwd,...
%           global_basin_2D, global_basin_bndry_list, global_basin_name_list, global_basin_centr_list,...
%           cond_TP, TP_bndry, TP_centr, dataset] = param_watersip('expmnt', expmnt, 'output_res', output_res, 'dt', dt, 'year_list', year_list, 'optimal_rr', optimal_rr,...
%                                                                  'RHc_dqc_scheme', RHc_dqc_scheme, 'from_basin', from_basin,...
%                                                                  'dataset', dataset, 'masterfolder', masterfolder);

%======================== Regime Classification ==========================%
regime_ver = 'V2';  %/   '': Original 
                    %/ 'V2': Integrate easterly in North India into Indian Monsoon (only JJA)
%==========================================================================

%/ Load the global basins first (save time for computing SR network if needed)
recompute = 0;
slct_reg  = 'hydrosheds_TPbasins_oceans';

[global_basin_2D, global_basin_bndry_list, global_basin_name_list, ~, global_basin_centr_list] = reg_extractor('slct_reg', slct_reg, 'lon', WSV_lon,...
                                           'lat', WSV_lat(2:end-1), 'data_folder', data_folder, 'savemat', 1, 'recompute', recompute, 'verbose', 1);

%/ Load the TP basin boundary
slct_reg = 'TP';
[cond_TP, TP_bndry, ~, ~, TP_centr] = reg_extractor('slct_reg', slct_reg, 'lon', WSV_lon, 'lat', WSV_lat(2:end-1),...
                                            'data_folder', data_folder, 'savemat', 1, 'recompute', recompute, 'verbose', 1);
cond_TP(isnan(cond_TP)) = 0;
cond_TP = logical(cond_TP);

%/ Load Land Cover Type from CERA-20C
if contains(expmnt, '_CERA_')
    LC_name = {'HiVeg', 'LoVeg'};
    for k = 1:length(LC_name)
        if isfield(dataset, LC_name{k})       %/ skip it if it has been loaded.
            disp(['!!! ', LC_name{k}, ' has already been loaded. Skip data retrieval. !!!'])
        else
            select_field = 'daily';
            domain_str = '_global';

            loadfile = strcat(data_folder, LC_name{k}, '_', select_field, domain_str, '.mat');
            disp(['Loading ', loadfile, ' ...']);
            load(loadfile, 'S');

            fld = fieldnames(S);
            for f = 1:length(fld)
                dataset.(LC_name{k}).(fld{f}) = S.(fld{f});
            end
            clear S;
        end
    end
end

%/ Other variables to load (for consistency)
labels_TP_basins        = {'TP-Indus',   'TP-Ganges', 'TP-Irrawaddy', 'TP-Mekong',...
                           'TP-Yangtze', 'TP-Yellow', 'TP-Gobi1', 'TP-Tarim', 'TP-Syr Darya', 'TP-Amu Darya',...
                           'TP-CapsianSea East Coast', 'TP-Helmand',...
                           'S. Inner TP', 'N. Inner TP', 'Qaidam'};

%/ Mid-lower basins (in connection w/ TP basins)
labels_ML_basins        = {'ML-Indus',   'ML-Ganges', 'ML-Irrawaddy', 'ML-Mekong',...
                           'ML-Yangtze', 'ML-Yellow', 'ML-Gobi1', 'ML-Tarim', 'ML-Syr Darya', 'ML-Amu Darya',...
                           'ML-CapsianSea East Coast', 'ML-Helmand'};
% labels_domain_wise       = {'Pan_IM_noTP', 'Pan_Westerly_noTP', 'Pan_EAM_noTP'};
% labels_domain_wise2      = {'Pan_IM_noTP', 'Pan_EAM_noTP'};

labels_domain_wise_exact = {'NH_MidLat_Westerly_noTP', 'NH_Polar_Westerly_noTP',...
                            'IM_noTP', 'EAM_noTP', 'Overturning_elsewhere_noTP', 'Easterly_noTP', 'SH_Westerly', 'Elsewhere_noTP'};
col_P             = [ 60  60  60]./255;
col_E             = [255 102   0]./255;
col_T2m           = [204  51 153]./255;
col_TP            = [186 215 233]./255;
col_local         = [72 152  177]./255;
col_nonlocalTP    = [153 235 245]./255;   % [218 240 243]./255;
col_Westerly      = [235 69   95]./255;
col_IM            = [43  52  103]./255;
col_EAM           = [242 242  46]./255;   % [252 255 231]./255;
col_Easterly      = [102 204 153]./255;
col_Overturning_elsewhere = [200 170 255]./255;
col_PolarWesterly = [255 188 229]./255;
col_SH_Westerly   = [248 173 158]./255;
col_Elsewhere     = [.7 .7 .7];
col_Unattributed  = [ 1  1  1];

colmap_domain_wise_exact = [col_local;
                            col_nonlocalTP;
                            col_Westerly;
                            col_PolarWesterly;
                            col_IM;
                            col_EAM;
                            col_Overturning_elsewhere;
                            col_Easterly;
                            col_SH_Westerly;
                            col_Elsewhere;
                            col_Unattributed];
license('inuse');


%% [Step 5.4] SR network (Pie Chart / Bar plots / Line plots / Table) (Figs. 4, 5)
close all; 
clc; plot_SR_pie = 0; plot_SR_chord = 0; plot_SR_bar = 0; plot_lines = 0; plot_legend_only = 0; plot_table = 0; SR = []; anom_base_period = []; slct_xdata = [];
select_field = 'daily'; x_intvl = []; model_list = []; exp_list = []; ens_list = []; DivByAreaOf = [];  %/default
map_xticks = []; map_xticklabels = []; show_y1_as_bar   = 0;  N_running_avg = []; SR = []; yshift = []; line_data_xlim = []; show_diagonal = 0;
st_month = []; st_day = []; ed_month = []; ed_day = []; bar_pve_col = []; bar_nve_col = []; titlefontsize = []; shift_title_y = []; regress_linewidth = []; all_basins_in_oneplot = 0; x_intvl = []; map_xticks = [];

savefig               = 1;
savemat               = 1;
savenc                = 0;     %/ whether to save SR time series into nc files
recompute_SR          = 0;     %<-- mind this!    
thres_agglomerate     = 3;     %/ group those contributions < x% into 'Others'
trend_mode            = 1;
trend_alpha           = 0.1; 
trend_test            = 'MK';   %/ [], 'MK'
asterisk              = {'(*)', '(**)', '(***)'}; 
alpha_level           = [  0.1,  0.05,   0.01]; 
show_regress_line     = 1;
trend_abs_or_rel      = 'abs';    %/ 'abs' or 'rel'
box_mode              = 'on';
plot_or_not           = 1; 
line_yerror_type      = 'percentile';   %/ 'sd', 'percentile'
%==========================================================================
Njob = 4; job_id = 4;   %/ r147, r146, r148, r149
%==========================================================================

% plot_SR_bar = 1;   %/ bar plots in Fig. S1
% plot_SR_pie = 1; 
% plot_lines = 1; plot_legend_only = 0; show_legend = 0; lgd_position = 'eastoutside'; %'eastoutside'; southoutside
% plot_lines = 1; plot_legend_only = 1; show_legend = 1; lgd_position = 'eastoutside'; %'eastoutside'; southoutside
plot_table = 1;  %/ get the numbers for Fig. 1 and Fig. S1

if plot_table
    %======================================================================
    line_setname = 'VarSet_Pm_frac_real'; %<- mind this!
    %======================================================================
    if isequal(line_setname, 'VarSet_Pm_frac_real')  %<- Get the estimates for Fig. 1 (1x1 TP grids)
        slct_data_list        = [ repmat({'Pm_frac_real'}, 1, 12)];   %/ compute total water budget
        unit_conv_list        = ones(1, length(slct_data_list));
        area_mean_or_sum_list = [repmat({'sum'}, 1, length(slct_data_list))];
        ins_or_acc_list       = [repmat({'acc'}, 1, length(slct_data_list))];
        grouping_mode_list    = [repmat({'domain-wise-exact'}, 1, 10), {'land-ocean'}, {'land-ocean'}];
        slct_src_list         = {'local',  'nonlocal TP',  'NH_MidLat_Westerly_noTP', 'NH_Polar_Westerly_noTP',  'IM_noTP', 'EAM_noTP', 'Overturning_elsewhere_noTP', 'Easterly_noTP', 'SH_Westerly', 'Elsewhere_noTP', 'land', 'ocean'};
        DivByAreaOf           = 'basin_itself';     
        year_range_bc         = repmat([year_list(1), year_list(end)], length(slct_data_list), 1);
        shared_year_list      = year_list;
        select_field_list     = repmat({'daily'}, 1, length(slct_data_list));
        anom_mode             = 0;   %/ whether to show absolute value or anomaly (i.e. abs - clim) 
        anom_base_period      = []; 
        % mth_list              = [13:16];      %/ Fig. S1
        mth_list              = [14, 20];  %/ Fig. 1
        basin_list            = job_distributor('array', 1:length(basin_catalog), 'Njob', Njob, 'job_id', job_id);
        % basin_list            = 1:length(basin_catalog);
    end
end

%----------------------------------------
if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
    shared_dates = date_array_gen('year_list', shared_year_list, 'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day);
end

%/ Double check if the input WSV var is valid
valid_WSV_dataname = {'Pm', 'Pm_BL', 'P_LA', 'contr_map', 'uptake', 'BL_uptake',...
                      'optimal_trajtime', 'CWRT', 'RH2', 'P_LA',  'rr_tot_L', 'rr_tot_NLL', 'rr_tot_NLO'};

% dataset = []; %/ clear out data to avoid bugs
dataset.placeholder = []; WSV = []; 

if any(ismember(grouping_mode_list, {'domain-wise-exact'}))
    str_regime_ver = strcat('_regime', regime_ver);
else
    str_regime_ver = '';
end

%/ Annual/Seasonal/Monthly Loop
for mth = mth_list
    if mth == 16 && length(shared_year_list) == 1   continue;  end   %/ skip DJF if only one year
    slct_data_old  = nan; area_mean_or_sum_old = nan; ins_or_acc_old = nan; grouping_mode_old = nan; %/ re-initialize them outside the variable loop
    
%     close all;
    %/ Variable Loop
    unit_list = cell(length(slct_data_list),1);
    for j = 1:length(slct_data_list)
        %/ Load queries
        slct_data           = slct_data_list{j};
        if ~isempty(model_list)  model = model_list{j};  else  model = '';  end
        if ~isempty(exp_list)    exp   = exp_list{j};    else  exp   = '';  end
        if ~isempty(ens_list)    ens   = ens_list{j};    else  ens   = '';  end
        area_mean_or_sum    = area_mean_or_sum_list{j};
        ins_or_acc          = ins_or_acc_list{j};
        grouping_mode       = grouping_mode_list{j};
        select_field        = select_field_list{j};
        if ~isempty(grouping_mode)
            str_grouping_mode = strcat('_', strrep(grouping_mode, '-', '_'));
        else
            str_grouping_mode = '';
        end
        
        %/ Dates
        str_years_meteo = sprintf('_%d-%d', year_range_bc(j,1), year_range_bc(j,end));
        if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
            str_timescale = '_interannual'; %/ To imply the mean/acc of the same specific period for different year.
            str_dates = sprintf('%d%02d%02d-%d%02d%02d', year_range_bc(1), st_month, st_day, year_range_bc(end), ed_month, ed_day);
        else
            if mth == 0        
                str_timescale = '_annual';    
                str_dates     = str_years_meteo; %/ [CAVEAT]: Do NOT modify it, otherwise will affect the reading of the processed WSV data.
            else
                str_timescale = strcat('_', str_mth{mth});
                str_dates     = strcat(str_years_meteo, str_timescale); 
            end
        end
        
        if anom_mode 
            if isempty(anom_base_period)
                str_anom_base = sprintf('_anom%d_%d', shared_year_list(1), shared_year_list(end));
            else
                str_anom_base = sprintf('_anom%d_%d', anom_base_period(1), anom_base_period(end));
            end
        else
            str_anom_base = '';
        end
        if trend_mode
            if ~isempty(trend_test)
                str_trend_test = strcat('_', trend_test);
            else
                str_trend_test = '';
            end
        else
            str_trend_test = '';
        end

        %/ Basin Loop
        cnt = 0; 
        for top = basin_list
            cnt  = cnt + 1;
            basin_name = ''; basin_bndry = []; str_basin = [];
            if top == -1      %/ global 
                basin_name  = 'global';
                basin_bndry = [];
            elseif top == -2  %/ global land 
                basin_name  = 'land';
                basin_bndry = [];
            else
                if from_basin
                    if from_basin == 4   
                        if top == 0             %/ then show all the basins 
                            basin_name = 'TP';  %/ Aggregated P_LA/Pm/Pm_BL from sub-basins were processed at [Step 5.1]
                            basin_bndry = TP_bndry;
                        elseif top == 16     %/ S. Inner TP + N. Inner TP
                            ind = [13, 14];
                            basin_name = 'Inner_TP';
                            bndry_data = cat(1,basin_catalog(ind).bndry);
                        elseif top == 17     %/ TP-Indus + TP-Ganges
                            ind = [1, 2];
                            basin_name = 'TP-Indus_TP-Ganges';
                            bndry_data = cat(1,basin_catalog(ind).bndry);
                        else
                            basin_name  = basin_catalog(top).name{:};
                            basin_bndry = basin_catalog(top).bndry;
                        end
                    else
                        basin_name  = basin_catalog(top).name{:};
                        basin_bndry = basin_catalog(top).bndry;
                    end
                    str_basin = strcat('_', basin_name);
                end
            end
            basin_name_strrep = strrep(strrep(strrep(basin_name, '-', '_'), ' ', '_'), '.', '_');

            if isequal(grouping_mode, 'domain-wise-exact')
                grouping_mode_labels = labels_domain_wise_exact;
            else
                grouping_mode_labels = [];  %/ no need to specify labels for other group modes; will be specified in get_SR()
            end
            
            year_list_bc    = year_range_bc(j,1):year_range_bc(j,end);
            savemat_prefix  = strcat('scheme', num2str(RHc_dqc_scheme), str_optimal);
            
            %/ Further post-processing 
            if isequal(slct_data, 'EP_ratio')
                raw_slct_data_list        = {'CERA_E', 'CERA_P'};           %/ Evaporation must be at the first index 
                raw_area_mean_or_sum_list = {'mean', 'mean'};
                raw_ins_or_acc_list       = {'acc',   'acc'};
            elseif isequal(slct_data, 'ERA5_EP_ratio')
                raw_slct_data_list        = {'ERA5_E', 'ERA5_P'}; %/ Evaporation must be at the first index  
                raw_area_mean_or_sum_list = {'mean', 'mean'};
                raw_ins_or_acc_list       = {'acc',   'acc'};
            elseif contains(slct_data, 'evspsbl_pr_ratio')
                raw_slct_data_list        = {'evspsbl', 'pr'};    %/ Evaporation must be at the first index  
                raw_area_mean_or_sum_list = {'mean', 'mean'};
                raw_ins_or_acc_list       = {'acc',   'acc'};
            else
                raw_slct_data_list        = {slct_data};
                raw_area_mean_or_sum_list = {area_mean_or_sum};
                raw_ins_or_acc_list       = {ins_or_acc};
            end

            %/ Due to my dumb coding, no need to recompute SR each time. Once suffices.
            if ~(isequal(slct_data, slct_data_old) && isequal(model, model_old) && isequal(exp, exp_old) && isequal(ens, ens_old) && ...
                 isequal(grouping_mode, grouping_mode_old) && isequal(area_mean_or_sum, area_mean_or_sum_old) && isequal(ins_or_acc, ins_or_acc_old))
                for jj = 1:length(raw_slct_data_list)
                    raw_slct_data        = raw_slct_data_list{jj};
                    raw_area_mean_or_sum = raw_area_mean_or_sum_list{jj};
                    raw_ins_or_acc       = raw_ins_or_acc_list{jj};

                    if isequal(slct_xdata, slct_data)
                        raw_anom_base_period = [];
                    else
                        raw_anom_base_period = anom_base_period;
                    end
                    %/ [IMPORTANT] Get SR Matrix
                    [SR, SR_filename_suffix, dataset] = get_SR('SR', SR, 'project_name', project_name, 'param', param, 'dataset', dataset, 'dataname', dataname, 'data_folder', data_folder, 'valid_WSV_dataname', valid_WSV_dataname,...
                                            'optimal_rr', optimal_rr, 'expmnt', expmnt, 'ldirect', ldirect, 'output_res', output_res, 'dt', dt, 'RHc_dqc_scheme', RHc_dqc_scheme, 'from_basin', from_basin, 'masterfolder', masterfolder,...
                                            'slct_data', raw_slct_data, 'model', model, 'exp', exp, 'ens', ens, 'area_mean_or_sum', raw_area_mean_or_sum, 'ins_or_acc', raw_ins_or_acc, 'grouping_mode', grouping_mode, 'grouping_mode_labels', grouping_mode_labels,...
                                            'regime_ver', regime_ver, 'basin_name', basin_name, 'basin_catalog', basin_catalog, 'thres_agglomerate', thres_agglomerate,...
                                            'select_field', select_field, 'DivByAreaOf', DivByAreaOf, 'mth', mth, 'str_mth', str_mth,...
                                            'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day, 'year_list', year_list_bc, 'shared_year_list', shared_year_list, 'anom_base_period', raw_anom_base_period,...
                                            'recompute_SR', recompute_SR, 'savemat', savemat, 'savemat_prefix', savemat_prefix, 'savenc', savenc, 'NumWorkers', NumWorkers);

                end

                %/ Further post-processing of the EP ratio
                if contains(slct_data, 'EP_ratio') || contains(slct_data, 'evspsbl_pr_ratio')
                    if contains(slct_data, 'evspsbl_pr_ratio')
                        str_cmip6 = strcat('_',model,'_',exp);
                    else
                        str_cmip6 = '';
                    end

                    SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale)) = SR.(basin_name_strrep).(strcat(raw_slct_data_list{1}, str_cmip6, str_timescale))./SR.(basin_name_strrep).(strcat(raw_slct_data_list{2}, str_cmip6, str_timescale))*100;
                    
                    %/ Compute anom of EP ratio
                    if ~isempty(anom_base_period)
                        str_cmip6_hist = strcat('_',model,'_historical');
                        hist_period = SR.(basin_name_strrep).([raw_slct_data_list{1}, str_cmip6_hist, '_date', str_timescale]);
                        ind_dates = findismember_loop(hist_period, anom_base_period);
                        if length(ind_dates) ~= length(anom_base_period)
                            error('Data period do not cover the ''anom_base_period''!');
                        end
                        
                        hist_data = SR.(basin_name_strrep).(strcat(raw_slct_data_list{1}, str_cmip6_hist, str_timescale))./SR.(basin_name_strrep).(strcat(raw_slct_data_list{2}, str_cmip6_hist, str_timescale))*100;
                        mean_hist_data = mean(hist_data(ind_dates), 'omitnan');
                        SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale, str_anom_base)) = SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale)) - mean_hist_data;
                    end
                    
                    SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale, '_clim')) = mean(SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale)), 'omitnan');
                    SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale, '_sd'))   = std(SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale)), 0, 2, 'omitnan'); 
                    SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, '_date', str_timescale)) = SR.(basin_name_strrep).(strcat(raw_slct_data_list{1}, str_cmip6, '_date', str_timescale)); 
                    SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale, '_unit')) = '%';
                end

                %/ Loop over individual MME models
                if contains(line_setname, 'CMIP6_MME')
                    slct_data_new      = strrep(strcat(raw_slct_data_list{1},'_',model,'_',exp), '-', '_');
                    year_list_bc       = year_range_bc(j,1):year_range_bc(j,end);
                    MME_model_list     = SR.(basin_name_strrep).(strcat(slct_data_new, '_model_list')); 
                    MME_model_ens_list = SR.(basin_name_strrep).(strcat(slct_data_new, '_ens_list')); 

                    %/ Initialization
                    if j == 1 && cnt == 1
                        MME_model_data = nan(length(shared_year_list), length(MME_model_list), length(slct_data_list), length(basin_list)); 
                        MME_model_sd   = nan(length(shared_year_list), length(slct_data_list), length(basin_list)); 
                        MME_model_25th = nan(length(shared_year_list), length(slct_data_list), length(basin_list)); 
                        MME_model_75th = nan(length(shared_year_list), length(slct_data_list), length(basin_list)); 
                    end

                    for mm = 1:length(MME_model_list)
                        each_model     = MME_model_list{mm};
                        each_model_ens = MME_model_ens_list{mm};

                        slct_data_new_individual = cell(length(raw_slct_data_list), 1);
                        for jj = 1:length(raw_slct_data_list)
                            raw_slct_data        = raw_slct_data_list{jj};
                            raw_area_mean_or_sum = raw_area_mean_or_sum_list{jj};
                            raw_ins_or_acc       = raw_ins_or_acc_list{jj};
                            slct_data_new_individual{jj} = strrep(strcat(raw_slct_data,'_',each_model,'_',exp), '-', '_');

                            %/ [IMPORTANT] Get SR Matrix
                            [SR, SR_filename_suffix, dataset] = get_SR('SR', SR, 'project_name', project_name, 'param', param, 'dataset', dataset, 'dataname', dataname, 'data_folder', data_folder, 'valid_WSV_dataname', valid_WSV_dataname,...
                                                    'optimal_rr', optimal_rr, 'expmnt', expmnt, 'output_res', output_res, 'dt', dt, 'RHc_dqc_scheme', RHc_dqc_scheme, 'from_basin', from_basin, 'masterfolder', masterfolder,...
                                                    'slct_data', raw_slct_data, 'model', each_model, 'exp', exp, 'ens', each_model_ens, 'area_mean_or_sum', raw_area_mean_or_sum, 'ins_or_acc', raw_ins_or_acc, 'grouping_mode', grouping_mode, 'grouping_mode_labels', grouping_mode_labels,...
                                                    'regime_ver', regime_ver, 'basin_name', basin_name, 'basin_catalog', basin_catalog, 'thres_agglomerate', thres_agglomerate,...
                                                    'select_field', select_field, 'DivByAreaOf', DivByAreaOf, 'mth', mth, 'str_mth', str_mth,...
                                                    'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day, 'year_list', year_list_bc, 'shared_year_list', shared_year_list, 'anom_base_period', anom_base_period,...
                                                    'recompute_SR', recompute_SR, 'savemat', savemat, 'savemat_prefix', savemat_prefix, 'savenc', savenc, 'NumWorkers', NumWorkers);
                        end

                        %/ Subset dates
                        ind_date = findismember_loop(shared_year_list, SR.(basin_name_strrep).([slct_data_new_individual{1}, str_grouping_mode, '_date', str_timescale]));

                        %/ Further post-processing 
                        if ~isempty(anom_base_period)  
                            if contains(slct_data, 'evspsbl_pr_ratio')
                                str_cmip6_hist = strcat('_',strrep(each_model, '-', '_'),'_historical');
                                str_cmip6      = strcat('_',strrep(strcat(each_model,'_',exp), '-', '_'));     
                                hist_period    = SR.(basin_name_strrep).([raw_slct_data_list{1}, str_cmip6_hist, '_date', str_timescale]);
                                ind_dates      = findismember_loop(hist_period, anom_base_period);
                                if length(ind_dates) ~= length(anom_base_period)
                                    error('Data period do not cover the ''anom_base_period''!');
                                end
                                
                                %/ Model time series in historical run
                                hist_data = SR.(basin_name_strrep).(strcat(raw_slct_data_list{1}, str_cmip6_hist, str_timescale))./SR.(basin_name_strrep).(strcat(raw_slct_data_list{2}, str_cmip6_hist, str_timescale))*100;
                                mean_hist_data = mean(hist_data(ind_dates), 'omitnan');

                                %/ Model time series in the queried experiment (historical/ssp245/ssp585) run
                                exp_data = SR.(basin_name_strrep).(strcat(raw_slct_data_list{1}, str_cmip6, str_timescale))./SR.(basin_name_strrep).(strcat(raw_slct_data_list{2}, str_cmip6, str_timescale))*100;
                                MME_model_data(ind_date,mm,j,cnt) = exp_data - mean_hist_data;
                            else
                                MME_model_data(ind_date,mm,j,cnt) = SR.(basin_name_strrep).(strcat(slct_data_new_individual{1}, str_grouping_mode, str_timescale, str_anom_base));
                            end
                        else
                            if contains(slct_data, 'evspsbl_pr_ratio')
                                MME_model_data(ind_date,mm,j,cnt) = SR.(basin_name_strrep).(strcat(slct_data_new_individual{1}, str_grouping_mode, str_timescale))./SR.(basin_name_strrep).(strcat(slct_data_new_individual{2}, str_grouping_mode, str_timescale))*100;
                            else
                                MME_model_data(ind_date,mm,j,cnt) = SR.(basin_name_strrep).(strcat(slct_data_new_individual{1}, str_grouping_mode, str_timescale));
                            end
                        end
                    end
                    
                    %/ Inter-model error
                    MME_model_sd(:,j,cnt)   = std(squeeze(MME_model_data(:,:,j,cnt)), [], 2);  
                    MME_model_75th(:,j,cnt) = prctile(squeeze(MME_model_data(:,:,j,cnt)), 75, 2); 
                    MME_model_25th(:,j,cnt) = prctile(squeeze(MME_model_data(:,:,j,cnt)), 25, 2); %/ along row (different from std)
                end
            end
        end
        slct_data_old           = slct_data;
        model_old               = model;
        exp_old                 = exp;
        ens_old                 = ens;
        grouping_mode_old       = grouping_mode;
        area_mean_or_sum_old    = area_mean_or_sum;
        ins_or_acc_old          = ins_or_acc;
        
        %/ Plot Pie/Chord Chart (after all basins loaded to get consistent colorings!)
        if (plot_SR_chord || plot_SR_pie || plot_SR_bar) && plot_or_not
            % close all;
            %/ Get the unique sources in the SR network of all basins
            sig_src_cell = [];
            all_basin_names = fieldnames(SR);
            for i = 1:length(all_basin_names)
                sig_src_cell = [sig_src_cell; SR.(all_basin_names{i}).([slct_data, str_grouping_mode, '_src_name'])];
            end
            sig_src = unique(sig_src_cell, 'stable');

            if isequal(grouping_mode, 'basin-wise')
                sort_mode       = [];   %/ no auto-sorting but manual sorting and re-ordering
                
                %/ Re-ordering of the sources
                ind_TPbasins    = findismember_loop(sig_src, [basin_catalog.name]);
                
                %/ Make sure a consistent ordering here and from reg_extractor 
                labels_ML_basins = {'ML-Indus', 'ML-Ganges', 'ML-Irrawaddy', 'ML-Mekong',...
                                    'ML-Yangtze', 'ML-Yellow', 'ML-Gobi1', 'ML-Tarim', 'ML-Syr Darya', 'ML-Amu Darya',...
                                    'ML-CapsianSea East Coast', 'ML-Helmand'};
%                 ind_MLTPbasins  = find(contains(sig_src, 'ML-'));
                ind_MLTPbasins  = findismember_loop(sig_src, labels_ML_basins);
                ind_Elsewhere   = findismember_loop(sig_src, {'Elsewhere'});
                ind_TP          = findismember_loop(sig_src, {'TP'});
                ind_remaining   = setdiff(1:length(sig_src), [ind_TPbasins', ind_MLTPbasins', ind_Elsewhere, ind_TP])';

                new_ordering    = [ind_TPbasins; ind_MLTPbasins; ind_remaining; ind_Elsewhere; ind_TP];
                sig_src_ordered = sig_src(new_ordering);

                col_TPbasins      = brewermap(length(ind_TPbasins  )+1, 'Blues');  col_TPbasins(1,:)   = []; %remove the lightest color.
                col_MLTPbasins    = brewermap(length(ind_MLTPbasins)+1, 'Greens'); col_MLTPbasins(1,:) = []; %remove the lightest color.
                col_remaining     = brewermap(length(ind_remaining )+1, 'Reds');   col_remaining(1,:)  = []; %remove the lightest color.
                if ~isempty(ind_TP)  col_TP_bc = [0 0 0];  else  col_TP_bc = [];  end
                sig_src_colorings = [col_TPbasins; col_MLTPbasins; col_remaining; col_Elsewhere; col_TP_bc];   
                
            elseif contains(grouping_mode, 'domain-wise')  %/ domain-wise or domain-wise-exact
                sort_mode   = [];   %/ no auto-sorting but manual sorting and re-ordering
                
                new_ordering = findismember_loop(sig_src, {'local', 'nonlocal TP', 'NH_MidLat_Westerly_noTP', 'IM_noTP', 'EAM_noTP',...
                                                           'Easterly_noTP', 'Overturning_elsewhere_noTP', 'NH_Polar_Westerly_noTP',...
                                                           'SH_Westerly', 'Elsewhere_noTP', 'Unattributed'});  
                sig_src_ordered = sig_src(new_ordering);

                if isequal(grouping_mode, 'domain-wise')  
                    sig_src_colorings = [ col_local;
                                          col_nonlocalTP;
                                          col_IM;
                                          col_Westerly;
                                          col_EAM;
                                          col_Elsewhere;];

                elseif isequal(grouping_mode, 'domain-wise-exact')
                    sig_src_colorings = [ col_local;
                                          col_nonlocalTP;
                                          col_Westerly;
                                          col_IM;
                                          col_EAM;
                                          col_Easterly;
                                          col_Overturning_elsewhere;
                                          col_PolarWesterly;
                                          col_SH_Westerly;
                                          col_Elsewhere;
                                          col_Unattributed];
                else
                    error('Code not set!')
                end
                                 
            elseif isequal(grouping_mode, 'land-ocean')
                sort_mode             = [];
                ind_ordered           = findismember_loop(sig_src, {'land', 'ocean'});
                sig_src_ordered       = sig_src(ind_ordered);

                sig_src_colorings = [ 143 227 207;
                                       43  52 103;]./255;
                                   
            elseif isequal(grouping_mode, 'local-nonlocal')
                sort_mode             = [];       %/ no sorting
                sig_src_ordered       = sig_src;  %/ no reordering is needed
                sig_src_colorings = [ 143 227 207;
                                       43  52 103;]./255;
            else
                error('code not set!')
            end
            sig_src_ordered_label = strcat(pad(num2roman(1:length(sig_src_ordered)), 'right')',{': '}, sig_src_ordered);
%             disp(sig_src_ordered_label)
            
            %/ Save the ordered unique sources and the colorings for plotting
            %/  the boundaries on map [Step 3.2]
            if from_basin == 4 && isequal(basin_list, 0:length(basin_catalog))
                SR_sig_src_filename = strcat(data_folder, 'SR_sig_src_', SR_filename_suffix, '_AllTPbasins', '.mat');
                if savemat
                    fprintf('*** Saving sig source data: %s *** \n', SR_sig_src_filename)
                    save(SR_sig_src_filename, 'sig_src_ordered', 'sig_src_colorings', '-v7.3');
                end
            end
            cnt = 0;
            for top = basin_list
                cnt = cnt + 1;
                basin_name = ''; basin_bndry = [];
                if from_basin
                    if from_basin == 4   
                        if top == 0             %/ then show all the basins 
                            basin_name = 'TP';  %/ Aggregated P_LA/Pm/Pm_BL from sub-basins were processed at [Step 5.1]
                            basin_bndry = TP_bndry;
                        elseif top == 16     %/ S. Inner TP + N. Inner TP
                            ind = [13, 14];
                            basin_name = 'Inner_TP';
                            bndry_data = cat(1,basin_catalog(ind).bndry);
                        elseif top == 17     %/ TP-Indus + TP-Ganges
                            ind = [1, 2];
                            basin_name = 'TP-Indus_TP-Ganges';
                            bndry_data = cat(1,basin_catalog(ind).bndry);
                        else
                            basin_name  = basin_catalog(top).name{:};
                            basin_bndry = basin_catalog(top).bndry;
                        end
                    else
                        basin_name  = basin_catalog(top).name{:};
                        basin_bndry = basin_catalog(top).bndry;
                    end
                    str_basin = strcat('_', basin_name);
                end
                basin_name_strrep = strrep(strrep(strrep(basin_name, '-', '_'), ' ', '_'), '.', '_');

                dataname4plot   = SR.(basin_name_strrep).([slct_data, str_grouping_mode, '_src_name']);
                data_ts         = SR.(basin_name_strrep).([slct_data, str_grouping_mode, str_timescale]);
                data_clim       = SR.(basin_name_strrep).([slct_data, str_grouping_mode, str_timescale, '_clim']); 
                errdata_clim    = SR.(basin_name_strrep).([slct_data, str_grouping_mode, str_timescale, '_sd']);

                %/ Rmb to update the ordering!
                ind             = findismember_loop(dataname4plot, sig_src_ordered);
                ind_local       = findismember_loop(dataname4plot, {basin_name, 'local'});
                data_ts         = data_ts(ind,:);         
                data_clim       = data_clim(ind);       
                errdata_clim    = errdata_clim(ind);   
                
                rowName         = sig_src_ordered_label;
                colName         = {basin_name};
                colmap          = sig_src_colorings;

                if plot_SR_chord        %/ Not recommended
                    chord_data = data_clim;
                    figure
                    set(gcf, 'color','w');
                    set(gcf,'position', [700 100 800 800]) % sets figure size

                    rowName = cellfun(@(x) char(strrep(x,{'-'}, '')), rowName, 'UniformOutput', false);
                    rowName = cellfun(@(x) char(strrep(x,{' '}, '')), rowName, 'UniformOutput', false);
                    rowName = cellfun(@(x) char(strrep(x,{'.'}, '')), rowName, 'UniformOutput', false);

                    colName = cellfun(@(x) char(strrep(x,{'-'}, '')), colName, 'UniformOutput', false);
                    colName = cellfun(@(x) char(strrep(x,{' '}, '')), colName, 'UniformOutput', false);
                    colName = cellfun(@(x) char(strrep(x,{'.'}, '')), colName, 'UniformOutput', false);

                    CC=chordChart(chord_data,'rowName',rowName,'colName',colName);
                    CC=CC.draw();
%                     disp(rowName_real)
                end     
                if plot_SR_pie          %/ Preferred (climatological mean)
                    fontsize   = 14;
                    linewidth  = 2;
                    show_lgd   = 0;  %<- mind this!
                    
                    %/ Put 'Unattributed' at the end, while the others are sorted in descending order
                    if isequal(grouping_mode, 'domain-wise-exact')
                        sort_mode       = [];
                        ind_unattr      = find(contains(rowName, {'Unattributed'}));
                        ind_others      = setdiff(1:length(rowName), ind_unattr);
                        [~, I]          = sort(data_clim(ind_others), 'descend');
                        ind_ordered     = [I; ind_unattr]; 
                        T_rowname_bc    = rowName(ind_ordered); 
                        pie_data        = data_clim(ind_ordered);
                        pie_colmap      = colmap(ind_ordered,:);
                    else
                        sort_mode       = 'descend';   
                        ind_ordered     = ind;
                        T_rowname_bc    = rowName;
                        pie_data        = data_clim;
                        pie_colmap      = colmap;
                    end
                    
                    T_varname_bc = colName;
                    str_numbering = num2roman(ind_ordered);
%                     str_numbering = string(ind_ordered);
%                     str_numbering(ind_local) = 'Local'; %/ sometimes a long string causes overlapping
                    pie_label = strcat(str_numbering, {': '}, string(round(pie_data, 1)), '%');
    %                 pie_label = strcat(string(round(pie_data, 0)), '%');

                    %/ Set the labels within large pies, otherwise outside the pie.
                    if isequal(line_setname, 'VarSet_domain-wise-exact')
                        thres_for_noshow = 7;      %/ no need to show all, as we want to see the dominant regime to the subbasin only.
%                         thres_for_noshow = 10;      %/ no need to show all, as we want to see the dominant regime to the subbasin only.
                    else
                        thres_for_noshow = thres_agglomerate; %/ [Viz] do not label those < xx%   
                    end
                    
%                     label_pos = rescale(pie_data, 0.5, 0.7); %/ no inverse
                    
                    label_pos = repmat(1.1, length(pie_data), 1);    %/ outside (default)
                    label_pos(pie_data >= thres_for_noshow & pie_data < 20) = 0.6; %/ inside
                    label_pos(pie_data >= 20) = 0.5;                 %/ more inside

                    %/ Use black (white) label (when it is inside the pie) for light (dark) pie.
                    label_col    = repmat([0 0 0], length(pie_data), 1);
                    ind_usewhite = find(mean(pie_colmap,2)<0.35 & label_pos < 1);
                    label_col(ind_usewhite,:) = repmat([1 1 1], length(ind_usewhite), 1);  

                    %/ Set label fontsize
                    label_fontsize = repmat(fontsize, length(pie_data), 1);
                    label_fontsize(pie_data >= thres_for_noshow & pie_data < 20) = fontsize*1.2;
                    label_fontsize(pie_data >= 20) = fontsize*1.4;

                    pie_legend = strcat({' '}, strrep(T_rowname_bc, '_', ' '));  %/ add a small placeholder
                    
                    pie_titlename = strrep(strcat(basin_name, str_timescale), '_', ' ');
                    if savefig
                        pie_figname = strrep(strcat('pie_', slct_data, str_grouping_mode, str_regime_ver, str_dates, '_', basin_name), ' ', '_');
                        disp(pie_figname)
                    else
                        pie_figname = [];
                    end

                    quickplot_pie('pie_data', pie_data, 'pie_label', pie_label, 'pie_legend', pie_legend, 'pie_colmap', pie_colmap,...
                              'thres_for_noshow', thres_for_noshow, 'sort_mode', sort_mode, 'label_pos', label_pos, 'label_col', label_col, 'label_fontsize', label_fontsize,...
                              'fontsize', fontsize, 'linewidth', linewidth, 'pie_titlename', pie_titlename,...
                              'pie_figname', pie_figname, 'plotting_folder', plotting_folder, 'show_lgd',  show_lgd, 'plot_lgd_only', 0);

                    %/ Plot legend
                    if cnt == length(basin_list) && mth == mth_list(end)
                        %/ make up a fake data to show all sources in a legend
                        pie_data_demo   = [1:length(sig_src_ordered)]/sum(1:length(sig_src_ordered))*100;
                        pie_label_demo  = strcat(string(round(pie_data_demo, 1)), '%');
                        pie_legend_demo = strcat({' '}, strrep(sig_src_ordered_label, '_', ' '));  %/ add a small placeholder
                        pie_colmap_demo = sig_src_colorings;
                        Orientation     = 'vertical';
                        if isequal(grouping_mode, 'basin-wise') 
                            NumColumns = 7;
                        else
                            NumColumns = 1;
                        end
                        fontsize_demo   = fontsize*0.75;
                        linewidth_demo  = linewidth*0.75;

                        quickplot_pie('pie_data', pie_data_demo, 'pie_label', pie_label_demo, 'pie_legend', pie_legend_demo, 'pie_colmap', pie_colmap_demo,...
                                      'thres_for_noshow', thres_for_noshow, 'sort_mode', sort_mode, 'label_pos', [], 'fontsize', fontsize_demo, 'linewidth', linewidth_demo,...
                                      'pie_titlename', pie_titlename, 'pie_figname', pie_figname, 'plotting_folder', plotting_folder,...
                                      'plot_lgd_only', 1, 'Orientation', Orientation, 'NumColumns', NumColumns);
                    end
                end
                if plot_SR_bar          
                    if isequal(bar_mode, 'grouped')      %/ Climatological mean
                        bar_data         = data_clim;     %/ group x bars
                        bar_error_upper  = errdata_clim;
                        bar_error_lower  = errdata_clim;
                        sort_mode        = 'descend';
                        bar_labels       = num2roman(1:length(rowName))'; 
                        bar_labels_intvl = 1;
                        stack_labels     = rowName;
                        bar_colmap       = sig_src_colorings;
                        bar_range        = [0, 59];
                        bar_rangeticks   = [];         %/ -100:10:100;
                        box_on           = 0;
                        %/ remove 'Unattributed'
                        if isequal(grouping_mode, 'domain-wise-exact')
                            ind_unattr = find(contains(stack_labels, {'Unattributed'}));
                            bar_data(ind_unattr,:)     = [];
                            bar_labels(ind_unattr,:)   = [];
                            stack_labels(ind_unattr,:) = [];
                            bar_colmap(ind_unattr,:)   = [];
                            % bar_data(ind_unattr,:)     = [];
                            % bar_labels(ind_unattr,:)   = [];
                            % stack_labels(ind_unattr,:) = [];
                            % bar_colmap(ind_unattr,:)   = [];
                        end

                    elseif isequal(bar_mode, 'stacked')   %/ Stacked bar of all source contr. over 40 years
                        bar_data         = data_ts';       %/ bars x stacks
                        bar_error_upper  = [];
                        bar_error_lower  = [];
                        sort_mode        = []; %'descend';
                        % bar_labels      = num2roman(1:length(rowName));
                        bar_labels       = year_list_bc;
                        bar_labels_intvl = 10;
                        stack_labels     = rowName';
                        bar_colmap       = sig_src_colorings;
                        bar_range        = [0, 120];
                        bar_rangeticks   = bar_range(1):20:bar_range(end);        %/ -100:10:100;
                        box_on           = 1;
                        %/ remove 'Unattributed'
                        if isequal(grouping_mode, 'domain-wise-exact')
                            ind_unattr = find(contains(stack_labels, {'Unattributed'}));
                            bar_data(:,ind_unattr)     = [];
                            stack_labels(:,ind_unattr) = [];
                            bar_colmap(ind_unattr,:)   = [];
                            % bar_data(ind_unattr,:)     = [];
                            % stack_labels(ind_unattr,:) = [];
                            % bar_colmap(ind_unattr,:)   = [];
                        end
                    else
                        error('''bar_mode'' can only take ''stacked'' or ''grouped''!');
                    end

                    show_value = 1; 
                    edgecolor = 'k';
                    y_labels    = [];
                    orientation = [];
                    panel_pos = [800 200 800 475];
                    fontsize  = 14.5;
                    linewidth = 2;
                    titlename = strrep(strcat(basin_name, str_timescale), '_', ' ');
                    
                    if savefig
                        FigName_underscore = strrep(strcat('bar_', slct_data, str_grouping_mode, str_regime_ver, str_dates, '_', basin_name), ' ', '_');
                        disp(FigName_underscore)
                    else
                        FigName_underscore = [];
                    end
                    plot_stackedbar('bar_mode', bar_mode, 'bar_data', bar_data, 'bar_error_upper', bar_error_upper, 'bar_error_lower', bar_error_lower, 'bar_range', bar_range, 'bar_rangeticks', bar_rangeticks,...
                                    'sort_mode', sort_mode, 'show_value', show_value, 'bar_labels', bar_labels, 'bar_labels_intvl', bar_labels_intvl, 'stack_labels', stack_labels,...
                                    'bar_colmap', bar_colmap, 'edgecolor', edgecolor, 'y_labels', y_labels, 'orientation', orientation,...
                                    'panel_pos', panel_pos, 'fontsize', fontsize, 'linewidth', linewidth, 'box_on', box_on,...
                                    'titlename', titlename, 'FigName_underscore', FigName_underscore, 'plotting_folder', plotting_folder, 'savefig', savefig);
                end
                
                if isequal(line_setname, 'VarSet_land-ocean')
                    fprintf('%s: traceability = %.1f %s %.1f%%\n', strrep(str_timescale, '_', ''), sum(data_clim), char(177), sqrt(sum(errdata_clim.^2)));
                end
            end
        end
    end
    
    %/ After loading all variables, plot lines if queried
    if plot_lines && plot_or_not
        pve_nve_bar_mode = 1; 
        show_y_max       = 0;              %/ whether to pinpoint the local max of the line.
        show_axislabel   = 1;              %/ whether to put labels on each y-axes.
        
        %/ Plot multiple lines for each basin
        cnt = 0; 
        for top = basin_list
            cnt = cnt + 1;
            basin_name = ''; basin_bndry = []; str_basin = [];
            if top == -1  %/ global 
                basin_name = 'global';
                basin_bndry = [];
            elseif top == -2  %/ global land 
                basin_name = 'land';
                basin_bndry = [];
            else
                if from_basin
                    if from_basin == 4   
                        if top == 0             %/ then show all the basins 
                            basin_name = 'TP';  %/ Aggregated P_LA/Pm/Pm_BL from sub-basins were processed at [Step 5.1]
                            basin_bndry = TP_bndry;
                        elseif top == 16     %/ S. Inner TP + N. Inner TP
                            ind = [13, 14];
                            basin_name = 'Inner_TP';
                            bndry_data = cat(1,basin_catalog(ind).bndry);
                        elseif top == 17     %/ TP-Indus + TP-Ganges
                            ind = [1, 2];
                            basin_name = 'TP-Indus_TP-Ganges';
                            bndry_data = cat(1,basin_catalog(ind).bndry);
                        else
                            basin_name  = basin_catalog(top).name{:};
                            basin_bndry = basin_catalog(top).bndry;
                        end
                    else
                        basin_name  = basin_catalog(top).name{:};
                        basin_bndry = basin_catalog(top).bndry;
                    end
                    str_basin = strcat('_', basin_name);
                end
            end
            basin_name_strrep = strrep(strrep(strrep(basin_name, '-', '_'), ' ', '_'), '.', '_');

            %===== line_ydata (ntime, nvar) =====%
            if all_basins_in_oneplot
                nvar = length(slct_data_list)*length(basin_list);
            else
                nvar = length(slct_data_list);
                line_yerror_upper = [];  line_yerror_lower = []; %/ then clear up line_yerror_upper, line_yerror_lower every time
            end

            %/ Initialization
            if cnt == 1
                if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
                    line_ydata        = nan(length(shared_dates), nvar);
                    line_ydata_ref    = nan(length(shared_dates), nvar);
                else
                    line_ydata        = nan(length(shared_year_list), nvar);
                    line_ydata_ref    = nan(length(shared_year_list), nvar);
                end
                line_lgd_str          = cell(nvar, 1);
                line_axislabel        = cell(nvar, 1);
                ind_ydata_firstNonNaN = nan(nvar,  1);
                str_period            = cell(nvar, 1);
            end

            for j = 1:length(slct_data_list)
                slct_data = slct_data_list{j};
                if ~isempty(model_list)  model = model_list{j};  else  model = '';  end
                if ~isempty(exp_list)    exp   = exp_list{j};    else  exp   = '';  end
                if ~isempty(ens_list)    ens   = ens_list{j};    else  ens   = '';  end
                grouping_mode = grouping_mode_list{j};

                %/ Whether to show all basins' time series in one plot
                if all_basins_in_oneplot
                    ind_var = j + (cnt-1)*length(slct_data_list); 
                else
                    ind_var = j;
                end
    
                %/ Set slct_data_callfile, slct_data_new to handle CMIP6 data    
                if ~isempty(model) && ~isempty(exp) && ~isempty(ens)
                    slct_data_new      = strrep(strcat(slct_data,'_',model,'_',exp), '-', '_');
                else
                    slct_data_new      = slct_data;
                end

                if ~isempty(grouping_mode)
                    str_grouping_mode = strcat('_', strrep(grouping_mode, '-', '_'));
                else
                    str_grouping_mode = '';
                end
                if contains(slct_data_new, 'Pm')
                    ind       = findismember_loop(SR.(basin_name_strrep).([slct_data_new, str_grouping_mode, '_src_name']), slct_src_list{j});
                    line_lgd_str{ind_var} = strcat(basin_name, '_', slct_src_list{j});
                else
                    ind       = 1;
                    line_lgd_str{ind_var} = strcat(basin_name, '_', slct_data_new);
                end
                if isequal(basin_name, 'TP') && isequal(slct_src_list{j}, 'nonlocal TP')
                    %/ No need to show 'nonlocal TP' contribution for 'TP'. 
                    continue;
                else
                    if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
                        %/ Plot the daily/monthly time series
                        ind_date = findismember_loop(shared_dates, SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, '_', 'date_yyyymmdd_AllYr')));
                        line_ydata(ind_date,ind_var) = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, '_', select_field))(ind, :);
                        line_axislabel{ind_var} = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, '_', select_field, '_unit'));
                    else
                        %/ Plot the annual/seasonal/monthly time series
                        ind_date = findismember_loop(shared_year_list, SR.(basin_name_strrep).([slct_data_new, str_grouping_mode, '_date', str_timescale]));

                        %/ [IMPORTANT]: if slct_data is not slct_xdata, then we load its anom value.
                        if ~isempty(anom_base_period) && ~isequal(slct_xdata, slct_data)
                            line_ydata(ind_date,ind_var) = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale, str_anom_base))(ind, :);
                        else
                            line_ydata(ind_date,ind_var) = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale))(ind, :);
                        end
                        %/ Compute line_ydata_ref as we need absolute reference to compute relative change.
                        line_ydata_ref(ind_date,ind_var) = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale))(ind, :);  
                        
                        line_axislabel{ind_var} = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale, '_unit'));
                    end
                end
                
                %/ Unit Conversion (if any)
                line_ydata(:,ind_var)     = line_ydata(:,ind_var)*unit_conv_list(j);
                if contains(line_setname, 'CMIP6_MME')
                    MME_model_sd(:,j,cnt)   = MME_model_sd(:,j,cnt)*unit_conv_list(j);
                    MME_model_75th(:,j,cnt) = MME_model_75th(:,j,cnt)*unit_conv_list(j);
                    MME_model_25th(:,j,cnt) = MME_model_25th(:,j,cnt)*unit_conv_list(j);
                end
                ind_ydata_firstNonNaN(ind_var) = find(~isnan(line_ydata(:,ind_var)), 1, 'first');
                str_period{ind_var} = convertStringsToChars( strjoin(string(year_range_bc(j,:)), '-') );
            end
            
            %/ Check whether to plot by each basin or by all
            if (all_basins_in_oneplot && cnt == length(basin_list)) || ~all_basins_in_oneplot 

                %/ Errors / Uncertainties
                if contains(line_setname, 'CMIP6_MME')
                    if isequal(line_yerror_type, 'sd')
                        %/ Reshape to concatenate the basin dimensions with var dimension
                        line_yerror_upper = reshape(MME_model_sd, size(MME_model_sd,1), []);   %/ take standard deviation along the model dimension
                        line_yerror_lower = reshape(MME_model_sd, size(MME_model_sd,1), []);   %/ take standard deviation along the model dimension
                    elseif isequal(line_yerror_type, 'percentile')
                        line_yerror_upper = reshape(MME_model_75th, size(MME_model_75th, 1), []); %/ 75th percentile
                        line_yerror_lower = reshape(MME_model_25th, size(MME_model_25th, 1), []); %/ 25th percentile
                    end
                end
                
                if ~isempty(line_yerror_upper) && ~isempty(line_yerror_lower)
                    lineerror_colmap = line_colmap; 
                    lineerror_alpha  = 0.15;
                    str_yerror       = strcat('_', line_yerror_type);
                else
                    lineerror_colmap = []; 
                    lineerror_alpha  = [];
                    str_yerror       = [];
                end
    
                if ~isempty(anom_base_period) || ~isempty(yshift)
                    show_zero_line = 1;
                else
                    show_zero_line = 0;
                end
    
                %/ Finally, apply moving averaging if required.
                if ~isempty(N_running_avg)
                    line_ydata = my_movmean('data', line_ydata,        'n', N_running_avg, 'dim', 1);
                    if ~isempty(line_yerror_upper) && ~isempty(line_yerror_lower)
                        line_yerror_upper = my_movmean('data', line_yerror_upper, 'n', N_running_avg, 'dim', 1);
                        line_yerror_lower = my_movmean('data', line_yerror_lower, 'n', N_running_avg, 'dim', 1);
                    end
                    str_N_moving_avg = sprintf('_%drunningavg', N_running_avg);
                else
                    str_N_moving_avg = '';
                end
                line_lgd_str = strrep(line_lgd_str, '_', ' ');
                
                %===== line_xdata (ntime, 1) ======%
                if ~isempty(slct_xdata)
                    %/ E.g. we can set xdata as 'T2m_global' to compute scaling
    
                    if all_basins_in_oneplot
                        %/ Then replicate slct_data_list by the # of basins first
                        ind = findismember_loop(repmat(slct_data_list, 1, length(basin_list)), slct_xdata); 
                    else
                        ind = findismember_loop(slct_data_list, slct_xdata); 
                    end
    
                    line_xdata = line_ydata(:,ind(1));        %/ get the first one (the other could be duplicated)
                    line_ydata(:,ind)           = [];         %/ remove the index for xdata from line_ydata!! [IMPORTANT]
                    line_ydata_ref(:,ind)       = [];         %/ remove [IMPORTANT]
                    if ~isempty(line_yerror_upper) && ~isempty(line_yerror_lower)
                        line_yerror_upper(:,ind)    = [];         %/ remove [IMPORTANT]
                        line_yerror_lower(:,ind)    = [];         %/ remove [IMPORTANT]
                    end
                    line_lgd_str(ind)           = [];         %/ remove [IMPORTANT]
                    line_axislabel(ind)         = [];         %/ remove [IMPORTANT]
                    ind_ydata_firstNonNaN(ind)  = [];         %/ remove [IMPORTANT]
                    str_period(ind)             = [];         %/ remove [IMPORTANT]
                    linestyle = 'none';
                    if contains(slct_xdata, {'T2m', 'tas'})
                        xgap_left  = 0.05;            %/ Avoid clipping or a gap at the 1st bar. 
                        xgap_right = 0.05;            %/ Avoid clipping or a gap at the last bar. 
                    else
                        xgap_left  = 0.1;            %/ Avoid clipping or a gap at the 1st bar. 
                        xgap_right = 0.1;            %/ Avoid clipping or a gap at the last bar. 
                    end
                    map_xtickangle = 0;
                else    
                    if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
                        line_xdata       = 1:length(shared_dates);
                        if isempty(x_intvl)    x_intvl    = 2;                             end  %/ interval: 2 days
                        if isempty(map_xticks) map_xticks = 1:x_intvl:length(line_xdata);  end
                        map_xticklabels  = string(shared_dates(map_xticks)); %/ Convert numeric dates into strings
                        map_xtickangle   = 90;    
                        xgap_left        = 1;            %/ Avoid clipping or a gap at the 1st bar. 
                        xgap_right       = 1;            %/ Avoid clipping or a gap at the last bar. 
                    else
                        line_xdata = 1:length(shared_year_list);
                        if contains(line_setname, 'CMIP6_MME')
                            if isempty(x_intvl)   x_intvl = 20;  end  %/ interval: 20 years
                        else
                            if isempty(x_intvl)   x_intvl = 10;  end  %/ interval: 10 years
                        end
                        if isempty(map_xticks) map_xticks = 1:x_intvl:length(line_xdata);  end
                        map_xticklabels  = shared_year_list(map_xticks); 
                        map_xtickangle   = 0;    
                        xgap_left        = 0.5;            %/ Avoid clipping or a gap at the 1st bar. 
                        xgap_right       = 0.5;            %/ Avoid clipping or a gap at the last bar.  
                    end
                    linestyle = '-';
                end
                
                %/ Obtain the reference ydata (the first non-NaN data) based on the mininum xdata
                ref = nan(1,size(line_ydata_ref,2));
                for j = 1:size(line_ydata_ref,2)
                    cond_nonNaN_y       = ~isnan(line_ydata_ref(:,j));
                    line_xdata_nonNaN_y = line_xdata;
                    line_xdata_nonNaN_y(~cond_nonNaN_y) = nan;
                    [B,I]  = sort(line_xdata_nonNaN_y, 'ascend');
                    ref(j) = line_ydata_ref(I(1),j);
                end
                line_ydata_ref = ref; 
    
                %/ Check and reverse the dim if line_xdata is not strictly increasing
                if diff(line_xdata(1:2)) < 0
                    line_xdata = flip(line_xdata, 1);
                    line_ydata = flip(line_ydata, 1);
                end
    
                if isequal(trend_abs_or_rel, 'rel')
                    str_abs_or_rel = strcat('_', trend_abs_or_rel);
                else
                    str_abs_or_rel = '';
                end
                    
                %/ Grid and label settings
                map_xlabel          = '';
                map_xticks_minor    = [];
                vert_lines          = []; 
                vert_lines_col      = 'r';
                str_dates_fig       = sprintf('_%d-%d%s', shared_year_list(1), shared_year_list(end), str_timescale);
                titlename           = strcat(basin_name, str_dates_fig, str_trend_test, str_N_moving_avg, str_abs_or_rel);
                FigName_underscore  = strrep(strcat('multilines_', line_setname, str_regime_ver, str_anom_base, str_trend_test, str_N_moving_avg, str_yerror, str_abs_or_rel, '_', basin_name, str_dates_fig), ' ', '_');
    
                [trend_data, pval_data, intercept_data, CI95_data, regress_x, regress_y] = plot_multi_ylines('line_ydata', line_ydata, 'line_ydata_ref', line_ydata_ref, 'line_xdata', line_xdata, 'line_lgd_str', line_lgd_str, 'line_axislabel', line_axislabel, 'line_colmap', line_colmap,...
                          'line_yerror_upper', line_yerror_upper, 'line_yerror_lower', line_yerror_lower, 'line_yerror_type', line_yerror_type, 'lineerror_colmap', lineerror_colmap, 'lineerror_alpha', lineerror_alpha, 'show_y1_as_bar', show_y1_as_bar, 'pve_nve_bar_mode', pve_nve_bar_mode, 'bar_pve_col', bar_pve_col, 'bar_nve_col', bar_nve_col, 'scatter_ydata', [], 'vert_lines', vert_lines, 'vert_lines_col', vert_lines_col,...
                          'yshift', yshift, 'show_y_max', show_y_max, 'show_zero_line', show_zero_line, 'show_axislabel', show_axislabel, 'map_xlabel', map_xlabel, 'map_xticks', map_xticks, 'map_xticklabels', map_xticklabels, 'map_xtickangle', map_xtickangle,...
                          'use_yyaxis', use_yyaxis, 'line_data_ylim', line_data_ylim, 'line_data_xlim', line_data_xlim, 'marker', marker, 'markersize', markersize, 'markeredgecolor', markeredgecolor, 'fontsize', fontsize, 'titlefontsize', titlefontsize, 'shift_title_y', shift_title_y, 'linewidth', linewidth, 'linestyle', linestyle,...
                          'titlename', titlename, 'FigName_underscore', FigName_underscore, 'show_legend', show_legend, 'lgd_position', lgd_position, 'plot_legend_only', plot_legend_only,...
                          'trend_mode', trend_mode, 'trend_abs_or_rel', trend_abs_or_rel, 'trend_alpha', trend_alpha, 'trend_test', trend_test, 'asterisk', asterisk, 'alpha_level', alpha_level, 'show_regress_line', show_regress_line, 'regress_linewidth', regress_linewidth, 'xgap_left', xgap_left, 'xgap_right', xgap_right,...
                          'fig_width', fig_width, 'fig_height', fig_height, 'fig_xpos', fig_xpos, 'plotting_folder', plotting_folder,...
                          'yyaxis1_col', yyaxis1_col, 'yyaxis2_col', yyaxis2_col, 'show_diagonal', show_diagonal, 'box_mode', box_mode, 'savefig', savefig);           
                
                %/ Plot multi-product ensemble P and E 
                if contains(line_setname, 'MultiProduct_P_E')
                    %======== P products ========%
                    if basin_list == -1 || basin_list == -2
                        EM_P_list   = {'ERA5_P', 'GPCC_P', 'GPCP_P', 'CRU_P'};  %/ If computing global values, we have to excluse the two regional dataset (TPR & HARv2)
                    else
                        EM_P_list   = {'ERA5_P', 'TPR_P', 'HARv2_P', 'GPCC_P', 'GPCP_P', 'CRU_P'};  %/ Excluding CERA
                    end
                    EM_P_data         = nan(size(line_ydata,1),length(EM_P_list)+1);
                    line_yerror_upper = nan(size(line_ydata,1),length(EM_P_list)+1);
                    line_yerror_lower = nan(size(line_ydata,1),length(EM_P_list)+1);
    
                    % ind = findismember_loop(slct_data_list, EM_P_list);
                    ind = findismember_loop(line_lgd_str',  strrep(strcat(basin_name, {' '}, EM_P_list), '_', ' '));
                    if length(ind) ~= length(EM_P_list)  
                        error('There are missing variables to compute ensemble mean!'); 
                    end
                    
                    EM_P_data(:,1:length(ind)) = line_ydata(:,ind); 
                    EM_P_data(:,end)           = mean(line_ydata(:,ind), 2); %/ use mean() but not nanmean() -> avoid artificial trend due to data scarcity
                    
                    if isequal(line_yerror_type, 'sd')
                        line_yerror_upper(:,end)   = std(line_ydata(:,ind), [], 2);
                        line_yerror_lower(:,end)   = std(line_ydata(:,ind), [], 2);
                    elseif isequal(line_yerror_type, 'percentile')
                        line_yerror_upper(:,end)   = prctile(line_ydata(:,ind), 75, 2);
                        line_yerror_lower(:,end)   = prctile(line_ydata(:,ind), 25, 2);
                    end
                    str_yerror = strcat('_', line_yerror_type);
    
                    %/ Finally, do moving averaging if required.
                    if ~isempty(N_running_avg)
                        EM_P_data  = my_movmean('data', EM_P_data,  'n', N_running_avg, 'dim', 1);
                        line_yerror_upper = my_movmean('data', line_yerror_upper, 'n', N_running_avg, 'dim', 1);
                        line_yerror_lower = my_movmean('data', line_yerror_lower, 'n', N_running_avg, 'dim', 1);
                        str_N_moving_avg = sprintf('_%drunningavg', N_running_avg);
                    else
                        str_N_moving_avg = '';
                    end
    
                    EM_start_from = findismember_loop(shared_year_list, 1983); % since then both EM_P and EM_E have values
                    
                    a = split(EM_P_list, '_');
                    str_EM_P_list = a(:,:,1);
                    EM_data_str   = strrep([str_EM_P_list, 'Ensemble'], '_', ' ');
                    EM_data_ylim  = repmat([0, 890], size(EM_P_data,2), 1);
                    EM_markersize = 12;
                    EM_fontsize   = 17.5;
                    EM_titlefontsize = 12;
                    EM_shift_title_y = 0.02;
                    EM_marker     = {'d', 's', '^', 'v', 'pentagram', 'hexagram', 'o'};
                    % EM_color      = [255 0 51]./255;
                    % EM_color      = [102  204  255]./255;
                    EM_color      = [0 102 204]./255;
                    line_colmap   = [repmat([.7 .7 .7], length(EM_P_list), 1); EM_color]; %/ append color for the multi-product
                    lineerror_colmap = line_colmap;
                    
                    yyaxis1_col = 'k'; 
                    yyaxis2_col = 'k';
    %                 yyaxis1_col = line_colmap(1,:);
    %                 yyaxis2_col = line_colmap(2,:);
                    lineerror_alpha  = 0.15;
                    titlename = strcat(basin_name, '_EM_P_and_products', str_dates_fig, str_trend_test);
                    FigName_underscore = strrep(strcat('multilines_EM_P_and_products', str_anom_base, str_trend_test, str_N_moving_avg, str_yerror, '_', basin_name, str_dates_fig), ' ', '_');
    
                    [EM_P_trend_data, EM_P_pval_data, EM_P_intercept_data, EM_P_CI95_data, EM_P_regress_x, EM_P_regress_y] = plot_multi_ylines('line_ydata', EM_P_data(EM_start_from:end,:), 'line_xdata', line_xdata(EM_start_from:end),...
                          'line_yerror_upper', line_yerror_upper(EM_start_from:end,:), 'line_yerror_lower', line_yerror_lower(EM_start_from:end,:), 'line_yerror_type', line_yerror_type, 'line_lgd_str', EM_data_str, 'line_axislabel', line_axislabel, 'line_colmap', line_colmap,...
                          'lineerror_colmap', lineerror_colmap, 'lineerror_alpha', lineerror_alpha, 'show_y1_as_bar', show_y1_as_bar, 'pve_nve_bar_mode', pve_nve_bar_mode, 'scatter_ydata', [], 'vert_lines', vert_lines, 'vert_lines_col', vert_lines_col,...
                          'show_y_max', show_y_max, 'show_zero_line', show_zero_line, 'show_axislabel', show_axislabel, 'map_xlabel', map_xlabel, 'map_xticks', map_xticks, 'map_xticklabels', map_xticklabels, 'map_xtickangle', map_xtickangle,...
                          'use_yyaxis', use_yyaxis, 'yyaxis1_col', yyaxis1_col,   'yyaxis2_col',  yyaxis2_col,   'line_data_ylim', EM_data_ylim, 'marker', EM_marker, 'markersize', EM_markersize, 'markeredgecolor', markeredgecolor, 'fontsize', EM_fontsize, 'titlefontsize', EM_titlefontsize,  'shift_title_y', EM_shift_title_y, 'linewidth', linewidth, 'linestyle', linestyle,...
                          'titlename', titlename, 'FigName_underscore', FigName_underscore, 'show_legend', show_legend, 'lgd_position', lgd_position, 'plot_legend_only', plot_legend_only,...
                          'trend_mode', trend_mode, 'trend_abs_or_rel', trend_abs_or_rel, 'trend_alpha', trend_alpha,  'trend_test', trend_test, 'asterisk', asterisk, 'alpha_level', alpha_level, 'show_regress_line', show_regress_line, 'regress_linewidth', regress_linewidth, 'xgap_left', xgap_left, 'xgap_right', xgap_right, 'fig_width', fig_width, 'fig_height', fig_height, 'fig_xpos', fig_xpos, 'plotting_folder', plotting_folder,...
                          'box_mode', box_mode, 'savefig', savefig);           
                    
                    %======== E products ========%
                    if basin_list == -1 || basin_list == -2
                        EM_E_list   = {'ERA5_E', 'GLEAM_E'};  %/ If computing global values, we have to excluse the two regional dataset (TPR & HARv2)
                    else
                        EM_E_list   = {'ERA5_E', 'TPR_E', 'HARv2_E', 'GLEAM_E'};  %/ Excluding CERA
                    end
                    EM_E_data   = nan(size(line_ydata,1),length(EM_E_list)+1);
                    line_yerror_upper = nan(size(line_ydata,1),length(EM_E_list)+1);
                    line_yerror_lower = nan(size(line_ydata,1),length(EM_E_list)+1);
                    
                    % ind = findismember_loop(slct_data_list, EM_E_list);

                    ind = findismember_loop(line_lgd_str',  strrep(strcat(basin_name, {' '}, EM_E_list), '_', ' '));
                    if length(ind) ~= length(EM_E_list)  
                        error('There are missing variables to compute ensemble mean!'); 
                    end
                    
                    EM_E_data(:,1:length(ind)) = line_ydata(:,ind); 
                    EM_E_data(:,end)           = mean(line_ydata(:,ind), 2); %/ use mean() but not nanmean() -> avoid artificial trend due to data scarcity
                    if isequal(line_yerror_type, 'sd')
                        line_yerror_upper(:,end)   = std(line_ydata(:,ind), [], 2);
                        line_yerror_lower(:,end)   = std(line_ydata(:,ind), [], 2);
                    elseif isequal(line_yerror_type, 'percentile')
                        line_yerror_upper(:,end)   = prctile(line_ydata(:,ind), 75, 2);
                        line_yerror_lower(:,end)   = prctile(line_ydata(:,ind), 25, 2);
                    end
                    
                    %/ Finally, do moving averaging if required.
                    if ~isempty(N_running_avg)
                        EM_E_data         = my_movmean('data', EM_E_data,  'n', N_running_avg, 'dim', 1);
                        line_yerror_upper = my_movmean('data', line_yerror_upper, 'n', N_running_avg, 'dim', 1);
                        line_yerror_lower = my_movmean('data', line_yerror_lower, 'n', N_running_avg, 'dim', 1);
                        str_N_moving_avg  = sprintf('_%drunningavg', N_running_avg);
                    else
                        str_N_moving_avg = '';
                    end
                    EM_start_from = findismember_loop(shared_year_list, 1983); % since then both EM_P and EM_E have values
                    
                    a = split(EM_E_list, '_');
                    str_EM_E_list = a(:,:,1);
                    EM_data_str   = strrep([str_EM_E_list, 'Ensemble'], '_', ' ');
                    EM_data_ylim  = repmat([105, 420], size(EM_E_data,2), 1);
                    EM_markersize = 12;
                    EM_fontsize   = 17.5;
                    EM_titlefontsize = 12;
                    EM_shift_title_y = 0.02;
                    EM_marker = {'d', 's', '^', 'v', 'o'};
                    % EM_color      = [0 204  0]./255;
                    EM_color      = [102  204  51]./255;
                    line_colmap   = [repmat([.7 .7 .7], length(EM_E_list), 1); EM_color]; %/ append color for the multi-product
                    lineerror_colmap = line_colmap;
                    titlename = strcat(basin_name, '_EM_E_and_products', str_dates_fig, str_trend_test);
                    FigName_underscore = strrep(strcat('multilines_EM_E_and_products', str_anom_base, str_trend_test, str_N_moving_avg, str_yerror, '_', basin_name, str_dates_fig), ' ', '_');
                    
                    [EM_E_trend_data, EM_E_pval_data, EM_E_intercept_data, EM_E_CI95_data, EM_E_regress_x, EM_E_regress_y] = plot_multi_ylines('line_ydata', EM_E_data(EM_start_from:end,:), 'line_xdata', line_xdata(EM_start_from:end),...
                          'line_yerror_upper', line_yerror_upper(EM_start_from:end,:), 'line_yerror_lower', line_yerror_lower(EM_start_from:end,:), 'line_yerror_type', line_yerror_type, 'line_lgd_str', EM_data_str, 'line_axislabel', line_axislabel, 'line_colmap', line_colmap,...
                          'lineerror_colmap', lineerror_colmap, 'lineerror_alpha', lineerror_alpha, 'show_y1_as_bar', show_y1_as_bar, 'pve_nve_bar_mode', pve_nve_bar_mode, 'scatter_ydata', [], 'vert_lines', vert_lines, 'vert_lines_col', vert_lines_col,...
                          'show_y_max', show_y_max, 'show_zero_line', show_zero_line, 'show_axislabel', show_axislabel, 'map_xlabel', map_xlabel, 'map_xticks', map_xticks, 'map_xticklabels', map_xticklabels, 'map_xtickangle', map_xtickangle,...
                          'use_yyaxis', use_yyaxis, 'yyaxis1_col', yyaxis1_col,   'yyaxis2_col',  yyaxis2_col,   'line_data_ylim', EM_data_ylim, 'marker', EM_marker, 'markersize', EM_markersize, 'markeredgecolor', markeredgecolor, 'fontsize', EM_fontsize, 'titlefontsize', EM_titlefontsize, 'shift_title_y', EM_shift_title_y, 'linewidth', linewidth, 'linestyle', linestyle,...
                          'titlename', titlename, 'FigName_underscore', FigName_underscore, 'show_legend', show_legend, 'lgd_position', lgd_position, 'plot_legend_only', plot_legend_only,...
                          'trend_mode', trend_mode, 'trend_abs_or_rel', trend_abs_or_rel, 'trend_alpha', trend_alpha,  'trend_test', trend_test, 'asterisk', asterisk, 'alpha_level', alpha_level, 'show_regress_line', show_regress_line, 'regress_linewidth', regress_linewidth, 'xgap_left', xgap_left, 'xgap_right', xgap_right, 'fig_width', fig_width, 'fig_height', fig_height, 'fig_xpos', fig_xpos, 'plotting_folder', plotting_folder,...
                          'box_mode', box_mode, 'savefig', savefig);               
                end
    
                %/ Plot multi-product ensemble EP ratio
                if contains(line_setname, 'MultiProduct_EP_ratio') 
                    %======== P products ========%
                    if basin_list == -1 || basin_list == -2
                        EM_E_list   = {'ERA5_E', 'GLEAM_E'}; 
                        EM_P_list   = {'ERA5_P', 'GPCC_P', 'GPCP_P', 'CRU_P'};  %/ If computing global values, we have to excluse the two regional dataset (TPR & HARv2)
                    else
                        EM_E_list   = {'ERA5_E', 'TPR_E', 'HARv2_E', 'GLEAM_E'};
                        EM_P_list   = {'ERA5_P', 'TPR_P', 'HARv2_P', 'GPCC_P',  'GPCP_P',  'CRU_P'};  %/ Excluding CERA
                    end
    
                    EM_EP_ratio_data  = nan(size(line_ydata,1),length(EM_P_list)+1);
                    line_yerror_upper = nan(size(line_ydata,1),length(EM_P_list)+1);
                    line_yerror_lower = nan(size(line_ydata,1),length(EM_P_list)+1);
                    EM_EP_ratio_unit  = repmat({'%'}, length(EM_P_list)+1, 1);
                    EM_EP_ratio_str_period = cell(length(EM_P_list)+1, 1);
                    
                    %/ Compute the EP ratio for each product 
                    for j = 1:length(EM_P_list)
                        ind_P = findismember_loop(line_lgd_str',  strrep(strcat(basin_name, {' '}, EM_P_list{j}), '_', ' '));
                        if isempty(ind_P)  error('%s is missing!', EM_P_list{j}); end
    
                        if isequal(EM_P_list{j}, 'ERA5_P')
                            ind_E = findismember_loop(line_lgd_str',  strrep(strcat(basin_name, {' '}, 'ERA5_E'), '_', ' '));
                        elseif isequal(EM_P_list{j}, 'TPR_P')
                            ind_E = findismember_loop(line_lgd_str',  strrep(strcat(basin_name, {' '}, 'TPR_E'), '_', ' '));
                        elseif isequal(EM_P_list{j}, 'HARv2_P')
                            ind_E = findismember_loop(line_lgd_str',  strrep(strcat(basin_name, {' '}, 'HARv2_E'), '_', ' '));
                        else
                            ind_E = findismember_loop(line_lgd_str',  strrep(strcat(basin_name, {' '}, 'GLEAM_E'), '_', ' '));
                        end
                        if isempty(ind_E)  
                            error('Missing evaporantion for %s to compute the EP ratio!', EM_P_list{j}); 
                        end
    
                        EM_EP_ratio_data(:,j) = line_ydata(:,ind_E)./line_ydata(:,ind_P)*100;
                        EM_EP_ratio_str_period{j} = str_period{ind_P};
                    end
    
                    %/ Finally, compute the multiproduct EP ratio
                    EM_EP_ratio_data(:,end)     = mean(EM_EP_ratio_data(:,1:end-1), 2); %/ use mean() but not nanmean() -> avoid artificial trend due to data scarcity
                    if isequal(line_yerror_type, 'sd')
                        line_yerror_upper(:,end)    = std(EM_EP_ratio_data(:,1:end-1), [], 2);
                        line_yerror_lower(:,end)    = std(EM_EP_ratio_data(:,1:end-1), [], 2);
                    elseif isequal(line_yerror_type, 'percentile')
                        line_yerror_upper(:,end)    = prctile(EM_EP_ratio_data(:,1:end-1), 75, 2);
                        line_yerror_lower(:,end)    = prctile(EM_EP_ratio_data(:,1:end-1), 25, 2);
                    end
                    str_yerror = strcat('_', line_yerror_type);
    
                    %/ Finally, do moving averaging if required.
                    if ~isempty(N_running_avg)
                        EM_EP_ratio_data  = my_movmean('data', EM_EP_ratio_data,  'n', N_running_avg, 'dim', 1);
                        line_yerror_upper = my_movmean('data', line_yerror_upper, 'n', N_running_avg, 'dim', 1);
                        line_yerror_lower = my_movmean('data', line_yerror_lower, 'n', N_running_avg, 'dim', 1);
                        str_N_moving_avg = sprintf('_%drunningavg', N_running_avg);
                    else
                        str_N_moving_avg = '';
                    end
    
                    EM_EP_ratio_str_period(end) = {'1983-2019'};
                    EM_start_from = findismember_loop(shared_year_list, 1983); % since then both EM_P and EM_E have values
                    
                    a = split(EM_P_list, '_');
                    str_EM_P_list = a(:,:,1);
                    EM_data_str   = strrep([str_EM_P_list, 'Ensemble'], '_', ' ')';
                    % EM_data_str   = strrep(strcat([str_EM_P_list, 'Ensemble'], {' EP ratio'}), '_', ' ')'; %/ append
                    % EM_data_ylim  = repmat([25, 75], size(EM_EP_ratio_data,2), 1);
                    EM_data_ylim  = repmat([17, 75], size(EM_EP_ratio_data,2), 1);
                    EM_markersize = 12;
                    EM_fontsize   = 17.5;
                    EM_titlefontsize = 12;
                    EM_shift_title_y = 0.02;
                    EM_marker     = {'d', 's', '^', 'v', 'pentagram', 'hexagram', 'o'};
                    EM_color      = [255 0 51]./255;
                    % EM_color      = [153   0  102;]./255;
                    line_colmap   = [repmat([.7 .7 .7], length(EM_P_list), 1); EM_color];  %/ append color for the multi-product
                    lineerror_colmap = line_colmap;
                    
                    yyaxis1_col = 'k'; 
                    yyaxis2_col = 'k';
                    lineerror_alpha  = 0.15;
                    titlename = strcat(basin_name, '_EM_EP_ratio_and_products', str_dates_fig, str_trend_test);
                    FigName_underscore = strrep(strcat('multilines_EM_EP_ratio_and_products', str_anom_base, str_trend_test, str_N_moving_avg, str_yerror, '_', basin_name, str_dates_fig), ' ', '_');
                    
                    [EM_EP_ratio_trend_data, EM_EP_ratio_pval_data, EM_EP_ratio_intercept_data, EM_EP_ratio_CI95_data, EM_EP_ratio_regress_x, EM_EP_ratio_regress_y] = plot_multi_ylines('line_ydata', EM_EP_ratio_data(EM_start_from:end,:), 'line_xdata', line_xdata(EM_start_from:end),...
                          'line_yerror_upper', line_yerror_upper(EM_start_from:end,:), 'line_yerror_lower', line_yerror_lower(EM_start_from:end,:), 'line_yerror_type', line_yerror_type, 'line_lgd_str', EM_data_str, 'line_axislabel', EM_EP_ratio_unit, 'line_colmap', line_colmap,...
                          'lineerror_colmap', lineerror_colmap, 'lineerror_alpha', lineerror_alpha, 'show_y1_as_bar', show_y1_as_bar, 'pve_nve_bar_mode', pve_nve_bar_mode, 'scatter_ydata', [], 'vert_lines', vert_lines, 'vert_lines_col', vert_lines_col,...
                          'show_y_max', show_y_max, 'show_zero_line', show_zero_line, 'show_axislabel', show_axislabel, 'map_xlabel', map_xlabel, 'map_xticks', map_xticks, 'map_xticklabels', map_xticklabels, 'map_xtickangle', map_xtickangle,...
                          'use_yyaxis', use_yyaxis, 'yyaxis1_col', yyaxis1_col,   'yyaxis2_col',  yyaxis2_col,   'line_data_ylim', EM_data_ylim, 'marker', EM_marker, 'markersize', EM_markersize, 'markeredgecolor', markeredgecolor, 'fontsize', EM_fontsize, 'titlefontsize', EM_titlefontsize, 'shift_title_y', EM_shift_title_y, 'linewidth', linewidth, 'linestyle', linestyle,...
                          'titlename', titlename, 'FigName_underscore', FigName_underscore, 'show_legend', show_legend, 'lgd_position', lgd_position, 'plot_legend_only', plot_legend_only,...
                          'trend_mode', trend_mode, 'trend_abs_or_rel', trend_abs_or_rel, 'trend_alpha', trend_alpha,  'trend_test', trend_test, 'asterisk', asterisk, 'alpha_level', alpha_level, 'show_regress_line', show_regress_line, 'regress_linewidth', regress_linewidth, 'xgap_left', xgap_left, 'xgap_right', xgap_right, 'fig_width', fig_width, 'fig_height', fig_height, 'fig_xpos', fig_xpos, 'plotting_folder', plotting_folder,...
                          'box_mode', box_mode, 'savefig', savefig);           
                end
    
                %/ Print the table for trend info
                if trend_mode
                    %/ Append the EM data (if needed)
                    if isequal(line_setname, 'MultiProduct_P_E')
                        trend_data     = [trend_data;   EM_P_trend_data(end);   EM_E_trend_data(end)];
                        pval_data      = [pval_data;    EM_P_pval_data(end);    EM_E_pval_data(end)];
                        CI95_data      = [CI95_data;    EM_P_CI95_data(end);    EM_E_CI95_data(end)];
                        line_ydata     = [line_ydata,   EM_P_data(:,end),       EM_E_data(:,end)];
                        line_lgd_str   = [line_lgd_str; {'Ensemble P'; 'Ensemble E'}];
                        line_axislabel = [line_axislabel; line_axislabel(end-1:end)];
                        str_period     = [str_period;   {'1983-2019'; '1983-2019'}];

                        %/ Get the first nonNaN ydata at the minimum xdata as the reference value for computing relative change
                        cond_nonNaN_y       = ~isnan(EM_P_data(:,end));
                        line_xdata_nonNaN_y = line_xdata;
                        line_xdata_nonNaN_y(~cond_nonNaN_y) = nan;
                        [~,I]    = sort(line_xdata_nonNaN_y, 'ascend');
                        EM_P_ref = EM_P_data(I(1),end);
                        
                        %/ Same, but for EM_E_ref
                        cond_nonNaN_y       = ~isnan(EM_E_data(:,end));
                        line_xdata_nonNaN_y = line_xdata;
                        line_xdata_nonNaN_y(~cond_nonNaN_y) = nan;
                        [~,I]    = sort(line_xdata_nonNaN_y, 'ascend');
                        EM_E_ref = EM_E_data(I(1),end);

                        line_ydata_ref = [line_ydata_ref, EM_P_ref, EM_E_ref]; 

                    elseif isequal(line_setname, 'MultiProduct_EP_ratio')
                        trend_data     = EM_EP_ratio_trend_data;
                        pval_data      = EM_EP_ratio_pval_data;
                        CI95_data      = EM_EP_ratio_CI95_data;
                        line_ydata     = EM_EP_ratio_data;
                        line_lgd_str   = EM_data_str;
                        % line_lgd_str   = [line_lgd_str; {'Ensemble EP ratio'}];
                        line_axislabel = EM_EP_ratio_unit;
                        str_period     = EM_EP_ratio_str_period;

                        %/ Recompute line_ydata_ref (since EP_ratio was post-processed) 
                        line_ydata_ref = nan(1, size(EM_EP_ratio_data,2));
                        for j = 1:length(line_ydata_ref)
                            cond_nonNaN_y       = ~isnan(EM_EP_ratio_data(:,j));
                            line_xdata_nonNaN_y = line_xdata;
                            line_xdata_nonNaN_y(~cond_nonNaN_y) = nan;
                            [~,I]    = sort(line_xdata_nonNaN_y, 'ascend');
                            line_ydata_ref(j) = EM_EP_ratio_data(I(1),j);
                        end
                    end
                    
                    if isempty(slct_xdata)              %/ that means we computed the trend of the time series
                        per_x_unit = {' decade^{-1}'};
                        trend_multiplier = 10;          %/ per yr -> per dec
                    elseif contains(slct_xdata, 'T2m') || contains(slct_xdata, 'tas') %/ that means we regress y onto surface temp.
                        per_x_unit = {' K^{-1}'};
                        trend_multiplier = 1;  
                    elseif isequal(slct_xdata, 'CERA_E')     %/ that means we regress y onto E
                        per_x_unit = {' mm yr^{-1}'};
                        trend_multiplier = 1;  
                    elseif contains(slct_xdata, {'EP_ratio', 'evspsbl_pr_ratio'})  %/ that means we regress y onto E
                        per_x_unit = {' %^{-1}'};
                        trend_multiplier = 1;  
                    else
                        error('code not set!')
                    end
                    
                    %/ Use asterisks to indicate significance level
                    str_sig = repmat({''}, length(pval_data), 1);
                    for i = 1:length(alpha_level)
                        ind = find(pval_data <= alpha_level(i));
                        str_sig(ind) = repmat(asterisk(i), length(ind), 1);
                    end
    
                    if isequal(trend_abs_or_rel, 'rel')
                        line_axislabel = repmat({'%'}, length(line_axislabel), 1);
                    end
    
                    trend_per_dec = strcat(line_lgd_str, {': '}, num2str(round(trend_data*trend_multiplier, 2, 'significant')),......
                                           char(177), num2str(round(CI95_data*trend_multiplier, 2, 'significant')),...
                                           line_axislabel, per_x_unit, {' (p = '}, num2str(round(pval_data, 3)), ')', str_sig, {' - '}, str_period);
                    fprintf('*** %s, %s, trend test = %s (%s)***\n', basin_name, strrep(str_timescale, '_', ' '), trend_test, strrep(str_N_moving_avg, '_', ''))
                    disp(trend_per_dec)
                    
                    %/ Get the start point of the regression line at which the true value is not a NaN
                    trend_per_dec_inpercentage = strcat(line_lgd_str, {': '}, num2str(round(trend_data*trend_multiplier./line_ydata_ref'*100, 2, 'significant')),...
                                                        char(177), num2str(round(CI95_data*trend_multiplier./line_ydata_ref'*100, 2, 'significant')), {'%'}, per_x_unit,...
                                                        {' (p = '}, num2str(round(pval_data, 3)), ')', str_sig, {' - '}, str_period);
                    disp(trend_per_dec_inpercentage)
                    
                    %/ Correlation matrix
                    %/ NOTE: 'Rows','Complete' will run into a bug when the
                    %/       matrix contains any column of all NaN.
                    %/       'Rows','pairwise' is more robust.
                    if size(line_xdata, 1) == 1
                        X = line_xdata';
                    else
                        X = line_xdata;
                    end
                    Y = line_ydata;
                    [R, P] = corrcoef([X, Y], 'Rows','pairwise');
    
                    str_sig = repmat({''}, size(R,1), size(R,2));
                    for i = 1:length(alpha_level)
                        cond = (P < alpha_level(i));
                        str_sig(cond) = asterisk(i);
                    end
                    
                    if size(line_lgd_str, 1) == 1
                        if ~isempty(slct_xdata)
                            VariableNames = strrep([{slct_xdata},line_lgd_str], ' ', '_');
                        else
                            VariableNames = strrep([{'x'},line_lgd_str], ' ', '_');
                        end
                    else
                        if ~isempty(slct_xdata)
                            VariableNames = strrep([{slct_xdata};line_lgd_str], ' ', '_');
                        else
                            VariableNames = strrep([{'x'};line_lgd_str], ' ', '_');
                        end
                    end
                    RowNames = VariableNames;
                    A = strcat(string(round(R,3)), str_sig);
                    T = array2table(A, 'RowNames', RowNames, 'VariableNames', VariableNames);
                    disp(T)
                end

            end
        end
    end
end

if plot_table && plot_or_not
    %/ Show the 'Unattributed' category in the table
    slct_data_list_bc     = slct_data_list;
    slct_src_list_bc      = slct_src_list;
    grouping_mode_list_bc = grouping_mode_list;
    ins_or_acc_list_bc    = ins_or_acc_list;
    if ~isempty([grouping_mode_list{:}])   %/ first check if it is an empty cell (otherwise will cause bugs)
        if any(ismember(grouping_mode_list, {'domain-wise-exact'}))
            slct_src_list_bc(end+1)      = {'Unattributed'};
            ind_local = findismember_loop(slct_src_list, 'local');
            slct_data_list_bc(end+1)     = slct_data_list(ind_local);
            grouping_mode_list_bc(end+1) = grouping_mode_list(ind_local); %/ copy the grouping mode of 'local' and paste it for 'Unattributed'
            ins_or_acc_list_bc(end+1)    = ins_or_acc_list(ind_local);    %/ copy the ins_or_acc    of 'local' and paste it for 'Unattributed'
        end
    end
    table_data    = nan(length(basin_list)*length(mth_list), length(slct_data_list_bc));
    table_data_sd = nan(length(basin_list)*length(mth_list), length(slct_data_list_bc));
    RowNames      = cell(length(basin_list)*length(mth_list), 1);
    
    
    cnt = 0; 
    for top = basin_list
        basin_name = ''; basin_bndry = []; str_basin = [];
        if top == -1  %/ global 
            basin_name = 'global';
            basin_bndry = [];
        elseif top == -2  %/ global land 
            basin_name = 'land';
            basin_bndry = [];
        else
            if from_basin 
                if from_basin == 4   
                    if top == 0             %/ then show all the basins 
                        basin_name = 'TP';  %/ Aggregated P_LA/Pm/Pm_BL from sub-basins were processed at [Step 5.1]
                        basin_bndry = TP_bndry;
                    elseif top == 16     %/ S. Inner TP + N. Inner TP
                        ind = [13, 14];
                        basin_name = 'Inner_TP';
                        bndry_data = cat(1,basin_catalog(ind).bndry);
                    elseif top == 17     %/ TP-Indus + TP-Ganges
                        ind = [1, 2];
                        basin_name = 'TP-Indus_TP-Ganges';
                        bndry_data = cat(1,basin_catalog(ind).bndry);
                    else
                        basin_name  = basin_catalog(top).name{:};
                        basin_bndry = basin_catalog(top).bndry;
                    end
                else
                    basin_name  = basin_catalog(top).name{:};
                    basin_bndry = basin_catalog(top).bndry;
                end
                str_basin = strcat('_', basin_name);
            end
        end
        basin_name_strrep = strrep(strrep(strrep(basin_name, '-', '_'), ' ', '_'), '.', '_');
        
        seasonal_ts = nan(4, length(shared_year_list), length(slct_data_list_bc));
        for mth = mth_list
            if mth == 16 && length(shared_year_list) == 1   continue;  end   %/ skip DJF if only one year
            cnt = cnt + 1;
            if mth == 0        
                str_timescale = '_annual';    
                str_dates     = str_years_meteo; %/ Do NOT modify it, otherwise will affect the reading of the processed WSV data.
            else
                str_timescale = strcat('_', str_mth{mth});
                str_dates     = strcat(str_years_meteo, str_timescale); 
            end
            RowNames{cnt} = strcat(basin_name,str_timescale);

            %===== line_ydata (ntime, nvar) =====%
            line_axislabel = cell(length(slct_data_list_bc), 1);
            for j = 1:length(slct_data_list_bc)
                slct_data     = slct_data_list_bc{j};
                if ~isempty(model_list)  model = model_list{j};  else  model = '';  end
                if ~isempty(exp_list)    exp   = exp_list{j};    else  exp   = '';  end
                if ~isempty(ens_list)    ens   = ens_list{j};    else  ens   = '';  end
                %/ Set slct_data_callfile, slct_data_new to handle CMIP6 data    
                if ~isempty(model) && ~isempty(exp) && ~isempty(ens)
                    slct_data_new      = strcat(slct_data,'_',model,'_',exp,'_',ens);
                else
                    slct_data_new      = slct_data;
                end

                grouping_mode = grouping_mode_list_bc{j};
                if ~isempty(grouping_mode)
                    str_grouping_mode = strcat('_', strrep(grouping_mode, '-', '_'));
                else
                    str_grouping_mode = '';
                end

                if contains(slct_data_new, {'Pm'})
                    %/ As we allow for seasonally varying regimes, we need to add up the seasonal values to the annual one
                    %/ CAVEAT: we also need to do it for Local Recycling!
                    ind = findismember_loop(SR.(basin_name_strrep).([slct_data_new, str_grouping_mode, '_src_name']), slct_src_list_bc{j});
                    table_data(cnt,j)       = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale, '_clim'))(ind, :);
                    table_data_sd(cnt,j)    = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale, '_sd'))(ind, :);
                    line_axislabel{j}       = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale, '_unit'));
                else
                    ind = 1;
                    table_data(cnt,j)       = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale, '_clim'))(ind, :);
                    table_data_sd(cnt,j)    = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale, '_sd'))(ind, :);
                    line_axislabel{j}       = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale, '_unit'));
                end
                %/ Store seasonal time series in case we compute annual value
                if ismember(mth , [13:16]) && ~isequal(slct_data_new, 'EP_ratio')
                    seasonal_ts(mth-12,:,j) = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale))(ind, :);
                end
            end
        end
    end
    VariableNames = strcat(slct_data_list_bc, {' '}, slct_src_list_bc);
    VariableNames = strrep(strrep(VariableNames, '/', 'over'), ' ', '_');
    VariableNames_withunit = strcat(VariableNames', {' ('}, line_axislabel, {')'});  % {'\x00B1'}, 
    
    %/ Contributions from 'Elsewhere'/'Unattributed (Error Propagation)
    if isequal(line_setname, 'VarSet_Pm_frac_adj')
        ind = findismember_loop(slct_src_list_bc, ['local', labels_domain_wise_exact]);
        elsewhere = 100 - sum(table_data(:,ind), 2);
        elsewhere_sd = sqrt(sum(table_data_sd(:,ind).^2, 2));
        table_data(:,end+1) = elsewhere;
        table_data_sd(:,end+1) = elsewhere_sd;
        VariableNames(end+1) = {'Elsewhere'};
%         fprintf('*** Annual contributions from ''elsewhere'': %.2f%s%.2f%% ***\n', round(annual_elsewhere, 1), char(177), round(annual_elsewhere_sd, 1));
    end

    T_cell = cellstr(strcat(string(round(table_data, 1)), {' '}, char(177), string(round(table_data_sd, 1)) ));
    T      = cell2table(T_cell, 'VariableNames', VariableNames, 'RowNames', RowNames);
    fprintf('*** mth = %d ***\n', mth)
    disp(T)
%     disp(VariableNames_withunit')

end

% a = [0.531, 0.2, 0.1, 0.1-0.031];
% a_adj = a/sum(a)
% 
% b = [0.531, 0.2, 0.15, 0.1-0.031];
% b_adj = b/sum(b)

%% Send email

send_email_to    = 'tfcheng@hawaii.edu'; %/ Set your email to send notification after the program completed. Do NOT use connect.ust.hk email. It will just go to '/var/spool/mail/<Your Account>' in the ust server.
server_name      = getenv('HOSTNAME');          %/ get the server name
server_name_cell = strsplit(server_name, '.');
server_name      = server_name_cell{1};

email_subject = sprintf('[%s] Job Progress: get_SR_program Completed', server_name);
first_sentence = sprintf('Job Progress: Your Matlab Program ''get_SR_program.m'' has just finished a job!\n');
email_content = [sprintf('Hello from %s!\n', server_name),...
     newline,...
     first_sentence];

% disp(email_subject)
% disp(email_content)
command = char(sprintf('echo "%s" | mail -s "%s" %s', email_content, email_subject, send_email_to));
[status, ~] = system(command);

if status == 0
    fprintf('=======================================================\n')
    fprintf('=== A progress notification sent to %s ===\n', send_email_to)
    fprintf('=======================================================\n')
else
    error('something wrong with sending a notification email.') 
end



%% Exit the program 
exit;  %/ [IMPORTANT]: This prevents 'Warning: Error reading character from command line' that causes the program to hang on forever..

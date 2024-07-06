%% [Step 1] Read data and FLEXPART-WaterSip experiments parameters (see also moisture_tracking.m)
% clearvars;

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

%/ Remark when saving into nc file
nc_remark       = 'Author: Franklin Cheng (franklin.cheng@ust.hk); Reference: Cheng et al. (2024), Increased sensitivity of Tibetan Plateau precipitation to local evapotranspiration under anthropogenic warming, Communications Earth and Environment';

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
% expmnt = 'domfill_CERA_MPI'; ldirect = 'bwd'; from_basin = 6; output_res = 1; dt = 3; year_list = 1971:2010; optimal_rr = 0.99; RHc_dqc_scheme = 23;  
% expmnt = 'domfill_CERA_MPI'; ldirect = 'bwd'; from_basin = 5; output_res = 1; dt = 3; year_list = 1971:2010; optimal_rr = 0.99; RHc_dqc_scheme = 23; 
expmnt = 'domfill_CERA_MPI'; ldirect = 'bwd'; from_basin = 4; output_res = 1; dt = 3; year_list = 1971:2010; optimal_rr = 0.99; RHc_dqc_scheme = 23; 
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

%/ Load Land Cover Type from CERA
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
disp(data_folder);





%% [Step 3] read_any (reanalysis / gauge / satellite / multi-product ensemble mean)
close all; time_dim = 3; stlevel = []; edlevel = []; slct_level = []; lon_range = []; lat_range = []; plot_to_check = 0; dataset = []; 
savefig = 0; 

savemat            = 1;
recompute          = 1;   %/ <- mind this!
set_conden_to_zero = 1;   %/ Default to set condensation to zeros for evap and prcp products.
NumWorkers         = 10;
load_or_not        = 0; 
if isequal(project_name, 'CMIP6_MJO')
    noleap = 1;    %/ Omitting leap days ensures a consistent length of days with CMIP6 models (most used a 365-day calendar)
else
    noleap = 0;  
end

select_field = {'monthly'}; %/ 'daily', 'subdaily', 'dailyAnom_PMC', 'monthly'
% slct_year    = 1971:2010;
slct_year    = 1983:2014;

slct_data = {'CERA_P', 'CERA_E'}; select_field = {'daily'}; slct_year = 1971:2010;
% slct_data = {'CERA_SRO'};
% slct_data = {'CERA_dTWS'};
% slct_data = {'ERA5_P'};
% slct_data = {'q500Omega300'}); slct_year = 1971:2010;
% slct_data = {'ERA5_u850', 'ERA5_v850', 'ERA5_u500', 'ERA5_v500', 'ERA5_u300', 'ERA5_v300'}); slct_year = 1971:2010;
% slct_data = {'ERA5_P', 'GPCC_P'}); slct_year = 1983:2019;
% slct_data = {'ERA5_u500', 'ERA5_v500'}); slct_year = 1983:2019;
% slct_data = {'ERA5_u700', 'ERA5_v700'}); slct_year = 1983:2019;
% slct_data = {'CRU_P'}); slct_year = 1971:2020; select_field = {'monthly'};
% slct_data = {'CRU_P', 'GPCP_P', 'HARv2_P', 'GPCC_P', 'ERA5_P', 'TPR_P', 'CERA_P'}); slct_year = 1983:2010;
% slct_data = {'CRU_P', 'GPCP_P', 'HARv2_P', 'GPCC_P', 'ERA5_P', 'TPR_P'}; slct_year = 1983:2014;
% slct_data = {'EM_P'}; slct_year = 1983:2014;
% slct_data = {'CRU_P', 'GPCP_P', 'HARv2_P', 'GPCC_P', 'ERA5_P', 'TPR_P', 'EM_P'}; slct_year = 1983:2019;
% slct_data = {'EM_P'}; slct_year = 1983:2019;
% slct_data = {'GLEAM_E', 'HARv2_E', 'TPR_E', 'ERA5_E', 'EM_E'}; slct_year = 1983:2019;
% slct_data = {'TPR_E', 'ERA5_E', 'EM_E'}; slct_year = 1983:2019;
% slct_data = {'TPR_E'}; slct_year = 1983:2019;

% slct_data = {'CRU_P', 'GPCC_P', 'GPCP_P', 'ERA5_P', 'TPR_P', 'HARv2_P',...
%              'GLEAM_E', 'HARv2_E', 'TPR_E', 'ERA5_E'}); slct_year = 1983:2014;
% slct_data = {'EM_P', 'EM_E'}; slct_year = 1983:2010; plot_to_check = 1; savefig = 0;

dataset = read_any('dataset', dataset,  'param', param, 'slct_data',  slct_data, 'select_field', select_field, 'slct_year', slct_year,...
                   'time_dim', time_dim, 'noleap', noleap, 'stlevel', stlevel, 'edlevel', edlevel, 'slct_level', slct_level, 'lon_range', lon_range, 'lat_range', lat_range,...
                   'data_folder', data_folder, 'savemat', savemat, 'recompute', recompute, 'set_conden_to_zero', set_conden_to_zero,...
                   'load_or_not', load_or_not, 'NumWorkers', NumWorkers, 'plot_to_check', plot_to_check, 'savefig', savefig);


%% ERA5 module

error('paramter module (ERA5) is deprecated. Load param_universal.m instead!')

% dataset.placeholder = [];
% 
% data_path_parts = [
%    {'/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/'}, {'ERA5_'}, {'_hrly_mth01-12_1x1_0_360E_-40_90N_'};
%    {'/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/'}, {'ERA5_'}, {'_hrly_mth01-12_1x1_0_360E_-90_90N_'};
%    {'/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/'}, {'ERA5_'}, {'_3hrly_mth01-12_1x1_0_360E_-90_90N_'};                         %/ 3 
%    {'/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/'}, {'ERA5_'}, {'_hrly_mth01-12_0.25x0.25_90_140E_10_50N_lv700-925_'};...
%    {'/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/'}, {'ERA5_'}, {'_6hrly_mth01-12_1x1_90_140E_10_70N_lv100-1000_'}; 
%    {'/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/'}, {'ERA5_'}, {'_6hrly_mth01-12_1.5x1.5_0_360E_-90_90N_lv100-1000_'};          %/ 6 
%    {'/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/'}, {'ERA5_'}, {'_6hrly_mth01-12_1x1_0_360E_-90_90N_lv100-1000_'}; 
%    {'/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/'}, {'ERA5_'}, {'_6hrly_mth01-12_1x1_40_180E_-30_50N_lv100-1000_'};];
% 
% % filename_part   = {'_hrly_mth01-12_1x1_0_360E_-40_90N_',...                  %/ Semi-Global        
% %                    '_hrly_mth01-12_1x1_0_360E_-90_90N_',...                  %/ Global
% %                    '_3hrly_mth01-12_1x1_0_360E_-90_90N_',...                 %/ Global
% %                    '_hrly_mth01-12_0.25x0.25_90_140E_10_50N_lv700-925_',...
% %                    '_6hrly_mth01-12_1x1_90_140E_10_70N_lv100-1000_',...
% %                    '_6hrly_mth01-12_1.5x1.5_0_360E_-90_90N_lv100-1000_',...
% %                    '_6hrly_mth01-12_1x1_0_360E_-90_90N_lv100-1000_',...
% %                    '_6hrly_mth01-12_1x1_40_180E_-30_50N_lv100-1000_'};       
% % datafolder      = repmat({'/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/'}, length(filename_part), 1);
% % filename_prefix = repmat({'ERA5_'}, length(filename_part), 1);
% 
% %/ Non-reanalysis data
% varname_nonrean           = {'NOAA_CDR_OLR', 'EM_P', 'EM_E', 'CRU_P', 'GPCP_P', 'TPR_P', 'TPR_E', 'TPR_SRO', 'TPR_SSRO', 'TPR_RO', 'HadCRUT5_T2m', 'GPCC_P', 'HARv2_P', 'HARv2_E', 'HARv2_T2m', 'GLEAM_SMsurf', 'GLEAM_E', 'GLASS_E', 'CMORPH_P', 'IMERG_P', 'traj_freq', 'total_init_contr', 'avg_init_contr', 'contr_map', 'Pm', 'BL_Pm', 'LAI', 'PDSI', 'NDVI', 'E_minus_P', 'EmPSST', 'P_minus_E_gte0', 'E_minus_P_gte0', 'waterbudget_residual', 'prop_crop', 'prop_pasture', 'prop_urbn', 'prop_primary',  'prop_secd', 'biomass_harv_prim', 'biomass_secd',        'prop_vc_transition', 'prop_vp_transition', 'prop_vu_transition'}; 
% dataname_callfile_nonrean = {'NOAA_CDR_OLR', 'EM_P', 'EM_E', 'CRU_P', 'GPCP_P', 'TPR_P', 'TPR_E', 'TPR_SRO', 'TPR_SSRO', 'TPR_RO', 'HadCRUT5_T2m', 'GPCC_P', 'HARv2_P', 'HARv2_E', 'HARv2_T2m', 'GLEAM_SMsurf', 'GLEAM_E', 'GLASS_E', 'CMORPH_P', 'IMERG_P', 'traj_freq', 'total_init_contr', 'avg_init_contr', 'contr_map', 'Pm', 'BL_Pm', 'LAI', 'PDSI', 'NDVI', 'E_minus_P', 'EmPSST', 'P_minus_E_gte0', 'E_minus_P_gte0', 'waterbudget_residual',     'gcrop',        'gpast',     'gurbn',        'gothr',      'gsecd',             'gvbh1', 'gssmb',                            'gflvc',              'gflvp',              'gflvu'};
% dataname_nonrean          = {'NOAA_CDR_OLR', 'EM_P', 'EM_E', 'CRU_P', 'GPCP_P', 'TPR_P', 'TPR_E', 'TPR_SRO', 'TPR_SSRO', 'TPR_RO', 'HadCRUT5_T2m', 'GPCC_P', 'HARv2_P', 'HARv2_E', 'HARv2_T2m', 'GLEAM_SMsurf', 'GLEAM_E', 'GLASS_E', 'CMORPH_P', 'IMERG_P', 'traj_freq', 'total_init_contr', 'avg_init_contr', 'contr_map', 'Pm', 'BL_Pm', 'LAI', 'PDSI', 'NDVI', 'E_minus_P', 'EmPSST', 'P_minus_E_gte0', 'E_minus_P_gte0', 'waterbudget_residual',      'crop',      'pasture',     'urban',      'primary',  'secondary', 'biomass_harv_prim', 'biomass_density_secd',      'primary2crop',    'primary2pasture',      'primary2urban'};
% dataunitconv_nonrean      = ones(1, length(varname_nonrean));           %/ dummy
% datatype_nonrean          = repmat({''}, 1, length(varname_nonrean));   %/ dummy
% whichfolder_nonrean       = ones(1, length(varname_nonrean));           %/ dummy
% fc_steps_nonrean          = cell(1, length(varname_nonrean));           %/ dummy
% StepsInADay_nonrean       = ones(1, length(varname_nonrean));           %/ dummy
% 
% %/ Check whether the ERA5 field is ins or acc.
% %/ https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Table2
% %/ https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Table3
% 
% %/ A set of 3D data (keep it updated!)/%
% %/ NOTE: In ERA5, a downward flux is positive. So -ve E = evaporation; +ve E = condensation
% %/       OLR = -TTR. For hrly OLR, we consider it as an ins field after
% %/       multiply hrly data with -1/3600 -> W/m^2.
% varname_3D           = { 'sp',  'vimd',   'tp',   'e', 'tcwv', 'p71.162', 'p72.162',  'msl',     'z',    'z',    'z', 'Z200sw',     'u',     'v', 'u200sw',     'u',     'v',     'u',     'v',         'w',   'u10',   'v10', 't2m', 'sst',   'ttr',   'str',   'tsr',  'tisr',  'ssr',  'slhf',  'sshf', 'cape', 'cin', 'blh'};
% dataname_callfile_3D = { 'sp',  'VIMD', 'prcp',   'E', 'tcwv',     'IVT',     'IVT',  'slp',  'Z850', 'Z500', 'Z200', 'Z200sw', 'UV200', 'UV200', 'u200sw', 'UV500', 'UV500', 'UV850', 'UV850',  'omega500', 'UV10m', 'UV10m', 'T2m', 'SST',   'TTR',   'STR',   'TSR',  'TISR',  'SSR',  'SLHF',  'SSHF', 'CAPE', 'CIN', 'BLH'};
% dataname_3D          = { 'sp',  'VIMD',    'P',   'E',   'pw',    'uIVT',    'vIVT',  'slp',  'Z850', 'Z500', 'Z200', 'Z200sw',  'u200',  'v200', 'u200sw',  'u500',  'v500',  'u850',  'v850',  'omega500',  'u10m',  'v10m', 'T2m', 'SST',   'OLR',   'STR',   'TSR',  'TISR',  'SSR',  'SLHF',  'SSHF', 'CAPE', 'CIN', 'BLH'};  
% dataunitconv_3D      = [    1,       1,   1000, -1000,      1,         1,         1,  1/100,  1/9.81, 1/9.81, 1/9.81,        1,       1,       1,        1,       1,       1,       1,       1, 24*3600/100,       1,       1,     1,     1, -1/3600, -1/3600,  1/3600,  1/3600, 1/3600, -1/3600, -1/3600,      1,     1,     1];                
% datatype_3D          = {'ins',   'ins',  'acc', 'acc',  'ins',     'ins',     'ins',  'ins',   'ins',  'ins',  'ins',    'ins',   'ins',   'ins',    'ins',   'ins',   'ins',   'ins',   'ins',       'ins',   'ins',   'ins', 'ins', 'ins',   'ins',   'ins',   'ins',   'ins',  'ins',   'ins',   'ins',  'ins', 'ins', 'ins'};    
% whichfolder_3D       = [    2,       2,      2,     2,      2,         2,         2,      1,       1,      1,      1,        1,       2,       2,        1,       1,       1,       3,       3,           1,       2,       2,     2,     3,       2,       2,       2,       2,      2,       2,       2,      2,     2,     2];
% StepsInADay_3D       = [   24,      24,     24,    24,     24,        24,        24,     24,      24,     24,     24,       24,      24,      24,       24,      24,      24,       8,       8,          24,      24,      24,    24,     8,      24,      24,      24,      24,     24,      24,      24,     24,    24,    24];
% 
% %/ a set of 4D data (appending)/%
% varname_4D           =        {  'q',   't',     'z',   'u',   'v',    'w',   'd',};
% dataname_callfile_4D =        {  'q',   'T',     'Z',   'U',   'V',    'W', 'div',};
% dataname_4D          =        {  'q',   'T',     'Z',   'U',   'V',    'W', 'div',};
% dataunitconv_4D      =        [    1,     1,  1/9.81,     1,     1,      1,     1,];
% datatype_4D          = repmat({'ins'}, 1, length(varname_4D));
% whichfolder_4D       =        [   7,     7,        7,     7,     7,      7,     8,];
% StepsInADay_4D       =        [   4,     4,        4,     4,     4,      4,     4,];
% 
% %/ pre-processed data
% varname_preprocc          = {'EPT', 'dEPTdy'};
% dataname_callfile_preproc = varname_preprocc;
% dataname_preproc          = varname_preprocc;
% dataunitconv_preproc      = ones(1, length(varname_preprocc));
% datatype_preproc          = repmat({'ins'}, 1, length(varname_preprocc));
% whichfolder_preproc       = ones(1, length(varname_preprocc))*-1;
% StepsInADay_preproc       = repmat(24, 1, length(varname_preprocc));
% 
% varname              = cat(2, varname_nonrean,           varname_3D,           varname_4D,           varname_preprocc         );
% dataname_callfile    = cat(2, dataname_callfile_nonrean, dataname_callfile_3D, dataname_callfile_4D, dataname_callfile_preproc);
% dataname             = cat(2, dataname_nonrean,          dataname_3D,          dataname_4D,          dataname_preproc         );
% dataunitconv         = cat(2, dataunitconv_nonrean,      dataunitconv_3D,      dataunitconv_4D,      dataunitconv_preproc     );
% datatype             = cat(2, datatype_nonrean,          datatype_3D,          datatype_4D,          datatype_preproc         );
% whichfolder          = cat(2, whichfolder_nonrean,       whichfolder_3D,       whichfolder_4D,       whichfolder_preproc      );
% StepsInADay          = cat(2, StepsInADay_nonrean,       StepsInADay_3D,       StepsInADay_4D,       StepsInADay_preproc      );
% 
% %/ double check if size is consistent.
% a = unique([length(varname), length(dataname_callfile), length(dataname), length(dataunitconv), length(datatype), length(whichfolder), length(StepsInADay)]);
% 
% if length(a) ~= 1
%    error('inconsistent size of dataname and others.'); 
% end
% 
% dataset.placeholder = [];
% 
% % filename = '/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/ERA5_div_6hrly_mth01-12_1x1_40_180E_-30_50N_lv100-1000_2013.nc';
% % filename = '/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/ERA5_W_6hrly_mth01-12_1x1_0_360E_-90_90N_lv100-1000_2013.nc';
% % filename = '/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/ERA5_SSHF_hrly_mth01-12_1x1_0_360E_-90_90N_2013.nc';
% % filename = '/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/ERA5_SLHF_hrly_mth01-12_1x1_0_360E_-90_90N_2013.nc';
% % filename = '/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/ERA5_VIMD_hrly_mth01-12_1x1_0_360E_-90_90N_2013.nc';
% filename = '/disk/r128/tfchengac/NOAA_CDR_data/www.ncei.noaa.gov/data/outgoing-longwave-radiation-daily/access/olr-daily_v01r02_20200101_20201231.nc';
% ncdisp(filename)
% % units = ncreadatt(filename,"vimd","units")
% 
% % ncread(filename, 'level')
% fprintf('*** paramter module (ERA5) has been loaded. ***\n')




















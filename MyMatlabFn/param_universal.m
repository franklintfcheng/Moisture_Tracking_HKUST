%% Universal module
opengl software  %/ This may help solve the freezing problem of graphics
dataset.placeholder = [];

str_mth = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec',....
           'MAM', 'JJA', 'SON', 'DJF', 'AMJJAS', 'ONDJFM', 'JFD', 'nonJJA', 'MJJASO', 'NDJFMA', 'MJJAS'};
stDate          = datetime(1980,   1,  1,    'Format','yyyy-MM-dd');
edDate          = datetime(1980,   12, 31,   'Format','yyyy-MM-dd');
date_OneYr      = (stDate:edDate)';
month_OneYr     = date_OneYr.Month;
day_OneYr       = date_OneYr.Day;
date_mmdd_OneYr = date_OneYr.Month*100 + date_OneYr.Day;

data_path_parts = [
   {'/disk/r128/tfchengac/fandyr128SOM/ERA5-Land_data_Annual/'}, {'ERA5-Land_'}, {'_24hrly_mth01-12_1x1_0_360E_-90_90N_'};
   {'/disk/v183.b/share/ERA5_data_Annual/'},                     {'ERA5_'},      {'_6hrly_mth01-12_1x1_0_360E_-90_90N_lv100-1000_'};
   {'/disk/v183.b/share/ERA5_data_Annual/'},                     {'ERA5_'},      {'_hrly_mth01-12_1x1_0_360E_-90_90N_'};         %/ 3 
   {'/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/'},   {'cera20c_'},   {'_step3_6_9_12_15_18_21_24_mth01-12_1x1_0_360E_-90_90N_'};...
   {'/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/'},   {'cera20c_'},   {'_3hrly_mth01-12_1x1_0_360E_-90_90N_'}; 
   {'/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/'},   {'cera20c_'},   {'_mthly_mth01-12_1x1_0_360E_-90_90N_'};        %/ 6 
   {'/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/'},   {'cera20c_'},   {'_step6_24_mth01-12_1x1_0_360E_-90_90N_'}; 
   {'/disk/r128/tfchengac/fandyr128SOM/Climate Indices/'},       {''},           {''}; 
   {'/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/'},   {'cera20c_'},   {'_6hrly_mth01-12_1x1_0_360E_-90_90N_'};        %/ 9 
   {'/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/'},   {'cera20c_'},   {''};
   {'/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/'},   {'cera20c_'},   {'_3hrly_mth01-12_1x1_20_150E_0_50N_'};         %/ 11
   {'/disk/v183.b/share/ERA5_data_Annual/'},                     {'ERA5_'},      {'_hrly_mth01-12_0.5x0.5_0_360E_-90_90N_'};];         

%/ CMIP6 data
cmip_varname           = ['pr', 'evspsbl', 'tas', 'rlut', 'tos', 'psl', 'od550aer',...
                          'evspsbl_pr_ratio', 'S400', 'SS', 'SS2', 'hus500wap300', 'uIVT', 'vIVT'...
                           strcat({'ua', 'va', 'wap', 'ta', 'hus', 'div'}, {'850'}),...
                           strcat({'ua', 'va', 'wap', 'ta', 'hus', 'div'}, {'700'}),...
                           strcat({'ua', 'va', 'wap', 'ta', 'hus', 'div'}, {'500'}),...
                           strcat({'ua', 'va', 'wap', 'ta', 'hus', 'div'}, {'300'}),...
                           strcat({'ua', 'va', 'wap', 'ta', 'hus', 'div'}, {'250'})];
cmip_dataname_callfile = cmip_varname;
cmip_dataname          = cmip_varname;
cmip_dataunitconv      = ones(1, length(cmip_varname));
cmip_whichfolder       = ones(1, length(cmip_varname));
cmip_StepsInADay       = ones(1, length(cmip_varname));
cmip_datatype          = repmat({''}, 1, length(cmip_varname));
cmip_dataunit          = repmat({''}, 1, length(cmip_varname));  %/ Update later
cmip_fc_steps          = cell(1, length(cmip_varname));
% cmip_dataname          = strcat('cmip6_', cmip_dataname); %/ Do not make things complicated!

%/ Gauge/Satellite/Model-based/Ensemble/Processed data
nonrean_varname           = {'NOAA_Interp_OLR', 'NOAA_CDR_OLR', 'EM_P', 'EM_E', 'CRU_P', 'GPCP_P', 'TPR_P', 'TPR_E', 'TPR_SRO', 'TPR_SSRO', 'TPR_RO', 'HadCRUT5_T2m', 'GPCC_P', 'HARv2_P', 'HARv2_E', 'HARv2_T2m', 'GLEAM_SMsurf', 'GLEAM_E', 'GLASS_E', 'CMORPH_P', 'IMERG_P', 'traj_freq', 'total_init_contr', 'avg_init_contr', 'contr_map', 'Pm', 'BL_Pm', 'LAI', 'PDSI', 'NDVI', 'E_minus_P', 'EmPSST', 'P_minus_E_gte0', 'E_minus_P_gte0', 'waterbudget_residual', 'prop_crop', 'prop_pasture', 'prop_urbn', 'prop_primary',  'prop_secd', 'biomass_harv_prim',         'biomass_secd', 'prop_vc_transition', 'prop_vp_transition', 'prop_vu_transition'}; 
nonrean_dataname_callfile = {'NOAA_Interp_OLR', 'NOAA_CDR_OLR', 'EM_P', 'EM_E', 'CRU_P', 'GPCP_P', 'TPR_P', 'TPR_E', 'TPR_SRO', 'TPR_SSRO', 'TPR_RO', 'HadCRUT5_T2m', 'GPCC_P', 'HARv2_P', 'HARv2_E', 'HARv2_T2m', 'GLEAM_SMsurf', 'GLEAM_E', 'GLASS_E', 'CMORPH_P', 'IMERG_P', 'traj_freq', 'total_init_contr', 'avg_init_contr', 'contr_map', 'Pm', 'BL_Pm', 'LAI', 'PDSI', 'NDVI', 'E_minus_P', 'EmPSST', 'P_minus_E_gte0', 'E_minus_P_gte0', 'waterbudget_residual',     'gcrop',        'gpast',     'gurbn',        'gothr',      'gsecd',             'gvbh1',                'gssmb',              'gflvc',              'gflvp',              'gflvu'};
nonrean_dataname          = {'NOAA_Interp_OLR', 'NOAA_CDR_OLR', 'EM_P', 'EM_E', 'CRU_P', 'GPCP_P', 'TPR_P', 'TPR_E', 'TPR_SRO', 'TPR_SSRO', 'TPR_RO', 'HadCRUT5_T2m', 'GPCC_P', 'HARv2_P', 'HARv2_E', 'HARv2_T2m', 'GLEAM_SMsurf', 'GLEAM_E', 'GLASS_E', 'CMORPH_P', 'IMERG_P', 'traj_freq', 'total_init_contr', 'avg_init_contr', 'contr_map', 'Pm', 'BL_Pm', 'LAI', 'PDSI', 'NDVI', 'E_minus_P', 'EmPSST', 'P_minus_E_gte0', 'E_minus_P_gte0', 'waterbudget_residual',      'crop',      'pasture',     'urban',      'primary',  'secondary', 'biomass_harv_prim', 'biomass_density_secd',       'primary2crop',    'primary2pasture',      'primary2urban'};
nonrean_dataunitconv      = ones(1, length(nonrean_varname));
nonrean_whichfolder       = ones(1, length(nonrean_varname));
nonrean_StepsInADay       = ones(1, length(nonrean_varname));
nonrean_datatype          = repmat({''}, 1, length(nonrean_varname));
nonrean_dataunit          = repmat({''}, 1, length(nonrean_varname));  %/ Update later
nonrean_fc_steps          = cell(1, length(nonrean_varname));

%/ Global ERA5 data (hourly)
ERA5_varname           = {'SxW',   'S',   't',     'q',    'w',  'sp', 'vimd',   'tp',   'tp',    'e', 'tcwv', 'p71.162', 'p72.162',  'msl',     'z',    'z',    'z', 'Z200sw',     'u',     'v', 'u200sw',     'u',     'v',     'u',     'v',    'w',   'u10',   'v10', 't2m', 'sst',    't',    't',    't',    't',   'ttr',   'str',   'tsr', 'tisr',  'ssr',  'slhf',  'sshf', 'cape', 'cin', 'blh', 'S500',    'q',     'q',  'S2_Wang1988', 'C0_Wang1988', 'alpha_Wang1988', 'q3_bar_Wang1988', 'q1_bar_Wang1988', 'dq_bar_Wang1988', 'I_Wang1988'};
ERA5_dataname_callfile = {'SxW',   'S',   'T',     'q',    'W',  'sp', 'VIMD', 'prcp', 'prcp',    'E', 'tcwv',     'IVT',     'IVT',  'slp',  'Z850', 'Z500', 'Z200', 'Z200sw', 'UV200', 'UV200', 'u200sw', 'UV500', 'UV500', 'UV850', 'UV850', 'W500', 'UV10m', 'UV10m', 'T2m', 'SST', 'T250', 'T300', 'T500', 'T700',   'TTR',   'STR',   'TSR', 'TISR',  'SSR',  'SLHF',  'SSHF', 'CAPE', 'CIN', 'BLH', 'S500', 'q850', 'q1000',  'S2_Wang1988', 'C0_Wang1988', 'alpha_Wang1988', 'q3_bar_Wang1988', 'q1_bar_Wang1988', 'dq_bar_Wang1988', 'I_Wang1988'};
ERA5_dataname          = {'SxW',   'S',   'T',     'q',    'W',  'sp', 'VIMD',    'P', 'P_05',    'E',   'pw',    'uIVT',    'vIVT',  'slp',  'Z850', 'Z500', 'Z200', 'Z200sw',  'u200',  'v200', 'u200sw',  'u500',  'v500',  'u850',  'v850', 'W500',  'u10m',  'v10m', 'T2m', 'SST', 'T250', 'T300', 'T500', 'T700',   'OLR',   'STR',   'TSR', 'TISR',  'SSR',  'SLHF',  'SSHF', 'CAPE', 'CIN', 'BLH', 'S500', 'q850', 'q1000',  'S2_Wang1988', 'C0_Wang1988', 'alpha_Wang1988', 'q3_bar_Wang1988', 'q1_bar_Wang1988', 'dq_bar_Wang1988', 'I_Wang1988'};
ERA5_dataunitconv      = [    1,     1,     1,       1,      1,     1,      1,   1000,   1000,  -1000,      1,         1,         1,  1/100,  1/9.81, 1/9.81, 1/9.81,        1,       1,       1,        1,       1,       1,       1,       1,      1,       1,       1,     1,     1,       1,     1,      1,      1, -1/3600, -1/3600,  1/3600, 1/3600, 1/3600, -1/3600, -1/3600,      1,     1,     1,      1,      1,       1,              1,             1,                1,                 1,                 1,                 1,            1];                
ERA5_whichfolder       = [    2,     2,     2,       2,      2,     3,      3,      3,     12,      3,      3,         3,         3,      3,       3,      3,      3,        3,       3,       3,        3,       3,       3,       3,       3,      3,       3,       3,     3,     3,       3,     3,      3,      3,       3,       3,       3,      3,      3,       3,       3,      3,     3,     3,      3,      3,       3,              3,             3,                3,                 3,                 3,                 3,            3];
ERA5_StepsInADay       = [    4,     4,     4,       4,      4,    24,     24,     24,     24,     24,     24,        24,        24,     24,      24,     24,     24,       24,      24,      24,       24,      24,      24,      24,      24,     24,      24,      24,    24,    24,      24,   nan,     24,     24,      24,      24,      24,      24,    24,      24,      24,     24,    24,    24,     24,     24,     nan,            nan,           nan,              nan,               nan,               nan,               nan,          nan];
ERA5_datatype          = {'ins', 'ins', 'ins',   'ins',  'ins', 'ins',  'ins',  'acc',  'acc',  'acc',  'ins',     'ins',     'ins',  'ins',   'ins',  'ins',  'ins',    'ins',   'ins',   'ins',    'ins',   'ins',   'ins',   'ins',   'ins',  'ins',   'ins',   'ins', 'ins', 'ins',   'ins', 'ins',  'ins',  'ins',   'ins',   'ins',   'ins',   'ins', 'ins',   'ins',   'ins',  'ins', 'ins', 'ins',  'ins',  'ins',   'ins',          'ins',         'ins',            'ins',             'ins',             'ins',             'ins',        'ins'};    
ERA5_dataunit          = {'K/s','K/Pa',   'K', 'kg/kg', 'Pa/s',  'Pa',  's-1', 'mm/d', 'mm/d', 'mm/d',   'mm',  'kg/m/s',  'kg/m/s',  'hPa',     'm',    'm',    'm',      'm',   'm/s',   'm/s',    'm/s',   'm/s',   'm/s',   'm/s',   'm/s', 'Pa/s',   'm/s',   'm/s',   'K',   'K',     'K',   'K',    'K',    'K',   'W/m2', 'W/m2',  'W/m2',  'W/m2','W/m2',  'W/m2',  'W/m2', 'J/kg','J/kg',   'm', 'K/Pa','kg/kg', 'kg/kg',  'm2 s-2 Pa-2',         'm/s',               '',           'kg/kg',           'kg/kg',           'kg/kg',           ''};
ERA5_fc_steps          = cell(1,     length(ERA5_varname));
ERA5_dataname          = strcat('ERA5_', ERA5_dataname);  %/ to distinguish them.

%/ Global CERA-20C data (ins 3hourly or fc monthly) 
CERA_an_varname           = {'sp',          'r',          'r', 'blh',    'z',    'z',     'u',      'v',     'u',      'v',     'u',      'v',          'w',          'w',       'div500', 'q500Omega300', 'dTWS',        'q',         'u',         'v',        'w',        'q',             'uIVT',             'vIVT',             'IVTdiv',             'VMT',             'dwdt', 't2m', 'sst', 'T2mSST', 'cape',     'swvl1',         'sd',        'lai_lv',         'lai_hv',    'tvl',   'tvh'};
CERA_an_dataname_callfile = {'sp', 'RH300-1000', 'RH300-1000', 'BLH', 'Z850', 'Z500', 'UV850',  'UV850', 'UV500',  'UV500', 'UV300',  'UV300',   'Omega500',   'Omega300',       'div500', 'q500Omega300', 'dTWS', 'q10-1000',  'U10-1000',  'V10-1000', 'W10-1000', 'q10-1000', 'CERA_uIVT10_1000', 'CERA_vIVT10_1000', 'CERA_IVTdiv10_1000', 'CERA_VMT10_1000', 'CERA_dwdt10_1000', 'T2m', 'SST', 'T2mSST', 'CAPE', 'SM_layer1',  'SnowDepth','LAI_lowhighveg', 'LAI_lowhighveg',  'LoVeg', 'HiVeg'};
CERA_an_dataname          = {'sp',      'RH850',      'RH500', 'BLH', 'Z850', 'Z500',  'u850',   'v850',  'u500',   'v500',  'u300',   'v300',   'Omega500',   'Omega300',       'div500', 'q500Omega300', 'dTWS',     'q500',         'U',         'V',        'W',        'q',             'uIVT',             'vIVT',             'IVTdiv',             'VMT',             'dwdt', 'T2m', 'SST', 'T2mSST', 'CAPE', 'SM_layer1',  'SnowDepth',        'LAI_lv',         'LAI_hv',  'LoVeg', 'HiVeg'};
CERA_an_dataunitconv      = [   1,            1,            1,     1, 1/9.81, 1/9.81,       1,        1,       1,        1,       1,        1,  24*3600/100,  24*3600/100,              1,              1,      1,          1,           1,           1,          1,          1,                  1,                  1,                    1,                 1,                  1,      1,     1,        1,      1,        1000,            1,               1,                1,        1,      1];
CERA_an_whichfolder       = [   5,            6,            6,     6,      5,      9,       5,        5,       5,        5,       5,        5,            5,            5,              1,              1,      1,         11,          11,          11,         11,         11,                nan,                nan,                  nan,               nan,                  5,      5,     5,      nan,      5,           5,            5,               5,                5,       10,     10];
CERA_an_StepsInADay       = [   8,          nan,          nan,   nan,      8,      4,       8,        8,       8,        8,       8,        8,            8,            8,              1,              1,      1,          8,           8,           8,          8,          8,                nan,                nan,                  nan,               nan,                  8,      8,     8,      nan,      8,           8,            8,               8,                8,      nan,    nan];
CERA_an_datatype          = repmat({'ins'}, 1, length(CERA_an_varname)); % instantaneous (mean) or accumulated (sum)   
CERA_an_dataunit          = repmat({''}, 1, length(CERA_an_varname));  %/ Update later
CERA_an_fc_steps          = cell(1, length(CERA_an_varname));
CERA_an_dataname          = strcat('CERA_', CERA_an_dataname);  %/ to distinguish them.

%/ Global CERA-20C data (fc, daily) - 2 time steps. (preferred)
CERA_fc_varname           = {  'tp',        'e',      'ro',       'sro',         'ssro',      'smlt',   'es'};
CERA_fc_dataname_callfile = { 'E_P',      'E_P',  'Runoff', 'SfcRunoff', 'SubsfcRunoff',  'SnowMelt',   'Es'};
CERA_fc_dataname          = {   'P',        'E',      'RO',       'SRO',         'SSRO',      'SMLT',   'Es'};
CERA_fc_dataunitconv      = [  1000,      -1000,      1000,        1000,           1000,        1000,  -1000];
CERA_fc_whichfolder       = repmat(7,                  1, length(CERA_fc_varname));
CERA_fc_StepsInADay       = repmat(2, 1, length(CERA_fc_varname));
CERA_fc_datatype          = repmat({'acc'},            1, length(CERA_fc_varname)); % instantaneous (mean) or accumulated (sum)
CERA_fc_dataunit          = repmat({''}, 1, length(CERA_fc_varname));  %/ Update later
CERA_fc_fc_steps          = repmat({[6, 24]},          1, length(CERA_fc_varname)); %/ assume the same steps for all fc vars.
CERA_fc_dataname          = strcat('CERA_', CERA_fc_dataname);  %/ to distinguish them.

%/ Global ERA5-land data (hourly)
ERA5Land_varname           = {'evabs', 'evavt', 'evaow', 'evatc',    'es',    'e'};
ERA5Land_dataname_callfile = {'Esoil', 'Etran', 'Eopwt', 'Ecano', 'Esnow', 'Etot'};  %/ Since they messed up with the varname.
ERA5Land_dataname          = {'Etran', 'Eopwt', 'Esoil', 'Ecano', 'Esnow', 'Etot'};
ERA5Land_dataunitconv      = [  -1000,   -1000,   -1000,   -1000,   -1000,  -1000];
ERA5Land_whichfolder       = ones(1, length(ERA5Land_varname));
ERA5Land_StepsInADay       = ones(1, length(ERA5Land_varname));
ERA5Land_datatype          = repmat({'acc'},            1, length(ERA5Land_varname)); % instantaneous (mean) or accumulated (sum)
ERA5Land_dataunit          = repmat({''}, 1, length(ERA5Land_varname));  %/ Update later
ERA5Land_fc_steps          = cell(1, length(ERA5Land_varname));
ERA5Land_dataname          = strcat('ERA5Land_', ERA5Land_dataname);  %/ to distinguish them.

varname           = cat(2,  CERA_fc_varname,           nonrean_varname,           cmip_varname,           CERA_an_varname,           ERA5Land_varname,           ERA5_varname);
dataname_callfile = cat(2,  CERA_fc_dataname_callfile, nonrean_dataname_callfile, cmip_dataname_callfile, CERA_an_dataname_callfile, ERA5Land_dataname_callfile, ERA5_dataname_callfile);
dataname          = unique(cat(2,  CERA_fc_dataname,   nonrean_dataname,          cmip_dataname,          CERA_an_dataname,          ERA5Land_dataname,          ERA5_dataname), 'stable');  %/ Only dataname should be unique, use 'stable' to avoid auto sorting!
dataunitconv      = cat(2,  CERA_fc_dataunitconv,      nonrean_dataunitconv,      cmip_dataunitconv,      CERA_an_dataunitconv,      ERA5Land_dataunitconv,      ERA5_dataunitconv);
whichfolder       = cat(2,  CERA_fc_whichfolder,       nonrean_whichfolder,       cmip_whichfolder,       CERA_an_whichfolder,       ERA5Land_whichfolder,       ERA5_whichfolder);
StepsInADay       = cat(2,  CERA_fc_StepsInADay,       nonrean_StepsInADay,       cmip_StepsInADay,       CERA_an_StepsInADay,       ERA5Land_StepsInADay,       ERA5_StepsInADay);
datatype          = cat(2,  CERA_fc_datatype,          nonrean_datatype,          cmip_datatype,          CERA_an_datatype,          ERA5Land_datatype,          ERA5_datatype);
dataunit          = cat(2,  CERA_fc_dataunit,          nonrean_dataunit,          cmip_dataunit,          CERA_an_dataunit,          ERA5Land_dataunit,          ERA5_dataunit);
fc_steps          = [       CERA_fc_fc_steps,          nonrean_fc_steps,          cmip_fc_steps,          CERA_an_fc_steps,          ERA5Land_fc_steps,          ERA5_fc_steps];

%/ Double check!!
if ~isequal(length(dataname), length(dataname_callfile))
    error('Inconsistent length between dataname (%d) and dataname_callfile (%d)!', length(dataname), length(dataname_callfile))
elseif ~isequal(length(dataname), length(varname))
    error('Inconsistent length between dataname (%d) and varname (%d)!', length(dataname), length(varname))
elseif ~isequal(length(dataname), length(dataunitconv))
    error('Inconsistent length between dataname (%d) and dataunitconv (%d)!', length(dataname), length(dataunitconv))
elseif ~isequal(length(dataname), length(datatype))
    error('Inconsistent length between dataname (%d) and datatype (%d)!', length(dataname), length(datatype))
elseif ~isequal(length(dataname), length(whichfolder))
    error('Inconsistent length between dataname (%d) and whichfolder (%d)!', length(dataname), length(whichfolder))
elseif ~isequal(length(dataname), length(dataunit))
    error('Inconsistent length between dataname (%d) and dataunit (%d)!', length(dataname), length(dataunit))
elseif ~isequal(length(dataname), length(fc_steps))
    error('Inconsistent length between dataname (%d) and fc_steps (%d)!', length(dataname), length(fc_steps))
elseif ~isequal(length(dataname), length(StepsInADay))
    error('Inconsistent length between dataname (%d) and StepsInADay (%d)!', length(dataname), length(StepsInADay))
end

% ind = find(ismember(ERA5_dataname, {'ERA5_P'}))
% ERA5_dataunitconv(ind)

% ind = find(ismember(dataname, {'ERA5_P'}))
% dataunitconv(ind)

% ncdisp('/disk/r128/tfchengac/GPCP_data/GPCPMON_L3_201912_V3.1.nc4')
% ncdisp('/disk/r128/tfchengac/HadCRUT5_Analysis_data/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc')
% ncdisp('/disk/r128/tfchengac/CRU_data/cru_ts4.07.1901.2022.pre.dat.nc')
% ncdisp('/disk/r128/tfchengac/TPR_data/biggeo.gvc.gu.se/TPReanalysis/TP9km/Daily/SubsurfaceRunoff/WRFOut_TP9km_Daily_Mean_SubsurfaceRunoff_2019.nc')
% ncdisp('/disk/r128/tfchengac/TPR_data/biggeo.gvc.gu.se/TPReanalysis/TP9km/Daily/SurfaceRunoff/WRFOut_TP9km_Daily_Mean_SurfaceRunoff_2018.nc')
% ncdisp('/disk/r128/tfchengac/TPR_data/biggeo.gvc.gu.se/TPReanalysis/TP9km/Hourly/SFCEVP/WRFOut_TP9km_HourlySFCEVP_2019_12.nc')
% ncdisp('/disk/r128/tfchengac/TPR_data/biggeo.gvc.gu.se/TPReanalysis/TP9km/Hourly/Precip/WRFOut_TP9km_HourlyP_1999_01.nc')
% ncdisp('/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/cera20c_E_P_step6_24_mth01-12_1x1_0_360E_-90_90N_2005.nc')
% ncdisp('/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/cera20c_sp_3hrly_mth01-12_1x1_0_360E_-90_90N_2004.nc')
% ncdisp('/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/cera20c_q10-1000_3hrly_mth01-12_1x1_20_150E_0_50N_2004.nc')

fprintf('*** param_universal.m loaded ***\n');


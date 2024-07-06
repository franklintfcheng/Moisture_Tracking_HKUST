%% [Step 1] ERA5/CERA data (module)
addpath(genpath('/home/tfchengac/MyMatlabPkg'));
addpath(genpath('/home/tfchengac/MyMatlabFn'));
addpath(genpath('/home/tfchengac/Flexpart_Codes'));
opengl('save','software'); %/ to solve slow figure issue (not very effective though).
masterfolder = '/disk/r059/tfchengac/FLEXPART/';

stDate = datetime(1980,   1,  1,    'Format','yyyy-MM-dd');
edDate = datetime(1980,   12, 31,   'Format','yyyy-MM-dd');
date_OneYr = [stDate:edDate]';
month_OneYr = date_OneYr.Month;
day_OneYr = date_OneYr.Day;
date_mmdd_OneYr = [date_OneYr.Month*100 + date_OneYr.Day];

fc_steps = [];
data_path_parts = [
   {'/disk/r128/tfchengac/fandyr128SOM/ERA5-Land_data_Annual/'}, {'ERA5-Land_'}, {'_24hrly_mth01-12_1x1_0_360E_-90_90N_'};...
   {'/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/'},      {'ERA5_'},      {'_3hrly_mth01-12_1x1_0_360E_-90_90N_'};...
   {'/disk/r128/tfchengac/fandyr128SOM/ERA5_data_Annual/'},      {'ERA5_'},      {'_hrly_mth01-12_1x1_0_360E_-90_90N_'};...
   {'/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/'},   {'cera20c_'},   {'_step3_6_9_12_15_18_21_24_mth01-12_1x1_0_360E_-90_90N_'};...
   {'/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/'},   {'cera20c_'},   {'_3hrly_mth01-12_1x1_0_360E_-90_90N_'};...
   {'/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/'},   {'cera20c_'},   {'_mthly_mth01-12_1x1_0_360E_-90_90N_'};...
   {'/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/'},   {'cera20c_'},   {'_step6_24_mth01-12_1x1_0_360E_-90_90N_'};...
   {'/disk/r128/tfchengac/fandyr128SOM/Climate Indices/'},       {''},           {''}];

%/ Global CERA-20C data (fc, daily) - for 2010 only, 8 time steps.
% varname_0           = {  'tp',        'e'};
% dataname_callfile_0 = { 'E_P',        'E_P'};
% dataname_0          = { 'P_d',        'E_d'};
% dataunitconv_0      = [  1000,     -1000];
% datatype_0          = repmat({'acc'}, 1, length(varname_0)); % instantaneous (mean) or accumulated (sum)
% whichfolder_0       = repmat(4,  1, length(varname_0));
% fc_steps            = [3:3:24]; %/ assume the same steps for all fc vars.
% StepsInADay_0       = repmat(length(fc_steps), 1, length(varname_0));

%/ Global CERA-20C data (fc, daily) - for all years, 2 time steps.
varname_0           = {  'tp',        'e'};
dataname_callfile_0 = { 'E_P',      'E_P'};
dataname_0          = {   'P',        'E'};
dataunitconv_0      = [  1000,      -1000];
datatype_0          = repmat({'acc'}, 1, length(varname_0)); % instantaneous (mean) or accumulated (sum)
whichfolder_0       = repmat(7,  1, length(varname_0));
fc_steps            = [6, 24]; %/ assume the same steps for all fc vars.
StepsInADay_0       = repmat(length(fc_steps), 1, length(varname_0));

%/ Non-reanalysis data
varname_1           = {'traj_freq', 'total_init_contr', 'avg_init_contr', 'contr_map', 'Pm', 'BL_Pm', 'LAI', 'PDSI', 'NDVI', 'E_minus_P', 'EmPSST', 'P_minus_E_gte0', 'E_minus_P_gte0', 'waterbudget_residual', 'prop_crop', 'prop_pasture', 'prop_urbn', 'prop_primary',  'prop_secd', 'biomass_harv_prim', 'biomass_secd',        'prop_vc_transition', 'prop_vp_transition', 'prop_vu_transition'}; 
dataname_callfile_1 = {'traj_freq', 'total_init_contr', 'avg_init_contr', 'contr_map', 'Pm', 'BL_Pm', 'LAI', 'PDSI', 'NDVI', 'E_minus_P', 'EmPSST', 'P_minus_E_gte0', 'E_minus_P_gte0', 'waterbudget_residual',     'gcrop',        'gpast',     'gurbn',        'gothr',      'gsecd',             'gvbh1', 'gssmb',                            'gflvc',              'gflvp',              'gflvu'};
dataname_1          = {'traj_freq', 'total_init_contr', 'avg_init_contr', 'contr_map', 'Pm', 'BL_Pm', 'LAI', 'PDSI', 'NDVI', 'E_minus_P', 'EmPSST', 'P_minus_E_gte0', 'E_minus_P_gte0', 'waterbudget_residual',      'crop',      'pasture',     'urban',      'primary',  'secondary', 'biomass_harv_prim', 'biomass_density_secd',      'primary2crop',    'primary2pasture',      'primary2urban'};
dataunitconv_1      = repmat(1,  1, length(varname_1), 1);
datatype_1          = repmat({''}, 1, length(varname_1), 1);
whichfolder_1       = repmat(1,  1, length(varname_1), 1);
StepsInADay_1       = repmat(1,  1, length(varname_1), 1);

%/ Global CERA-20C data (ins 3hourly or fc monthly)
varname_2           = {'HadISST', 'ersstv5',          'r',          'r', 'blh',    'z',     'u',      'v',         'u',         'v',         'q',    'uIVT',  'vIVT',  'IVTdiv',       'u',        'v',       'q', 'uIVT1-1000', 'vIVT1-1000', 'IVTdiv1-1000', 't2m', 'sst', 'T2mSST',     'swvl1',         'lai_lv',         'lai_hv'};
dataname_callfile_2 = {'HadISST', 'ersstv5', 'RH300-1000', 'RH300-1000', 'BLH', 'Z850', 'UV850',  'UV850', 'U300-1000', 'V300-1000', 'q300-1000',    'uIVT',  'vIVT',  'IVTdiv', 'U1-1000',  'V1-1000', 'q1-1000', 'uIVT1-1000', 'vIVT1-1000', 'IVTdiv1-1000', 'T2m', 'SST', 'T2mSST', 'SM_layer1', 'LAI_lowhighveg', 'LAI_lowhighveg'};
dataname_2          = {'HadISST', 'ersstv5',      'RH850',      'RH500', 'BLH', 'Z850',  'u850',   'v850', 'U300_1000', 'V300_1000', 'q300_1000',    'uIVT',  'vIVT',  'IVTdiv', 'U1_1000',  'V1_1000', 'q1_1000', 'uIVT1-1000', 'vIVT1-1000', 'IVTdiv1-1000', 'T2m', 'SST', 'T2mSST', 'SM_layer1',         'LAI_lv',         'LAI_hv'};
dataunitconv_2      = [        1,         1,            1,            1,     1, 1/9.81,       1,        1,           1,           1,           1,         1,       1,         1,          1,         1,         1,            1,            1,              1,     1,     1,        1,           1,                1,                1];
datatype_2          = repmat({'ins'}, 1, length(varname_2)); % instantaneous (mean) or accumulated (sum)
whichfolder_2       = [        8,         8,            6,            6,     6,      5,       5,        5,           5,           5,           5,         5,       5,         5,          5,         5,         5,            5,            5,              5,     6,     6,        6,           6,                6,                6];
StepsInADay_2       = repmat(8, 1, length(varname_2));

%/ Global ERA5-land data (hourly)
varname_3           = {'evabs', 'evavt', 'evaow', 'evatc',    'es',    'e'};
dataname_callfile_3 = {'Esoil', 'Etran', 'Eopwt', 'Ecano', 'Esnow', 'Etot'};  %/ Since they messed up with the varname.
dataname_3          = {'Etran', 'Eopwt', 'Esoil', 'Ecano', 'Esnow', 'Etot'};
dataunitconv_3      = [  -1000,   -1000,   -1000,   -1000,   -1000,  -1000];
datatype_3          = {  'acc',   'acc',   'acc',   'acc',   'acc',  'acc'};
whichfolder_3       = repmat(1,  1, length(varname_3));
StepsInADay_3       = repmat(1,  1, length(varname_3));

varname           = cat(2, varname_0,           varname_1,           varname_2,           varname_3);
dataname_callfile = cat(2, dataname_callfile_0, dataname_callfile_1, dataname_callfile_2, dataname_callfile_3);
dataname          = cat(2, dataname_0,          dataname_1,          dataname_2,          dataname_3);
dataunitconv      = cat(2, dataunitconv_0,      dataunitconv_1,      dataunitconv_2,      dataunitconv_3);
datatype          = cat(2, datatype_0,          datatype_1,          datatype_2,          datatype_3);
whichfolder       = cat(2, whichfolder_0,       whichfolder_1,       whichfolder_2,       whichfolder_3);
StepsInADay       = cat(2, StepsInADay_0,       StepsInADay_1,       StepsInADay_2,       StepsInADay_3);

%NOTE:
%   - Every time step in ERA5 is reanalysis for both ins and acc
%   - Time step 0 in ERA5-Land is the daily acc value

%/ ERA5-Land daily 1x1 data
%/IMPORTANT: Esoil, Etran and Eopwt are swapped. Check the known issues. https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation
% varname_1st           = {'evabs', 'evavt', 'evaow', 'evatc',    'es',    'e'};
% dataname_callfile_1st = {'Esoil', 'Etran', 'Eopwt', 'Ecano', 'Esnow', 'Etot'};
% dataname_1st          = {'Etran', 'Eopwt', 'Esoil', 'Ecano', 'Esnow', 'Etot'};  
% dataunitconv_1st      = [  -1000,   -1000,   -1000,   -1000,   -1000,  -1000];
% datatype_1st          = {  'acc',   'acc',   'acc',   'acc',   'acc',  'acc'};
% whichfolder_1st       = repmat(1,  1, length(varname_1st));
% StepsInADay_1st       = repmat(1,  1, length(varname_1st));

%/ Global ERA5 *3hrly* 1x1 data (ins)
% varname_2nd           = {'z',    'z', 'sst',   'u',     'v',        'u',     'v',      [],      [],   [],  [],  [],     [],       [],       [],       []};
% dataname_callfile_2nd = {'Z850', 'Z500', 'SST', 'UV200', 'UV200', 'UV850', 'UV850',    [],      [],   [],  [],  [],     [],       [],       [],       []};
% dataname_2nd          = {'Z850', 'Z500', 'SST', 'u200',  'v200',   'u850',  'v850',   'T2mSST','S','S1','S2','uWAF','vWAF','psi200','u200chi','v200chi'};
% dataunitconv_2nd      = [1/9.81, 1/9.81,     1,   1,       1,           1,       1,     1,        1,    1,   1,   1,      1,        1,       1,        1];
% datatype_2nd          = repmat({'ins'}, 1, length(varname_2nd)); % instantaneous (mean) or accumulated (sum)
% whichfolder_2nd       = repmat(2, 1, length(varname_2nd));
% StepsInADay_2nd       = repmat(8, 1, length(varname_2nd));

%/ Global ERA5 *hrly* 1x1 data (acc)
% varname_3rd           = {  'tp'};
% dataname_callfile_3rd = {'prcp'};
% dataname_3rd          = {'prcp'};
% dataunitconv_3rd      = [  1000];
% datatype_3rd          = repmat({'acc'}, 1, length(varname_3rd)); % instantaneous (mean) or accumulated (sum)
% whichfolder_3rd       = repmat(3,  1, length(varname_3rd));
% StepsInADay_3rd       = repmat(24, 1, length(varname_3rd));

%/ Global CERA-20C data
% varname_4th           = {       'tp',       'e'};
% dataname_callfile_4th = {      'E_P',     'E_P'};
% dataname_4th          = {'prcp_CERA',  'E_CERA'};
% dataunitconv_4th      = [       1000,    -1000];
% datatype_4th          = repmat({'acc'}, 1, length(varname_4th)); % instantaneous (mean) or accumulated (sum)
% whichfolder_4th       = repmat(4,  1, length(varname_4th));
% StepsInADay_4th       = repmat(2, 1, length(varname_4th));
% fc_steps               = [6, 24]; %/ assume the same steps for all fc vars.

% varname           = cat(2, varname_1st,           varname_2nd,           varname_3rd,           varname_4th);
% dataname_callfile = cat(2, dataname_callfile_1st, dataname_callfile_2nd, dataname_callfile_3rd, dataname_callfile_4th);
% dataname          = cat(2, dataname_1st,          dataname_2nd,          dataname_3rd,          dataname_4th);
% dataunitconv      = cat(2, dataunitconv_1st,      dataunitconv_2nd,      dataunitconv_3rd,      dataunitconv_4th);
% datatype          = cat(2, datatype_1st,          datatype_2nd,          datatype_3rd,          datatype_4th);
% whichfolder       = cat(2, whichfolder_1st,       whichfolder_2nd,       whichfolder_3rd,       whichfolder_4th);
% StepsInADay       = cat(2, StepsInADay_1st,       StepsInADay_2nd,       StepsInADay_3rd,       StepsInADay_4th);

%% [Step 2] Read Flexpart Experiment Parameters
WSV = [];
%/ Derive forwad trajectories from domain-fill and the uptake map
fprintf('*** See Matlab script "flex_save_traj_uptake.m" for deriving trajs and uptake map *** \n')

%===========================
hotspot_list_ver = 3;        
%===========================

%/ partial data available
% year_list = 1992:2010; ldirect = 'fwd'; from_hs = 1; str_remark = '_RH2_correct';  RHc = 85; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = [];   BLH_factor = 1; str_src_z = '_nozbound';

%/ full data available
% year_list = 1971:2010; ldirect = 'bwd'; from_hs = 1; str_remark = '_RH2';  RHc = 85; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = []; BLH_factor = 1; str_src_z = '_nozbound';
% year_list = 1971:2010; ldirect = 'bwd'; from_hs = 0; str_remark = '_RH2';  RHc = 85; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = []; BLH_factor = 1; str_src_z = '_nozbound';
year_list = 1971:2010; ldirect = 'fwd'; from_hs = 1; str_remark = '_RH2_correct';  RHc = 85; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = [];   BLH_factor = 1; str_src_z = '_nozbound';


% ldirect = 'bwd';  RHc = 85; dqc = 0.5; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = []; BLH_factor = 1; str_remark = '_RH2'; str_src_z = '_nozbound';
% ldirect = 'bwd';  RHc = 85; dqc = 0.2; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = []; BLH_factor = 1; str_remark = '_RH2'; str_src_z = '_nozbound';

% ldirect = 'bwd';  RHc = 90; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = []; BLH_factor = 1; str_remark = '_RH2'; str_src_z = '_nozbound';
% ldirect = 'bwd';  RHc = 87.5; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = []; BLH_factor = 1; str_remark = '_RH2'; str_src_z = '_nozbound';
% ldirect = 'bwd';  RHc = 80; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = []; BLH_factor = 1; str_remark = '_RH2'; str_src_z = '_nozbound';
% ldirect = 'bwd'; RHc = 80; dqc = 0.2; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = []; BLH_factor = 1; str_remark = '_RH2'; str_src_z = '_nozbound';
% ldirect = 'bwd'; RHc = 80; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 1; rm_when_dz = 5000; BLH_factor = 1; str_remark = '_new'; str_src_z = '_nozbound';
% ldirect = 'bwd'; RHc = 80; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 1; rm_when_dz = 5000; BLH_factor = 1.5; str_remark = '_new'; str_src_z = [];
% ldirect = 'bwd'; RHc = 85; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = [];   BLH_factor = 1.5; str_remark = []; str_src_z = [];  
% ldirect = 'bwd'; RHc = 90; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = [];   BLH_factor = 1.5; str_remark = []; str_src_z = [];    
% ldirect = 'bwd'; RHc = 80; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = [];   BLH_factor = 1.5; str_remark = []; str_src_z = [];    
% ldirect = 'bwd'; RHc = 80; dqc = 0.2; optimal_tracking = 1; traj_rm_jump = 0; rm_when_dz = [];   BLH_factor = 1.5; str_remark = []; str_src_z = [];  
% ldirect = 'bwd'; RHc = 80; dqc = 0.1; optimal_tracking = 1; traj_rm_jump = 1; rm_when_dz = 5000; BLH_factor = 1;   str_remark = []; str_src_z = [];    %/ full data in 1971-2010
% ldirect = 'bwd'; RHc = 80; dqc = 0.2; optimal_tracking = 0; traj_rm_jump = 1; rm_when_dz = 5000; BLH_factor = 1;   str_remark = []; str_src_z = [];  

if hotspot_list_ver == 1       %/ 17 hotspots (NOTE: since functions like box_region has changed. 
    slct_var       = 'BL_Pm';
    include_oceans = 1;
    str_figver     = ''; 
        
elseif hotspot_list_ver == 2   %/ 16 hotspots (only change AYR_hotspot, not AYR_hotspot_fwd!)
    slct_var       = 'BL_Pm';
    include_oceans = 1;  
    str_figver     = '_v2';                          
    
elseif hotspot_list_ver == 3   %/ 17 terrestrial hotspots based on Pm.
    slct_var       = 'Pm'; 
    include_oceans = 0;  
    str_figver     = '_v3';
else
    error('wrong input of hotspot_list_ver!')
end

if include_oceans == 0     str_noOceans = '_noOcns';  else  str_noOceans = '';  end


if optimal_tracking
    if ismember(ldirect, {'bwd'})
        maxtraj_day = 20;                        %/ max traj days
        optimal_rr  = 0.9;                       %/ if not [], then track recursively until rr >= optimal value
        str_optimal = '_optimal';
        
    elseif ismember(ldirect, {'fwd'})
        maxtraj_day   = 15;                      %/ max traj days
        optimal_prcnt = 0.1;                     %/ if not [], then track recursively until rr >= optimal value
        str_optimal   = '_optimal';
        
    end
else
    optimal_rr = [];
    optimal_prcnt  = [];                     %/ if not [], then track recursively until rr >= optimal value
    maxtraj_day = 10;                        %/ max traj days
    str_optimal = [];
end

%/ Load the parameters                 
expmnt = 'domfill_CERA_MPI';
dt_slct = 3;
dt  = 3;                                 %/ FLEXPART time interval (in hr)
leap = dt_slct/dt;                       %/ Leap from FLEXPART time interval to user-defined interval
maxtraj_ts = 24*maxtraj_day/dt_slct + 1; %/ traj time steps, +1 to include the very last time step
res = 1;
lon = 0:res:360-res;                     %/ this is the lon that I input to calcuate uptake map, although in FLexpart it is in [-179, 180].
lat = -90:res:90;

%/ set the strings
seq = sequence('data', year_list);
str_years = '';
for i = 1:length(seq)
    a = strjoin({num2str(seq{i}(1)), num2str(seq{i}(end))}, '-');
    str_years = strjoin({str_years, a}, '_');
end

if traj_rm_jump == 0      str_traj_rm_jump = '_intact';                        else str_traj_rm_jump = [];   end
if BLH_factor ~= 1        str_BLH_factor   = sprintf('_%.1fBLH', BLH_factor);  else  str_BLH_factor   = [];  end

%/ indicate where to load the WSV data from
WSV_dir     = cell(length(year_list),1);
WSV_dir_hs3 = cell(length(year_list),1);
flexpart_date_dt = cell(length(year_list),1);
for y = 1:length(year_list)
    year = year_list(y);
    
    %========================= Set WSV folder =============================
    if ismember(ldirect, {'fwd'})
        WSV_folder     = '/disk/r034/tfchengac/FLEXPART/';   
        WSV_folder_hs3 = '/disk/r037/tfchengac/FLEXPART/';  
        
    elseif ismember(ldirect, {'bwd'})
        if year >= 1991
            WSV_folder = '/disk/r034/tfchengac/FLEXPART/'; 
        else
            WSV_folder = '/disk/r059/tfchengac/FLEXPART/'; 
        end
        WSV_folder_hs3 = WSV_folder;  %/ just the same
    end
    WSV_dir{y}     = strcat(WSV_folder,     expmnt, '/prcssd_data_', num2str(year),'/');
    WSV_dir_hs3{y} = strcat(WSV_folder_hs3, expmnt, '/prcssd_data_', num2str(year),'/');
    %========================= Load the header ============================
    if year >= 2001 && year <= 2010
        flexpart_output_dir = strcat('/disk/r037/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
    else if year >= 1991 && year <= 2000
        flexpart_output_dir = strcat('/disk/r034/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
    else if year >= 1981 && year <= 1990
        flexpart_output_dir = strcat('/disk/r012/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
    else if year >= 1971 && year <= 1980
        flexpart_output_dir = strcat('/disk/r014/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
    end
    end
    end
    end
    
    readp = 1; %/ read release points (0/1)
    [header, fail] = flex_header(flexpart_output_dir, 0, readp, 0);
    header.dates_str = string(importdata([flexpart_output_dir 'dates'])); %/ read dates file in string
    header.dates_dt = datetime(header.dates_str,'InputFormat','yyyyMMddHHmmss', 'Format', 'yyyy-MM-dd HH:mm:ss');
    area = calc_grid_area_header(header); %/ my function to calculating gird area.
    
    %/ Since we added a buffer month before each year's expmnt, 
    %/ we have to save the datetime in cells without concatenating them!!
    if      leap == 1    flexpart_date_dt{y} = header.dates_dt(1:leap:end);
    else if leap == 2    flexpart_date_dt{y} = header.dates_dt(2:leap:end); %/ since it starts from year(t-1)-12-01 03:00
    end
    end
end

%/ load AYR_hotspot (based on glbland)
AYR_hs_data_filename = strcat(masterfolder, expmnt, '/prcssd_data_4plotting/',...
                             'AYR_hs_p95_area_sum_', slct_var, str_noOceans,'_RHc85_dqc0.1_bwd_20d_optimal_3h_glbland_nozbound_intact_RH2_1971-2010', str_figver,'.mat');
fprintf('*** The queried data are found. *** \n*** Loading: %s *** \n', AYR_hs_data_filename)
load(AYR_hs_data_filename);                                            %/ output is AYR_hotspot (based on glbland)


%/ set the experiment name and folder name.
AYR_hotspot_fwdbwd_bc = [];
if ismember(ldirect, {'fwd'})
    
    str_domain     = '';   %/ use the old list of AYR_hotspot_fwd.
    str_expmntinfo = strcat({' RHc'}, num2str(RHc), {' dqc'}, num2str(dqc),  {' '}, ldirect, {' '}, num2str(maxtraj_day),...
                            {'d'}, str_optimal, {' '}, num2str(dt_slct), 'h',  str_domain, str_src_z, str_traj_rm_jump, str_BLH_factor, str_remark);

    foldername = strcat({' RHc'}, num2str(RHc), {' dqc'}, num2str(dqc),  {' '}, ldirect, {' '}, num2str(maxtraj_day),...
                            {'d'}, str_optimal, {' '}, num2str(dt_slct), 'h', str_src_z, str_traj_rm_jump, str_BLH_factor, str_remark);

    %/ load AYR_hotspot_fwd
    AYR_hs_fwd_data_filename = string(strcat(masterfolder, expmnt, '/prcssd_data_4plotting/', 'AYR_hs_fwd', str_expmntinfo, str_years, str_figver,  '.mat'));
    AYR_hs_fwd_data_filename = strrep(AYR_hs_fwd_data_filename, ' ', '_');
    if isfile(AYR_hs_fwd_data_filename)
        fprintf('*** The queried data are found. *** \n*** Loading: %s *** \n', AYR_hs_fwd_data_filename)
        tic; load(AYR_hs_fwd_data_filename); toc;                          %/ output is AYR_hotspot_fwd
        AYR_hotspot_fwdbwd_bc = AYR_hotspot_fwd;                           %/ make a broadcast var for AYR_hotspot_fwd/AYR_hotspot_bwd for easy coding.
    else
        warning('!!! %s Data Not Found. Skip loading. !!!\n', AYR_hs_fwd_data_filename)
    end
    
elseif ismember(ldirect, {'bwd'})
    str_domain = '_glbland';
    str_expmntinfo = strcat({' RHc'}, num2str(RHc), {' dqc'}, num2str(dqc),  {' '}, ldirect, {' '}, num2str(maxtraj_day),...
                            {'d'}, str_optimal, {' '}, num2str(dt_slct), 'h',  str_domain, str_src_z, str_traj_rm_jump, str_BLH_factor, str_remark);
    foldername = str_expmntinfo;

    %/ load AYR_hotspot_bwd
    AYR_hs_bwd_data_filename = string(strcat(masterfolder, expmnt, '/prcssd_data_4plotting/', 'AYR_hs_bwd', str_expmntinfo, str_years, str_figver, '.mat')); %/ since not yet complete 40-yr running.
    AYR_hs_bwd_data_filename = strrep(AYR_hs_bwd_data_filename, ' ', '_');
    if isfile(AYR_hs_bwd_data_filename)
        fprintf('*** The queried data are found. *** \n*** Loading: %s *** \n', AYR_hs_bwd_data_filename)
        tic; load(AYR_hs_bwd_data_filename); toc;                          %/ output is AYR_hotspot_bwd
        AYR_hotspot_fwdbwd_bc = AYR_hotspot_bwd;                           %/ make a broadcast var for AYR_hotspot_fwd/AYR_hotspot_bwd for easy coding.
    else
        warning('!!! %s Data Not Found. Skip loading. !!!\n', AYR_hs_bwd_data_filename)
    end
end

plotting_folder = strcat(masterfolder, expmnt, '/Plottings', strrep(foldername, ' ', '_'), '/');
mkdir(char(plotting_folder));

%/ how would dqc = 0.2 g/kg means in terms of mm/day?
A = mean(mean(area));
mass = 5.095e18 / 4995000;

P_LA = dqc * 10e-3 * mass * 24/dt / A;
% fprintf('*** dqc = %.3f g/kg corresponds to %.3f mm/day ***\n', dqc, P_LA);

%% [Step 3] Climate Indices
clc; close all;
ClimIndices = [];
datafolder = '/disk/r128/tfchengac/fandyr128SOM/Climate Indices/';

indexname =  {'MJO',           'GMT',     'SAM',   'IOD_HadISST',  'EP_ENSO', 'CP_ENSO', 'Mixed_ENSO',    'ONI', 'NINO34_HadISST',  'NINO3_HadISST',  'NINO4_HadISST',    'AMO',    'PDO',     'NP',    'PNA',     'WP',     'AO',    'NAO',  'SolarFlux'};
dataformat = {'%f',             '%f',      '%f',            '%f',       '%f',      '%f',         '%f',     '%f',             '%f',             '%f',             '%f',     '%f',     '%f',     '%f',     '%f',     '%f',     '%f',    '%f,',     '%f'};
type =       {'MJO_format', 'Stndrd',  'Stndrd',        'Stndrd',   'Stndrd',  'Stndrd',     'Stndrd', 'Stndrd',         'Stndrd',         'Stndrd',         'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'TwoCol', 'Stndrd', 'OneCol',   'OneCol'};
skipheadings=[   2,                0,        0,                0,          0,         0,            0,        0,                0,                0,                0,        0,        0,        0,        0,        9,        0,        0,        0];

% indexname =  {'SAM',    'ONI',    'NINO34', 'NINO4',  'AMO',    'PDO',    'NP',     'PNA',    'WP',     'AO',     'NAO',    'SolarFlux',  'EuraSnowCover', 'NAGLSnowCover', 'NHSnowCover'};
% dataformat = {'%f',     '%f',      '%f',     '%f',     '%f',     '%f',     '%f',     '%f',     '%f',     '%f',     '%f,',    '%f',         '%f,',           '%f,',           '%f,'};
% type =       {'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'TwoCol', 'Stndrd', 'OneCol', 'OneCol',     'OneCol',        'OneCol',        'OneCol'};
% skipheadings=[   0,        0,        0,        0,        0,        0,        0,        0,        9,        0,        0,        0,            4,               4,               4];
for i = 1:length(indexname)
    a = read_ClimIndex('type', type{i}, 'datafolder', datafolder, 'skipheadings', skipheadings(i), 'indexname', indexname{i},...
                       'select_year', year_list, 'dataformat', dataformat{i});
    fld = fieldnames(a);
    for f= 1:length(fld)
        ClimIndices.(fld{f}) = a.(fld{f});
    end
end

%/ EP ENSO (ONI)
%  Warm (red) and cold (blue) periods based on a threshold of +/- 0.5oC for the Oceanic Nino Index (ONI)
%  for at least 5 consecutive overlapping seasons
%  https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php

% NOTE: my result is slightly different from that in the website which
%       is based on a ONI index rounded to 1 digit only.

ClimIndices.ENSO_event = CI2events('CI_dates', ClimIndices.ONI(:,1), 'CI_data', ClimIndices.ONI(:,2), 'thres', 0.5, 'n', 5);
% ClimIndices.EP_ENSO_event = CI2events('CI_dates', ClimIndices.EP_ENSO(:,1), 'CI_data', ClimIndices.EP_ENSO(:,2), 'thres', 0.7, 'n', 5); %/ my guess
% ClimIndices.CP_ENSO_event = CI2events('CI_dates', ClimIndices.CP_ENSO(:,1), 'CI_data', ClimIndices.CP_ENSO(:,2), 'thres', 0.7, 'n', 5); %/ my guess

%/ IOD
%  For monitoring the IOD, Australian climatologists consider sustained values above +0.4 °C 
%  as typical of a positive IOD, and values below −0.4 °C as typical of a negative IOD.
%  http://www.bom.gov.au/climate/enso/indices/about.shtml

%/ NOTE: Since no statement on how 'sustain' it should be. Assume n = 5 months.

ClimIndices.IOD_event = CI2events('CI_dates', ClimIndices.IOD_HadISST(:,1), 'CI_data', ClimIndices.IOD_HadISST(:,2), 'thres', 0.4, 'n', 5);

ClimIndices.AMO_event = CI2events('CI_dates', ClimIndices.AMO(:,1), 'CI_data', ClimIndices.AMO(:,2), 'thres', 0.1, 'n', 5); %/ here 0.1 refers to 0.5sd

%/ MJO events 
%/ Amp > 1 and persistent for >= 3 days (my choice)
MJO_dates = ClimIndices.MJO(:,1);
MJO_amp   = ClimIndices.MJO(:,2);
MJO_phase = ClimIndices.MJO(:,3);

thres = 1;
% n = 1;     
n = 3;
M_hankel_amp   = hankel(MJO_amp(1:n),   MJO_amp(n:end));
M_hankel_phase = hankel(MJO_phase(1:n), MJO_phase(n:end));

MJO_mature_phase = zeros(size(MJO_phase));
for k = 1:8
    %/ find each moving n length of data with the same phase and strong enough amp
    mjo_phase_ind = find(all(M_hankel_amp > thres, 1) == 1   &...
                         all(M_hankel_phase == k, 1) == 1);  
                     
    for i = 1:length(mjo_phase_ind)
        MJO_mature_phase(mjo_phase_ind(i):(mjo_phase_ind(i) + n - 1), 1) = k; %/ NOTE: whatever a column of n data met the criteria, we mark *all the n data* instead of one.                                                       
    end
end
ClimIndices.MJO_event(:,1) = MJO_dates;
ClimIndices.MJO_event(:,2) = MJO_mature_phase;

% for k = 1:8
%     disp(length(find( MJO_mature_phase  == k)))
% end

%/ Zhang et al 2019's method to classify EP, CP and Mixed ENSO (Assume N3 and N4 have the same date)
N3_dates = ClimIndices.NINO3_HadISST(:,1);
N4_dates = ClimIndices.NINO4_HadISST(:,1);
if ~isequal(N3_dates, N4_dates)  error('Make sure N3 and N4 have the same date!!'); end

N3 = ClimIndices.NINO3_HadISST(:,2);
N4 = ClimIndices.NINO4_HadISST(:,2);

%/ 3-month running mean
N3 = movmean(N3, 3);
N4 = movmean(N4, 3);

ClimIndices.UCEI(:,1) = N3_dates;
ClimIndices.UCEI(:,2) = sqrt(2*(N3.^2 + N4.^2));                           %/ r
ClimIndices.UCEI(:,3) = atand((N3-N4)./(N3+N4));                           %/ theta, atand in degree

ind_La = find(N3+N4 < 0);
ClimIndices.UCEI(ind_La,3) = ClimIndices.UCEI(ind_La,3) - 180;
                              
complex_ENSO_event      = zeros(length(N3_dates), 2);                      
complex_ENSO_event(:,1) = N3_dates;

r        = ClimIndices.UCEI(:,2);
thres    = 0.5;   n = 5;                                                   %/ Condition: r > 0.5 for at least 5 months
M_hankel = hankel(r(1:n), r(n:end));
enso_ind = find(all(M_hankel > thres, 1) == 1);  
for i = 1:length(enso_ind)
    complex_ENSO_event(enso_ind(i):(enso_ind(i) + n - 1), 2) = 999;        %/ NOTE: whatever a column of 5 months met the criteria, we mark *all the 5 months* instead of one.
                                                                           %/       mark as 999 since we will further classify the type based on theta
end
theta_range =     [  15   90    1;        %/ EP    El (1)
                    -15   15    2;        %/ Mixed El (2)
                    -90  -15    3;        %/ CP    El (3)
                   -165  -90   -1;        %/ EP    La (-1)
                   -195  -165  -2;        %/ Mixed La (-2)
                   -270  -195  -3;];      %/ CP    La (-3)
                   
theta = ClimIndices.UCEI(:,3);
for i = 1:length(theta_range)
    
    ind_phase = find(complex_ENSO_event(:,2) == 999   &...
                     theta > theta_range(i,1) &...
                     theta < theta_range(i,2));
                 
    complex_ENSO_event(ind_phase,2) = theta_range(i,3);
end

ClimIndices.complex_ENSO_event = complex_ENSO_event;

%% Plot MJO weather maps




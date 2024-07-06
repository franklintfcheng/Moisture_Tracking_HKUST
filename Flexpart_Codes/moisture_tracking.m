%% Main
%================================= NOTE ===================================
%/      Author: Franklin Cheng (fandycheng@ust.hk)
%/ Last update: 28 Jun 2024
%/
%/ Description: This program generates backward (bwd) or forward (fwd) 
%/              trajectories and the corresponding moisture contributions
%/              using WaterSip or WaterDrip algorithm.
%/
%/ 28 Jun 2024: Enabled WaterSip to compute source contributions from
%/              Atmospheric Rivers (ARs)
%/
%/ Run this program in Linux (for example): 
%/              tmux
%/              cd /disk/r059/tfchengac/FLEXPART/domfill_EA_0.5deg_MPI_202309-10_ScotCase
%/              nohup matlab2023a -nodisplay -nosplash -nodesktop -r "run('/home/tfchengac/Flexpart_Codes/moisture_tracking.m');" > tracking.log </dev/null &
%==========================================================================
clearvars;

% matlab_package_path = '/disk/r059/tfchengac/codes4flexpart_tutorial/'; %/ Set the path of the folder where you store the matlab packages
% masterfolder        = '/disk/r059/tfchengac/FLEXPART/';                %/ Set your master path (i.e., where you want your data to be stored)
% send_email_to       = [];                                              %/ Set your email to send notification after the program completed. Do NOT use connect.ust.hk email. It will just go to '/var/spool/mail/<Your Account>' in the ust server.

matlab_package_path = '/home/tfchengac/';                     %/ Set the path of the folder where you store the matlab packages
masterfolder        = '/disk/r059/tfchengac/FLEXPART/';       %/ Set your master path (i.e., where you want your data to be stored)
send_email_to       = 'franklintfcheng@gmail.com';            %/ Set your email to send notification after the program completed. Do NOT use connect.ust.hk email. It will just go to '/var/spool/mail/<Your Account>' in the ust server.

%==========================================================================
addpath(genpath(fullfile(matlab_package_path,'MyMatlabPkg')));
addpath(genpath(fullfile(matlab_package_path,'MyMatlabFn')));
addpath(genpath(fullfile(matlab_package_path,'m_map1.4o')));
addpath(genpath(fullfile(matlab_package_path,'Flexpart_Codes')));

%==========================================================================
test_mode         = 0;                                    %/ See the test_mode setting below
Njob              = 1;                                    %/ Whether to split the job into N jobs (1 by default)
job_id            = 1;                                    %/ Which splitted job to run (1 by default)
reverse_loop      = 0;                                    %/ Whether to reverse the t_list loop

save_catalog      = 1;                                    %/ Whether to save basin catalog
recompute_catalog = 0;                                    %/ Whether to re-generate basin catalog
compute_AR_contr  = 1;                                    %/ Whether to compute source contribution by Atmospheric Rivers

%========= Set the following paths and variables to suit your needs! =======================================================
expmnt           = 'domfill_EA_1.0deg_MPI_201806_ALA';
year_list        = 2018;                                 %/ The year period, vector
stdate           = 20180615;                             %/ Start date (if set [], then will tracking from the first date available)
eddate           = 20180619;                             %/ End date (if set [], then will tracking till the last date available)
trajdirect       = 'bwd';  from_basin = 11;              %/ See the description below    
store_WSV_yearly = 0;                                    %/ Whether to save data into prcssd_data_2001, prcssd_data_2002, etc. folders (useful when server space is limited)
output_res       = 0.25;                                 %/ Controls the resolution of output source map (regardless of the native resolution of FLEXPART)
dt_slct          = 3;                                    %/ traj time interval (in hr; user-defined)
dt               = 3;                                    %/ FLEXPART time interval (in hr)
optimal_rr       = 0.99;                                 %/ Track recursively until rr >= optimal_rr; for 'fwd', track the contribution until f < (1 - optimal_rr)
RHc_dqc_scheme   = 1;                                    %/ RHc & dqc scheme (for WaterSip/WaterDrip), set 1 for the default scheme (RHc=85, dqc=0.1)
FLEXPART_folder  = '/disk/r148/ysongbn/flexpart_v10.4/output/domfill_EA_1.0deg_MPI_201806_ALA/';  %/ Set your own folder where your FLEXPART model output was stored (Set [] will lead to the default folder)
WSV_folder       = strcat(masterfolder, expmnt, '/prcssd_data/');                                   %/ Set your own folder where you would like to store the moisture tracking output (Set [] will lead to the default folder)

% expmnt           = 'domfill_EA_0.5deg_MPI_202309-10_ScotCase';  %/ Input the folder name of your FLEXPART experiment
% year_list        = 2023;                                 %/ The year period, vector
% stdate           = 20231006;                             %/ Start date (if set [], then will tracking from the first date available)
% eddate           = 20231008;                             %/ End date (if set [], then will tracking till the last date available)
% trajdirect       = 'bwd';  from_basin = 10;              %/ See the description below    
% store_WSV_yearly = 0;                                    %/ Whether to save data into prcssd_data_2001, prcssd_data_2002, etc. folders (useful when server space is limited)
% output_res       = 0.25;                                 %/ Controls the resolution of output source map (regardless of the native resolution of FLEXPART)
% dt_slct          = 3;                                    %/ traj time interval (in hr; user-defined)
% dt               = 3;                                    %/ FLEXPART time interval (in hr)
% optimal_rr       = 0.99;                                 %/ Track recursively until rr >= optimal_rr; for 'fwd', track the contribution until f < (1 - optimal_rr)
% RHc_dqc_scheme   = 2;                                    %/ RHc & dqc scheme (for WaterSip/WaterDrip), set 1 for the default scheme (RHc=85, dqc=0.1)
% % RHc_dqc_scheme   = 2;                                    %/ RHc & dqc scheme (for WaterSip/WaterDrip), set 1 for the default scheme (RHc=85, dqc=0.1)
% FLEXPART_folder  = [];                                   %/ Set your own folder where your FLEXPART model output was stored (Set [] will lead to the default folder)
% WSV_folder       = [];                                   %/ Set your own folder where you would like to store the moisture tracking output (Set [] will lead to the default folder)

% expmnt           = 'domfill_EA_0.5deg_MPI_202202_AusCase';  %/ Input the folder name of your FLEXPART experiment
% year_list        = 2022;                                 %/ The year period, vector
% stdate           = 20220222;                             %/ Start date (if set [], then will tracking from the first date available)
% eddate           = 20220228;                             %/ End date (if set [], then will tracking till the last date available)
% trajdirect       = 'bwd';  from_basin = 9;               %/ See the description below    
% store_WSV_yearly = 0;                                    %/ Whether to save data into prcssd_data_2001, prcssd_data_2002, etc. folders (useful when server space is limited)
% output_res       = 0.25;                                 %/ Controls the resolution of output source map (regardless of the native resolution of FLEXPART)
% dt_slct          = 3;                                    %/ traj time interval (in hr; user-defined)
% dt               = 3;                                    %/ FLEXPART time interval (in hr)
% optimal_rr       = 0.99;                                 %/ Track recursively until rr >= optimal_rr; for 'fwd', track the contribution until f < (1 - optimal_rr)
% RHc_dqc_scheme   = 0;                                    %/ RHc & dqc scheme (for WaterSip/WaterDrip), set 1 for the default scheme (RHc=85, dqc=0.1)
% % RHc_dqc_scheme   = 1;                                    %/ RHc & dqc scheme (for WaterSip/WaterDrip), set 1 for the default scheme (RHc=85, dqc=0.1)
% FLEXPART_folder  = [];                                   %/ Set your own folder where your FLEXPART model output was stored (Set [] will lead to the default folder)
% WSV_folder       = [];                                   %/ Set your own folder where you would like to store the moisture tracking output (Set [] will lead to the default folder)

% expmnt           = 'domfill_EA_0.5deg_MPI_202207-08_PakCase'; %/ Input the folder name of your FLEXPART experiment
% year_list        = 2022;                                 %/ The year period, vector
% stdate           = 20220810;                             %/ Start date (if set [], then will tracking from the first date available)
% eddate           = 20220824;                             %/ End date (if set [], then will tracking till the last date available)
% trajdirect      = 'bwd';  from_basin = 7;                %/ See the description below   
% % trajdirect      = 'fwd';  from_basin = 8;             %/ See the description below   
% store_WSV_yearly = 0;                                    %/ Whether to save data into prcssd_data_2001, prcssd_data_2002, etc. folders (useful when server space is limited)
% output_res       = 0.25;                                 %/ Controls the resolution of output source map (regardless of the native resolution of FLEXPART)
% dt_slct          = 3;                                    %/ traj time interval (in hr; user-defined)
% dt               = 3;                                    %/ FLEXPART time interval (in hr)
% optimal_rr       = 0.99;                                 %/ Track recursively until rr >= optimal_rr; for 'fwd', track the contribution until f < (1 - optimal_rr)
% RHc_dqc_scheme   = 0;                                    %/ RHc & dqc scheme (for WaterSip/WaterDrip), set 1 for the default scheme (RHc=85, dqc=0.1)
% % RHc_dqc_scheme   = 1;                                    %/ RHc & dqc scheme (for WaterSip/WaterDrip), set 1 for the default scheme (RHc=85, dqc=0.1)
% FLEXPART_folder  = [];                                   %/ Set your own folder where your FLEXPART model output was stored (Set [] will lead to the default folder)
% WSV_folder       = [];                                   %/ Set your own folder where you would like to store the moisture tracking output (Set [] will lead to the default folder)

%/ Tracking from from TP baisns (Cheng et al. 2024)
% expmnt           = 'domfill_CERA_MPI';                   %/ Input the folder name of your FLEXPART experiment
% year_list        = 1971:2010;                            %/ The year period, vector
% stdate           = 101;                                  %/ Start date (mmdd or yyyymmdd) (if set [], then will tracking from the first date available)
% eddate           = 1230;                                 %/ End date   (mmdd or yyyymmdd) (if set [], then will tracking till the last date available)
% trajdirect       = 'bwd';  from_basin = 4;               %/ See the description below    
% store_WSV_yearly = 1;                                    %/ Whether to save data into prcssd_data_2001, prcssd_data_2002, etc. folders (useful when server space is limited)
% output_res       = 1;                                    %/ Controls the resolution of output source map (regardless of the native resolution of FLEXPART)
% dt_slct          = 3;                                    %/ traj time interval (in hr; user-defined)
% dt               = 3;                                    %/ FLEXPART time interval (in hr)
% optimal_rr       = 0.99;                                 %/ Track recursively until rr >= optimal_rr; for 'fwd', track the contribution until f < (1 - optimal_rr)
% RHc_dqc_scheme   = 23;                                   %/ RHc & dqc scheme (for WaterSip/WaterDrip), set 1 for the default scheme (RHc=85, dqc=0.1)
% FLEXPART_folder  = [];                                   %/ Set your own folder where your FLEXPART model output was stored (Set [] will lead to the default folder)
% WSV_folder       = [];                                   %/ Set your own folder where you would like to store the moisture tracking output (Set [] will lead to the default folder)

%/ Tracking from global land (Cheng and Lu 2023)
% expmnt           = 'domfill_CERA_MPI';                   %/ Input the folder name of your FLEXPART experiment
% year_list        = 2010;                                 %/ The year period, vector
% stdate           = 101;                                  %/ Start date (mmdd or yyyymmdd) (if set [], then will tracking from the first date available)
% eddate           = 1230;                                 %/ End date   (mmdd or yyyymmdd) (if set [], then will tracking till the last date available)
% trajdirect       = 'bwd';  from_basin = 0;               %/ See the description below    
% store_WSV_yearly = 1;                                    %/ Whether to save data into prcssd_data_2001, prcssd_data_2002, etc. folders (useful when server space is limited)
% output_res       = 1;                                    %/ Controls the resolution of output source map (regardless of the native resolution of FLEXPART)
% dt_slct          = 3;                                    %/ traj time interval (in hr; user-defined)
% dt               = 3;                                    %/ FLEXPART time interval (in hr)
% optimal_rr       = 0.9;                                  %/ Track recursively until rr >= optimal_rr; for 'fwd', track the contribution until f < (1 - optimal_rr)
% RHc_dqc_scheme   = 10;                                   %/ RHc & dqc scheme (for WaterSip/WaterDrip), set 1 for the default scheme (RHc=85, dqc=0.1)
% FLEXPART_folder  = [];                                   %/ Set your own folder where your FLEXPART model output was stored (Set [] will lead to the default folder)
% WSV_folder       = [];                                   %/ Set your own folder where you would like to store the moisture tracking output (Set [] will lead to the default folder)

%========= Set the backward (bwd) or forward (fwd) tracking setting to suit your need! ==========================
if isequal(trajdirect, 'bwd')
    output_rr_map = 0;              %/ For 'bwd' only, [1]: Output continental precip. recycling ratio map (for from_basin = 0)
    output_LA_3D  = 0;              %/ For 'bwd' only, [1]: Output 3D traj matrix centered at the basin (for from_basin ~= -1 or 0)
    maxtraj_day   = 20;             %/ max tracking days
    sharpcut      = 0;             %/ For 'fwd' only, [1]: Sharp cut the forward tracking by the end day 
    str_remark    = 'RH2_';         %/ Do not change it
else
    output_rr_map = 0;             %/ For 'bwd' only, [1]: Output continental precip. recycling ratio map (for from_basin = 0)
    output_LA_3D  = 0;             %/ For 'bwd' only, [1]: Output 3D traj matrix centered at the basin (for from_basin ~= -1 or 0)
    maxtraj_day   = 20;            %/ max tracking days         
    sharpcut      = 1;             %/ For 'fwd' only, [1]: Sharp cut the forward tracking by the end day 
    str_remark    = 'RH2_';        %/ Do not change it
end
%========= Default Setting =============================================================================
traj_rm_jump     = 0; rm_when_dz = [];          %/ Toggle on to slice traj of a particle with a *big* jump
BLH_factor       = 1;                           %/ Multiplier on BLH to account for small-scale variability
flag_successful  = 1;
server_name      = getenv('HOSTNAME');          %/ get the server name
server_name_cell = strsplit(server_name, '.');
server_name      = server_name_cell{1};
data_folder      = strcat(masterfolder, expmnt, '/prcssd_data_4plotting/');
maxtraj_ts       = 24*maxtraj_day/dt_slct + 1;      %/ Maximum traj time steps, +1 to include the very last time step
lon              = 0:output_res:360-output_res;     %/ The lon of output source map 
lat              = -90:output_res:90;               %/ The lat of output source map 
lon_m179_180     = conv_to_lon_m179_180(lon);       %/ obtain a full lon array in [-179, 180)
mkdir(data_folder)

%========= Generate basin_catalog (from_basin) =========================================
%/       -1]: global,                   0]: from global land,               1]: from 17 original hotspots, 
%/        2]: from 3 new AYR hotspot,   3]: 17 most-updated AYR hotspots,   4]: from the TP basins,
%/        5]: from TP main LoVeg area,  6]: from TP grid,                   7]: from Pakistan (a box region)    
%/        8]: from IPCC regions (for Pakistan event)
%/        9]: from Australia (a box region)   
%/       10]: from Scotland (a box region)   
%==========================================================================

%/ Below are the regions related to from_basin index (for my own projects)
if     from_basin == -1  slct_reg = {'global'};  
elseif from_basin == 0   slct_reg = {'land'};    %/ Global land
elseif from_basin == 4   slct_reg = {'TP_basins'};  
elseif from_basin == 5   slct_reg = {'TP_Grass', 'TP_Semidesert', 'TP_Tundra'};  
elseif from_basin == 6   slct_reg = {'TP_grids_1x1'};  
elseif from_basin == 7   slct_reg = {'Pakistan_box'};  
elseif from_basin == 8   slct_reg = {'IPCC-PAK'};  
elseif from_basin == 9   slct_reg = {'Australia_box'};  
elseif from_basin == 10  slct_reg = {'Scotland_box'};  
elseif from_basin == 11  slct_reg = {'ALA'};
else
    error('code not set!');
end

%/ Create/Load the basin catalog
%/ NOTE: For any new slct_reg, remember to adapt the function 'reg_extractor.m' 
[basin_catalog, str_domain_trajfile, str_domain] = create_basin_catalog('slct_reg', slct_reg, 'lon', lon, 'lat', lat,...
                                                            'data_folder', data_folder, 'save_catalog', save_catalog, 'recompute_catalog', recompute_catalog);

%/ Load AR logical matrix 
AR = [];
if compute_AR_contr
    AR_filename = fullfile(data_folder, 'atmospheric_river_logical_matrix.nc');
    ncdisp(AR_filename)
    AR.logical  = logical(ncread(AR_filename, 'atmospheric_river_region'));
    AR.lon      = ncread(AR_filename, 'longitude');
    AR.lat      = ncread(AR_filename, 'latitude');
    AR.time     = ncread(AR_filename, 'time');
    date_format = 'yyyyMMddHHmm';
    AR.dates    = timesince2date('filename', AR_filename, 'date_format', date_format);
    AR.dates_dt = int2datetime(AR.dates, date_format);
end

%=============================== Deprecated Code ==========================
%/ Generate basin_catalog for each slct_reg
% [~, reg_bndry_list, reg_name_list, reg_id_list, ~] = reg_extractor('lon', lon, 'lat', lat, 'slct_reg', slct_reg,...
                                                                   % 'data_folder', data_folder, 'savemat', 1, 'recompute_catalog', recompute_catalog);
%/ Create the basin catalog
% basin_catalog = create_basin_catalog('slct_reg', slct_reg, 'lon', lon, 'lat', lat, 'id', reg_id_list, 'name', reg_name_list);

% %/ Save the basin catalog
% catalog_filename = strcat(data_folder, 'catalog_', strjoin(slct_reg,'_'), '.mat');
% if save_catalog
%     fprintf('*** Saving basin_catalog into %s *** \n', catalog_filename)
%     save(catalog_filename, 'basin_catalog', '-v7.3');
% end

%/ Load the one based on from_basin
% [basin_catalog, str_from_basin_traj, str_from_basin] = load_from_basin('from_basin', from_basin);
% str_domain_trajfile = str_from_basin_traj;
% str_domain          = str_from_basin;     
%==========================================================================

NBasin = length(basin_catalog);  %/ IMPORTANT! Otherwise some parfor loops may not work!

%/ Set different no. of workers for paralell I/O & computation
if ismember(server_name, {'hqlx145', 'hqlx146', 'hqlx147', 'hqlx148', 'hqlx149'})
    if contains(expmnt, '1.0deg')
        NumWorkers_IO  = 20; NumWorkers_HPC = 20; 
    else
        NumWorkers_IO  = 10; NumWorkers_HPC = 10; %/ This is almost the *fastest* setting for 0.5deg case study!
    end
else
    NumWorkers_IO  = 5;  NumWorkers_HPC = 5;  %/ NOTE: the idled workers will still occupy a large amount of memory as they inherit the data from the main worker
end

if test_mode == 1
    savemat                       = 1;  %/ [1]: save mat data
    save_traj                     = 1;  %/ [1]: save trajectories   
    % savemat                       = 0;  %/ [1]: save mat data
    % save_traj                     = 0;  %/ [1]: save trajectories
    reload_traj                   = 0;  %/ [0]: skip if file exists.  [1]: always reload trajectories 
    overwrite_watersip_waterdrip  = 1;  %/ [0]: skip if file exists.  [1]: always overwrite  [2]: always skip 
    plot_traj                     = 0;  %/ [1]: plot trajectories (for demo)
    save_tracking_info            = 1;  %/ [1]: save moisture tracking info for analysis

elseif test_mode == 2                   %/ To plot trajectories (demo)
    savemat                       = 0;  %/ [1]: save mat data
    save_traj                     = 0;  %/ [1]: save trajectories
    reload_traj                   = 0;  %/ [0]: skip if file exists.  [1]: always reload trajectories 
    overwrite_watersip_waterdrip  = 1;  %/ [0]: skip if file exists.  [1]: always overwrite  [2]: always skip
    plot_traj                     = 1;  %/ [1]: plot trajectories (for demo)
    save_tracking_info            = 1;  %/ [1]: save moisture tracking info for analysis

else
    savemat                       = 1;  %/ [1]: save mat data
    save_traj                     = 1;  %/ [1]: save trajectories
    reload_traj                   = 0;  %/ [0]: skip if file exists.  [1]: always reload trajectories 
    overwrite_watersip_waterdrip  = 0;  %/ [0]: skip if file exists.  [1]: always overwrite  [2]: always skip 
    plot_traj                     = 0;  %/ [1]: plot trajectories (for demo)
    save_tracking_info            = 1;  %/ [1]: save moisture tracking info for analysis
end

if job_id > Njob   error('job_id cannot be larger than the total number of jobs (Njob)!');  end %/ debug
if size(year_list, 2) == 1  year_list = year_list'; end  %/ year_list must be a row vector.

if contains(expmnt, '_EA_') 
    forcing = 'EA';
elseif contains(expmnt, '_CERA_') 
    forcing = 'CERA';
else
    error('forcing not set!')
end

%/ Read RHc & dqc Scheme
[dqc, RHc, RHc_map_lon, RHc_map_lat, str_RHc_dqc] = read_RHc_dqc_scheme(RHc_dqc_scheme, forcing, output_res);

%/ Pre-process cond_land & cond_ocean [Save time and Map_Toolbox license!]
lon_grids = lon;
lat_grids = lat;
filename_cond_landocean = strcat(data_folder, sprintf('cond_landocean_lon%.2f_%.2f_from%.2fto%.2f_lat%.2f_%.2f_from%.2fto%.2f_res%.2fdeg',...
                                 min(lon_grids), max(lon_grids), lon_grids(1), lon_grids(end), min(lat_grids), max(lat_grids), lat_grids(1), lat_grids(end), output_res), '.mat');
if ~isfile(filename_cond_landocean)
    [~, cond_land, cond_ocean] = show_land_or_ocean_hydrosheds('lon_grids', lon_grids, 'lat_grids', lat_grids);
    fprintf('\n*** To avoid using Map_Toolbox, saving cond_land & cond_ocean into %s ...***\n', filename_cond_landocean)
    save(filename_cond_landocean, 'cond_land', 'cond_ocean', '-v7.3');
else
    fprintf('\n*** To avoid using Map_Toolbox, loading cond_land & cond_ocean from %s ...***\n', filename_cond_landocean)
    load(filename_cond_landocean, 'cond_land', 'cond_ocean');
end

src_lon = []; src_lat = []; src_z = [];
if traj_rm_jump == 0  str_traj_rm_jump = 'intact_';                        else  str_traj_rm_jump = [];  end
if BLH_factor ~= 1    str_BLH_factor   = sprintf('%.1fBLH_', BLH_factor);  else  str_BLH_factor   = [];  end
if isempty(src_z)     str_src_z        = 'nozbound_';                      else  str_src_z        = [];  end

if isequal(forcing, 'EA')
    if ~isempty(optimal_rr)
        str_optimal_traj = '_optimal';
        str_optimal      = sprintf('_optimal%.2f', optimal_rr);
    else
        str_optimal      = [];
        str_optimal_traj = [];
    end

elseif isequal(forcing, 'CERA')
    if ~isempty(optimal_rr)
        str_optimal_traj = '_optimal';     
        if optimal_rr ~= 0.9
            str_optimal = sprintf('_optimal%.2f', optimal_rr);   
        else
            str_optimal = '_optimal';
        end
    else
        str_optimal_traj = [];  
    end
end
str_expmntinfo = strcat({' '}, str_RHc_dqc, {' '}, trajdirect, {' '}, num2str(maxtraj_day),...
                        {'d'}, str_optimal, {' '}, num2str(dt_slct), 'h_', str_src_z, str_traj_rm_jump, str_BLH_factor, str_remark);

plotting_folder = strcat(masterfolder, expmnt, '/Plottings', strrep(str_expmntinfo, ' ', '_'), '/');

%/ Set WSV Directory (for output WaterSip / WaterDrip variables)
if ~isempty(WSV_folder)
    if ischar(WSV_folder)
        WSV_folder = {WSV_folder};
    end
    WSV_dir = repmat(WSV_folder, length(year_list),1);

else
    WSV_dir = cell(length(year_list),1);
    for y = 1:length(year_list)
        year = year_list(y);
        %/ Set the directory where the post-processed data (traj + WSV) are stored
        if ismember(trajdirect, {'bwd'})
            if isequal(forcing, 'EA')
                WSV_folder = '/disk/r149/tfchengac/FLEXPART/';
                
            elseif isequal(forcing, 'CERA')
                if RHc_dqc_scheme == 1
                    if year >= 1991
                        WSV_folder = '/disk/r034/tfchengac/FLEXPART/';  %/ occupies 10 TB
                    else
                        WSV_folder = '/disk/r059/tfchengac/FLEXPART/';  %/ occupies 10 TB
                    end
                else
                    if year >= 2001
                        WSV_folder = '/disk/r149/tfchengac/FLEXPART/';  %/ 10 yrs -> 5 TB (~500 Gb a year)
                    elseif year >= 1981 && year < 2001
                        WSV_folder = '/disk/r012/tfchengac/FLEXPART/';  %/ 20 yrs -> 10 TB (~500 Gb a year)
                    elseif year >= 1971 && year < 1981
                        WSV_folder = '/disk/r148/tfchengac/FLEXPART/';  %/ 10 yrs -> 5 TB (~500 Gb a year)
                    end
                end
            end

        elseif ismember(trajdirect, {'fwd'})
            if isequal(forcing, 'EA')
                WSV_folder = '/disk/r149/tfchengac/FLEXPART/';

            elseif isequal(forcing, 'CERA')
                if RHc_dqc_scheme == 1
                    if from_basin == 1
                        WSV_folder = '/disk/r034/tfchengac/FLEXPART/';   
                    else
                        WSV_folder = '/disk/r037/tfchengac/FLEXPART/';  
                    end
                else
                    WSV_folder = '/disk/r149/tfchengac/FLEXPART/'; 
                end
            end
        end

        if store_WSV_yearly
            suffix_yearly = strcat('_', num2str(year));
        else
            suffix_yearly = '';
        end
        WSV_dir{y}     = strcat(WSV_folder, expmnt, '/prcssd_data', suffix_yearly,'/');
    end
end

%/ Set the FLEXPART Model Output Directory
if ~isempty(FLEXPART_folder)
    if ischar(FLEXPART_folder)
        FLEXPART_folder = {FLEXPART_folder};
    end
    FLEXPART_dir = repmat(FLEXPART_folder, length(year_list),1);
else
    FLEXPART_dir = cell(length(year_list),1);
    for y = 1:length(year_list)
        year = year_list(y);
        if isequal(forcing, 'EA')
            % FLEXPART_folder = strcat('/disk/r149/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
            FLEXPART_folder = strcat('/disk/r149/tfchengac/flexpart_output/', expmnt, '/');
        elseif isequal(forcing, 'CERA')
            if year >= 2001 && year <= 2010
                FLEXPART_folder = strcat('/disk/r037/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
            elseif year >= 1991 && year <= 2000
                FLEXPART_folder = strcat('/disk/r034/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
            elseif year >= 1981 && year <= 1990
                FLEXPART_folder = strcat('/disk/r012/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
            elseif year >= 1971 && year <= 1980
                FLEXPART_folder = strcat('/disk/r014/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
            end
        end
        FLEXPART_dir{y} = FLEXPART_folder;
    end
end
str_sharpcut = '';

%/ Save moisture tracking info for analysis (Recommended)
if save_tracking_info
    S                     = [];
    S.expmnt              = expmnt;
    S.stdate              = stdate;                             
    S.eddate              = eddate;                             
    S.output_res          = output_res;
    S.dt                  = dt;
    S.year_list           = year_list;
    S.optimal_rr          = optimal_rr;
    S.str_optimal         = str_optimal;
    S.RHc_dqc_scheme      = RHc_dqc_scheme;
    S.trajdirect          = trajdirect;
    S.maxtraj_day         = maxtraj_day;
    S.str_remark          = str_remark;
    S.str_traj_rm_jump    = str_traj_rm_jump;
    S.str_BLH_factor      = str_BLH_factor;
    S.str_src_z           = str_src_z;
    S.dt_slct             = dt_slct;
    S.data_folder         = data_folder;
    S.plotting_folder     = plotting_folder;
    S.WSV_dir             = WSV_dir;
    S.str_expmntinfo      = str_expmntinfo;
    S.str_sharpcut        = str_sharpcut;
    S.basin_catalog       = basin_catalog;
    S.str_domain_trajfile = str_domain_trajfile;
    S.str_domain          = str_domain;

    tracking_info_filename = sprintf('tracking_info_%s_%s_%.2f_%d_%d-%d_%.2f_%d_%d.mat',...
                                     expmnt, trajdirect, output_res, dt, year_list(1),year_list(end),optimal_rr,from_basin,RHc_dqc_scheme);
    tracking_info_fullpath = fullfile(data_folder, tracking_info_filename);
    save(tracking_info_fullpath, 'S', '-v7.3');
    fprintf('*** Tracking info is saved: %s ***\n', tracking_info_fullpath);
end

%==========================================================================
%==================== Main program starts from here =======================
%==========================================================================
% try 
    for y = 1:length(year_list)
        year             = year_list(y);
        st_time          = datetime('now');
        FLEXPART_folder  = FLEXPART_dir{y};

        %==== Always read the header file ====%
        nest             = 0;
        readp            = 1; %/ read release points (0/1)
        calcarea         = 1;
        [header, ~]   = flex_header(FLEXPART_folder, nest, readp, calcarea); %/ Postprocessing tool from flexpart
        header.dates_str = unique(string(importdata(fullfile(FLEXPART_folder, 'dates'))));    %/ Read dates file in string, use unique to remove redundant dates (possibly due to MPI parallel)
        header.dates_dt  = datetime(header.dates_str,'InputFormat','yyyyMMddHHmmss', 'Format', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');

        ind = find(ismember(string(header.dates), string(header.ireleasestart)));
        header.ireleasestart_dt = header.dates_dt(ind, :);

        ind = find(ismember(string(header.dates), string(header.ireleaseend)));
        header.ireleaseend_dt = header.dates_dt(ind, :);             
        area = calc_grid_area('lon', lon, 'lat', lat(2:end-1));
        
        %/ Load the FLEXPART run dates
        run_dates_dt = header.dates_dt;
        run_dates_yyyymmddHHMM = datetime2int(run_dates_dt, 'yyyymmddHHMM');
        run_dates_yyyymmdd = yyyymmdd(run_dates_dt);

        %/ Subset the dates for the year(t) when tracking is requested
        if ~isempty(stdate) && ~isempty(eddate)
            if ismember(numel(num2str(stdate)), [3,4]) && ismember(numel(num2str(eddate)), [3,4])  %/ If given in mmdd
                %/ Ensuring bwd tracking always start from a date of the year (e.g., Jan 1), 
                %/ since FLEXPART run dates began from Dec 1 the prev year.
                stdate_dt = int2datetime(year*1e4+stdate, 'yyyyMMdd'); 
                eddate_dt = int2datetime(year*1e4+eddate, 'yyyyMMdd');          
            else
                stdate_dt = int2datetime(stdate, 'yyyyMMdd');
                eddate_dt = int2datetime(eddate, 'yyyyMMdd');
            end
            event_dates = double(date_array_gen('st_year', stdate_dt.Year, 'st_month', stdate_dt.Month, 'st_day', stdate_dt.Day,...
                                                'ed_year', eddate_dt.Year, 'ed_month', eddate_dt.Month, 'ed_day', eddate_dt.Day));

            t_ind = findismember_loop(run_dates_yyyymmdd, event_dates);
            t_st = t_ind(1);
            t_ed = t_ind(end);
        else
            %/ If stdate & eddate not given, tracking for all FLEXPART run dates available
            t_st = maxtraj_ts + 1;
            t_ed = length(run_dates_dt); 
        end

        if sharpcut && ismember(trajdirect, {'fwd'})
            %/ [IMPORTANT] for forward tracking, start from XX days before the event
            date_st              = run_dates_dt(t_st) - days(maxtraj_day); 
            date_begins_sharpcut = run_dates_dt(t_ed) - days(maxtraj_day); 

            %/ Update t_st 
            t_st       = find(run_dates_yyyymmdd == date_st.Year*1e4 + date_st.Month*1e2 + date_st.Day, 1);  %/ Find the first index of the date
            t_sharpcut = find(run_dates_yyyymmdd == date_begins_sharpcut.Year*1e4 + date_begins_sharpcut.Month*1e2 + date_begins_sharpcut.Day, 1, 'last');  %/ Find the last index of the date since it is from the last index of the end date
        end

        if isempty(t_st)  error('t_st returns empty! Check your code!'); end
        if isempty(t_ed)  error('t_ed returns empty! Check your code!'); end

        %/ split the job
        t_split = floor((t_ed - t_st)/Njob);
        if job_id == Njob
            t_list = (t_st + t_split*(job_id-1)):t_ed; %/ to workaround the uneven job split.
        else
            t_list = (1:t_split) + t_st - 1 + t_split*(job_id-1);
        end

        if test_mode      
            fprintf('=============================\n');
            fprintf('===== Test Mode %d is On =====\n', test_mode);
            fprintf('=============================\n');
            t_list = t_list(1);                         
        end

        if reverse_loop
            t_list = flip(t_list);
        end

        disp(run_dates_dt(t_list,:))
        fprintf('============================================\n');
        fprintf('=== NOTE: RHc_dqc_scheme %02d is selected. ===\n', RHc_dqc_scheme)
        fprintf('===       Njob = %d, job_id = %d           ===\n', Njob, job_id)
        fprintf('============================================\n');

        if ismember(trajdirect, {'bwd'})
            %/ Backward Tracking Timesteps
            trajtime_loop_direct = 0:-1:-(maxtraj_ts-1);                       %/ backward traj time loop
            trajtime             = trajtime_loop_direct*dt_slct;
            for t = t_list   
                %/ [IMPORTANT]
                if t == length(run_dates_dt)
                    fprintf('*** NOTE: The last time step of the experiment is skipped cos the output vars contain NaNs only. ***\n');
                    continue;
                end
                slct_one_date_dt = run_dates_dt(t);
                clear domfill; %/ release memory
                fprintf('\n*** reload_traj   = %d *** \n',   reload_traj);
                fprintf('*** overwrite_watersip_waterdrip = %d *** \n', overwrite_watersip_waterdrip);
                
                if from_basin == -1
                    domain = 1;
                elseif from_basin == 0
                    domain = 2;
                else
                    domain = 3;
                end
                
                %===== Set filepaths =====%
                %/ NOTE: different str_domain for traj file and the other files
                traj_suffix = strcat(str_RHc_dqc, '_', trajdirect,'_',...
                                num2str(maxtraj_day), 'd', str_sharpcut, str_optimal_traj, '_', num2str(dt_slct), 'h_', str_domain_trajfile,'_', str_src_z, str_BLH_factor, str_remark,...
                                datestr(run_dates_dt(t,:), 'yyyymmddHHMM'), '.mat');
            
                suffix = strcat(str_RHc_dqc, '_', trajdirect,'_',...
                                num2str(maxtraj_day), 'd', str_sharpcut, str_optimal, '_', num2str(dt_slct), 'h_', str_domain,'_', str_traj_rm_jump, str_src_z, str_BLH_factor, str_remark,...
                                datestr(run_dates_dt(t,:), 'yyyymmddHHMM'), '.mat');

                % prcssd_dir = strcat(WSV_folder, expmnt, '/prcssd_data_', num2str(year),'/');
                prcssd_dir = WSV_dir{y};
                
                [status, msg] = mkdir(char(prcssd_dir));  
                if status == 0   error(msg);   end
                traj_fullpath             = strcat(prcssd_dir,  'traj_',         traj_suffix);

                %/ for bwd mode
    %             BL_uptake_fullpath        = strcat(prcssd_dir,  'BL_uptake_',         suffix);
    %             uptake_fullpath           = strcat(prcssd_dir,  'uptake_',            suffix);
    %             EminusP_LA_fullpath       = strcat(prcssd_dir,  'EminusP_LA_',        suffix);
                BL_Pm_fullpath            = strcat(prcssd_dir,  'BL_Pm_',             suffix);
                Pm_fullpath               = strcat(prcssd_dir,  'Pm_',                suffix);
                Pm_AR_fullpath            = strcat(prcssd_dir,  'Pm_AR_',             suffix);
                watersip_rr_fullpath      = strcat(prcssd_dir,  'watersip_',          suffix);
                optimal_trajtime_fullpath = strcat(prcssd_dir,  'optimal_trajtime_',  suffix);
                CWRT_fullpath             = strcat(prcssd_dir,  'CWRT_',              suffix);
                P_LA_fullpath             = strcat(prcssd_dir,  'P_LA_',              suffix);
                rr_tot_L_fullpath         = strcat(prcssd_dir,  'rr_tot_L_',          suffix);
                rr_tot_NLL_fullpath       = strcat(prcssd_dir,  'rr_tot_NLL_',        suffix);
                rr_tot_NLO_fullpath       = strcat(prcssd_dir,  'rr_tot_NLO_',        suffix);
                RH2_fullpath              = strcat(prcssd_dir,  'RH2_',               suffix);

                if from_basin == 0 || from_basin == -1  LA_3D_fullpath = '';
                else                                    LA_3D_fullpath = strcat(prcssd_dir,  'LA_3D_', [basin_catalog.name], '_', suffix); %/ each hostpot has one file.
                end
                %===== Skip if all data exist (save time without loading traj) =======%
    %             && isfile(BL_uptake_fullpath) && isfile(BL_Pm_fullpath) && isfile(uptake_fullpath) ...
                all_file_exist = (isfile(traj_fullpath) && isfile(Pm_fullpath) && isfile(BL_Pm_fullpath) && (~isempty(AR) && isfile(Pm_AR_fullpath))...
                                   && isfile(watersip_rr_fullpath) && isfile(optimal_trajtime_fullpath) && isfile(CWRT_fullpath) && isfile(P_LA_fullpath) ...
                                   && (output_rr_map == 0 || isfile(rr_tot_L_fullpath) && isfile(rr_tot_NLL_fullpath) && isfile(rr_tot_NLO_fullpath)) ...
                                   && isfile(RH2_fullpath) && (output_LA_3D  == 0 || all(isfile(LA_3D_fullpath))) ...
                                   && reload_traj == 0 && overwrite_watersip_waterdrip == 0);
    
                if all_file_exist
                    fprintf('!!! [y = %d, t = %d/%d] All traj, uptake and WaterSip files exist, and no overwriting is needed. !!! \n', year, t, max(t_list))
                    continue;
                end

                %===== Check if traj data exists, load from it if overwrite = 0 =======%
                traj_exist = 0;
                if reload_traj == 0
                    if isfile(traj_fullpath)
                        tic
                        fprintf('!!! [y = %d, t = %d/%d] traj data are found: %s. Loading traj data... !!! \n', year, t, max(t_list), traj_fullpath)
                        load(traj_fullpath, 'domfill');   
                        traj_exist = 1;
                        toc;
                    end
                end

                %===== Derive Traj (slct_traj_parfor) =====%
                if traj_exist == 0
                    tic
                    clear partoutput; %/ release memory
                    numpart_esti = 'unknown';

                    %/ preallocate struct field for parfor loop later.
                    domfill_parfor = struct('slct_id_list', cell(1,1),...
                                            'traj',     cell(1,length(trajtime_loop_direct)),...
                                            'q',        cell(1,length(trajtime_loop_direct)),...
                                            'BLH',      cell(1,length(trajtime_loop_direct)),...
                                            'T',        cell(1,length(trajtime_loop_direct)),...
                                            'topo',     cell(1,length(trajtime_loop_direct)),...
                                            'mass',     cell(1,length(trajtime_loop_direct)),...
                                            'RH2',      cell(1,length(trajtime_loop_direct)),...
                                            'dates_dt', cell(1,length(trajtime_loop_direct)));

                    trajtime_loop_direct_init = trajtime_loop_direct(1:2);
                    for cnt = 1:length(trajtime_loop_direct_init)                  %/ run the first 2 timesteps
                        j = trajtime_loop_direct_init(cnt);

                        fprintf('*** Loading traj ... ***\n')
                        disp(numpart_esti)
                        date_flag = datestr(run_dates_dt(t+j,:), 'yyyymmddHHMMSS');
                        if isequal(numpart_esti, 'unknown')
                            numpart_esti = readpart10(FLEXPART_folder, date_flag, numpart_esti); %/ If numpart is unknown, input 'unknown' to estimate it.
                        end

                        partoutput = readpart10(FLEXPART_folder, date_flag, numpart_esti);           
                        fprintf('*** [y = %d, t = %d/%d] [%.1f days %s from %s]: %d trajs on %s (%d/%d) are loaded. ***\n',...
                                year, t, max(t_list), maxtraj_day, trajdirect, run_dates_dt(t,:), numpart_esti, run_dates_dt(t+j,:), cnt, length(trajtime_loop_direct))

                        %======================== allocation ======================
                        if cnt == 1
                            %/ ====================================================================================
                            %/ Selection 1: *ONLY* handle those particles in the given domain
                            %               (e.g. from a box region or from the land grids)
                            %/ ====================================================================================
                            if domain == 3        %/ Select trajs over all the target regions (from_basin)
                                ind = [];
                                for top = 1:NBasin
                                    %/ Always convert into [-179, 180] as the only discontinuity will be in the ocean (dateline)
                                    traj_x = partoutput.xyz(:,1);
                                    traj_y = partoutput.xyz(:,2);
                                    traj_x(traj_x > 180) = traj_x(traj_x > 180) - 360;
                                    
                                    bndry_data = basin_catalog(top).bndry;
                                    bndry_x = bndry_data(:,1);
                                    bndry_x(bndry_x > 180) = bndry_x(bndry_x > 180) - 360;
                                    bndry_data(:,1) = bndry_x;

                                    try
                                        [in, ~] = inpoly2([traj_x, traj_y], bndry_data);                          %/ inpoly2 is 600xx faster than inpolygon!! 
                                    catch 
%                                         warning(E.identifier, 'Error caught!: %s. Proceed with inpolygon().', E.message)
                                        [in, ~] = inpolygon(traj_x, traj_y, bndry_data(:,1), bndry_data(:,2));    %/ sometimes inpoly2 could run into a bug, then we use inpolygon.
                                    end
                                    ind = [ind; find(in == 1)];
                                end
                                ind = unique(ind);  %/ remove redundant traj indices
                                
                                if isempty(ind)
                                    %/ NOTE: do NOT skip the loop! Keep it running to save empty files! Or the program will always rerun..
                                    fprintf('!!! No trajs from all basins. !!!\n'); 
                                end
                                
                            elseif domain == 2     %/ Select trajs over global land
                                ind = show_land_or_ocean_hydrosheds('pos_lon_array', partoutput.xyz(:,1), 'pos_lat_array', partoutput.xyz(:,2),...
                                                         'lon_grids', lon, 'lat_grids', lat, 'land_or_ocean', 1,...
                                                         'cond_land', cond_land, 'cond_ocean', cond_ocean);

                            elseif domain == 1     %/ Select trajs globally
                                ind = 1:length(partoutput.xyz(:,1)); %/ simply all
%                                 ind = find(partoutput.xyz(:,1) >= src_lon(1) & partoutput.xyz(:,1) <= src_lon(2) & ...
%                                            partoutput.xyz(:,2) >= src_lat(1) & partoutput.xyz(:,2) <= src_lat(2));
                                       
                            else
                                error('Code not set for ''domain'' == %d!', domain)
                            end
                            domfill_parfor(1).slct_id_list = partoutput.npoint(ind);                

                            %/ Preallocation (will help the codes run faster)
                            domfill_parfor(cnt).traj     = partoutput.xyz(ind, :);
                            domfill_parfor(cnt).q        = partoutput.vars(ind, 1)*1000;
                            domfill_parfor(cnt).BLH      = partoutput.vars(ind, 2);
                            domfill_parfor(cnt).T        = partoutput.vars(ind, 3);
                            domfill_parfor(cnt).topo     = partoutput.vars(ind, 4);
                            domfill_parfor(1).mass       = partoutput.xmass(ind,1); 
                            domfill_parfor(cnt).dates_dt = run_dates_dt(t+j);

                        elseif cnt == 2 
                            ind = findismember_loop(partoutput.npoint, domfill_parfor(1).slct_id_list); %/ get indices based on id list (no auto sorting)
                            domfill_parfor(cnt).traj     = partoutput.xyz(ind, :);
                            domfill_parfor(cnt).q        = partoutput.vars(ind, 1)*1000;
                            domfill_parfor(cnt).BLH      = partoutput.vars(ind, 2);
                            domfill_parfor(cnt).T        = partoutput.vars(ind, 3);
                            domfill_parfor(cnt).topo     = partoutput.vars(ind, 4);
                            domfill_parfor(cnt).dates_dt = run_dates_dt(t+j);

                            %/ ====================================================================================
                            %/ Selection 2: *ONLY* handle trajs ended up with RH > RHc and dq > dqc (much less than the total no. of parcels)
                            %/              --> update domfill.slct_id_list (speedup by > 1 s for each loop!)
                            %/ ====================================================================================
                            %/ This updated fn has a more consistent coding.
                            domfill_parfor = slct_traj_parfor('trajdirect', trajdirect, 'domfill', domfill_parfor,...
                                                              'RHc', RHc, 'dqc', dqc, 'str_RHc_dqc', str_RHc_dqc, 'str_remark', str_remark,...
                                                              'RHc_map_lon', RHc_map_lon, 'RHc_map_lat', RHc_map_lat, 'slct_one_date_dt', slct_one_date_dt, 'BLH_factor', BLH_factor); 

                        end
                        clear partoutput; %/ release memory 
                    end

                    %======== Parfor loop over the rest traj time !!! ============%
                    if isempty(gcp('nocreate')) && ~isempty(NumWorkers_IO) %/ if set worker number
                        % parpool('Threads', NumWorkers_IO)  %/ It cannot work with the reading
                        parpool('Processes', max([NumWorkers_IO,NumWorkers_HPC]))
                    end

                    slct_id_list = domfill_parfor(1).slct_id_list;
                    parfor (k = cnt+1:length(trajtime_loop_direct), NumWorkers_IO)    %/ parfor for the rest (starting from cnt+1) - much faster.
                    % for k = cnt+1:length(trajtime_loop_direct)    %/ for testing 
                        %/ Broadcast variable
                        slct_date_dt_bc         = run_dates_dt;
                        t_list_bc               = t_list;
                        trajtime_loop_direct_bc = trajtime_loop_direct;
                        j                       = trajtime_loop_direct_bc(k);
                        
                        fprintf('*** Loading traj ... ***\n')
                        date_flag = datestr(slct_date_dt_bc(t+j,:), 'yyyymmddHHMMSS');

                        numpart_esti = 'unknown';
                        numpart_esti = readpart10(FLEXPART_folder, date_flag, numpart_esti);    
                        partoutput = readpart10(FLEXPART_folder, date_flag, numpart_esti);            
                        fprintf('*** [y = %d, t = %d/%d] [%.1f days %s from %s]: %d trajs on %s (%d/%d) are loaded. *** \n',...
                                year,  t, max(t_list_bc), maxtraj_day, trajdirect, slct_date_dt_bc(t,:), numpart_esti, slct_date_dt_bc(t+j,:), k, length(trajtime_loop_direct_bc))

                        %======================== allocation ======================
                        ind = findismember_loop(partoutput.npoint, slct_id_list);  %/ get indices based on id list (no auto sorting)
                            
                        if length(ind) ~= length(slct_id_list)
                            error('[k == %d] Not all of the queried particle IDs are found in partoutput.npoint! Check FLEXPART output!', k);
                        end

                        domfill_parfor(k).traj     = partoutput.xyz(ind, :);  
                        domfill_parfor(k).q        = partoutput.vars(ind,1) * 1000; %/ from kg/kg to g/kg.
                        domfill_parfor(k).BLH      = partoutput.vars(ind,2);
                        domfill_parfor(k).T        = partoutput.vars(ind,3);
                        domfill_parfor(k).topo     = partoutput.vars(ind,4);
                        domfill_parfor(k).dates_dt = run_dates_dt(t+j);
                    end
                    domfill.traj     = cat(3, domfill_parfor(:).traj);             %/ change nonscalar domfill_parfor to scalar domfill (compatible for other codes) 
                    domfill.q        = cat(2, domfill_parfor(:).q);                %/ and concatenate the numbered fields 
                    domfill.BLH      = cat(2, domfill_parfor(:).BLH);
                    domfill.topo     = cat(2, domfill_parfor(:).topo);
                    domfill.dates_dt = cat(2, domfill_parfor(:).dates_dt);
                    domfill.T        = cat(2, domfill_parfor(:).T);   
                    domfill.mass     = domfill_parfor(1).mass;
                    domfill.RH2      = domfill_parfor(1).RH2;
%                     domfill_parfor   = [];   %/ release memory

                    domfill.trajdirect   = trajdirect;
                    domfill.src_lon      = src_lon;
                    domfill.src_lat      = src_lat;
                    domfill.src_z        = src_z;
                    domfill.dt_slct      = dt_slct;
                    domfill.maxtraj_ts   = maxtraj_ts;
                    domfill.maxtraj_day  = maxtraj_day;
                    domfill.trajtime     = trajtime;

                    %===== Save traj data =====%
                    if save_traj
                        fprintf('*** [y = %d, t = %d/%d] Saving traj into: %s ... *** \n', year, t, max(t_list), traj_fullpath)
                        save(traj_fullpath, 'domfill', '-v7.3');             %/ save/overwrite
                    end
                    toc
                end

                %===== WaterSip Diagnostic =====%
                if overwrite_watersip_waterdrip == 2
                    fprintf('!!! [y = %d, t = %d/%d] uptake and watersip data are skipped as per user''s request. !!!\n', year, t, max(t_list))
                else
                    tic
                    %===== Check if the uptake_map data exists, load from it if overwrite = 0 =======%
    %                 uptake_exist = 0; BL_uptake_exist = 0; 
                    BL_Pm_exist = 0; Pm_exist = 0; Pm_AR_exist = 0; watersip_exist = 0; optimal_trajtime_exist = 0; CWRT_exist = 0;
                    P_LA_exist = 0; rr_tot_L_exist = 0; rr_tot_NLL_exist = 0; rr_tot_NLO_exist = 0; RH2_exist = 0;
                    LA_3D_exist = zeros(1, NBasin);
                    if overwrite_watersip_waterdrip == 0
                        %/ bwd mode
    %                     uptake_exist            = isfile(uptake_fullpath);                                                       
    %                     BL_uptake_exist         = isfile(BL_uptake_fullpath);                         
                        BL_Pm_exist             = isfile(BL_Pm_fullpath);       
                        Pm_exist                = isfile(Pm_fullpath);     
                        Pm_AR_exist             = isfile(Pm_AR_fullpath);     
                        watersip_exist          = isfile(watersip_rr_fullpath);                      
                        optimal_trajtime_exist  = isfile(optimal_trajtime_fullpath);  
                        CWRT_exist              = isfile(CWRT_fullpath);  
                        P_LA_exist              = isfile(P_LA_fullpath);                             
                        rr_tot_L_exist          = isfile(rr_tot_L_fullpath);                         
                        rr_tot_NLL_exist        = isfile(rr_tot_NLL_fullpath);                        
                        rr_tot_NLO_exist        = isfile(rr_tot_NLO_fullpath); 
                        RH2_exist               = isfile(RH2_fullpath); 
                        LA_3D_exist             = isfile(LA_3D_fullpath);                               
                    end

                    if output_rr_map
                        if ~(from_basin == 0 || from_basin == -1)
                            error('Current code is not set for basin-wise tracking.');
                        else
                            all_WSV_file_exist = (BL_Pm_exist && Pm_exist && watersip_exist && optimal_trajtime_exist && CWRT_exist && ...
                                                 P_LA_exist && rr_tot_L_exist && rr_tot_NLL_exist && rr_tot_NLO_exist && ...
                                                 RH2_exist && all(LA_3D_exist));
                        end
                    else
                        all_WSV_file_exist = (BL_Pm_exist && Pm_exist && watersip_exist && optimal_trajtime_exist && CWRT_exist && ...
                                                P_LA_exist && RH2_exist && all(LA_3D_exist));
                    end
                    
                    %/ Either overwrite is required or any watersip files not exists, we will run the algorithm.
                    if overwrite_watersip_waterdrip || all_WSV_file_exist == 0
                        %/ [perform a hotspot loop here if from_basin == 1, 2 or 3]
                        if from_basin ~= 0
                            fprintf('*** [y = %d, t = %d/%d] Start WaterSip algorithm for each hotspot... *** \n', year, t, max(t_list))
                            BL_Pm_map                 = cell(NBasin, 1);    %/ Since in WaterSip, it is multipled by area without the two poles.
                            Pm_map                    = cell(NBasin, 1);    %/ Since in WaterSip, it is multipled by area without the two poles.
                            P_LA_map                  = cell(NBasin, 1);    %/ Since in WaterSip, it is multipled by area without the two poles.
                            mean_optimal_trajtime_map = cell(NBasin, 1);
                            mean_CWRT_map             = cell(NBasin, 1);
                            rr_tot_L_map              = cell(NBasin, 1);    
                            rr_tot_NLL_map            = cell(NBasin, 1);    
                            rr_tot_NLO_map            = cell(NBasin, 1);  
                            RH2_map                   = cell(NBasin, 1);
                            watersip_rr               = cell(NBasin, 1);
%                             traj_x = domfill.traj(:,1,1);
%                             traj_y = domfill.traj(:,2,1);
%                             traj_x(traj_x > 180) = traj_x(traj_x > 180) - 360;
                                                        
                            for top = 1:NBasin
                                basin_catalog_bc = basin_catalog;
                                domfill_basin    = domfill;
                                traj_x = domfill_basin.traj(:,1,1);
                                traj_y = domfill_basin.traj(:,2,1);

                                %/ Always convert into [-179, 180] as the only discontinuity will be in the ocean (dateline)
                                traj_x(traj_x > 180) = traj_x(traj_x > 180) - 360;

                                bndry_data       = basin_catalog_bc(top).bndry;   %<<-- important!
                                basin_name       = basin_catalog_bc(top).name;
                                bndry_x          = bndry_data(:,1);
                                bndry_x(bndry_x > 180) = bndry_x(bndry_x > 180) - 360;
                                bndry_data(:,1) = bndry_x;

                                lon_m179_180_bc  = lon_m179_180;
                                lon_bc           = lon;
                                lat_bc           = lat;
                                t_list_bc        = t_list;
                                try
                                    [in, ~] = inpoly2([traj_x, traj_y], bndry_data);                          %/ inpoly2 is 600xx faster than inpolygon!! 
                                catch 
%                                     warning(E.identifier, 'Error caught!: %s. Proceed with inpolygon().', E.message)
                                    [in, ~] = inpolygon(traj_x, traj_y, bndry_data(:,1), bndry_data(:,2));    %/ sometimes inpoly2 could run into a bug, then we use inpolygon.
                                end
                                ind = find(in == 1);

                                if isempty(ind)
                                    %/ NOTE: do NOT skip the loop! Keep it running to save empty files! Or the program will always rerun..
                                    fprintf('!!! No trajs from basin #%d. Skipped it.!!!\n', top); 
                                end

                                domfill_basin.traj  = domfill_basin.traj(ind,:,:);
                                domfill_basin.q     = domfill_basin.q   (ind,:);
                                domfill_basin.BLH   = domfill_basin.BLH (ind,:);
                                domfill_basin.topo  = domfill_basin.topo(ind,:);
                                domfill_basin.T     = domfill_basin.T   (ind,:);
                                domfill_basin.mass  = domfill_basin.mass(ind,:);
                                domfill_basin.RH2   = domfill_basin.RH2 (ind,:);

%                                 %/ Check if there are points outside TP
%                                 %/ (could be due to the non-closure of bndry_data)
%                                 if ~isempty(ind) && from_basin == 4
%                                     x_bc = squeeze(domfill_basin.traj(:,1,1));
%                                     y_bc = squeeze(domfill_basin.traj(:,2,1));
%                                     check = find(x_bc > 82 & x_bc < 85 & y_bc > 38 & y_bc < 42);
%                                     if ~isempty(check)
%                                         error('Detected %d points outside TP! Check inpoly2!', length(check));
%                                     end
%                                 end

                                %/ An extented area covering the target basin (+-40 in lon, +-30 in lat)
                                basin_lon_min = min(bndry_data(:,1)) - 40;
                                basin_lon_max = max(bndry_data(:,1)) + 40;
                                basin_lat_min = min(bndry_data(:,2)) - 30;
                                basin_lat_max = max(bndry_data(:,2)) + 30;

                                if basin_lon_max > 360  error('code not set up for this!');   end
                                if basin_lon_min < 0
                                    basin_lon_range = lon_m179_180_bc(lon_m179_180_bc >= basin_lon_min & lon_m179_180_bc <= basin_lon_max);
                                else
                                    basin_lon_range = lon_bc(lon_bc >= basin_lon_min & lon_bc <= basin_lon_max);
                                end
                                basin_lat_range     = lat_bc(lat_bc >= basin_lat_min & lat_bc <= basin_lat_max);
                                if isempty(basin_lon_range)    error('basin_lon_range is empty! Check if lon is from -180 to 180!');  end
                                if isempty(basin_lat_range)    error('basin_lat_range is empty!');  end
                                z_intval         = 200;
                                basin_z_range    = (0:z_intval:10000)+z_intval/2;   %/ 51 blocks

                                fprintf('*** [y = %d, t = %d/%d, basin = %d/%d, %s] Start WaterSip algorithm... *** \n', year, t, max(t_list_bc), top, length(basin_catalog_bc), basin_name{:})
                                [~, Pm_map_basin, ~, BL_Pm_map_basin, watersip_rr_basin, watersip_stacktable_basin,...
                                 mean_optimal_trajtime_map_basin, mean_CWRT_map_basin, P_LA_map_basin, ~, ~, ~, RH2_map_basin, LA_3D, Pm_AR_map] ...
                                        = WaterSip('domfill', domfill_basin, 'lon', lon, 'lat', lat, 'dqc', dqc, 'BLH_factor', BLH_factor, 'area', area,...
                                                   'cond_land', cond_land, 'cond_ocean', cond_ocean, 'traj_rm_jump', traj_rm_jump, 'rm_when_dz', rm_when_dz, 'NumWorkers', NumWorkers_HPC, 'optimal_rr', optimal_rr,...
                                                   'output_rr_map', output_rr_map, 'output_LA_3D', output_LA_3D, 'basin_lon_range', basin_lon_range, 'basin_lat_range', basin_lat_range, 'basin_z_range', basin_z_range, 'basin_name', basin_name,...
                                                   'AR', AR);

                                BL_Pm_map{top}                 = BL_Pm_map_basin;
                                Pm_map{top}                    = Pm_map_basin;
                                mean_optimal_trajtime_map{top} = mean_optimal_trajtime_map_basin;
                                mean_CWRT_map{top}             = mean_CWRT_map_basin;
                                P_LA_map{top}                  = P_LA_map_basin;
                                RH2_map{top}                   = RH2_map_basin;
                                watersip_rr{top}               = watersip_rr_basin;            %/ final recycling ratios for all trajs from the hotspot
                                watersip_stacktable            = watersip_stacktable_basin';   %/ a full table just for double check.
                                
                                %/ save LA_3D for each basin to save loading time afterwards!
                                if savemat
                                    if LA_3D_exist(top) == 0 && output_LA_3D   par_save(LA_3D_fullpath{top},  'LA_3D');  end  
                                end
                            end
                            
                            %/ Concatenate cell arrays along the 3rd dim
                            BL_Pm_map = cat(3, BL_Pm_map{:});
                            Pm_map = cat(3, Pm_map{:});
                            mean_optimal_trajtime_map = cat(3, mean_optimal_trajtime_map{:});
                            mean_CWRT_map = cat(3, mean_CWRT_map{:});
                            P_LA_map = cat(3, P_LA_map{:});
                            RH2_map = cat(3, RH2_map{:});
                            
                        else
                            fprintf('*** [y = %d, t = %d/%d] Start WaterSip algorithm... *** \n', year, t, max(t_list))
                            [~, Pm_map, ~, BL_Pm_map, watersip_rr, watersip_stacktable,...
                                 mean_optimal_trajtime_map, mean_CWRT_map, P_LA_map, rr_tot_L_map, rr_tot_NLL_map, rr_tot_NLO_map, RH2_map, ~, Pm_AR_map] ...
                                        = WaterSip('domfill', domfill, 'lon', lon, 'lat', lat, 'dqc', dqc, 'BLH_factor', BLH_factor, 'area', area,...
                                                   'cond_land', cond_land, 'cond_ocean', cond_ocean, 'traj_rm_jump', traj_rm_jump, 'rm_when_dz', rm_when_dz, 'NumWorkers', NumWorkers_HPC, 'optimal_rr', optimal_rr,...
                                                   'output_rr_map', output_rr_map, 'output_LA_3D', output_LA_3D,...
                                                   'basin_lon_range', [], 'basin_lat_range', [], 'basin_z_range', []);
                        end

                        if savemat
                            %/ WARNING: when overwrite_watersip_waterdrip = 1, all the below WSV will be overwritten.
                            if BL_Pm_exist == 0                                save(BL_Pm_fullpath,              'BL_Pm_map',                 '-v7.3'); end  %/ only for general validation.
                            if Pm_exist == 0                                   save(Pm_fullpath,                 'Pm_map',                    '-v7.3'); end
                            if Pm_AR_exist == 0 && ~isempty(AR)                save(Pm_AR_fullpath,              'Pm_AR_map',                 '-v7.3'); end
                            if optimal_trajtime_exist == 0                     save(optimal_trajtime_fullpath,   'mean_optimal_trajtime_map', '-v7.3'); end
                            if CWRT_exist == 0                                 save(CWRT_fullpath,               'mean_CWRT_map',             '-v7.3'); end
                            if P_LA_exist == 0                                 save(P_LA_fullpath,               'P_LA_map',                  '-v7.3'); end
                            if watersip_exist == 0                             save(watersip_rr_fullpath,        'watersip_rr',               '-v7.3'); end
                            if rr_tot_L_exist == 0    && output_rr_map         save(rr_tot_L_fullpath,           'rr_tot_L_map',              '-v7.3'); end
                            if rr_tot_NLL_exist == 0  && output_rr_map         save(rr_tot_NLL_fullpath,         'rr_tot_NLL_map',            '-v7.3'); end
                            if rr_tot_NLO_exist == 0  && output_rr_map         save(rr_tot_NLO_fullpath,         'rr_tot_NLO_map',            '-v7.3'); end
                            if RH2_exist == 0                                  save(RH2_fullpath,                'RH2_map',                   '-v7.3'); end
                            fprintf('*** [y = %d, t = %d/%d] WaterSip products saved under %s *** \n', year, t, max(t_list), prcssd_dir)
                        end
                        % plot_contfmap('contf_data', Pm_AR_map, 'contf_lon', lon, 'contf_lat', lat(2:end-1))
                    else
                        fprintf('*** [y = %d, t = %d/%d] All the WaterSip output files already exist. Skipped the loop. *** \n', year, t, max(t_list))
                    end
                    toc
                end
                
                %===== Plot trajectories [For demo ONLY] =====%
                if plot_traj
                    %=================================================
                    fprintf('*** Plotting trajs... ***\n');
                    close all
                    savefig             = 0;               %/ <- mind this!
                    cbar_location       = 'eastoutside'; 

                    %=================================================
                    slct_intvl          = 10;               %/ Select the interval to subset the trajs
                    slct_max_trajtime   = 24*10;  %/ Select the backtrack time (in hours) of the trajs
                    % slct_max_trajtime   = 24*maxtraj_day;  %/ Select the backtrack time (in hours) of the trajs
                    %=================================================

                    ind_slct_traj    = 1:slct_intvl:size(domfill.traj,1);  
                    traj_data        = permute(domfill.traj(ind_slct_traj,:,:), [1,3,2]);  %/ permute to [NoOfTraj, TrajTime, XYZ]
                    traj_data(:,:,3) = traj_data(:,:,3) + domfill.topo(ind_slct_traj,:);   %/ z-pos + topo = altitude

                    %/ Select the max backtracking time of trajs
                    traj_time             = 1:size(traj_data, 2); 
                    ind_slct_max_trajtime = 1:floor(slct_max_trajtime/dt);

                    traj_time = traj_time(ind_slct_max_trajtime);
                    traj_data = traj_data(:,ind_slct_max_trajtime,:);
                    size(traj_data)
                    min(traj_data(:,:,3), [], 'all')  %/ the minimum altitude (sometimes be -ve, e.g., -10 m)
                    max(traj_data(:,:,3), [], 'all')

                    %/ Set the colorbar
                %     traj_levels = -500:500:12000;  %/ altitude of the particles.
                    traj_levels = -500:500:7000;  %/ altitude of the particles.
                    traj_colmap      = viridis(length(traj_levels)-1);
                    traj_unit   = 'm'; 

                    %/ Map setting
                    grid_mode       = 0;
                    fontsize        = 14;
                    cbar_fontsize   = fontsize;
                    linewi          = 1.5;
                    coast_wi        = 2;
                %     coast_col       = 'k';
                    coast_col       = 'none';
                    coast_patch_col = [240 239 221]./255;
                    backcolor       = [151 182 226]./255;
                    fig_fmt         = 'pdf';

                    if ismember(from_basin, [7,8])
                        map_lon_lower = 20; map_lon_upper = 120; map_lat_lower = -25; map_lat_upper = 40; markersize = 1; %/ Avoid south / north pole -> problematic in drawing trajs.
                    elseif ismember(from_basin, 9)
                        map_lon_lower = 40; map_lon_upper = 280; map_lat_lower = -70; map_lat_upper = 10; markersize = 1; %/ Avoid south / north pole -> problematic in drawing trajs.
                    elseif ismember(from_basin, 10)
                        map_lon_lower = -40; map_lon_upper = 40; map_lat_lower = 30; map_lat_upper = 80; markersize = 1; %/ Avoid south / north pole -> problematic in drawing trajs.
                    else
                        map_lon_lower = -179; map_lon_upper = 180; map_lat_lower = -80; map_lat_upper = 80; markersize = 1; %/ Avoid south / north pole -> problematic in drawing trajs.
                    end

                    savepath = [];
                    if savefig
                        date_yyyymmddHHMM = datetime2int(domfill.dates_dt(1), 'yyyymmddHHMM');
                        savepath = strcat(plotting_folder, sprintf('traj_%d_intvl%d_maxtrajtime%dh', date_yyyymmddHHMM, slct_intvl, slct_max_trajtime));
                    end
                    
                    draw_cbar_only = 0;
                    plot_contfmap('map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper,...
                                  'coast_wi', coast_wi, 'fontsize', fontsize, 'create_fig', 1, 'grid_mode', grid_mode,...
                                  'traj_data', traj_data, 'traj_time', traj_time, 'traj_levels', traj_levels, 'traj_colmap', traj_colmap, 'traj_unit', traj_unit,...
                                  'markersize', markersize, 'linewi', linewi, 'coast_col', coast_col, 'coast_patch_col', coast_patch_col, 'backcolor', backcolor,...
                                  'savepath', savepath, 'cbar_location', cbar_location, 'draw_cbar_only', draw_cbar_only, 'cbar_fontsize', cbar_fontsize,  'fig_fmt', fig_fmt);

                    % if t == t_list(1)
                    %     draw_cbar_only = 1;
                    %     plot_contfmap('map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper,...
                    %               'coast_wi', coast_wi, 'fontsize', fontsize, 'create_fig', 1, 'grid_mode', grid_mode,...
                    %               'traj_data', traj_data, 'traj_time', traj_time, 'traj_levels', traj_levels, 'traj_colmap', traj_colmap, 'traj_unit', traj_unit,...
                    %               'markersize', markersize, 'linewi', linewi, 'coast_col', coast_col, 'coast_patch_col', coast_patch_col, 'backcolor', backcolor,...
                    %               'savepath', savepath, 'cbar_location', cbar_location, 'draw_cbar_only', draw_cbar_only, 'cbar_fontsize', cbar_fontsize, 'fig_fmt', fig_fmt);       
                    % end
                end
            end
        end

        if ismember(trajdirect, {'fwd'}) 
            %/ Backward Tracking Timesteps
            trajtime_loop_direct = -1:1:maxtraj_ts-1;                          %/ forward traj time loop (starting from -1)
            trajtime             = trajtime_loop_direct*dt_slct;
            for t = t_list
                %/ Broadcast var
                trajtime_loop_direct_sharpcut = trajtime_loop_direct;
                trajtime_sharpcut             = trajtime;
                
                if sharpcut
                    if t > t_sharpcut 
                        t_cut = t - t_sharpcut - 1; %/ - 1 such that there will be three time steps to compute initial gain and Cf map.
                        str_sharpcut = sprintf('_sharpcut%d', run_dates_yyyymmddHHMM(t_ed));
                    else
                        t_cut = 0;
                    end
                    trajtime_loop_direct_sharpcut(end-t_cut+1:end) = []; %/ sharp cut the tracking length when approaching the end day
                    trajtime_sharpcut(end-t_cut+1:end) = [];
                end

                %/ [IMPORTANT]
                if t == length(run_dates_dt)
                    fprintf('*** NOTE: The last time step of the experiment is skipped cos the output vars contain NaNs only. ***\n');
                    continue;
                end
                clear domfill; %/ release memory
                fprintf('\n*** reload_traj   = %d *** \n',   reload_traj);
                fprintf('*** overwrite_watersip_waterdrip = %d *** \n', overwrite_watersip_waterdrip);
                slct_one_date_dt = run_dates_dt(t);

                %===== Set filepaths =====%
                %/ NOTE: different str_domain for traj file and the other files
                traj_suffix = strcat(str_RHc_dqc, '_', trajdirect,'_',...
                                num2str(maxtraj_day), 'd', str_sharpcut, str_optimal_traj, '_', num2str(dt_slct), 'h_', str_domain_trajfile,'_', str_src_z, str_BLH_factor, str_remark,...
                                datestr(run_dates_dt(t,:), 'yyyymmddHHMM'), '.mat');

                suffix = strcat(str_RHc_dqc, '_', trajdirect,'_',...
                                num2str(maxtraj_day), 'd', str_sharpcut, str_optimal, '_', num2str(dt_slct), 'h_', str_domain,'_', str_traj_rm_jump, str_src_z, str_BLH_factor, str_remark,...
                                datestr(run_dates_dt(t,:), 'yyyymmddHHMM'), '.mat');

                basin_name = [basin_catalog.name];
                basin_suffix = strcat(str_RHc_dqc, '_', trajdirect,'_',...
                                num2str(maxtraj_day), 'd', str_sharpcut, str_optimal, '_', num2str(dt_slct), 'h_', basin_name,'_', str_traj_rm_jump, str_src_z, str_BLH_factor, str_remark,...
                                datestr(run_dates_dt(t,:), 'yyyymmddHHMM'), '.mat');

                %/ dir to store processed data
                if yyyymmdd(run_dates_dt(t,:)) < year*1e4     %/ then modify the output dir for year (-1)
                    prcssd_dir = WSV_dir{y-1};
                else
                    prcssd_dir = WSV_dir{y};
                end
                [status, msg] = mkdir(char(prcssd_dir));
                if status == 0   error(msg);   end
                traj_fullpath             = strcat(prcssd_dir,  'traj_',         traj_suffix);

                %/ for fwd mode
                Cf_map_fullpath           = strcat(prcssd_dir,  'Cf_map_',            basin_suffix);
                Cf_map_dates_fullpath     = strcat(prcssd_dir,  'Cf_map_dates_',      suffix);
                optimal_trajtime_fullpath = strcat(prcssd_dir,  'optimal_trajtime_',  basin_suffix);

                all_fullpath = [Cf_map_fullpath, Cf_map_dates_fullpath, optimal_trajtime_fullpath];  %/ Bundle all the filepaths

                %===== Skip if all data exist and do not overwrite them =======%
                if all(isfile(all_fullpath)) &&  reload_traj == 0 && overwrite_watersip_waterdrip == 0
                    fprintf('!!! [y = %d, t = %d/%d] All traj, uptake and WaterDrip files exist, and no overwriting is needed. !!! \n', year, t, max(t_list))
                    continue;
                end

                %===== Check if traj data exists, load from it if overwrite = 0 =======%
                traj_exist = 0;
                if reload_traj == 0
                    if isfile(traj_fullpath)
                        fprintf('!!! [y = %d, t = %d/%d] traj data are found: %s. Loading traj data... !!! \n', year, t, max(t_list), traj_fullpath)
                        load(traj_fullpath);
                        traj_exist = 1;
                    end
                end

                if traj_exist == 0
                    tic
                    clear partoutput;   numpart_esti = 'unknown';
                    fwd = cell(NBasin, 1);

                    trajtime_loop_direct_init = trajtime_loop_direct_sharpcut(1:2);
                    for cnt = 1:length(trajtime_loop_direct_init)                  %/ run the first 2 timesteps
                        j = trajtime_loop_direct_init(cnt);

                        fprintf('*** Loading traj ... ***\n')
                        disp(numpart_esti)
                        date_flag = datestr(run_dates_dt(t+j,:), 'yyyymmddHHMMSS');
                        if ismember(numpart_esti, 'unknown')
                            numpart_esti = readpart10(FLEXPART_folder, date_flag, numpart_esti); %/ If numpart is unknown, input 'unknown' to estimate it.
                        end

                        partoutput = readpart10(FLEXPART_folder, date_flag, numpart_esti);           
                        fprintf('*** [y = %d, t = %d/%d] [%.1f days %s from %s]: %d trajs on %s (%d/%d) are loaded. *** \n',...
                                year, t, max(t_list), maxtraj_day, trajdirect, run_dates_dt(t,:), numpart_esti, run_dates_dt(t+j,:), cnt, length(trajtime_loop_direct_sharpcut))

                        %/ loop over basins
                        for top = 1:NBasin
                            basin_name = [basin_catalog.name];
                            bndry_data = basin_catalog(top).bndry;   %<<-- important!
                            fprintf('*** [cnt = %d, basin %d/%d] %s ... ***\n', cnt, top, NBasin, basin_name{top})

                            if cnt == 1
                                %/ preallocate struct field for EACH basin!
                                domfill_parfor = struct('slct_id_list', cell(1,1),...
                                                         'traj',cell(1,length(trajtime_loop_direct_sharpcut)),...
                                                         'q',   cell(1,length(trajtime_loop_direct_sharpcut)),...
                                                         'BLH', cell(1,length(trajtime_loop_direct_sharpcut)),...
                                                         'T',    cell(1,length(trajtime_loop_direct_sharpcut)),...
                                                         'topo', cell(1,length(trajtime_loop_direct_sharpcut)),...
                                                         'dates_dt', cell(1,length(trajtime_loop_direct_sharpcut)));
                            else
                                domfill_parfor = fwd{top};  %/retrieve the struct of the basin in the cell from cnt == 1
                            end

                            %======================== allocation ======================
                            %/ store q and id one timestep before the start of forward traj.
                            id = partoutput.npoint;
                            domfill_parfor(cnt).traj     = partoutput.xyz;
                            domfill_parfor(cnt).q        = partoutput.vars(:,1) * 1000; %/ from kg/kg to g/kg.
                            domfill_parfor(cnt).BLH      = partoutput.vars(:,2);
                            domfill_parfor(cnt).T        = partoutput.vars(:,3);
                            domfill_parfor(cnt).topo     = partoutput.vars(:,4);
                            domfill_parfor(cnt).dates_dt = run_dates_dt(t+j);

                            if cnt == 2                                        %/ subset trajs from basin at cnt == 2 since it is the *real* t = 0.
                                if ~isempty(bndry_data(:,1) < 0)               %/ convert to [0, 360) if not so.
                                    ind = find(bndry_data(:,1) < 0);
                                    bndry_data(ind,1) = bndry_data(ind,1) + 360;  
                                end 

                                traj_x = partoutput.xyz(:,1);
                                traj_y = partoutput.xyz(:,2);
                                if ~isempty(traj_x < 0)                        %/ convert to [0, 360) if not so.
                                    ind = find(traj_x < 0);
                                    traj_x(ind) = traj_x(ind) + 360;  
                                end

                                %/ Select the traj within the basin boundary (do NOT use bndry_data for 'rest'!)
                                if from_basin == 8 && contains(basin_name{top}, 'rest')
                                    if top ~= NBasin
                                        error('Make sure the rest of area (''rest'') is placed at last');
                                    end
                                    in = out;  %/ The remaining unassigned trajs will be for the rest of area
                                else
                                    try
                                        [in, ~] = inpoly2([traj_x, traj_y], bndry_data);                          %/ inpoly2 is 600xx faster than inpolygon!! 
                                    catch 
    %                                     warning(E.identifier,'Error caught!: %s. Proceed with inpolygon().', E.message)
                                        [in, ~] = inpolygon(traj_x, traj_y, bndry_data(:,1), bndry_data(:,2));    %/ sometimes inpoly2 could run into a bug, then we use inpolygon.
                                    end
                                    if top == 1
                                        out = ~in;
                                    else
                                        out = (out & ~in); 
                                    end
                                end
                                ind = find(in == 1);

                                %/ check if ind is empty.
                                if isempty(ind)
                                    warning('No forward trajs found from the basin')
                                    domfill_parfor(1).slct_id_list = [];
                                    domfill_parfor(1).mass         = [];
                                    for k = 1:2
                                        domfill_parfor(k).traj   = [];
                                        domfill_parfor(k).q      = [];
                                        domfill_parfor(k).BLH    = [];
                                        domfill_parfor(k).T      = [];
                                        domfill_parfor(k).topo   = [];
                                    end
                                else
                                    domfill_parfor(1).slct_id_list = id(ind);
                                    domfill_parfor(1).mass         = partoutput.xmass(ind,1); 
                                    for k = 1:2
                                        domfill_parfor(k).traj   = domfill_parfor(k).traj(ind,:);
                                        domfill_parfor(k).q      = domfill_parfor(k).q   (ind);
                                        domfill_parfor(k).BLH    = domfill_parfor(k).BLH (ind);
                                        domfill_parfor(k).T      = domfill_parfor(k).T   (ind);
                                        domfill_parfor(k).topo   = domfill_parfor(k).topo(ind);
                                    end
                                end

                                domfill_parfor = slct_traj_parfor('trajdirect', trajdirect, 'domfill', domfill_parfor,...
                                                              'RHc', RHc, 'dqc', dqc, 'str_RHc_dqc', str_RHc_dqc, 'str_remark', str_remark,...
                                                              'RHc_map_lon', RHc_map_lon, 'RHc_map_lat', RHc_map_lat, 'slct_one_date_dt', slct_one_date_dt, 'BLH_factor', BLH_factor); 

                            end
                            fwd{top} = domfill_parfor;                         %/update the struct in the cell
                        end
                        clear partoutput; %/ release memory 
                    end
                    toc

                    %/ First replicate the first time step of fwd{top}
                    %/ to the rest time steps to avoid bugs due to empty cells.
                    for top = 1:NBasin
                        for k = cnt+1:length(trajtime_loop_direct_sharpcut)
                            fwd{top}(:,k) = fwd{top}(:,1);
                        end
                    end
                    fwd = cat(1, fwd{:});  %/ remove nests (so that it can be looped by parfor)

                    %======== Parfor loop over the rest traj time !!! ============%
                    if isempty(gcp('nocreate')) && ~isempty(NumWorkers_IO) %/ if set worker number
                        % parpool('Threads', NumWorkers_IO)  %/ It cannot work with the reading
                        parpool('Processes', max([NumWorkers_IO,NumWorkers_HPC]))
                    end
                    tic        
                    parfor (k = cnt+1:length(trajtime_loop_direct_sharpcut), NumWorkers_IO)  %/ parfor for the rest (starting from cnt+1)
                        %/ Broadcast vars
                        % fwd_bc = fwd;
                        t_list_bc               = t_list;
                        slct_date_dt_bc         = run_dates_dt;
                        trajtime_loop_direct_bc = trajtime_loop_direct_sharpcut;
                        j                       = trajtime_loop_direct_bc(k);

                        fprintf('*** Loading traj ... ***\n')
                        date_flag = datestr(slct_date_dt_bc(t+j,:), 'yyyymmddHHMMSS');

                        partoutput = readpart10(FLEXPART_folder, date_flag, numpart_esti);            
                        fprintf('*** [y = %d, t = %d/%d] [%.1f days %s from %s]: %d trajs on %s (%d/%d) are loaded. *** \n',...
                                 year, t, max(t_list_bc), maxtraj_day, trajdirect, slct_date_dt_bc(t,:), numpart_esti, slct_date_dt_bc(t+j,:), k, length(trajtime_loop_direct_sharpcut))

                        %======================== allocation ======================
                        for top = 1:NBasin
                            domfill_parfor_k = fwd(top,k);
                            slct_id_list = domfill_parfor_k(1).slct_id_list;
    
                            if isempty(slct_id_list)
                                continue;                   %/ skip if no trajs from the basin
                            end
    
                            ind = findismember_loop(partoutput.npoint, slct_id_list);  %/ get indices based on id list (no auto sorting)
    
                            domfill_parfor_k.traj     = partoutput.xyz(ind, :);  
                            domfill_parfor_k.q        = partoutput.vars(ind,1) * 1000; %/ from kg/kg to g/kg.
                            domfill_parfor_k.BLH      = partoutput.vars(ind,2);
                            domfill_parfor_k.T        = partoutput.vars(ind,3);
                            domfill_parfor_k.topo     = partoutput.vars(ind,4);
                            domfill_parfor_k.dates_dt = run_dates_dt(t+j);
   
                            fwd(top,k) = domfill_parfor_k;  %/ UPDATE
                        end
                    end
                    toc

                    %/ Assemble the processed struct data into one, for each basin
                    tic
                    fprintf('*** Re-structure the data (fwd -> fwd_final) ... ***\n')
                    fwd_final = cell(NBasin, 1);
                    for top = 1:NBasin
                        domfill_parfor   = fwd(top,:);                             %/ now use this complete domfill_parfor.
                        domfill.traj     = cat(3, domfill_parfor(:).traj);         %/ change nonscalar domfill_parfor to scalar domfill (compatible for other codes) 
                        domfill.q        = cat(2, domfill_parfor(:).q);            %/ and concatenate the numbered fields 
                        domfill.BLH      = cat(2, domfill_parfor(:).BLH);
                        domfill.dates_dt = cat(2, domfill_parfor(:).dates_dt);
                        domfill.T        = cat(2, domfill_parfor(:).T);
                        domfill.topo     = cat(2, domfill_parfor(:).topo);
    %                     domfill.T        = domfill_parfor(1).T;
                        domfill.mass     = domfill_parfor(1).mass;
            %             domfill_parfor = [];   %/ release memory

                        domfill.trajdirect   = trajdirect;
                        domfill.src_lon      = src_lon;
                        domfill.src_lat      = src_lat;
                        domfill.src_z        = src_z;
                        domfill.dt_slct      = dt_slct;
                        domfill.maxtraj_ts   = maxtraj_ts;
                        domfill.maxtraj_day  = maxtraj_day;
                        domfill.trajtime     = trajtime_sharpcut;

                        fwd_final{top} = domfill;
                    end
                    toc

                    %===== Save traj data =====%
                    if save_traj
                        fprintf('*** [y = %d, t = %d/%d] Saving traj into: %s ... *** \n', year, t, max(t_list), traj_fullpath)
                        save(traj_fullpath, 'fwd_final', '-v7.3');             %/ save/overwrite
                    end
                end

                %===== WaterDrip algorithm =====%
                if overwrite_watersip_waterdrip == 2
                    fprintf('!!! [y = %d, t = %d/%d] uptake and watersip data are skipped as per user''s request. !!!', year, t, max(t_list))
                else
                    tic
                    %===== Check if the uptake_map data exists, load from it if overwrite = 0 =======%
                    for top = 1:NBasin
                        basin_name = [basin_catalog.name];
                        domfill    = fwd_final{top};

                        if overwrite_watersip_waterdrip == 0
                            if isfile(Cf_map_fullpath{top})            Cf_map_exist = 1;             else Cf_map_exist = 0;           end
                            if isfile(Cf_map_dates_fullpath)           Cf_map_dates_exst = 1;        else Cf_map_dates_exst = 0;      end
                            if isfile(optimal_trajtime_fullpath{top})  optimal_trajtime_exist = 1;   else optimal_trajtime_exist = 0; end
                        end
                        
                        %/ Either overwrite is required or any WaterDrip files not exists, we will run the algorithm.
                        if overwrite_watersip_waterdrip == 1 || (Cf_map_exist == 0 || Cf_map_dates_exst == 0 || optimal_trajtime_exist == 0)
    
                            Cf_map_dates  = run_dates_dt(t+trajtime_loop_direct_sharpcut(1:end-1)); %/ minus 1 time step to correctly assign the contribution

                            fprintf('*** [basin %d/%d %s, y = %d, t = %d/%d] Start WaterDrip (forward) algorithm... *** \n', top, NBasin, basin_name{top}, year, t, max(t_list))
                            [Cf_map, watersip_fwd_stacktable, mean_optimal_trajtime_map] ...
                                    = WaterDrip('domfill', domfill, 'lon', lon, 'lat', lat, 'dqc', dqc, 'RHc', RHc, 'BLH_factor', BLH_factor, 'area', area,...
                                                   'NumWorkers', NumWorkers_HPC, 'optimal_rr', optimal_rr);

                            % max(abs(Cf_map), [], 'all')
                            % figure
                            % contourf(Cf_map(:,:,2))
                            if savemat
                                if Cf_map_exist == 0             save(Cf_map_fullpath{top},             'Cf_map',                    '-v7.3'); end
                                if Cf_map_dates_exst == 0        save(Cf_map_dates_fullpath,            'Cf_map_dates',              '-v7.3'); end
                                if optimal_trajtime_exist == 0   save(optimal_trajtime_fullpath{top},   'mean_optimal_trajtime_map', '-v7.3'); end
                                fprintf('*** [y = %d, t = %d/%d] WaterDrip products saved under %s *** \n', year, t, max(t_list), prcssd_dir)
                            end
                        end
                    end
                    toc
                end
            end
        end
        
        ed_time = datetime('now');
        %/ Job Progress Notification (No need if only one year is to be processed.)
        if year ~= year_list(end) && ~isempty(send_email_to)
            email_subject = sprintf('[%s] Job Progress: Year %d Completed', server_name, year);
            first_sentence = sprintf('Job Progress: Your Matlab Program ''moisture_tracking.m'' has just finished the job for the year %d!\n', year);
            email_content = [sprintf('Hello from %s!\n', server_name),...
                 newline,...
                 first_sentence,...
                 newline,...
                 sprintf('======= Basic Details ======= \n'),...
                 sprintf('Start Time: %s\n', st_time),...
                 sprintf('End Time: %s\n',   ed_time),...
                 sprintf('job_id:  %d\n',     job_id),...
                 sprintf('Njob:  %d\n',     Njob),...
                 sprintf('year_list:  %d-%d\n',  year_list(1), year_list(end)),...
                 sprintf('RHc_dqc_scheme:  %d\n', RHc_dqc_scheme),...
                 newline,...
                 sprintf('test_mode:  %d\n',     test_mode),...
                 sprintf('savemat:  %d\n',  savemat),...
                 sprintf('save_traj:  %d\n', save_traj),...
                 sprintf('reload_traj:  %d\n', reload_traj),...
                 sprintf('overwrite_watersip_waterdrip:  %d\n', overwrite_watersip_waterdrip),...
                 newline,...
                 sprintf('======= Further Details ======= \n'),...
                 sprintf('trajdirect:  %s\n', trajdirect),...
                 sprintf('from_basin:  %d\n', from_basin),...
                 sprintf('output_rr_map:  %d\n', output_rr_map),...
                 sprintf('output_LA_3D:  %d\n', output_LA_3D),...
                 sprintf('maxtraj_day:  %d\n', maxtraj_day),...
                 sprintf('optimal_rr:  %.3G\n', optimal_rr),...
                 sprintf('str_remark:  %s\n', str_remark),...
                 sprintf('NumWorkers:  %d\n', NumWorkers_IO),...
                 ];
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
        end
    end
    fprintf('\n !!!!!!!!!!!!!!!!!!!!!!')
    fprintf('\n !!! Job completed. !!!\n')
    fprintf(' !!!!!!!!!!!!!!!!!!!!!!\n')

    %/ Shut down parpool
    % poolobj = gcp('nocreate');
    % delete(poolobj);

% catch E  %/ e is an MException struct
%     warning('========================================================')
%     warning(E.identifier, 'There was an error! The message was:\n%s', E.message);
%     warning('========================================================')
%     flag_successful = 0;
% end

%% Send a notification email when job is completed or fails
if flag_successful
    email_subject = sprintf('[%s] Job Completed', server_name);
    first_sentence = sprintf('Your Matlab Program ''flexpart_save_traj_update.m'' has just successfully completed a job!\n');
else
    email_subject = sprintf('[%s] Job Failed!', server_name);
    first_sentence = [sprintf('An error just occurred in your Matlab Program ''flexpart_save_traj_update.m''!\n'),...
                      newline,...
                      sprintf('*** Error message ***\n'),...
                      sprintf('%s\n', E.message),...
                      newline,...
                      sprintf('Please have a check!\n')];
end

email_content = [sprintf('Hello from %s!\n', server_name),...
                 newline,...
                 first_sentence,...
                 newline,...
                 sprintf('======= Basic Details ======= \n'),...
                 sprintf('job_id:  %d\n',     job_id),...
                 sprintf('Njob:  %d\n',     Njob),...
                 sprintf('year_list:  %d-%d\n',  year_list(1), year_list(end)),...
                 sprintf('RHc_dqc_scheme:  %d\n', RHc_dqc_scheme),...
                 newline,...
                 sprintf('test_mode:  %d\n',     test_mode),...
                 sprintf('savemat:  %d\n',  savemat),...
                 sprintf('save_traj:  %d\n', save_traj),...
                 sprintf('reload_traj:  %d\n', reload_traj),...
                 sprintf('overwrite_watersip_waterdrip:  %d\n', overwrite_watersip_waterdrip),...
                 newline,...
                 sprintf('======= Further Details ======= \n'),...
                 sprintf('trajdirect:  %s\n', trajdirect),...
                 sprintf('from_basin:  %d\n', from_basin),...
                 sprintf('output_rr_map:  %d\n', output_rr_map),...
                 sprintf('output_LA_3D:  %d\n', output_LA_3D),...
                 sprintf('maxtraj_day:  %d\n', maxtraj_day),...
                 sprintf('optimal_rr:  %.3G\n', optimal_rr),...
                 sprintf('str_remark:  %s\n', str_remark),...
                 sprintf('NumWorkers:  %d\n', NumWorkers_IO),...
                 ];
% disp(email_subject)
disp(email_content)
command = char(sprintf('echo "%s" | mail -s "%s" %s', email_content, email_subject, send_email_to));
[status, ~] = system(command);

if status == 0
    fprintf('=======================================================\n')
    fprintf('=== A notification sent to %s ===\n', send_email_to)
    fprintf('=======================================================\n')
else
    error('something wrong with sending a notification email.') 
end

%% Exit the program 
exit;  %/ [IMPORTANT]: This prevents 'Warning: Error reading character from command line' that causes the program to hang on forever..

%%
%/ Load AR logical matrix (from yurong)
AR_filename = fullfile(data_folder, 'atmospheric_river_logical_matrix.nc');
ncdisp(AR_filename)
AR      = ncread(AR_filename, 'atmospheric_river_region');
AR_lon  = ncread(AR_filename, 'longitude');
AR_lat  = ncread(AR_filename, 'latitude');
AR_time = ncread(AR_filename, 'time');
date_format = 'yyyymmddHHMM';
AR_date = timesince2date('filename', AR_filename, 'date_format', date_format);

close all
for t = 5:10
    contf_data = AR(:,:,t);
    plot_contfmap('contf_data', contf_data, 'contf_lon', AR_lon, 'contf_lat', AR_lat);
end
% dn = datenum(1900, 1, 1, 0, 0, 0);  %/ in days
% date_yyyymmdd     = int64(str2num(datestr(dn + time_ori/24+timeshift,'yyyymmdd')));

%%

T = watersip_stacktable_basin{1,1};
T.loss(1)
sum(T.dq_BL, 'omitnan') + sum(T.dq_FT, 'omitnan')
sum(T.f, 'omitnan')
sum(T.e, 'omitnan')

%%
% for t = t_list
% 
%     %/ Broadcast var
%     trajtime_loop_direct_sharpcut = trajtime_loop_direct;
%     trajtime_sharpcut             = trajtime;
% 
%     str_sharpcut = '';
%     if sharpcut
%         % if t >= t_sharpcut  
%         %     t_cut = t - t_sharpcut; 
%         %     str_sharpcut = '_sharpcut'; %/ To distinguish the sharp-cut traj file from the complete ones
%         % end
%         if t > t_sharpcut  %/ At the final date, there will be three time steps to compute initial gain and Cf map.
%             t_cut = t - t_sharpcut - 1; 
%         else
%             t_cut = 0;
%         end
%         trajtime_loop_direct_sharpcut(end-t_cut+1:end) = []; %/ sharp cut the tracking length when approaching the end day
%     end
%     fprintf('length(trajtime_loop_direct_cut) = %d \n', length(trajtime_loop_direct_sharpcut))
% end

%% Distribution of recylcing ratio 
% close all
% figname = '/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/Plottings/fixed_vs_optimal';
% plot_hist_cdf('x', 0:0.05:1, 'y', watersip_rr(:,5),  'y2', watersip_rr_optimal(:,5),...
%               'yname', '10-d', 'y2name', 'optimal till 20-d', 'y_lim', [0, 12000], 'savepath', figname);

%%
% a_cell = cell(1000, 1);
% a = ones(1, 360).*ones(181,1);
% 
% for i = 1:length(a_cell)
%     disp(i)
%     a_cell{i} = a;
% end
% b = cell2mat(a_cell);
% b = sum(b,1);
% 
% 
% b = reshape(b, size(a, 1), [])';
% 
% %%
% b = sum(cat(3, a_cell{:}),3);

%% Rename files starting with xxx
% cd('/disk/r034/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_2010/')
% files = dir('EminusP_LA_RHc80_dqc0.2_bwd_20d_optimal_3h_glbland_*.mat'); % Get all files starting with xxx
% % Loop through each
% for id = 1:length(files)
% %     % Get the file name (minus the extension)
%     
%     [~, f] = fileparts(files(id).name);
%     
%     filename_parts = strsplit(f, '_');
%     
%     filename_parts_new = {filename_parts{1:9}, 'intact', filename_parts{10:end}};
%     
%     f_new = strjoin(filename_parts_new, '_');
%     disp(f_new)
% 
%     movefile(files(id).name, sprintf('%s.mat', f_new));
% 
% end





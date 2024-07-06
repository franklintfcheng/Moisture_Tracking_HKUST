%% [Step 1] Read data and FLEXPART-WaterSip experiments parameters
clearvars;

%==========================================================================
matlab_package_path = '/disk/r059/tfchengac/Moisture_Tracking_HKUST/';  %/ Set the path of the folder where you store the matlab packages
% matlab_package_path = '/home/tfchengac/';               %/ Set the path of the folder where you store the matlab packages
masterfolder        = '/disk/r059/tfchengac/FLEXPART/'; %/ Set your own master folder

%/ Set your moisture tracking details (based on moisture_tracking.m)
expmnt = 'domfill_EA_1.0deg_MPI_201806_ALA';         ldirect = 'bwd'; from_basin = 11; RHc_dqc_scheme = 1;  output_res = 0.25; dt = 3; year_list = 2018; optimal_rr = 0.99;    
% expmnt = 'domfill_EA_0.5deg_MPI_202309-10_ScotCase'; ldirect = 'bwd'; from_basin = 10; RHc_dqc_scheme = 0;  output_res = 0.25; dt = 3; year_list = 2023; optimal_rr = 0.99;    
% expmnt = 'domfill_EA_0.5deg_MPI_202202_AusCase';     ldirect = 'bwd'; from_basin = 9;  RHc_dqc_scheme = 0;  output_res = 0.25; dt = 3; year_list = 2022; optimal_rr = 0.99; 
% expmnt = 'domfill_EA_0.5deg_MPI_202207-08_PakCase';  ldirect = 'bwd'; from_basin = 7;  RHc_dqc_scheme = 2;  output_res = 0.25; dt = 3; year_list = 2022; optimal_rr = 0.99; 
% expmnt = 'domfill_EA_0.5deg_MPI_202207-08_PakCase';  ldirect = 'fwd'; from_basin = 8;  RHc_dqc_scheme = 1;  output_res = 0.25; dt = 3; year_list = 2022; optimal_rr = 0.99; 
% expmnt = 'domfill_CERA_MPI';                         ldirect = 'bwd'; from_basin = 6;  RHc_dqc_scheme = 23; output_res = 1;    dt = 3; year_list = 1971:2010; optimal_rr = 0.99; %/ *LATEST*

%==========================================================================
%/ Explanation
%/ 'from_basin':   0]: from global land
%/               1-3]: from global source hotspots
%/                 4]: from TP basins 
%/                 5]: from TP low-veg covers
%/                 6]: from TP grids (0.5x0.5)
%/                 7]: from Pakistan box region
%/                 8]: from IPCC regions (for Pakistan case)
%/                 9]: from Australian box region 
%/                10]: from Scotland box region 
%/                11]: from Alaska (ALA) box region
%/     'expmnt': The name of your FLEXPART experiment (i.e., the folder name)
%/ 'output_res': The output WSV spatial resolution (in deg)
%/         'dt': FLEXPART time interval (in hr)
%/ See also 'moisture.tracking.m'

addpath(genpath(fullfile(matlab_package_path,'MyMatlabPkg')));
addpath(genpath(fullfile(matlab_package_path,'MyMatlabFn')));
addpath(genpath(fullfile(matlab_package_path,'m_map1.4o')));
addpath(genpath(fullfile(matlab_package_path,'Flexpart_Codes')));
addpath(genpath(fullfile(matlab_package_path,'MCS_Codes')));

param = 'param_universal';
project_name = 'flexpart';
run(param);

[str_RHc_dqc, str_remark, str_src_z, str_years, WSV_lon, WSV_lat,...
 maxtraj_day, str_BLH_factor, str_optimal, str_traj_rm_jump, dt_slct, data_folder,...
 plotting_folder, WSV_dir, str_expmntinfo, basin_catalog, str_domain, str_domain_trajfile, str_sharpcut, forcing, stdate, eddate] = ...
    load_tracking_param('expmnt', expmnt, 'ldirect', ldirect, 'output_res', output_res, 'dt', dt, 'year_list', year_list,...
                        'optimal_rr', optimal_rr, 'RHc_dqc_scheme', RHc_dqc_scheme, 'from_basin', from_basin, 'masterfolder', masterfolder);

valid_WSV_dataname = {'Pm', 'Pm_BL', 'Pm_AR', 'P_LA', 'Cf_map', 'uptake', 'BL_uptake',...
                      'Pm_frac_real', 'Pm_frac_adj', 'Pm_frac',...
                      'optimal_trajtime', 'CWRT', 'RH2', 'rr_tot_L', 'rr_tot_NLL', 'rr_tot_NLO'};

labels_IPCC_PAK = {'IPCC-PAK-NEAF', 'IPCC-PAK-SEAF', 'IPCC-PAK-WCA', 'IPCC-PAK-TIB',...
                   'IPCC-PAK-ARP', 'IPCC-PAK-SAS', 'IPCC-PAK-ARS', 'IPCC-PAK-BOB',...
                   'IPCC-PAK-EIO', 'IPCC-PAK-SIO', 'IPCC-PAK-rest'}; 
license('inuse');

%% IPCC AR6 regions

mask_filename = '/disk/r059/tfchengac/FLEXPART/domfill_EA_0.5deg_MPI_202207-08_PakCase/Masks/IPCCregions_Pakistancase.nc';
ncdisp(mask_filename)
mask       = ncread(mask_filename, 'mask');
mask_lon   = ncread(mask_filename, 'lon');
mask_lat       = ncread(mask_filename, 'lat');
mask_region_id = ncread(mask_filename, 'region');
mask_region_abbrevs = ncread(mask_filename, 'abbrevs');
mask_region_fullnames = ncread(mask_filename, 'names');

%% [Step 2.2] Plot topography / global basins / labels_domain_wise_exact (Figs. 1-3)
sig_src_ordered = []; sig_src_colorings = []; text_data = []; map_center = [];  SR = []; 
outline_by_grid = 0;   draw_rings = 0;  check_grouping = 0; trans_bg = 0; 
show_topo = 0; draw_river = 0; show_target_basins = 0; show_bndry_in_patch = 0; show_name = 0; mth_list = 0;
marker = 'none';  markersize = 5; markerfacecolor = [255 51 204]./255; markeredgecolor = 'none';
pcolor_mode = 0; shadedrelief_mode = 0; ColorByDom = 0; slct_reg = [];  slct_point = []; point_data = []; str_slct_reg = []; str_slct_point = [];

close all;
savefig         = 1;
savemat         = 1;
recompute       = 0;     %<- Mind this! It controls whether to recompute bndry/centr of terrestrial/oceanic basins.
fig_fmt         = 'pdf'; %   pdf png
draw_cbar_only  = 0;     %<- mind this!
cbar_mode       = 0;
cbar_location   = 'southoutside';  cbar_position = [100 100 1100 800]; cbar_fontsize = 10; %/ x y w h
% cbar_location   = 'eastoutside';  cbar_position = [100 100 400 400]; cbar_fontsize = 10; %/ x y w h

% show_topo = 0; draw_river = 0; show_bndry_in_patch = 0; patch_alpha = 0.8; show_name = 0; slct_region = 0;  fig_fmt = 'png'; %/ 
show_topo = 0; draw_river = 0; show_bndry_in_patch = 0; patch_alpha = 0.8; show_name = 1; slct_region = 2;  fig_fmt = 'png'; %/ 
% show_topo = 0; show_bndry_in_patch = 1; patch_alpha = 1; show_name = 0; slct_region = 0;  %/ for viz & name checking
% show_topo = 0; show_bndry_in_patch = 1; patch_alpha = 1; show_name = 0; slct_region = 1;  %/ for Pakistan case

%/================= Load the sig sources for TP baisns =================%
% SR_sig_src_filename = strcat(data_folder, 'SR_sig_src_Pm_sum_ins_basin-wise_agglo3.0perc_1971-2010_anom1971_2010_AllTPbasins.mat');
% load(SR_sig_src_filename, 'sig_src_ordered', 'sig_src_colorings');
% sig_src_ordered(ismember(sig_src_ordered, {'Elsewhere'})) = [];  %/ remove 'Elsewhere' before inputting into 'reg_extractor'
% slct_reg = sig_src_ordered;

%/================= Otherwise, Load other regions =================%
% slct_reg = {'Pakistan_box'}; recompute = 0;
slct_reg = {'IPCC-PAK'}; recompute = 0;
% slct_reg = {'IPCC-PAK-rest'}; recompute = 0;
% slct_reg = {'IPCC-PAK-NEAF', 'IPCC-PAK-SEAF', 'IPCC-PAK-WCA', 'IPCC-PAK-TIB',...
%             'IPCC-PAK-ARP', 'IPCC-PAK-SAS', 'IPCC-PAK-ARS', 'IPCC-PAK-BOB',...
%             'IPCC-PAK-EIO', 'IPCC-PAK-SIO', 'IPCC-PAK-rest', }; recompute = 0;

%/=========================================================================

%/ Retrieve bndry_data

for mth = mth_list
    if mth == 0 
        str_timescale = '_annual';   
        str_mth_bc    = '';
    else
        str_timescale = strcat('_', str_mth{mth});
        str_mth_bc    = sprintf('_%s', str_mth{mth});
    end
    
    %/ 'slct_reg' as contf_data/bndry_data
    if ~isempty(slct_reg)
        str_slct_reg = strcat('_', slct_reg);
        if show_bndry_in_patch
            output_HR_ocean  = 1;
        else
            output_HR_ocean  = 0;
        end

        lon = WSV_lon;
        lat = WSV_lat(2:end-1);

        [reg_2D, reg_bndry_list, reg_name_list, ~, reg_centr_list] = reg_extractor('lon', lon, 'lat', lat, 'mth', mth, 'slct_reg', slct_reg,...
                                                                            'output_HR_ocean', output_HR_ocean, 'outline_by_grid', outline_by_grid, 'draw_rings', draw_rings,...
                                                                            'regime_ver', [], 'data_folder', data_folder, 'savemat', savemat, 'recompute', recompute);
        % a = reg_2D;
        % a(isnan(a)) = 0;
        % unique(reshape(a, 1, []))
        if show_bndry_in_patch
            bndry_patch_mode = ones(length(reg_bndry_list), 1);
        else
            bndry_patch_mode = zeros(length(reg_bndry_list), 1);
        end
        
        %/ Colormap / patch_color / Color
        bndry_data   = reg_bndry_list;
        color        = repmat([0 0 0]./255, length(reg_name_list), 1); %/ by default, all outlined in black
        if isequal(slct_reg, {'TP_basins'})  %/ reverse the position to avoid bug when slct_reg is not a single element
            col_exdor = [255 188 200]./255; 
            col_endor = [0 204 255]./255;
        %     col_exdor = [0 204 255]./255;
        %     col_endor = [255 188 200]./255;
            colmap = [repmat(col_exdor, 12, 1);...
                      repmat(col_endor,  3, 1);];
            color = repmat([1 1 1], length(bndry_data), 1);
    %         color(16:end,:) = repmat([1 1 1], 12, 1);

            bndry_data(end+1)       = TP_bndry;
            color(end+1,:)          = [0 0 0];
            colmap(end+1,:)         = [0 0 0];
            bndry_patch_mode(end+1) = 0;  
            
        elseif isequal(slct_reg, {'hydrosheds'})
        %     patch_color        = repmat(col_land,  NoOfColors, 1);
            colmap        = brewermap(NoOfColors, 'spectral');
            rng(129)
            a = 1:size(colmap,1);
            a_rand = a(randperm(length(a)));
            colmap = colmap(a_rand,:);

        elseif isequal(slct_reg, {'IPCC-PAK'})
            colmap = nclCM('amwg_blueyellowred', length(reg_name_list)); colmap(6,:) = [0 204 153]./255;
        else
            rng('default')
            a = 1:length(reg_name_list);
            a_rand = a(randperm(length(a)));
            colmap = nclCM('MPL_Blues', length(a));
            colmap = colmap(a_rand,:);
        end
        patch_color = colmap;
    end

    %/ 'slct_point' as point_data
    if ~isempty(slct_point)
        str_slct_point = strcat('_', slct_point);
        output_HR_ocean = 1;

        if isequal(slct_point, {'TP_grids_0.5x0.5'})
            lon = [0:0.5:360-0.5];
            lat = [-90:0.5:90];
        else
            lon = WSV_lon;
            lat = WSV_lat;
        end

        [~, ~, point_name_list, ~, point_centr_list] = reg_extractor('lon', lon, 'lat', lat, 'mth', mth, 'slct_reg', slct_point,...
                                                                       'output_HR_ocean', output_HR_ocean, 'outline_by_grid', outline_by_grid, 'draw_rings', draw_rings,...
                                                                       'data_folder', data_folder, 'savemat', savemat, 'recompute', recompute);
        point_data = point_centr_list;
                                                                   
        %/ Colormap / patch_color / Color
        if ismember(slct_point, {'TP_and_ML_basins'}) || contains(slct_point, {'TP_grids'}) 
            marker          = 'o';  
            markersize      = 7;
            
            %/ Color them based on their dominant moisture conveyor (ignoring 'unattributed') (Fig. 3)
            if ColorByDom   
                %/ 1. Create a set of colmap (to be adapted)
                if ismember(slct_point, {'TP_and_ML_basins'})
                    markerfacecolor  = [repmat([1 1 1], 15, 1);
                                        repmat([200 200 200]./255, 12, 1);];
                    markeredgecolor  = repmat([1 1 1], length(point_name_list), 1);
                else
                    markerfacecolor  = repmat([1 1 1], length(point_name_list), 1);
                    markeredgecolor  = repmat([1 1 1], length(point_name_list), 1);
                end
                
                %/ 2. Loop over the regions to get their SR matrices [If not computed, go to Step 5.3] 
                SR = [];
                for k = 1:length(point_name_list)
                    basin_name          = point_name_list{k};
                    slct_data           = 'Pm_frac_real';
                    str_grouping_mode   = strcat('_', strrep(grouping_mode, '-', '_'));
                    savemat_prefix      = strcat('scheme', num2str(RHc_dqc_scheme), str_optimal);
                    [SR, SR_filename_suffix] = get_SR('SR', SR, 'project_name', project_name, 'param', param, 'dataset', dataset, 'dataname', dataname, 'data_folder', data_folder, 'valid_WSV_dataname', valid_WSV_dataname,...
                                                'optimal_rr', optimal_rr, 'expmnt', expmnt, 'output_res', output_res, 'dt', dt, 'RHc_dqc_scheme', RHc_dqc_scheme, 'from_basin', from_basin, 'masterfolder', masterfolder,...
                                                'slct_data', slct_data, 'area_mean_or_sum', 'sum', 'ins_or_acc', 'acc', 'grouping_mode', grouping_mode, 'grouping_mode_labels', labels_domain_wise_exact,...
                                                'regime_ver', regime_ver, 'basin_name', basin_name, 'basin_catalog', basin_catalog, 'thres_agglomerate', [],...
                                                'select_field', 'daily', 'DivByAreaOf', [], 'mth', mth, 'str_mth', str_mth, 'year_list', year_list, 'shared_year_list', year_list, 'anom_base_period', [],...
                                                'derive_SR', 0, 'recompute_SR', 0, 'savemat', 0, 'savemat_prefix', savemat_prefix);
                    if ~isempty(SR)
                        basin_name_strrep = strrep(strrep(strrep(basin_name, '-', '_'), ' ', '_'), '.', '_');
                        data4plot         = SR.(basin_name_strrep).([slct_data, str_grouping_mode, str_timescale, '_clim']);
                        data4plot         = data4plot(1:end-1,:);  
                        dataname4plot     = SR.(basin_name_strrep).([slct_data, str_grouping_mode, '_src_name']);
                        [~, I] = max(data4plot);
                        markerfacecolor(k,:) = colmap_domain_wise_exact(I,:);
                    end
                end
            else
                error('code not ready yet!');
            end
        else
            marker          = 'o';  
            markersize      = 5;
            markerfacecolor = [255 51 204]./255; 
            markeredgecolor = 'none';
        end
    end
    
    %/ Append target basins (e.g., TP basins) to the list of boundary data
    if show_target_basins
        reg_name_list = [reg_name_list', basin_catalog.name];
        bndry_data = [bndry_data; {basin_catalog.bndry}'];
        bndry_patch_mode = [bndry_patch_mode; zeros(length(basin_catalog), 1)];
        color = [color; repmat([0 0 0]./255, length(basin_catalog), 1)];
    end
                                               
    %/ Map Domain
    %/      0: global; 1: MC; 2: Caspian Sea; 3: Eastern Hemisphere;  4: TP,  5: Eurasia, 
    %/      6: TP (zoom-in); 7: for SR network 
    backcolor       = [1 1 1]; %[.7 .7 .7];
    coast_patch_col = [];
    col_land        = [131 255 190]./255;
    col_ocean       = [51 153 204]./255;  % [0 43 66]./255;
    if slct_region == 0        %/ global 
        glb_data_mode = 1;
        fontsize      = 18;
        markersize    = 2; 
        grid_mode     = 0;
        coast_col     = 'k';
        linewi        = 1.5;
        coast_wi      = 2;
        map_proj      = 'equidistant cylindrical';
        % map_proj      = 'robin';
        map_lon_lower = -179;
        map_lon_upper = 180;
        map_lat_lower = -89;
        map_lat_upper = 89;
        str_regional = '_global';
        text_fontsize = fontsize*0.5;

    elseif slct_region == 1   %/ Pakistan
        glb_data_mode = 0;
        fontsize      = 18;
        linewi        = 1;
        coast_wi      = 1.5;
        grid_mode     = 6;     %/ no need to draw lon lat
        map_proj      = 'Miller Cylindrical';
    %     map_proj      = 'lambert';
        map_lon_lower = 46;
        map_lon_upper = 90;
        map_lat_lower = 0;
        map_lat_upper = 35; 
        str_regional = '_Pakistan';
        trans_bg = 1; %/ output fig with transparent bg = 1
        coast_col       = 'none';
        coast_patch_col = [240 239 221]./255;
        backcolor       = [151 182 226]./255;
        text_fontsize = fontsize*0.5;

    elseif slct_region == 2   %/ Pakistan (larger)
        glb_data_mode = 0;
        fontsize      = 14;
        linewi        = 1;
        coast_wi      = 1.5;
        grid_mode     = 3;     %/ no need to draw lon lat
        map_proj      = 'Miller Cylindrical';
    %     map_proj      = 'lambert';
        map_lat_lower = -33;
        map_lat_upper = 48; 
        map_lon_lower = 15;
        map_lon_upper = 130;
        str_regional = '_Pakistan';
        trans_bg = 1; %/ output fig with transparent bg = 1
        coast_col       = [255 255 51]./255;
        coast_patch_col = 'none';
        backcolor       = 'none';
        text_fontsize = 10;
    else
        error('slct_region not set!');
    end
    title_pos = [0.12,1.01,1.0,0];
    
    glb_plateau_mode        = 0;
    plateau_hgt             = 3000;
    plateau_col             = [255 51 204]./255;
 
    if show_topo
        m_proj('Miller Cylindrical','longitudes',[0 360], 'latitudes', [-90 90]);
        if isequal(map_proj, 'ortho')
            topo_domain = [0, 360, -90, 90];
        else
            topo_domain = [map_lon_lower, map_lon_upper, map_lat_lower, map_lat_upper];
        end
        %/ Extract 1-minute topography from etopo1
        [contf_data, lon_2D_topo, lat_2D_topo]=m_etopo2(topo_domain);

        %/ Extract 5-minute topography from tbase
%             [contf_data, lon_2D_topo, lat_2D_topo]=m_tbase(topo_domain); %REGION =[west east south north];
%             
%             topo_intvl = 225;
%             contf_data(contf_data < 0) = -100;
%             contf_levels = [-topo_intvl:topo_intvl:topo_intvl*26];
%             colmap = my_colormap(length(contf_levels)-1,'topo_land_27lev');
%             colmap = copper(length(contf_levels)-1);
%             colmap = nclCM('OceanLakeLandSnow', length(contf_levels)-1); colmap(1:2,:) = [30 160 244;28 89 30]./255;

        topo_intvl = 100; topo_max = 6000 + topo_intvl;
        contf_data(contf_data < 0) = -100;
        contf_levels = [-topo_intvl:topo_intvl:topo_max];
%             colmap = [[30  160 244]./255;m_colmap('gland', topo_max/topo_intvl)];
%             colmap = copper(length(contf_levels)-1);
%             colmap = brewermap(length(contf_levels)-1, 'Greys');
        colmap = nclCM('OceanLakeLandSnow', length(contf_levels)-1); colmap(1:2,:) = [30 160 244;28 89 30]./255;


%             contf_levels = [-topo_max:topo_intvl:topo_max];
%             colmap = nclCM('GMT_relief', length(contf_levels)-1);
%             colmap(1:length(contf_levels)/2, :) = [];     %/ discard half of the ocean colors.
%             colmap(1,:) = [30  160 244]./255;             %/
%             contf_levels(1:length(contf_levels)/2) = [];  %/ discard half of the ocean colors.


        contf_lon       = lon_2D_topo(1,:);
        contf_lat       = lat_2D_topo(:,1);
        contf_unit      = '';

        pcolor_mode     = 0;  shadedrelief_mode = 1; 
%             pcolor_mode     = 1;  shadedrelief_mode = 0; %/ pcolor generates a larger fig
%             pcolor_mode     = 0;  shadedrelief_mode = 0;
        cbar_interval   = 2;
        cbar_YTick      = contf_levels(2:10:end-1);
    %     cbar_YTick      = [-4000:1000:7000];
        cbar_YTickLabel = cbar_YTick;
        coast_col       = 'none';
        str_show_topo   = 'topo_';
    else
        if show_bndry_in_patch
            contf_data    = [];
        else
            contf_data    = reg_2D;
        end
        contf_lon     = lon;
        contf_lat     = lat;
        contf_levels  = 1:length(reg_name_list)+1;   %/ since id in reg_2D starts from 1. 
        pcolor_mode   = 1;
        contf_unit    = '';
        NoOfColors    = length(contf_levels)-1;
        cbar_interval   = 1;
        cbar_YTick      = contf_levels(1:end-1) + 0.5;
        cbar_YTickLabel = reg_name_list; 
        str_show_topo   = [];
    end

    %/ Create text_data to label the basins
    if show_name
        text_data = cell(length(reg_name_list),4);
        for i = 1:length(reg_name_list)
            if isequal(slct_reg, {'TP_basins_with_ds'}) && contains(reg_name_list{i}, 'TP-')
                % do not label it.
            else
                %/ Skip labeling if it is outside the map domain.
                if reg_centr_list(i,1) < map_lon_lower || reg_centr_list(i,1) > map_lon_upper || ...
                   reg_centr_list(i,2) < map_lat_lower || reg_centr_list(i,2) > map_lat_upper
                    continue;
                end

                text_data{i,1} = reg_centr_list(i,1);
                text_data{i,2} = reg_centr_list(i,2);
                if ~isempty(sig_src_ordered)
                    text_data{i,3} = num2str(i);
                else
                    text_data(i,3) = strrep(reg_name_list(i), '_', ' ');
                end
            end
        end
        text_data(:,4) = {[0 0 0]};   %/ set text color to be black
        %/ Set text to be white if color is too dark
        ind = find(mean(patch_color,2)<0.35);
        for i = 1:length(ind)
            text_data(ind(i),4) = {[1 1 1]};
        end
    end

    create_fig = 1;
    
    if ~isempty(str_mth_bc)
        titlename = strrep(str_mth_bc, '_', ' ');
    else
        titlename = [];
    end
    if savefig
        if ~isempty(sig_src_colorings)
            savepath = strcat(plotting_folder, 'basin_bndry_SR_network_AllTPbasins', str_slct_point, str_regional, str_mth_bc);
        else
            savepath = strcat(plotting_folder, str_show_topo, 'basin_bndry', str_slct_reg, str_slct_point, str_regional, str_mth_bc);
        end
    else
        savepath = [];
    end
    
    % close all;
    plot_contfmap('contf_data', contf_data, 'contf_lon', contf_lon, 'contf_lat', contf_lat, 'contf_levels', contf_levels,...
                  'contf_unit', contf_unit, 'colmap', colmap, 'cbar_interval', cbar_interval, 'pcolor_mode', pcolor_mode, 'shadedrelief_mode', shadedrelief_mode,...
                  'point_data', point_data, 'marker', marker, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'markeredgecolor', markeredgecolor, 'linewi', linewi, 'color', color, ...
                  'bndry_data', bndry_data, 'bndry_patch_mode', bndry_patch_mode, 'patch_color', patch_color, 'patch_alpha', patch_alpha, 'text_data', text_data, 'text_fontsize', text_fontsize,...
                  'titlename', titlename, 'title_pos', title_pos, 'savepath', savepath, 'fig_fmt', fig_fmt, 'trans_bg', trans_bg,...
                  'map_proj', map_proj, 'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper, 'map_center', map_center,...
                  'coast_col', coast_col, 'coast_wi', coast_wi, 'coast_patch_col', coast_patch_col, 'backcolor', backcolor,...
                  'glb_data_mode', glb_data_mode, 'glb_plateau_mode', glb_plateau_mode, 'plateau_hgt', plateau_hgt, 'plateau_col', plateau_col,...
                  'draw_river', draw_river, 'cbar_location', cbar_location, 'draw_cbar_only', draw_cbar_only, 'cbar_fontsize', cbar_fontsize, 'cbar_position', cbar_position, 'fontsize', fontsize,...
                  'create_fig', create_fig, 'grid_mode', grid_mode, 'cbar_mode', cbar_mode, 'cbar_YTick', cbar_YTick, 'cbar_YTickLabel', cbar_YTickLabel);
    fprintf('done \n')
    
    %/ Checking
%     a = bndry_data{:};
%     bndry_lon = a(:,1);
%     bndry_lat = a(:,2);
%     [ind_nan, ~] = find(isnan(bndry_lon));
%     m_line(bndry_lon(ind_nan-1), bndry_lat(ind_nan-1), 'Marker', 'o', 'Markersize', 5, 'Color', 'w', 'markerfacecolor', 'g', 'linest', 'none');

end

%% [Step 3] read_any (reanalysis / gauge / satellite / multi-product ensemble mean)
close all; time_dim = 3; stlevel = []; edlevel = []; lon_range = []; lat_range = []; plot_to_check = 0; dataset = []; 
savefig = 0; 

savemat            = 1;
recompute          = 0;   %/ <- mind this!
set_conden_to_zero = 1;   %/ Default to set condensation to zeros for evap and prcp products.
NumWorkers         = 10;

if isequal(project_name, 'CMIP6_MJO')
    noleap = 1;    %/ Omitting leap days ensures a consistent length of days with CMIP6 models (most used a 365-day calendar)
else
    noleap = 0;  
end

select_field = {'daily'}; %/ 'daily', 'subdaily', 'dailyAnom_PMC', 'monthly'
slct_year    = year_list;

% slct_data = {'ERA5_E'};
slct_data = {'ERA5_P'};
% slct_data = {'ERA5_P_05'};

dataset = read_any('dataset', dataset,  'param', param, 'slct_data',  slct_data, 'select_field', select_field, 'slct_year', slct_year,...
                   'time_dim', time_dim, 'noleap', noleap, 'stlevel', stlevel, 'edlevel', edlevel, 'lon_range', lon_range, 'lat_range', lat_range,...
                   'data_folder', data_folder, 'savemat', savemat, 'recompute', recompute, 'set_conden_to_zero', set_conden_to_zero,...
                   'NumWorkers', NumWorkers, 'plot_to_check', plot_to_check, 'savefig', savefig);

%% [Step 4] Geographic Map (Pm, Pm_AR, Cf_map, etc.)
slct_contf_str = []; slct_cont_str = []; slct_vector_str = []; slct_vector2_str = []; slct_vector_str2 = ''; slct_vector2_str2 = ''; slct_hatch_str = []; get_TP_LoVeg_bndry = 0;
contf_data = []; pcolor_mode = 0; SSIM = []; str_trend = []; 
Udata = []; Vdata = []; Udata_dates = []; Vdata_dates = []; uv_lon = []; uv_lat = []; vec_step_lon = []; vec_step_lat = []; vector_color = []; vector_edgecolor = []; vecscale = []; vecscale2 = []; shaftwidth = []; headlength = []; vec_lbs = []; vec_mag_ref = []; vec_ref_fontsize = [];
U2data = []; V2data = []; U2data_dates = []; V2data_dates = []; uv2_lon = []; uv2_lat = []; vec2_step_lon = []; vec2_step_lat = []; vector2_color = []; vector2_edgecolor = []; vec2scale = []; vec2scale2 = []; shaft2width = []; head2length = []; vec2_lbs = []; vec2_mag_ref = []; vec2_ref_fontsize = [];
point_data = []; marker = []; markersize = []; markerfacecolor = []; linewi = []; color = []; bndry_data = [];
text_data = []; text_backgroundcolor = []; text_edgecolor = []; text_fontsize = []; 
hatch_data = []; hatch_data_dates = []; hatch_lon = []; hatch_lat = []; hatch_thres_pve = []; hatch_thres_nve = []; color_hatch_pve = []; color_hatch_nve =[]; hatch_mode = 0; hatch_linewi = []; hatch_intvl = []; 
cont_data = []; cont_data_raw = []; cont_data_dates = []; cont_lon = []; cont_lat = []; cont_levels = []; cont_colmap = []; cont_linewi = []; cont_labelsize = []; cont_label_col = []; skip_zero_cont = [];
contf_with_stipples = 0; contf_with_hatch = 0; draw_country = 0; compute_anombyAM = 0; alpha = []; drywet_mode = 0; str_drywet = []; landocean_text_mode = 0; bndry_patch_mode = 0;
slct_list = [];  str_land_or_ocean = []; RHc_map_scheme = []; landocean_text_mode = 0; show_TP_only = 0; glb_plateau_mode = 0; all_in_one = 0;
map_center = []; coast_patch_col = []; patch_color = []; patch_alpha = []; topx = []; prct = []; criteria = []; str_TP_only = ''; fig_fmt = 'pdf'; save_diffmap = 0; plot_sc_TP_RHc_error = 0; slct_date = []; str_hatch = '';
title_fontsize = []; season = []; AWM_lon_extent = []; AWM_lat_extent = [];

dates = date_array_gen('year_list', year_list, 'st_month', 1, 'st_day', 1, 'ed_month', 12, 'ed_day', 31);
if year_list(end) == 2010       dates(end) = [];     end  %/ remove 2010-12-31 since no CERA data

%================================================  
savefig             = 1;
savemat             = 1;
savenc              = 0;  %/ [IMPORTANT] Save WSV products as nc files (send to collaborators)
recompute           = 0;  %/ recompute the WSV product
cbar_mode           = 1;
trans_bg            = 1;
draw_cbar_only      = 0;  
% cbar_location = 'southoutside'; cbar_fontsize = 11; cbar_position = [100 100 500 250];  %/ x y w h
cbar_location = 'eastoutside'; cbar_fontsize = 14; cbar_position = [100 100 400 800];  %/ x y w h
draw_refvec_only    = 0;  %/ draw reference vec only
curvature           = 0;  %/ curved vectors [Not recommended]
sig_contf_in_hatch  = 0;  %/ <-- mind this! 1]: denote sig contf values in hatching (not robust), otherwise in stippling (more accurate)
plot_land_or_ocean  = 0;  %/ <-- mind this! Only for reanalysis data. Whether to filter land [1] or ocean [2] grids, or all [0].
plot_map            = 1;     
plot_pie_mode       = 0;   savefig_pie         = 0;  %/ draw pie chart instead of showing text on map. (Fig. 3)
plot_zonal_mean     = 0;   savefig_zonal_mean  = 0;   
sig_mode            = 2;                             %/ 1]: with hatching  2]: no hatching, only plot sig as contf
trend_test          = 'MK';  trend_alpha = 0.05;  save_trend = 1;  recompute_trend = 0;   %/ compute trends? 
select_field        = 'daily';
set_PE_mm_per_yr    = 0;
plateau_hgt         = 3000;
plateau_col         = [.1 .1 .1]; %/ [255 51 204]./255;
show_TP             = 1;
NumWorkers          = [];     
show_sink_region    = 1;
compute_contf_AWM   = 1;  
fig_fmt = 'png'; png_dpi = 200;   

%/ NOTE: 'slct_region':   0]: global (center at 0E)    1]: global (center at dateline)
%/                        2]: North Pacific rim        3]: true-global                  4]: TP 
%/                        5]: Eurasia                  6]: 3D sphere
%/          'mth_list':   0]: Annual                1-12]: Monthly                  13-16]: Seasonal
%/                       17]: AMJJAS                  18]: ONDJFM                      19]: JFD
%/                       20]: nonJJA                  21]: MJJASO                      22]: NDJFMA
%/                       23]: MJJAS

%/=========================== Region & Time ===============================
mth_list = 0; %/ dummy
basin_list = 1:length(basin_catalog); 
if contains(expmnt, 'PakCase') 
    st_month = 8; st_day = 10; ed_month = 8; ed_day = 24;   %/ Pakistan case (whole period mean)
    slct_region = 5; 

elseif contains(expmnt, 'AusCase')
    % st_month = 2; st_day = 22; ed_month = 2; ed_day = 28;   %/ Australian case (whole period mean)
    st_month = 2; st_day = 25; ed_month = 2; ed_day = 25;   %/ Australian case (the most extreme flood day)
    slct_region = 6;  
    
elseif contains(expmnt, 'ScotCase')
    if isequal(expmnt, 'domfill_EA_0.5deg_MPI_202208_ScotCase_test')
        st_month = 8; st_day = 10; ed_month = 8;  ed_day = 12;   %/ testing
    else
        st_month = 10; st_day = 6; ed_month = 10; ed_day = 8;    
    end
    slct_region = 7;
else
    %/ [By default] Use the tracking period 
    if numel(num2str(stdate)) == 8
        year_list = floor(stdate./1e4):floor(eddate./1e4); 
        
        st_month = floor(mod(stdate, 1e4)./1e2); st_day = mod(stdate, 1e2); 
        ed_month = floor(mod(eddate, 1e4)./1e2); ed_day = mod(eddate, 1e2); 
    end
    if contains(expmnt, 'ALA')
        slct_region = 2;
    else
        slct_region = 0;
    end
end
year_list_bc = year_list;

%==========================================================================
for top = basin_list
    close all; 
    basin_name           = basin_catalog(top).name{:};   
    basin_bndry          = basin_catalog(top).bndry;   %/ input to 'get_drywet_days' fn. 
    
    %/=========================== Variables ===================================
                     %contf          %cont      %vector  %vector2   %hatch   %trend_mode
    % slct_list = { 'Pm',               '',        '',     '',     '',   0;}; mean_or_sum = 'sum'; 
    % slct_list = { 'Pm_AR',            '',        '',     '',     '',   0;}; mean_or_sum = 'sum';  

    slct_list = { 'P_LA',             '',        '',     '',     '',   0;   
                  'ERA5_P',           '',        '',     '',     '',   0;}; mean_or_sum = 'sum'; AWM_lon_extent = [min(basin_bndry(:,1)), max(basin_bndry(:,1))]; AWM_lat_extent = [min(basin_bndry(:,2)), max(basin_bndry(:,2))]; 
        
    % slct_list = { 'P_LA',             '',        '',     '',     '',   0;    
                  % 'ERA5_P_05',        '',        '',     '',     '',   0;}; mean_or_sum = 'sum';slct_region = slct_region + 0.5; %/ Zoom-in
    
    % slct_list = { 'ERA5_P',           '',        '',     '',     '',   0;}; slct_region = slct_region + 0.5; %/ Zoom-in
    % slct_list = { 'Cf_map',           '',        '',     '',     '',   0;}; 
    % slct_list = {'ERA5_IVTdiv',       '',        '',     '',     '',   1;}; 
    % slct_list = {'CWRT',              '',        '',     '',     '',   0;};
    % slct_list = {'CWRT',              '',        '',     '',     '',   1;};  %/ Contribution-weighted Lagrangian residence time 
    
    %/=========================== Other Settings ==============================
    str_years_bc = sprintf('_%d-%d', year_list_bc(1), year_list_bc(end));
    if contf_with_stipples                            str_marker   = sprintf('_p%dstipple', prct); 
    elseif contf_with_hatch                           str_marker   = sprintf('_p%dhatch', prct);    
    else                                              str_marker   = [];                                        end
    if draw_cbar_only                                 str_cbar     = '_cbar';         else    str_cbar = [];          end
    if compute_anombyAM                               str_anombyAM    = '_anombyAM';  else    str_anombyAM = [];      end
    if ~isempty(alpha)    
        if sig_mode == 1        str_alpha = strcat(num2str(alpha), 'hatch');  
        else                    str_alpha = num2str(alpha);                       end
    else
        str_alpha = [];
    end
    if compute_anombyAM && isempty(alpha)             error('Input alpha if compute_anombyAM is on!');     end
    if ismember(ldirect, {'fwd'}) && from_basin == 0  error('from_basin must be 1 for fwd mode!');         end
    if drywet_mode == 4 && mth_list ~= 0              error('mth_list must set 0 for drywet_mode == 4!');  end
    if plot_land_or_ocean ~= 0                      
        if plot_land_or_ocean == 1
            str_land_or_ocean = '_LandOnly';
        elseif plot_land_or_ocean == 2
            str_land_or_ocean = '_OceanOnly';
        end
    end
    if show_TP_only         str_TP_only = '_TP_only';     end
    if isempty(trend_test)
        str_trend_test = '_t-test';
    elseif isequal(trend_test, 'MK')
        str_trend_test = '_MK_test';
    else
        error('Invalid input of ''trend_test''!');
    end
    if sig_contf_in_hatch == 1   
        str_hatch = '_hatch';
        warning('The hatching bug in plot_contfmap.m (fails to plot holes) has not yet been resolved!');  
    end

    %/ loop over multiple sets of vars (if exist)
    th_l_mth = [];  cnt_l_mth = 0; all_titlename = {};
    for l = 1:size(slct_list, 1)
        cnt_l_mth        = cnt_l_mth + 1;
        slct_contf_str   = slct_list(l,1);
        slct_cont_str    = slct_list(l,2);
        slct_vector_str  = slct_list(l,3);
        slct_vector2_str = slct_list(l,4);
        slct_hatch_str   = slct_list(l,5);
        trend_mode       = slct_list{l,6};
        
        %/ determine what to load given the slct_contf_str
        if ismember(slct_contf_str, {'Pm_frac'})
            slct_contf_str_WSV = {'Pm'};
            slct_contf_str_an = {''};
            
        elseif ismember(slct_contf_str, {'Pm_frac_adj'})
            slct_contf_str_WSV = {'Pm'};
            slct_contf_str_an = {'P'};
            
        elseif ismember(slct_contf_str, {'Cf_map_frac'})
            slct_contf_str_WSV = {'Cf_map'};
            slct_contf_str_an = {''};
            
        elseif contains(slct_contf_str, {'_diff'})
            slct_contf_str_WSV = {''};
            str_parts          = split(slct_contf_str{:}(1:end-5), '_');
            slct_contf_str_an  = [str_parts(1), strjoin(str_parts(2:end),'_')];

        elseif ismember(slct_contf_str, valid_WSV_dataname)
            slct_contf_str_WSV = slct_contf_str;
            slct_contf_str_an  = {''};
        else
            slct_contf_str_WSV = {''};
            slct_contf_str_an  = slct_contf_str;
        end

        %/ determine what to load given the slct_cont_str
        if ismember(slct_cont_str, valid_WSV_dataname)
            slct_cont_str_WSV = slct_cont_str;
            slct_cont_str_an  = {''};
        else
            slct_cont_str_WSV = {''};
            slct_cont_str_an  = slct_cont_str;
        end

        if ismember(slct_hatch_str, {'BL_Pm_ratio'})
            slct_hatch_str_WSV = {'Pm', 'BL_Pm'};
            slct_hatch_str_an  = {''};
        else
            slct_hatch_str_WSV = {''};
            slct_hatch_str_an  = slct_hatch_str;  %/ by default.
        end
        
        th_mth = []; cnt_mth = 0; str_trend = [];
        for mth = mth_list
            cnt_mth = cnt_mth + 1;
                                
            %/ Dates
            if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
                str_dates = sprintf('%d%02d%02d-%d%02d%02d', year_list_bc(1), st_month, st_day, year_list_bc(end), ed_month, ed_day);
            else
                if mth == 0        str_dates = strcat(str_years_bc);
                else               str_dates = strcat(str_years_bc,'_', str_mth{mth});    end
            end
            %/ Define bndry_data, color, str_basin, etc. 
            %/ (Place the code block here to avoid redundant appending)
            if from_basin
                bndry_data    = {basin_bndry};
                color         = [255 51 204]./255;    %/ source region boundary color      
                str_basin     = sprintf('_%s', basin_name);
                linewi = 2;
            else
                linewi = 2;
                basin_name = '';
                str_basin = '';
            end
            if show_TP
                lon = WSV_lon;
                lat = WSV_lat;
                slct_reg = 'TP';
                [~, TP_bndry, ~, ~, ~] = reg_extractor('lon', lon, 'lat', lat, 'slct_reg', slct_reg,...
                                                        'data_folder', data_folder, 'savemat', 1, 'recompute', 0);
                bndry_data = [bndry_data; TP_bndry];
                color      = [color; [0 0 0]./255];
            end
            if show_sink_region
                bndry_data = [bndry_data; basin_bndry];
                color      = [color; [255 0 0]./255];
            end
    
            %============= Initialization ===================%
            contf_data = []; cont_data = []; contf_data_daily = []; contf_data_sig = []; cont_data_daily = [];...
            Udata = []; Vdata = []; U2data = []; V2data = []; hatch_data = []; point_data = [];
            Udata_daily = []; Vdata_daily = []; U2data_daily = []; V2data_daily = []; hatch_data_daily = [];
            
            %============= Read Reanalysis Fields ===================%
            slct_an_data     = [slct_contf_str,    slct_cont_str,    slct_hatch_str,    strcat('u',slct_vector_str), strcat('v',slct_vector_str), strcat('u',slct_vector2_str), strcat('v',slct_vector2_str)]; 
            slct_an_data_raw = {slct_contf_str_an, slct_cont_str_an, slct_hatch_str_an, strcat('u',slct_vector_str), strcat('v',slct_vector_str), strcat('u',slct_vector2_str), strcat('v',slct_vector2_str)};     
            for j = 1:length(slct_an_data)
                if findismember_loop(slct_an_data_raw{j}(:), {'P_LA', 'Pm', 'BL_Pm', 'Cf_map'})  continue;  end
                
                if isempty([slct_an_data_raw{j}{:}])   continue;    end
                if mth == 0
                    mean_field_bc   = 'annual_mean';
                else
                    mean_field_bc   = strcat(str_mth{mth},'_mean');
                    subset_field_bc = strcat(str_mth{mth},'_', select_field);
                end
                select_data = findismember_loop(dataname, unique(slct_an_data_raw{j}, 'stable')); %/ unique() to avoid bug.
                for k = select_data
                    str_years_meteo = str_years_bc;
                    if isfield(dataset, dataname{k})                                   %/ skip it if it has been loaded.
                        disp(['!!! ', dataname{k}, ' has already been loaded. Skip data retrieval. !!!'])
                    else
                        if contains(dataname{k}, {'IVT', 'VMT', 'dwdt'}) 
                            str_lnP_coor    = '_lnP';
                            str_linear_mode = '_linear';
                            dataname_bc     = dataname_callfile{k};
                        else
                            str_lnP_coor    = '';
                            str_linear_mode = '';
                            dataname_bc     = dataname{k};
                        end
                        domain_str = '_global';
                        loadfile = strcat(data_folder, dataname_bc, '_', select_field, domain_str, str_years_meteo, str_lnP_coor, str_linear_mode, '.mat');
                        disp(['Loading ', loadfile, ' ...']);
                        load(loadfile, 'S');

                        fld = fieldnames(S);
                        for f = 1:length(fld)
                            dataset.(dataname{k}).(fld{f}) = S.(fld{f});
                        end
                        clear S;
                    end
                    
                    %/ Compute_anom outputs data on the correct dates, and seasonal anomalies if queried
                    %/ NOTE: In CERA-20C, the 1st day of each year is missing for acc var.
                    [X, X_mean, X_mean_sig, date_subset] = compute_anom('dataset', dataset, 'var', dataname{k}, 'select_field', select_field,...
                                                                        'mth', mth, 'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day,...
                                                                        'year_list', year_list_bc, 'compute_anombyAM', compute_anombyAM, 'alpha', alpha);
                    

                    %/ You may select some dates for data checking                                     
                    if ~isempty(slct_date)
                        ind_date = findismember_loop(date_subset, slct_date);
                        X_mean = X(:,:,ind_date);
                        str_dates = strcat('_', join(num2str(date_subset(ind_date)), '_'));  %/ update str_dates
                    end
                    
                    %/ For GPCC/GLEAM data, if the mean is strictly 0, that means missing value 
                    if contains(dataname{k}, {'GPCC', 'GLEAM'}) 
                        cond = (X_mean == 0);
                        X_mean(cond)     = nan;
                        X_mean_sig(cond) = nan;
                    end
                    
                    %/ if want to show land or ocean grid only
                    if plot_land_or_ocean ~= 0
                        X_mean = show_land_or_ocean_hydrosheds('matrix', X_mean, 'lon_grids', dataset.(dataname{k}).lon, 'lat_grids', dataset.(dataname{k}).lat,...
                                                               'land_or_ocean', plot_land_or_ocean);
                    end
                    dataset.(dataname{k}).(mean_field_bc) = X_mean;
                    dataset.(dataname{k}).date_subset = date_subset;
                end
                
                %/ Save the advanced var into dataset and also as contf_data/cont_data/Udata/Vdata/hatch_data
                if ~isempty(select_data) && ~contains(slct_an_data(j), {'P_LA', 'Pm'})
                    %/ Post-processing (if relevant)
                    if contains(slct_an_data(j), {'_diff'}) %/ Difference of the mean
                        data1     = dataset.(dataname{select_data(1)}).(mean_field_bc);
                        data1_lon = dataset.(dataname{select_data(1)}).lon;
                        data1_lat = dataset.(dataname{select_data(1)}).lat;
                        
                        data2     = dataset.(dataname{select_data(2)}).(mean_field_bc);
                        data2_lon = dataset.(dataname{select_data(2)}).lon;
                        data2_lat = dataset.(dataname{select_data(2)}).lat;
                        
                        if contains(slct_an_data(j), {'_prct_diff'})
                            show_prct_diff = 1;
                        else
                            show_prct_diff = 0;
                        end

                        %/ 'diff_2D_landmask' compute the pcorr, MSE and
                        %/ SSIM of data1 and data2 in a given region (omitting NaN)
                        domain_lon_range = [min(data1_lon), max(data1_lon)];
                        domain_lat_range = [min(data1_lat), max(data1_lat)];
                        [X_mean, SSIM, r, pval, RMSE, X_lon, X_lat]  = diff_2D_landmask('data1', data1, 'data1_lon', data1_lon, 'data1_lat', data1_lat,...
                                                                         'data2', data2, 'data2_lon', data2_lon, 'data2_lat', data2_lat,...
                                                                         'show_prct_diff', show_prct_diff, 'interp_method', 'linear',...
                                                                         'domain_lon_range', domain_lon_range, 'domain_lat_range', domain_lat_range,...
                                                                         'show_TP_only', show_TP_only, 'data_folder', data_folder);
                        X = []; %/ After getting the difference in mean, no need to keep the daily difference.
                    
                    elseif contains(slct_an_data{j}, {'uVWS','vVWS'})
                        X_lon  = dataset.(dataname{select_data(1)}).lon;
                        X_lat  = dataset.(dataname{select_data(1)}).lat;
                        if contains(slct_an_data{j}, '_sig')  %/ then we show only the significant values
                            if length(dataset.(dataname{select_data(2)}).date_subset) ...
                                < length(dataset.(dataname{select_data(1)}).date_subset) 
                                ii = 2;
                            else
                                ii = 1;
                            end
                            ind_date_1  = findismember_loop(dataset.(dataname{select_data(1)}).date_yyyymmdd_AllYr,...
                                                            dataset.(dataname{select_data(ii)}).date_subset);
                            ind_date_2  = findismember_loop(dataset.(dataname{select_data(2)}).date_yyyymmdd_AllYr,...
                                                            dataset.(dataname{select_data(ii)}).date_subset);
                            %/ Compute daily VWS
                            data = dataset.(dataname{select_data(1)}).(select_field)(:,:,ind_date_1)...
                                    - dataset.(dataname{select_data(2)}).(select_field)(:,:,ind_date_2);

                            X_mean = ttest_sig_fn(data,0.01,3);  %/ alpha = 0.01 (dy default)
                        else
                            data1  = dataset.(dataname{select_data(1)}).(mean_field_bc);
                            data2  = dataset.(dataname{select_data(2)}).(mean_field_bc);
                            X_mean = data1 - data2;
                        end
                        
                    elseif contains(slct_an_data{j}, {'totVWS', 'uv300x850'})
                        X_lon  = dataset.(dataname{select_data(1)}).lon;
                        X_lat  = dataset.(dataname{select_data(1)}).lat;
                        
                        %/ Assume the order to be in 'u300', 'u850', 'v300', 'v850'
                        if contains(slct_an_data{j}, {'_sig', '_regime', '_ttest_ci'})  %/ then we show only the significant values
                            if length(dataset.(dataname{select_data(2)}).date_subset) ...
                                < length(dataset.(dataname{select_data(1)}).date_subset) 
                                ii = 2;
                            else
                                ii = 1;
                            end
                            ind_date_1  = findismember_loop(dataset.(dataname{select_data(1)}).date_yyyymmdd_AllYr,...
                                                            dataset.(dataname{select_data(ii)}).date_subset);
                            ind_date_2  = findismember_loop(dataset.(dataname{select_data(2)}).date_yyyymmdd_AllYr,...
                                                            dataset.(dataname{select_data(ii)}).date_subset);
                            ind_date_3  = findismember_loop(dataset.(dataname{select_data(3)}).date_yyyymmdd_AllYr,...
                                                            dataset.(dataname{select_data(ii)}).date_subset);
                            ind_date_4  = findismember_loop(dataset.(dataname{select_data(4)}).date_yyyymmdd_AllYr,...
                                                            dataset.(dataname{select_data(ii)}).date_subset);
           
                            if contains(slct_an_data{j}, {'totVWS'})        %/ Compute daily VWS
                                data   = (dataset.(dataname{select_data(1)}).(select_field)(:,:,ind_date_1)...
                                        - dataset.(dataname{select_data(2)}).(select_field)(:,:,ind_date_2))...
                                        + (dataset.(dataname{select_data(3)}).(select_field)(:,:,ind_date_1)...
                                        - dataset.(dataname{select_data(4)}).(select_field)(:,:,ind_date_2));
                                [X_mean, X_ttest_ci] = ttest_sig_fn(data,0.01,3);  %/ alpha = 0.01 (dy default)
                                
                            elseif contains(slct_an_data{j}, {'uv300x850'}) %/ Compute daily dot product of uv300 and uv850
                                data   = (dataset.(dataname{select_data(1)}).(select_field)(:,:,ind_date_1).*dataset.(dataname{select_data(2)}).(select_field)(:,:,ind_date_2))...
                                         + (dataset.(dataname{select_data(3)}).(select_field)(:,:,ind_date_1).*dataset.(dataname{select_data(4)}).(select_field)(:,:,ind_date_2));
                                [X_mean, X_ttest_ci] = ttest_sig_fn(data,0.01,3);  %/ alpha = 0.01 (dy default)
                                
                                %/ Determine regimes
                                if isequal(slct_an_data{j}, 'uv300x850_regime')
                                    u300_mean_sig = ttest_sig_fn(dataset.(dataname{select_data(1)}).(select_field)(:,:,ind_date_1),0.01,3);
                                    cond_I_pve = (X_mean > 0); 
                                    cond_I_nve = (X_mean < 0); 
                                    cond_u300_pve = (u300_mean_sig > 0);
                                    cond_u300_nve = (u300_mean_sig < 0);
                                    regime = nan(size(X_mean));
                                    regime(cond_I_pve & cond_u300_pve) = 1; %/ Westerly 
                                    regime(cond_I_pve & cond_u300_nve) = 2; %/ Easterly 
                                    regime(cond_I_nve)                 = 3; %/ Overturning 
                                    regime_name = {'Westerly', 'Easterly', 'Overturning'};
                                    regime_lon  = X_lon;
                                    regime_lat  = X_lat;
                                    X_mean      = regime;
                                    
                                    regime_filename = strcat(data_folder, slct_an_data{j}, str_years_bc, '_', str_mth{mth}, '.mat');
                                    if savemat
                                        save(regime_filename, 'regime', 'regime_name', 'regime_lon', 'regime_lat', '-v7.3');
                                        fprintf('!!! Regime matrix saved into %s. !!!\n', regime_filename);
                                    end
                                elseif isequal(slct_an_data{j}, 'uv300x850_ttest_ci')
                                    X_mean = X_ttest_ci;
                                end
                            end
                        else
                            data1  = dataset.(dataname{select_data(1)}).(mean_field_bc);
                            data2  = dataset.(dataname{select_data(2)}).(mean_field_bc);
                            
                            if ismember(slct_an_data{j}, {'totVWS'})        %/ Compute daily VWS
                                X_mean = data1 - data2;
                            elseif ismember(slct_an_data{j}, {'uv300x850'}) %/ Compute daily product of u300 and u850
                                X_mean = data1.*data2;
                            end
                        end
                    elseif contains(slct_an_data{j}, {'u300x850'})
                        X_lon  = dataset.(dataname{select_data(1)}).lon;
                        X_lat  = dataset.(dataname{select_data(1)}).lat;
                        
                        %/ Assume the order to be in 'u300', 'u850', 
                        if contains(slct_an_data{j}, {'_sig', '_regime', '_ttest_ci'})  %/ then we show only the significant values
                            if length(dataset.(dataname{select_data(2)}).date_subset) ...
                                < length(dataset.(dataname{select_data(1)}).date_subset) 
                                ii = 2;
                            else
                                ii = 1;
                            end
                            ind_date_1  = findismember_loop(dataset.(dataname{select_data(1)}).date_yyyymmdd_AllYr,...
                                                            dataset.(dataname{select_data(ii)}).date_subset);
                            ind_date_2  = findismember_loop(dataset.(dataname{select_data(2)}).date_yyyymmdd_AllYr,...
                                                            dataset.(dataname{select_data(ii)}).date_subset);
           
                            data   = (dataset.(dataname{select_data(1)}).(select_field)(:,:,ind_date_1).*dataset.(dataname{select_data(2)}).(select_field)(:,:,ind_date_2));
                            [X_mean, X_ttest_ci] = ttest_sig_fn(data,0.01,3);  %/ alpha = 0.01 (dy default)

                            %/ Determine regimes
                            if isequal(slct_an_data{j}, 'u300x850_regime')
                                u300_mean_sig = ttest_sig_fn(dataset.(dataname{select_data(1)}).(select_field)(:,:,ind_date_1),0.01,3);
                                cond_I_pve = (X_mean > 0); 
                                cond_I_nve = (X_mean < 0); 
                                cond_u300_pve = (u300_mean_sig > 0);
                                cond_u300_nve = (u300_mean_sig < 0);
                                regime = nan(size(X_mean));
                                regime(cond_I_pve & cond_u300_pve) = 1; %/ Westerly 
                                regime(cond_I_pve & cond_u300_nve) = 2; %/ Easterly 
                                regime(cond_I_nve)                 = 3; %/ Overturning 
                                regime_name = {'Westerly', 'Easterly', 'Overturning'};
                                regime_lon  = X_lon;
                                regime_lat  = X_lat;
                                X_mean      = regime;

                                regime_filename = strcat(data_folder, slct_an_data{j}, str_years_bc, '_', str_mth{mth}, '.mat');
                                if savemat
                                    save(regime_filename, 'regime', 'regime_name', 'regime_lon', 'regime_lat', '-v7.3');
                                    fprintf('!!! Regime matrix saved into %s. !!!\n', regime_filename);
                                end
                            elseif isequal(slct_an_data{j}, 'u300x850_ttest_ci')
                                X_mean = X_ttest_ci;
                            end
                        else
                            data1  = dataset.(dataname{select_data(1)}).(mean_field_bc);
                            data2  = dataset.(dataname{select_data(2)}).(mean_field_bc);
                            
                            if ismember(slct_an_data{j}, {'totVWS'})        %/ Compute daily VWS
                                X_mean = data1 -data2;
                            elseif ismember(slct_an_data{j}, {'uv300x850'}) %/ Compute daily product of u300 and u850
                                X_mean = data1.*data2;
                            end
                        end
                    else
                        X_lon = dataset.(dataname{select_data(1)}).lon;
                        X_lat = dataset.(dataname{select_data(1)}).lat;
                    end
                    
                    %/ Save the advanced var into dataset.
                    dataset.(slct_an_data{j}).(mean_field_bc)     = X_mean;
                    dataset.(slct_an_data{j}).lon                 = X_lon;
                    dataset.(slct_an_data{j}).lat                 = X_lat;
                    dataset.(slct_an_data{j}).date_yyyymmdd_AllYr = dataset.(dataname{select_data(1)}).date_yyyymmdd_AllYr;
                    
                    %/ Only assign contf_data if slct_contf_str = slct_contf_str_an
                    if isequal(slct_contf_str{:}, slct_an_data{j}) && isempty(contf_data)
                        contf_data       = X_mean;
                        contf_data_daily = X;                   
                        contf_data_dates = date_subset;
                        contf_lon        = dataset.(slct_an_data{j}).lon;
                        contf_lat        = dataset.(slct_an_data{j}).lat;
                    end

                    if isequal(slct_cont_str{:}, slct_an_data{j}) && isempty(cont_data)
                        cont_data       = X_mean;
                        cont_data_daily = X;                   
                        cont_data_dates = date_subset;
                        cont_lon        = dataset.(slct_an_data{j}).lon;
                        cont_lat        = dataset.(slct_an_data{j}).lat;
                    end

                    if isequal(char(strcat('u',slct_vector_str)), slct_an_data{j}) && isempty(Udata)
                        if compute_anombyAM && ~isempty(alpha)     Udata = X_mean_sig;     
                        else                                       Udata = X_mean;              end
                        Udata_daily = X;
                        Udata_dates = date_subset;
                        uv_lon  = dataset.(slct_an_data{j}).lon;
                        uv_lat  = dataset.(slct_an_data{j}).lat;
                    end

                    if isequal(char(strcat('v',slct_vector_str)), slct_an_data{j}) && isempty(Vdata)
                        if compute_anombyAM && ~isempty(alpha)     Vdata = X_mean_sig;     
                        else                                       Vdata = X_mean;              end
                        Vdata_daily = X;
                        Vdata_dates = date_subset;
                        uv_lon  = dataset.(slct_an_data{j}).lon;
                        uv_lat  = dataset.(slct_an_data{j}).lat;
                    end

                    if isequal(char(strcat('u',slct_vector2_str)), slct_an_data{j}) && isempty(U2data)
                        if compute_anombyAM && ~isempty(alpha)     U2data = X_mean_sig;     
                        else                                       U2data = X_mean;              end
                        U2data_daily = X;
                        U2data_dates = date_subset;
                        uv2_lon  = dataset.(slct_an_data{j}).lon;
                        uv2_lat  = dataset.(slct_an_data{j}).lat;
                    end

                    if isequal(char(strcat('v',slct_vector2_str)), slct_an_data{j}) && isempty(V2data)
                        if compute_anombyAM && ~isempty(alpha)     V2data = X_mean_sig;     
                        else                                       V2data = X_mean;              end
                        V2data_daily = X;
                        V2data_dates = date_subset;
                        uv2_lon  = dataset.(slct_an_data{j}).lon;
                        uv2_lat  = dataset.(slct_an_data{j}).lat;
                    end
                    
                    if isequal(slct_hatch_str{:}, slct_an_data{j}) && isempty(hatch_data)
                        hatch_data       = X_mean;
                        hatch_data_daily = X;                   
                        hatch_data_dates = date_subset;
                        hatch_lon        = dataset.(slct_an_data{j}).lon;
                        hatch_lat        = dataset.(slct_an_data{j}).lat;
                    end
                end
                X = [];
            end
            
            %============= Read WaterSip Vars (retrieve_WSV) =================%
            WSV = []; %/ clear up WSV to spare memory first.
            flag_retrieve_WSV = 0; %/ check if any WSV product is retrieved. If so, emphasizing the RHc_dqc_scheme.
            slct_WSV_data        = [slct_contf_str, slct_cont_str, slct_hatch_str];         %/ advanced WSVs
            slct_raw_WSV_data    = {slct_contf_str_WSV, slct_cont_str_WSV, slct_hatch_str_WSV};     %/ raw WSVs in *nested* cells
            for j = 1:length(slct_WSV_data)
                if isempty([slct_raw_WSV_data{j}{:}])   continue;    end
                if ~isequal(select_field, 'daily')   
                    warning('Detected ''select_field'' is not ''daily''. retrieve_WSV may not be able to handle that!');
                end
                
                %/ Since some advanced WSVs require multiple raw WSVs
                ind_WSV = findismember_loop(valid_WSV_dataname, slct_raw_WSV_data{j});
                final_data_daily = []; final_data_dates = [];
                for i = 1:length(ind_WSV)
                    flag_retrieve_WSV = 1;
                    WSV_name = valid_WSV_dataname{ind_WSV(i)};
                    
                    WSV_data_filename = string(strcat(data_folder, 'WSV_', WSV_name, str_basin,...
                                              strrep(str_expmntinfo, ' ', '_'), str_dates, '.mat'));

                    %/ Load the WSV if exists and need not to recompute
                    if isfile(WSV_data_filename) && recompute == 0
                        if from_basin ~= 0
                            %/ ALWAYS load data since Pm field is different for different hs!
                            fprintf('!!! The queried %s for the %s region is found. Loading... !!! \n', WSV_name, basin_name);
                            L = load(WSV_data_filename, 'WSV');
                            flds = fieldnames(L.WSV);
                            for ii = 1:length(flds)
                                WSV.(flds{ii}) = L.WSV.(flds{ii});  %/ recursively load into WSV struct. Otherwsie there will always be one WSV field.
                            end
                        else
                            %/ No need to always load data if from_basin == 0;
                            if isfield(WSV, WSV_name)                                   %/ skip it if it has been loaded.
                                disp(['!!! ', WSV_name, ' has already been loaded. Skip data retrieval. !!!'])
                            else
                                fprintf('!!! The queried %s is found. Loading... !!! \n', WSV_name);
                                L = load(WSV_data_filename, 'WSV');
                                flds = fieldnames(L.WSV);
                                for ii = 1:length(flds)
                                    WSV.(flds{ii}) = L.WSV.(flds{ii});  %/ recursively load into WSV struct. Otherwsie there will always be one WSV field.
                                end
                            end
                        end
                    else
                        fprintf('*** %s Not Found ***\n', WSV_name)
                        fprintf('*** Now implementing retrieve_WSV... ***\n')
                        
                        [WSV_data_daily, WSV_dates, ~, ~] = retrieve_WSV('WSV_name', WSV_name,  'WSV_dir',  WSV_dir, 'year_list', year_list_bc, 'mth', mth,...
                                                                    'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day, 'str_RHc_dqc', str_RHc_dqc,...
                                                                    'ldirect', ldirect, 'from_basin', from_basin, 'maxtraj_day', maxtraj_day, 'str_optimal', str_optimal,...
                                                                    'dt_slct', dt_slct, 'str_traj_rm_jump', str_traj_rm_jump, 'str_BLH_factor', str_BLH_factor, 'str_remark', str_remark,...
                                                                    'str_src_z', str_src_z, 'str_domain', str_domain, 'str_domain_trajfile', str_domain_trajfile, 'str_sharpcut', str_sharpcut, 'forcing', forcing,...
                                                                    'NumWorkers', NumWorkers, 'basin_name', basin_name, 'basin_catalog', basin_catalog);

                        %/ WaterSip Var (WSV)
                        WSV.(WSV_name).(select_field) = WSV_data_daily;   
                        WSV.(WSV_name).slct_dates = WSV_dates;

                        if size(WSV_data_daily, 2) == length(WSV_lat(2:end-1))
                            WSV.(WSV_name).lon  = WSV_lon;
                            WSV.(WSV_name).lat  = WSV_lat(2:end-1);
                        elseif size(WSV_data_daily, 2) == length(WSV_lat)
                            WSV.(WSV_name).lon  = WSV_lon;
                            WSV.(WSV_name).lat  = WSV_lat;
                        else
                            error('lat dimension of WSV_data_daily does not fit WSV_lat or WSV_lat(2:end-1)!');
                        end
                        if from_basin  WSV.(WSV_name).basin_name = basin_name;  
                        else           WSV.(WSV_name).basin_name = '';              end
                        
                        %/ Save only the whole year or seasonal WSV.
                        if savemat
                            fprintf('*** Saving WSV data: %s *** \n', WSV_data_filename{:})
                            save(WSV_data_filename, 'WSV', '-v7.3');
                        end
                    end

                    if savenc
                        ncfilename = strrep(WSV_data_filename, '.mat', '.nc');
                        fprintf('*** Saving WSV data as NetCDF file: %s *** \n', ncfilename{:});

                        %/ Accumulated (as requested by Chris/Emme/Jessica
                        data = squeeze(sum(WSV.(WSV_name).(select_field), 3, 'omitnan'));
                        date = WSV.(WSV_name).slct_dates;

                        if isequal(WSV_name, 'Pm')
                            data_shortname    = 'Cb';
                            data_standardname = 'Cb';
                            data_units        = 'mm/day';
                            data_longname     = sprintf('Accumulated backward contributions to %s over %d-%d', basin_name, date(1), date(end));
                        elseif isequal(WSV_name, 'Cf_map')
                            data_shortname    = 'Cf';
                            data_standardname = 'Cf';
                            data_units        = 'mm/day';
                            data_longname     = sprintf('Accumulated forward contributions from %s; accumulating to the period %d-%d', basin_name, date(1), date(end));
                        elseif isequal(WSV_name, 'P_LA')
                            data_shortname    = 'P_LA';
                            data_standardname = 'P_LA';
                            data_units        = 'mm/day';
                            data_longname     = sprintf('Precipitation estimated by Lagrangian WaterSip in %s; averaged over the period %d-%d', basin_name, date(1), date(end));
                        else
                            error('Set the shortname, standardname, etc.!');
                        end
                        date_format = 'yyyymmdd';
                        remark      = strrep(str_expmntinfo, '_', ' ');  %/ Moisture tracking info
                        write_nc('ncfilename', ncfilename, 'data', data, 'data_shortname', data_shortname,...
                                 'data_standardname', data_standardname, 'data_longname', data_longname, 'data_units', data_units,...
                                 'lon', WSV.(WSV_name).lon, 'lat', WSV.(WSV_name).lat, 'plev', [], 'time', [], 'time_unit', [], ...
                                 'date', [], 'date_format', [], 'remark', remark, 'othervar', []);
                        ncdisp(ncfilename);
                    end
                    final_data_daily = WSV.(WSV_name).(select_field);
                    final_data       = mean(final_data_daily, 3, 'omitnan'); 
                end
                
                %/ Post-Processing of WSVs 
                if ~isempty(ind_WSV)
                    land_or_ocean = [];
                    if ismember(slct_WSV_data(j), {'Pm_frac'})
                        Pm               = WSV.('Pm').(select_field);
                        Area             = calc_grid_area('lon', WSV.('Pm').lon, 'lat', WSV.('Pm').lat);
                        final_data_daily = (Pm.*Area)./sum(Pm.*Area, [1, 2], 'omitnan')*100;  %/ Unit: %; Upscaling; Area-weighted
                        final_data       = mean(final_data_daily, 3, 'omitnan'); 
                        Pm = [];
                    elseif ismember(slct_WSV_data(j), {'Pm_frac_adj'})
                        Pm      = WSV.('Pm').(select_field);
                        Area    = calc_grid_area('lon', WSV.('Pm').lon, 'lat', WSV.('Pm').lat);
                        Pm_frac = (Pm.*Area)./sum(Pm.*Area, [1, 2], 'omitnan')*100;  %/ Unit: %; Upscaling; Area-weighted
                        
                        %/ 1. Get the regional sum of P (times grid area)
                        P_gridtime = reshape(dataset.('P').(select_field), [], size(dataset.('P').(select_field),3));
                        P_lon = dataset.('P').lon;
                        P_lat = dataset.('P').lat;
                        P_Area        = calc_grid_area('lon', P_lon, 'lat', P_lat);
                        P_Area_2Dto1D = reshape(P_Area, [], 1);
                        [P_lon_2D, P_lat_2D] = meshgrid(P_lon, P_lat);
                        P_lon_2D = P_lon_2D'; P_lat_2D = P_lat_2D';
                        P_lon_2Dto1D = reshape(P_lon_2D, [], 1);
                        P_lat_2Dto1D = reshape(P_lat_2D, [], 1);
                        if from_basin == 4 
                            if top == 0     %/ then show all the basins 
                                target_bndry = TP_bndry{:};
                            else
                                target_bndry = basin_catalog(top).bndry;
                            end
                        else
                            error('code not set!');
                        end
                        [in,~]   = inpoly2([P_lon_2Dto1D, P_lat_2Dto1D], target_bndry); %/ inpoly2 is 600xx faster than inpolygon!!
                        ind_in   = find(in == 1);
                        P_regsum = sum(P_gridtime(ind_in,:).*P_Area_2Dto1D(ind_in,:), 1)*1e-12; %/ mm/day m^2 -> km^3/day 
                        
                        %/ 2. Subset the overlapped dates betn Pm and P 
                        common_dates     = intersect(WSV.('Pm').slct_dates, dataset.('P').date_yyyymmdd_AllYr);
                        ind_Pm_frac_date = findismember_loop(WSV.('Pm').slct_dates,             common_dates);
                        ind_P_date       = findismember_loop(dataset.('P').date_yyyymmdd_AllYr, common_dates);
                        
                        Pm_frac     = Pm_frac(:,:,ind_Pm_frac_date);
                        P_regsum    = P_regsum(:,ind_P_date);
                        
                        %/ 4.3 Check no-rain or no-Pm days!! [IMPORTANT to avoid BUG!]
                        %/      1. Due to no rain on some days, 
                        %/          the total ratio contributions = 0. 
                        %/          This might cause a misleadingly large 'Elsewhere'!
                        %/          Set nans for those days!
                        %/      2. Since P_LA does NOT perfectly match P, for some small
                        %/          areas we may find non-zero Pm but zero P, or zero Pm but non-zero P!
                        %/          Workaround: remove all zero-P or zero-Pm days.
                        ind_no_rain = find(P_regsum == 0 | squeeze(nansum(Pm_frac, [1,2]))' == 0); %<- Important!
                        if ~isempty(ind_no_rain)
                            fprintf('*** [mth: %d] [basin: %s]: Detected %d days when P(basin) = 0 OR sum(Pm) = 0! Setting values on those days to be NaN. ***\n',...
                                    mth, basin_name, length(ind_no_rain));
                        end
                        Pm_frac(:,:,ind_no_rain) = nan; %/ Note that here Pm_frac is 3D
                        P_regsum(:,ind_no_rain) = nan;  %<- Don't forget to modify P data also!

                        P_regsum    = reshape(P_regsum, 1, 1, length(P_regsum)); %/ reshape to (1, 1, ndate)
                        Pm_frac_x_P = Pm_frac.*P_regsum;
                        
                        %/ 3. Accumulate over time
                        [Pm_frac_x_P_ts, ~, ~, ~] = daily2any('data_daily', Pm_frac_x_P, 'dates', common_dates, 'mth', mth, 'ins_or_acc', 'acc');                                       
                        [P_regsum_ts, ~, ~, ~]    = daily2any('data_daily', P_regsum,    'dates', common_dates, 'mth', mth, 'ins_or_acc', 'acc');
                        
                        final_data_daily = [];      
                        final_data       = nanmean(Pm_frac_x_P_ts./P_regsum_ts, 3);   %/ mean sum(P.*Pm_frac)./sum(P)*100%
                        fprintf('*** Global sum of Pm_frac_adj = %.2f%% (Should be 100%%) ***\n', nansum(final_data, 'all'));
                        Pm = []; Pm_frac = []; P_gridtime = [];     %/ Spare memory
                        
                    elseif ismember(slct_WSV_data(j), {'Cf_map_frac'})
                        Cf_map        = WSV.('Cf_map').(select_field);
                        Area             = calc_grid_area('lon', WSV.('Cf_map').lon, 'lat', WSV.('Cf_map').lat);
                        final_data_daily = (Cf_map.*Area)./nansum(Cf_map.*Area, [1, 2])*100; %/ consider grid area also
                        final_data       = nanmean(final_data_daily, 3); %/ NOTE; the global sum will still be 100%.
                        Cf_map = [];

                    elseif ismember(slct_WSV_data(j), {'BL_Pm_ratio'})
                        %/ NOTE: The ratio of BL_Pm to Pm.
                        final_data_daily = WSV.('BL_Pm').(select_field)./WSV.('Pm').(select_field)*100;
                        final_data       = nanmean(final_data_daily, 3);
                        
                    elseif ismember(slct_WSV_data(j), {'P_LA_P_diff', 'P_LA_P_prct_diff'})
                        P_LA = nanmean(WSV.('P_LA').(select_field), 3);
                        P_LA(P_LA == 0) = nan;
                        if mth == 0
                            P = dataset.('P').(strcat('annual_mean'));
                        else
                            P = dataset.('P').(strcat(str_mth{mth},'_mean'));
                        end
    %                         P    = nanmean(dataset.('P').(select_field), 3);
                        P    = P(:, 2:end-1);
                        if ismember(slct_WSV_data(j), {'P_LA_P_prct_diff'})
                            show_prct_diff = 1;
                        else
                            show_prct_diff = 0;
                        end

                        %/ 'diff_2D_landmask' helps correct the different lat ordering.
                        %/ So, use it to compute percentage bias!!
                        domain_lon_range = [min(WSV.('P_LA').lon), max(WSV.('P_LA').lon)];
                        domain_lat_range = [min(WSV.('P_LA').lat), max(WSV.('P_LA').lat)];
                        [final_data, SSIM, r, pval, RMSE]  = diff_2D_landmask('data1', P_LA, 'data1_lon', WSV.('P_LA').lon,  'data1_lat', WSV.('P_LA').lat,...
                                                                             'data2', P,    'data2_lon', dataset.('P').lon, 'data2_lat', dataset.('P').lat(2:end-1),...
                                                                             'show_prct_diff', show_prct_diff,...
                                                                             'domain_lon_range', domain_lon_range, 'domain_lat_range', domain_lat_range,...
                                                                             'show_TP_only', show_TP_only, 'data_folder', data_folder);

                    elseif ismember(slct_WSV_data(j), {'Pm_E_diff', 'Pm_E_prct_diff'})
                        Pm = nanmean(WSV.('Pm').(select_field), 3);
%                         Pm(Pm == 0) = nan;

                        %/ Since Pm for continental prcp can have non-trivial value over oceans
                        %/  we will consider only the land values using land-sea mask
                        land_or_ocean = 1;

                        if mth == 0
                            E = dataset.('E').(strcat('annual_mean'));
                        else
                            E = dataset.('E').(strcat(str_mth{mth},'_mean'));
                        end
    %                         E = nanmean(dataset.('E').(select_field), 3); 
                        E = E(:, 2:end-1);
                        if ismember(slct_WSV_data(j), {'Pm_E_prct_diff'})
                            show_prct_diff = 1;
                        else
                            show_prct_diff = 0;
                        end
                        
                        domain_lon_range = [min(WSV.('Pm').lon), max(WSV.('Pm').lon)];
                        domain_lat_range = [min(WSV.('Pm').lat), max(WSV.('Pm').lat)];

                        [final_data, SSIM, r, pval, RMSE]  = diff_2D_landmask('data1', Pm, 'data1_lon', WSV.('Pm').lon, 'data1_lat', WSV.('Pm').lat,...
                                                                              'data2', E,  'data2_lon', dataset.('E').lon, 'data2_lat', dataset.('E').lat(2:end-1),...
                                                                              'show_prct_diff', show_prct_diff,...
                                                                              'domain_lon_range', domain_lon_range, 'domain_lat_range', domain_lat_range,...
                                                                              'show_TP_only', show_TP_only, 'data_folder', data_folder);

                    elseif ismember(slct_WSV_data(j), {'BL_Pm_E_diff'})
                        BL_Pm = nanmean(WSV.('BL_Pm').(select_field), 3);
                        if mth == 0
                            E = dataset.('E').(strcat('annual_mean'));
                        else
                            E = dataset.('E').(strcat(str_mth{mth},'_mean'));
                        end
    %                         E = nanmean(dataset.('E').(select_field), 3);
                        E = E(:, 2:end-1);
                        domain_lon_range = [min(WSV.('Pm').lon), max(WSV.('Pm').lon)];
                        domain_lat_range = [min(WSV.('Pm').lat), max(WSV.('Pm').lat)];
                        [final_data, SSIM, r, pval, RMSE]  = diff_2D_landmask('data1', BL_Pm, 'data1_lon', WSV.('Pm').lon, 'data1_lat', WSV.('Pm').lat,...
                                                                              'data2', E,     'data2_lon', dataset.('E').lon, 'data2_lat', dataset.('E').lat(2:end-1),...
                                                                              'show_prct_diff', show_prct_diff,...
                                                                              'domain_lon_range', domain_lon_range, 'domain_lat_range', domain_lat_range,...
                                                                              'show_TP_only', show_TP_only, 'data_folder', data_folder);

                    end

                    %/ Save the advanced var into WSV.
                    WSV.(slct_WSV_data{j}).(select_field)    = final_data_daily;
                    WSV.(slct_WSV_data{j}).('data_for_plot') = final_data;
                    WSV.(slct_WSV_data{j}).lon               = WSV.(slct_contf_str_WSV{1}).lon;
                    WSV.(slct_WSV_data{j}).lat               = WSV.(slct_contf_str_WSV{1}).lat;
                    WSV.(slct_WSV_data{j}).slct_dates        = WSV.(slct_contf_str_WSV{1}).slct_dates;

                    %/ Assign to either contf_data or cont_data
                    if ismember(slct_WSV_data(j), slct_contf_str) && isempty(contf_data)
                        contf_data_daily = WSV.(slct_contf_str{:}).(select_field);
                        contf_data       = WSV.(slct_contf_str{:}).('data_for_plot');
                        contf_data_dates = WSV.(slct_contf_str{:}).slct_dates;

                        %/ Retrieve slct_contf_str_WSV{:} field from WSV,
                        %since some advanced var (e.g. P_LA_P_prct_diff) 
                        if isempty(alpha) && ~isempty(prct) && ismember(slct_contf_str, {'Pm', 'BL_Pm'})
                            prctile_val = prctile(reshape(contf_data, 1,[]), prct); %/ NOTE: universal prctile_val applied to each grid.
                            if contf_with_stipples
                                a = contf_data;                                         %/ dummy var
                                a(a <= prctile_val) = nan;                              %/ filtered 

                                [ind_lon, ind_lat] = find(~isnan(a));
                                point_data = [WSV.(slct_contf_str{:}).lon(ind_lon); WSV.(slct_contf_str{:}).lat(ind_lat)]';

                            elseif contf_with_hatch
                                hatch_data = contf_data;
                                hatch_lon = WSV.(slct_contf_str{:}).lon;
                                hatch_lat = WSV.(slct_contf_str{:}).lat; 
                                hatch_thres_pve =    prctile_val;
                                hatch_thres_nve = -1*prctile_val;
                                color_hatch_pve = [.4 .4 .4]; 
                                color_hatch_nve = [1 1 1];                     %/ not shown anyway.
                                hatch_mode = 1;
                            end
                        end

                    elseif ismember(slct_WSV_data(j), slct_cont_str) && isempty(cont_data)
                        cont_data_daily  = WSV.(slct_cont_str{:}).(select_field);
                        cont_data        = WSV.(slct_cont_str{:}).('data_for_plot');
                        cont_data_dates  = WSV.(slct_cont_str{:}).slct_dates;  %/ for trend_mode
                        
                    elseif ismember(slct_WSV_data(j), slct_hatch_str) && isempty(hatch_data)
                        hatch_data       = WSV.(slct_hatch_str{:}).('data_for_plot');
                        hatch_data_dates = WSV.(slct_hatch_str{:}).slct_dates;  %/ for trend_mode
                        hatch_lon        = WSV.(slct_hatch_str{:}).lon;
                        hatch_lat        = WSV.(slct_hatch_str{:}).lat;
                    else
                        error('what WSV data it is?  Check the code.')
                    end

                    if ~isempty(SSIM)
                        str_diff_metrics = sprintf(' (SSIM=%.2f r=%.2f RMSE=%.2f)', SSIM, r, RMSE); %/ This str will be used.
                        disp(str_diff_metrics)
                    end
                end
            end
            
            str_scheme = '';
            if flag_retrieve_WSV
                str_scheme = sprintf(' scheme%d', RHc_dqc_scheme);
            end
            
            %======== contf setting =======%
            contf_plus     = 0; contf_multiply = 1;   %/ default
            cbar_interval  = 2;
            cbar_YTick     = []; cbar_YTickLabel = []; 
            backcolor      = 'none';
            coast_col = [.1 .1 .1];  
            if show_TP_only   
                pcolor_mode = 1;    
            else
                pcolor_mode = 0;
            end
            if ismember(slct_contf_str, {'P_LA', 'CERA_P',  'ERA5_P', 'ERA5_P_05', 'CMORPH_P', 'IMERG_P', 'GPCC_P', 'GPCP_P', 'CRU_P', 'HARv2_P', 'TPR_P', 'EM_P',...
                                         'CERA_E', 'ERA5_E', 'GLEAM_E', 'HARv2_E', 'TPR_E', 'EM_E'})
                if ismember(slct_contf_str, {'P_LA'})
                    contf_lon  = WSV.(slct_contf_str{:}).lon;
                    contf_lat  = WSV.(slct_contf_str{:}).lat;
                    contf_data(contf_data == 0) = nan;  %/ testing
                else
                    contf_lon  = dataset.(slct_contf_str{:}).lon;
                    contf_lat  = dataset.(slct_contf_str{:}).lat;
                end
                if trend_mode
                    contf_multiply = 10;       %/ mm/yr -> mm/decade
                    if mth == 0
                        contf_unit = 'mm/yr/decade';
                        contf_levels = -60:10:60;       %/ for trend acc: mm/yr
                    elseif any(ismember(13:18, mth))
                        contf_unit = 'mm/season/decade';  
                        contf_levels = (-60:10:60)/2;     %/ for trend acc: mm/yr
                    else
                        error('code not set');
                    end
                    colmap = my_colormap(length(contf_levels)-1, 'BuYlRd_12lev', 'flip');
%                     colmap = brewermap(length(contf_levels)+1, 'RdBu');
%                     colmap([1,end],:) = [];
                else
                    if set_PE_mm_per_yr
                        contf_multiply = 365.25;
                        contf_levels = 0:300:3300;
                        colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                        contf_unit = 'mm/yr';
                    else
                        if isequal(select_field, 'monthly')
                            contf_multiply = 1/30;  %/ get the average daily prcp rate from the monthly rate
                            contf_levels = [0, 0.5, 1, 1.5, 2, 3, 5, 8, 12, 17, 23, 30];
                            colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                            contf_unit = 'mm/day';
                        else
                            contf_levels = 0:2:22;
                            colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                            contf_unit = 'mm/day';
                        end
                    end
                    % colmap(1,:) = [.7 .7 .7];
                    cbar_interval = 3;
                end

            elseif contains(slct_contf_str, {'IVTdiv', 'VMT'})  
                if trend_mode
                    contf_multiply = 10*24*3600;   %/ for trend acc: mm/yr/yr to mm/yr/decade
                    contf_levels = [-30:5:30];      %/ for trend acc: mm/yr
                    contf_unit = 'mm/yr/decade';
                    colmap = flip(nclCM('precip_diff_12lev', length(contf_levels)-1), 1); %/ Bluish -> -ve (convergence)
                else
                    contf_multiply = 24*3600;
                    contf_levels = [-6:1:6]; %/ mm/day
                    contf_unit = 'mm/day';
                    colmap = brewermap(length(contf_levels)+1, '*BrBG');
                    colmap([1,end],:) = [];  %/ remove the darkest colors at two ends
                end

            elseif ismember(slct_contf_str, {'BL_Pm', 'Pm', 'Pm_AR'})
                contf_lon  = WSV.(slct_contf_str{:}).lon;
                contf_lat  = WSV.(slct_contf_str{:}).lat;
%                 contf_data(contf_data == 0) = nan; %/ testing
%                 coast_col  = 'w';
                coast_col  = 'k';
%                 coast_col = [0.15 0.15 0.15];
                if trend_mode
                    contf_multiply = 10;       %/ unit per yr -> unit per decade
                    if mth == 0
                        contf_unit   = 'mm/yr/decade';
                        contf_levels = -12:2:12;

                    elseif any(ismember(13:18, mth))
                        contf_unit   = 'mm/season/decade';  %/ trend_mode will be based on 'acc' -> trend of the annual/seasonal sum 
                        contf_levels = (-12:2:12)/2;
                    else
                        error('code not set');
                    end
                    colmap = my_colormap(length(contf_levels)-1, 'BuYlRd_12lev', 'flip');
                    cbar_interval = 2;
                    coast_col   = [51 255 0]./255;
                    
                elseif from_basin
                    contf_data(contf_data == 0) = nan;
                    contf_levels = [0, 0.005, 0.02, 0.05, 0.1, 0.25, 0.5, 1, 1.5, 2, 3, 4]; %/ uneven level
                    colmap = my_colormap(length(contf_levels)-1,'precip3_11lev');
%                     colmap(1,:) = [.85 .85 .85];
                    colmap(1,:) = [1 1 1];
%                     cbar_YTick      = contf_levels(2:2:end-1); 
%                     cbar_YTickLabel = cbar_YTick;
                    contf_unit = 'mm/day';
                else
                    if isequal(select_field, 'monthly')
                        contf_levels = (0:1:11)*30;
                        contf_unit = 'mm/month';
                    else
                        contf_levels = 0:1:11;
                        contf_unit = 'mm/day';
                    end
                    colmap = brewermap(length(contf_levels)-1, '*RdYlBu');
                    colmap(1,:) = [.7 .7 .7];
                    cbar_interval = 3;
                end
                pcolor_mode = 1;

            elseif ismember(slct_contf_str, {'Cf_map'})
                contf_lon  = WSV.(slct_contf_str{:}).lon;
                contf_lat  = WSV.(slct_contf_str{:}).lat;

                contf_levels = [0, 0.005, 0.02, 0.05, 0.1, 0.25, 0.5, 1, 1.5, 2, 3, 4]; %/ uneven level
                contf_unit = 'mm/day';
                cbar_interval = 3;
                colmap = brewermap(length(contf_levels)-1, '*RdYlBu');
                colmap(1,:) = [1 1 1];
                
            elseif ismember(slct_contf_str, {'Pm_frac', 'Pm_frac_adj'})
                contf_lon  = WSV.(slct_contf_str{:}).lon;
                contf_lat  = WSV.(slct_contf_str{:}).lat;
                
                if trend_mode
                    contf_multiply  = 10;  %/ unit per yr -> unit per decade
                    contf_unit      = '%/decade';  %/ trend_mode will be based on 'acc' -> trend of the annual/seasonal sum
                    contf_levels    = [-0.012:0.002:0.012];
                    colmap          = my_colormap(length(contf_levels)-1, 'BuYlRd_12lev', 'flip');
                    cbar_YTick      = contf_levels(2:2:end-1); 
                    cbar_YTickLabel = cbar_YTick;
                    
                elseif from_basin 
                    contf_data(contf_data == 0) = nan;
                    contf_levels = [0, 0.001, 0.005, 0.05:0.05:0.4, 1, 1.5]; %/ uneven level
                    contf_unit = '%';
                    colmap = my_colormap(length(contf_levels)-1,'radar_12lev');
                    colmap(1,:) = [.85 .85 .85];
                    cbar_YTick      = contf_levels(2:2:end-1); 
                    cbar_YTickLabel = cbar_YTick;
                else
                    contf_levels = [0:0.25:7.25]/10;  %/ in %
                    contf_unit = '';
                    colmap = brewermap(length(contf_levels)-1, '*RdYlBu');
                    colmap(1,:)   = [1 1 1];
                    cbar_YTick      = [0:1:7]/10; 
                    cbar_YTickLabel = cbar_YTick;
                end
                cbar_interval = 3;
%                 backcolor     = [.7 .7 .7];

                %/ How much fraction do all grid cells > xx% explain?
                %/    Analogous to Precipitationshed
                a = nansum(contf_data(contf_data > contf_levels(2)));
                fprintf('*** The sum of %s in all grid cells > %.3f%% = %.2f%%  ***\n', slct_contf_str{:}, contf_levels(2), a);

            elseif ismember(slct_contf_str, {'Cf_map_frac'})
                contf_lon  = WSV.(slct_contf_str{:}).lon;
                contf_lat  = WSV.(slct_contf_str{:}).lat;

                contf_levels = [0:0.25:7.25]/10;  %/ in %
                contf_unit = '%';
                colmap = brewermap((length(contf_levels)-1), '*RdYlGn');
                colmap(1,:)   = [.7 .7 .7];
                cbar_YTick      = [0:1:7]/10; 
                cbar_YTickLabel = cbar_YTick;
                cbar_interval = 3;
%                 backcolor     = [.7 .7 .7];

            elseif ismember(slct_contf_str, {'optimal_trajtime'})
                contf_lon  = WSV.(slct_contf_str{:}).lon;
                contf_lat  = WSV.(slct_contf_str{:}).lat;
                if trend_mode
                    contf_plus     = 0;          %/ no need
                    contf_multiply = dt_slct*10; %/ hour -> hour/decade    
                    contf_unit     = 'hour/decade';
                    contf_levels   = [-12:2:12];
%                     contf_multiply = dt_slct*10; %/ day -> day/decade    
%                     contf_unit     = 'day/decade';
%                     contf_levels   = [-0.6:0.1:0.6];
                    cbar_interval  = 2;
                    colmap         = my_colormap(length(contf_levels)-1, 'sunshine_12lev_white');
                    pcolor_mode    = 1;
                else
                    contf_plus      = -1;  %/ -1 means -3h, since two time steps give one duration.
                    contf_multiply  = dt_slct/24;  %/ day
                    contf_unit      = 'day';
                    contf_levels    = [0:1:12]*2;
                    cbar_YTick      = 0:2:contf_levels(end-1);
%                     contf_levels    = 0:1:12;
%                     cbar_YTick      = [0, 3, 7, 11];  %/ corresponds to different colorings in 'radar_12lev'
                    cbar_YTickLabel = cbar_YTick;
                    colmap = my_colormap(length(contf_levels)-1, 'radar_12lev');
%                     colmap = m_colmap('jet', length(contf_levels)-1);
                    pcolor_mode = 1;
                end
                
            elseif contains(slct_contf_str, '_diff')
                if contains(slct_contf_str, {'Pm', 'P_LA'})
                    contf_lon  = WSV.(slct_contf_str{:}).lon;
                    contf_lat  = WSV.(slct_contf_str{:}).lat;
                end
                if contains(slct_contf_str, '_prct_diff') 
                    contf_unit = '%';
                    contf_levels = -120:20:120;
                elseif isequal(select_field, 'daily')
                    contf_unit = 'mm/day';
                    if show_TP_only
                        contf_levels = (-6:1:6)./2;
                    else
                        contf_levels = -6:1:6;
                    end
                elseif isequal(select_field, 'monthly')
                    contf_unit = 'mm/month';
                    contf_levels = (-6:1:6)*30;
                else
                    error('contf_levels not set!')
                end
                
                colmap = my_colormap(length(contf_levels)-1, 'BuYlRd_12lev');
                pcolor_mode = 1;
            else
                error('contf_levels and contf_unit have not set!!')
            end
            
            %======= trend mode ======% 
            if trend_mode  
                %/ If 'acc', then compute annual *total*. 
                %/ If 'ins', then compute annual *average*. 
                ListOfVar4acc = {'P', 'CMORPH_P', 'IMERG_P', 'ERA5_P', 'ERA5_P_05', 'HARv2_P', 'GPCC_P', 'GPCP_P', 'CRU_P', 'TPR_P', 'EM_P',...
                                 'E', 'ERA5_E', 'HARv2_E', 'GLEAM_E', 'TPR_E', 'EM_E', 'IVTdiv', 'VMT', 'dwdt', 'Pm', 'BL_Pm', 'P_LA', 'RO', 'SRO', 'SSRO',...
                                 'ERA5_RO', 'ERA5_SRO', 'ERA5_SSRO', 'TPR_RO', 'TPR_SRO', 'TPR_SSRO', 'SMLT', 'Es'};
                if ismember(slct_contf_str, ListOfVar4acc)  ins_or_acc_contf = 'acc'; else ins_or_acc_contf = 'ins';  end
                if ismember(slct_cont_str,  ListOfVar4acc)  ins_or_acc_cont  = 'acc'; else ins_or_acc_cont  = 'ins';  end
                ins_or_acc_vec = 'ins';  %/ Assume vectors are always ins field.
                if ~isempty(contf_data_daily)
                    if ~isempty(findismember_loop(slct_contf_str, {'Pm', 'BL_Pm', 'P_LA'}))
                        data_name = strcat(slct_contf_str,str_basin);
                    else
                        data_name = slct_contf_str;    %/ no need to specify basin 
                    end
                    [contf_data, contf_data_sig, ins_or_acc] = compute_trend('data_daily', contf_data_daily, 'ins_or_acc', ins_or_acc_contf, 'lon', contf_lon, 'lat', contf_lat, 'data_dates', contf_data_dates, 'data_name', data_name,...
                                                                             'year_list', year_list_bc, 'mth', mth, 'trend_test', trend_test, 'trend_alpha', trend_alpha, 'save_trend', save_trend, 'recompute_trend', recompute_trend, 'data_folder', data_folder);
                end
                if ~isempty(cont_data_daily)
                    if ~isempty(findismember_loop(slct_cont_str, {'Pm', 'BL_Pm', 'P_LA'}))
                        data_name = strcat(slct_cont_str,str_basin);
                    else
                        data_name = slct_cont_str;    %/ no need to specify basin 
                    end
                    [     ~, cont_data, ins_or_acc] = compute_trend('data_daily', cont_data_daily, 'ins_or_acc', ins_or_acc_cont,  'lon', cont_lon, 'lat', cont_lat, 'data_dates', cont_data_dates, 'data_name', data_name,...
                                                                    'year_list', year_list_bc, 'mth', mth, 'trend_test', trend_test, 'trend_alpha', trend_alpha, 'save_trend', save_trend, 'recompute_trend', recompute_trend, 'data_folder', data_folder);
                end
                if ~isempty(Udata_daily) && ~isempty(Vdata_daily)
                    data_name = strcat('u',slct_vector_str);    %/ no need to specify basin if it is not WSV data.
                    [ Udata, Udata_sig,          ~] = compute_trend('data_daily', Udata_daily, 'ins_or_acc', ins_or_acc_vec,  'lon', uv_lon, 'lat', uv_lat, 'data_dates', Udata_dates, 'data_name', data_name,...
                                                                    'year_list', year_list_bc, 'mth', mth, 'trend_test', trend_test, 'trend_alpha', trend_alpha, 'save_trend', save_trend, 'recompute_trend', recompute_trend, 'data_folder', data_folder);
                    
                    data_name = strcat('v',slct_vector_str);    %/ no need to specify basin if it is not WSV data.
                    [ Vdata, Vdata_sig, ins_or_acc] = compute_trend('data_daily', Vdata_daily, 'ins_or_acc', ins_or_acc_vec,  'lon', uv_lon, 'lat', uv_lat, 'data_dates', Vdata_dates, 'data_name', data_name,...
                                                                    'year_list', year_list_bc, 'mth', mth, 'trend_test', trend_test, 'trend_alpha', trend_alpha, 'save_trend', save_trend, 'recompute_trend', recompute_trend, 'data_folder', data_folder);
                
                    %/ Restore the magnitude if either component of the vector is sig.
                    cond_U_sig = (~isnan(Udata_sig));
                    cond_V_sig = (~isnan(Vdata_sig));
                    Udata(~cond_U_sig & ~cond_V_sig) = nan;
                    Vdata(~cond_U_sig & ~cond_V_sig) = nan;
                end
                if ~isempty(U2data_daily) && ~isempty(V2data_daily)
                    data_name = strcat('u',slct_vector2_str);    %/ no need to specify basin if it is not WSV data.
                    [ U2data, U2data_sig,          ~] = compute_trend('data_daily', U2data_daily, 'ins_or_acc', ins_or_acc_vec, 'lon', uv2_lon, 'lat', uv2_lat, 'data_dates', U2data_dates, 'data_name', data_name,...
                                                                    'year_list', year_list_bc, 'mth', mth, 'trend_test', trend_test, 'trend_alpha', trend_alpha, 'save_trend', save_trend, 'recompute_trend', recompute_trend, 'data_folder', data_folder);
                    
                    data_name = strcat('v',slct_vector2_str);    %/ no need to specify basin if it is not WSV data.
                    [ V2data, V2data_sig, ins_or_acc] = compute_trend('data_daily', V2data_daily, 'ins_or_acc', ins_or_acc_vec, 'lon', uv2_lon, 'lat', uv2_lat, 'data_dates', V2data_dates, 'data_name', data_name,...
                                                                    'year_list', year_list_bc, 'mth', mth, 'trend_test', trend_test, 'trend_alpha', trend_alpha, 'save_trend', save_trend, 'recompute_trend', recompute_trend, 'data_folder', data_folder);
                
                    %/ Restore the magnitude if either component of the vector is sig.
                    cond_U2_sig = (~isnan(U2data_sig));
                    cond_V2_sig = (~isnan(V2data_sig));
                    U2data(~cond_U2_sig & ~cond_V2_sig) = nan;
                    V2data(~cond_U2_sig & ~cond_V2_sig) = nan;
                end
                str_trend = strcat('_trend_', ins_or_acc, str_trend_test);
            end
           
            %======= cont setting ======%
            cont_plus = 0;  cont_multiply = 1;   %/ default
            if ~isempty(slct_cont_str{:})
                cont_linewi     = 2;
                cont_labelsize  = 14;
                cont_label_col  = 'k';
                skip_zero_cont  = 1;    %/ NOTE: if drawing sig data, then there is no zeros!!
                if ismember(slct_cont_str, {'Cf_map_frac', 'Pm_frac'}) 
                    cont_lon    = WSV.(slct_cont_str{:}).lon;
                    cont_lat    = WSV.(slct_cont_str{:}).lat;
                    cont_levels = [0.25:0.75:7.25]/10;   %/ in %
        %             cont_colmap = jet(length(cont_levels));
                    cont_colmap = brewermap(length(cont_levels), '*Reds');
        %             cont_colmap    = repmat([0 153 0]./255, length(cont_levels), 1); 
        %             cont_colmap    = repmat([0 0 255]./255, length(cont_levels), 1); 

                elseif ismember(slct_cont_str, {'Z850'})
                    %/ output an 1-deg averaged topo (from 5-min high topo data)
                    topo = interp_topo('lon_grids', cont_lon, 'lat_grids', cont_lat);
                    cond_gt1500 = (topo > 1500);
                    cont_data(cond_gt1500) = nan;
                    if trend_mode
                        cont_levels = [-50:5:50];      %/ in m
                        cont_colmap = create_div_colmap([33 121 255]./255, [255 102 51]./255, cont_levels);
                    else
                        cont_levels = [1400:10:1530];   %/ in m
%                         cont_label_col  = [75 75 75]./255;
                        cont_label_col  = [153 0 102]./255;
%                         cont_colmap = brewermap(length(cont_levels), '*PuOr');
                        cont_colmap = repmat(cont_label_col, length(cont_levels), 1);
%                         cont_colmap = repmat([204 0 102]./255, length(cont_levels), 1);
%                         cont_colmap = repmat([180 40 225]./255, length(cont_levels), 1); 
%                         cont_colmap = jet(length(cont_levels)); 
                    end

                elseif ismember(slct_cont_str, {'Z500'})
                    if trend_mode
                        cont_levels = [-50:10:50];      %/ in m
                    else
                        cont_levels = [5350:50:5850];   %/ in m
%                         cont_colmap = brewermap(length(cont_levels), '*PuOr');
                        cont_colmap = repmat([180 40 225]./255, length(cont_levels), 1);
%                         cont_colmap = repmat([180 40 225]./255, length(cont_levels), 1); 
%                         cont_colmap = jet(length(cont_levels)); 
                    end
                    
                elseif ismember(slct_cont_str, {'P'})
                    if drywet_mode ~= 0
                        cont_levels = [-20:2:20];   %/ in mm/day
                        cont_colmap = create_div_colmap([33 121 255]./255, [255 102 51]./255, cont_levels);
                    else
                        cont_levels = [0:1:22];     %/ in mm/day
                        cont_colmap = jet(length(cont_levels)); 
                    end
                    
                elseif ismember(slct_cont_str, {'E'})
                    if drywet_mode ~= 0
                        cont_levels = [-10:1:10];   %/ in mm/day
                        cont_colmap = create_div_colmap([33 121 255]./255, [255 102 51]./255, cont_levels);
                    else
                        cont_levels = [0:1:22];     %/ in mm/day
                        cont_colmap = jet(length(cont_levels)); 
                    end
                    
                elseif ismember(slct_cont_str, {'Omega500'})
                    if drywet_mode ~= 0
                        cont_levels = [-50:5:50];      %/ in hPa/day
                        cont_colmap = create_div_colmap([33 121 255]./255, [255 102 51]./255, cont_levels);
                    else
                        cont_levels = [-100:10:100];   %/ in hPa/day
                        cont_colmap = jet(length(cont_levels)); 
                    end
                    
                elseif contains(slct_cont_str, {'IVTdiv', 'VMT'})
                    if drywet_mode ~= 0
                        cont_levels = [-flip(2:4:20), 0, 2:4:20]*1e-5;    %/ in kg/m^{2}/s; 2e-5 aligns with hatch_thres for double-checking.
                        cont_colmap = create_div_colmap([33 121 255]./255, [255 102 51]./255, cont_levels);
                    else
                        cont_levels = [-2:.4:2]*1e-4;    %/ in kg/m^{2}/s
                        cont_colmap = create_div_colmap([33 121 255]./255, [255 102 51]./255, cont_levels);
                    end
                elseif ismember(slct_cont_str, {'P_LA_P_diff', 'Pm_E_diff'})  
                    cont_lon    = WSV.(slct_cont_str{:}).lon;
                    cont_lat    = WSV.(slct_cont_str{:}).lat;
                    cont_levels = [-6:1:6];
                    cont_colmap = create_div_colmap([33 121 255]./255, [255 102 51]./255, cont_levels);
                elseif ismember(slct_cont_str, {'P_LA_P_prct_diff', 'Pm_E_prct_diff'})  
                    cont_lon    = WSV.(slct_cont_str{:}).lon;
                    cont_lat    = WSV.(slct_cont_str{:}).lat;
                    cont_levels = [-120:20:120]*2;
                    cont_colmap = create_div_colmap([33 121 255]./255, [255 102 51]./255, cont_levels);
                elseif contains(slct_cont_str, {'VWS'})
                    cont_levels = [-50:5:50];      %/ in m/s
                    cont_colmap = create_div_colmap([175 175 175]./255, [0 0 0]./255, cont_levels);
                elseif contains(slct_cont_str, {'corr_UL'})
                    cont_levels = [-1:.1:1];     
                    cont_colmap = create_div_colmap([175 175 175]./255, [0 0 0]./255, cont_levels);
                else
                    error('cont_levels not pre-set for the input cont data!')
                end
            end

            %======== vector setting =======%
            U_plus = 0; U_multiply = 1;
            V_plus = 0; V_multiply = 1;
            if ~isempty(slct_vector_str{:})
                if contains(slct_vector_str, {'IVT'})
%                     vector_color = [.7 .7 .7];  
                    vector_color = [255 102 177]./255;  
%                     vector_color = [0 57 107]./255; % navy 
%                     vector_edgecolor = 'none';
                    vector_edgecolor = 'k';
                    vec_mag_ref = 400;
                    vec_unit = {' kg/m/s'};
                    if trend_mode
                        U_multiply   = 10;
                        V_multiply   = 10;
                        vec_step_lon = 4;
                        vec_step_lat = 2;
                        vecscale     = 15;      % the smaller value the bigger vector. for winds
                        vecscale2    = 1.5;       % control shaft length.
                        shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                            
                    elseif all_in_one
                        vec_step_lon = 4;
                        vec_step_lat = 2;
                        vecscale     = 1600;     % the smaller value the bigger vector. for winds
                        vecscale2    = 1;       % control shaft length.
                        shaftwidth   = 1;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*3.5; % control length of the head, the larger value the longer

                    elseif slct_region ~= 0
                        vec_step_lon = 5;
                        vec_step_lat = 3;
                        vecscale     = 600;     % the smaller value the bigger vector. for winds
                        vecscale2    = 1;       % control shaft length.
                        shaftwidth   = 2;       % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*7; % control length of the head, the larger value the longer
%                         vec_step_lon = 3;
%                         vec_step_lat = 1;
%                         vecscale     = 200;     % the smaller value the bigger vector. for winds
%                         vecscale2    = 1;       % control shaft length.
%                         shaftwidth   = 1.8;    % control width of the shaft, the larger value the thicker
%                         headlength   = shaftwidth*4; % control length of the head, the larger value the longer

                    else
                        vec_step_lon = 15;
                        vec_step_lat = 6;
                        vecscale     = 500;     % the smaller value the bigger vector. for winds
                        vecscale2    = 1;       % control shaft length.
                        shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer

                    end
                    
                elseif ismember(slct_vector_str, {'850', '500', '300'})
                    if trend_mode
                        U_multiply   = 10;  %/ unit per dec
                        V_multiply   = 10;  %/ unit per dec
                        vec_unit = {' m/s/dec'};
                        if ismember(slct_vector_str, {'850'})
                            vector_color = col_IM;  % [102 204 255]./255; 
                            vec_mag_ref  = 1;
                        elseif ismember(slct_vector_str, {'500'})
                            vector_color = [20 20 20]./255;  % [102 204 255]./255; 
                            vec_mag_ref  = 1;
                        else
                            vector_color = col_Westerly; %[.2 .2 .2]; %dark 
                            vec_mag_ref  = 3;
                        end
                    else
                        vec_unit = {' m/s'};
                        if ismember(slct_vector_str, {'850'})
                            vector_color = col_IM;  % [102 204 255]./255; 
                            vec_mag_ref  = 10;
                        elseif ismember(slct_vector_str, {'500'})
                            vector_color = [20 20 20]./255;  % [102 204 255]./255; 
                            vec_mag_ref  = 20;
                        else
                            vector_color = col_Westerly; %[.2 .2 .2]; %dark 
                            vec_mag_ref  = 50;
                        end
                    end
%                     vector_edgecolor = 'w';
                    vector_edgecolor = 'none';
                    
                    if trend_mode
                        if slct_region == -1
                            vec_step_lon = 6;
                            vec_step_lat = 3;
                            if mth == 0
                                vecscale     = 0.5;     % the smaller value the bigger vector. for winds
                            else
                                vecscale     = 1;     % the smaller value the bigger vector. for winds
                            end
                            vecscale2    = 1;       % control shaft length.
                            shaftwidth   = 2.5;    % control width of the shaft, the larger value the thicker
                            headlength   = shaftwidth*3.5; % control length of the head, the larger value the longer
                        elseif slct_region == 5
                            vec_step_lon = 2;
                            vec_step_lat = 2;
                            vecscale     = 3;      % the smaller value the bigger vector. for winds
                            vecscale2    = 3;       % control shaft length.
                            shaftwidth   = 2.5;    % control width of the shaft, the larger value the thicker
                            headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                        else
                            vec_step_lon = 4;
                            vec_step_lat = 2;
                            vecscale     = 3;      % the smaller value the bigger vector. for winds
                            vecscale2    = 3;       % control shaft length.
                            shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                            headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                        end
                    elseif all_in_one
                        vec_step_lon = 4;
                        vec_step_lat = 2;
                        vecscale     = 160;     % the smaller value the bigger vector. for winds
                        vecscale2    = 1;       % control shaft length.
                        shaftwidth   = 1;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*3.5; % control length of the head, the larger value the longer

                    elseif slct_region == 5
                        vec_step_lon = 4;
                        vec_step_lat = 2;
                        vecscale     = 12;     % the smaller value the bigger vector. for winds
                        vecscale2    = 1;       % control shaft length.
                        shaftwidth   = 2;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*7; % control length of the head, the larger value the longer
                    
                    elseif slct_region == 0
                        vec_step_lon = 8;
                        vec_step_lat = 4;
                        vecscale     = 25;     % the smaller value the bigger vector. for winds
                        vecscale2    = 1;       % control shaft length.
                        shaftwidth   = 2;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    
                    elseif slct_region == -1
                        vec_step_lon = 6;
                        vec_step_lat = 3;
                        vecscale     = 25;     % the smaller value the bigger vector. for winds
                        vecscale2    = 1;       % control shaft length.
                        shaftwidth   = 2.5;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*3.5; % control length of the head, the larger value the longer
                    else
                        vec_step_lon = 7;
                        vec_step_lat = 3;
                        vecscale     = 50;     % the smaller value the bigger vector. for winds
                        vecscale2    = 1;       % control shaft length.
                        shaftwidth   = 2.5;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*3.5; % control length of the head, the larger value the longer
                    end
                    if ismember(slct_vector_str, {'300'}) %/ uv300 is usually stronger than uv850
                        vecscale = vecscale*3;
                    elseif ismember(slct_vector_str, {'500'}) %/ uv500 is usually stronger than uv850
                        vecscale = vecscale*2;
                    end
                end

                if compute_anombyAM || drywet_mode ~= 0
                    vec_mag_ref = vec_mag_ref*0.5;
                end
                vec_lbs = strcat(num2str(vec_mag_ref), vec_unit);
                slct_vector_str2 = strcat('uv', slct_vector_str);
            end
            
            %======== vector2 setting =======%
            U2_plus = 0; U2_multiply = 1;
            V2_plus = 0; V2_multiply = 1;
            if ~isempty(slct_vector2_str{:})
                if contains(slct_vector2_str, {'IVT'})
%                     vector2_color = [.7 .7 .7];  
%                     vector2_color = [35 117 206]./255;  
                    vector2_color = [5 57 112]./255; % navy 
                    vector2_edgecolor = 'none';
%                     vector2_edgecolor = 'w';
                    vec2_mag_ref = 500;
                    vec2_unit = {' kg/m/s'};
                    if trend_mode
                        U2_multiply   = 10;
                        V2_multiply   = 10;
                        vec2_step_lon = 4;
                        vec2_step_lat = 2;
                        vec2scale     = 15;      % the smaller value the bigger vector. for winds
                        vec2scale2    = 1.5;       % control shaft length.
                        shaft2width   = 3.5;    % control width of the shaft, the larger value the thicker
                        head2length   = shaft2width*4; % control length of the head, the larger value the longer
                            
                    elseif all_in_one
                        vec2_step_lon = 4;
                        vec2_step_lat = 2;
                        vec2scale     = 1600;     % the smaller value the bigger vector. for winds
                        vec2scale2    = 1;       % control shaft length.
                        shaft2width   = 1;    % control width of the shaft, the larger value the thicker
                        head2length   = shaft2width*3.5; % control length of the head, the larger value the longer

                    elseif slct_region == 5
                        vec2_step_lon = 3;
                        vec2_step_lat = 2;
                        vec2scale     = 400;     % the smaller value the bigger vector. for winds
                        vec2scale2    = 1;       % control shaft length.
                        shaft2width   = 2;    % control width of the shaft, the larger value the thicker
                        head2length   = shaft2width*7; % control length of the head, the larger value the longer

                    else
                        vec2_step_lon = 15;
                        vec2_step_lat = 6;
                        vec2scale     = 500;     % the smaller value the bigger vector. for winds
                        vec2scale2    = 1;       % control shaft length.
                        shaft2width   = 3.5;    % control width of the shaft, the larger value the thicker
                        head2length   = shaft2width*4; % control length of the head, the larger value the longer

                    end
                    
                elseif ismember(slct_vector2_str, {'850', '300'})
                    if trend_mode
                        U2_multiply   = 10;  %/ unit per dec
                        V2_multiply   = 10;  %/ unit per dec
                        vec2_unit = {' m/s/dec'};
                        if ismember(slct_vector2_str, {'850'})
                            vector2_color = col_IM;  % [102 204 255]./255; 
                            vec2_mag_ref  = 1;
                        else
                            vector2_color = col_Westerly; %[.2 .2 .2]; %dark 
                            vec2_mag_ref  = 3;
                        end
                    else
                        vec2_unit = {' m/s'};
                        if ismember(slct_vector2_str, {'850'})
                            vector2_color = col_IM;  % [102 204 255]./255; 
                            vec2_mag_ref  = 10;
                        else
                            vector2_color = col_Westerly; %[.2 .2 .2]; %dark 
                            vec2_mag_ref  = 30;
                        end
                    end
%                     vector2_edgecolor = 'w';
                    vector2_edgecolor = 'none';
                    if trend_mode
                        U2_multiply   = 10;
                        V2_multiply   = 10;
                        if slct_region == -1
                            vec2_step_lon = 6;
                            vec2_step_lat = 3;
                            if mth == 0
                                vec2scale     = 0.5;     % the smaller value the bigger vector. for winds
                            else
                                vec2scale     = 1;     % the smaller value the bigger vector. for winds
                            end
                            vec2scale2    = 1;       % control shaft length.
                            shaft2width   = 2.5;    % control width of the shaft, the larger value the thicker
                            head2length   = shaftwidth*3.5; % control length of the head, the larger value the longer
                        elseif slct_region == 5
                            vec2_step_lon = 2;
                            vec2_step_lat = 2;
                            vec2scale     = 3;      % the smaller value the bigger vector. for winds
                            vec2scale2    = 3;       % control shaft length.
                            shaft2width   = 2.5;    % control width of the shaft, the larger value the thicker
                            head2length   = shaft2width*4; % control length of the head, the larger value the longer
                        else
                            vec2_step_lon = 4;
                            vec2_step_lat = 2;
                            vec2scale     = 3;      % the smaller value the bigger vector. for winds
                            vec2scale2    = 3;       % control shaft length.
                            shaft2width   = 3.5;    % control width of the shaft, the larger value the thicker
                            head2length   = shaft2width*4; % control length of the head, the larger value the longer
                        end
                    elseif all_in_one
                        vec2_step_lon = 4;
                        vec2_step_lat = 2;
                        vec2scale     = 160;     % the smaller value the bigger vector. for winds
                        vec2scale2    = 1;       % control shaft length.
                        shaft2width   = 1;    % control width of the shaft, the larger value the thicker
                        head2length   = shaft2width*3.5; % control length of the head, the larger value the longer
                        
                    elseif slct_region == 5
                        vec2_step_lon = 4;
                        vec2_step_lat = 2;
                        vec2scale     = 12;     % the smaller value the bigger vector. for winds
                        vec2scale2    = 1;       % control shaft length.
                        shaft2width   = 2;    % control width of the shaft, the larger value the thicker
                        head2length   = shaft2width*7; % control length of the head, the larger value the longer

                    elseif slct_region == 0
                        vec2_step_lon = 8;
                        vec2_step_lat = 4;
                        vec2scale     = 25;     % the smaller value the bigger vector. for winds
                        vec2scale2    = 1;       % control shaft length.
                        shaft2width   = 2;    % control width of the shaft, the larger value the thicker
                        head2length   = shaft2width*4; % control length of the head, the larger value the longer
                    
                    elseif slct_region == -1
                        vec2_step_lon = 6;
                        vec2_step_lat = 3;
                        vec2scale     = 25;     % the smaller value the bigger vector. for winds
                        vec2scale2    = 1;       % control shaft length.
                        shaft2width   = 2.5;    % control width of the shaft, the larger value the thicker
                        head2length   = shaft2width*3.5; % control length of the head, the larger value the longer
                    
                    else
                        vec2_step_lon = 7;
                        vec2_step_lat = 3;
                        vec2scale     = 50;     % the smaller value the bigger vector. for winds
                        vec2scale2    = 1;       % control shaft length.
                        shaft2width   = 2.5;    % control width of the shaft, the larger value the thicker
                        head2length   = shaft2width*3.5; % control length of the head, the larger value the longer
                    end
                    if ismember(slct_vector2_str, {'300'}) %/ uv300 is usually stronger than uv850
                        vec2scale = vec2scale*3;
                    end
                end

%                 if compute_anombyAM || drywet_mode ~= 0
%                     vec2_mag_ref = vec2_mag_ref*0.5;
%                 end
                vec2_lbs = strcat(num2str(vec2_mag_ref), vec2_unit);
                slct_vector2_str2 = strcat('uv', slct_vector2_str);
            end
            
            %======== point setting (to denote sig contf values if hatch is not set) =======%
            if ~isempty(contf_data_sig) && sig_contf_in_hatch == 0
                marker          = 'o'; 
%             markerfacecolor = [255 0 204]./255;  %/ magenta
%             markerfacecolor = [255 255 51]./255;
                if ismember(slct_contf_str, {'CWRT'})
                    markerfacecolor = [204 255 153]./255;
                else
                    markerfacecolor = [0 0 0]./255;
                end
                cond_sig_2D     = ~isnan(contf_data_sig);
                cond_sig_2Dto1D = reshape(cond_sig_2D, [], 1);
                
                if isvector(contf_lon) && isvector(contf_lat)
                    [contf_lon_2D, contf_lat_2D] = meshgrid(contf_lon, contf_lat);
                    contf_lon_2D = contf_lon_2D'; contf_lat_2D = contf_lat_2D';
                else
                    contf_lon_2D = contf_lon;
                    contf_lat_2D = contf_lat;
                end
                contf_lon_2Dto1D = reshape(contf_lon_2D, [], 1);
                contf_lat_2Dto1D = reshape(contf_lat_2D, [], 1);
                
                point_data = [contf_lon_2Dto1D(cond_sig_2Dto1D), contf_lat_2Dto1D(cond_sig_2Dto1D)];
            end
            
            %======= hatch setting (to denote sig contf values if trend_mode = 1) [overwriting] ======%
            if trend_mode && sig_contf_in_hatch
                hatch_data = contf_data_sig;
                hatch_lon = contf_lon;
                hatch_lat = contf_lat; 
                hatch_thres_pve = [];  %/ since the sign is not important, set them [] to avoid bug
                hatch_thres_nve = [];
                color_hatch_pve = [.1 .1 .1]; 
                color_hatch_nve = [.1 .1 .1];
                hatch_mode = 1;
                if contains(slct_contf_str, {'TPR', 'HARv2'})
                    hatch_linewi = 1;
                    hatch_intvl  = 6;
                else
                    hatch_linewi = 2;
                    hatch_intvl  = 6;
                end
            elseif compute_anombyAM && sig_contf_in_hatch   %/ check if contf_data_daily and alpha co-exist
                %/ Then, we show sig data as hatching.
                hatch_data = ttest_sig_fn(contf_data_daily, alpha, 3); 
                hatch_lon = contf_lon;
                hatch_lat = contf_lat; 
                hatch_thres_pve = 0;
                hatch_thres_nve = 0;
                color_hatch_pve = [.4 .4 .4]; 
                color_hatch_nve = [.4 .4 .4];                    
                hatch_mode = 1;
                hatch_linewi = 1.5;
                hatch_intvl = 6;
            elseif ismember(slct_hatch_str, {'Omega500'})
                hatch_thres_pve =  5;   %/ hPa/day
                hatch_thres_nve = -5;   %/ hPa/day
                color_hatch_pve = [255 102 51]./255; 
                color_hatch_nve = [33 121 255]./255;                     %/ not shown anyway.
                hatch_mode = 1;
                hatch_linewi = 2;
                hatch_intvl = 12;

            elseif ismember(slct_hatch_str, {'P'})
                hatch_thres_pve =  1;   %/ mm/day
                hatch_thres_nve = -1;   %/ mm/day
                color_hatch_pve = [255 102 51]./255; 
                color_hatch_nve = [33 121 255]./255;                     %/ not shown anyway.
                hatch_mode = 1;
                hatch_linewi = 2;
                hatch_intvl = 6;

            elseif contains(slct_hatch_str, {'IVTdiv', 'VMT'})
                hatch_thres_pve =  2e-5;   %/ kg m**-2 s**-1
                hatch_thres_nve = -2e-5;   %/ kg m**-2 s**-1
                
%                     color_hatch_pve = [246 224 155]./255;    % yellow
%                     color_hatch_nve = [ 97 130 180]./255;    % blue
                color_hatch_pve = [255 102 51]./255;    % orange
                color_hatch_nve = [33 121 255]./255;    % blue
%                     color_hatch_pve = [.2 .2 .2];      % dark grey
%                     color_hatch_nve = [1 1 1];      % white         
                hatch_mode = 1;
                hatch_linewi = 2;
                hatch_intvl = 12;

            elseif ismember(slct_hatch_str, {'BL_Pm_ratio'})
                hatch_thres_pve = 66;   %/ %
                hatch_thres_nve = [];   %/ %
                color_hatch_pve = [1 1 1]; 
%                     color_hatch_pve = [255 102 51]./255;    
                color_hatch_nve = [];    % blue       
                hatch_mode = 1;

                if all_in_one == 1
                    hatch_linewi = 0.7;
                    hatch_intvl = 8;
                elseif all_in_one == 3
                    hatch_linewi = 0.7;
                    hatch_intvl = 2;
                else
                    hatch_linewi = 2;
                    hatch_intvl = 7;
                end
            end

            %/ Map setting
            title_pos = [];
            cell_all_vars = [slct_contf_str, slct_cont_str, slct_vector_str2, slct_vector2_str2, slct_hatch_str];
            cell_all_vars = cell_all_vars(~cellfun('isempty',cell_all_vars));
            str_all_vars  = strjoin(cell_all_vars, '_');
            
            %/ Select region to plot
            title_fontsize = 12;
            if slct_region == -1   %/ 3D Sphere
                glb_data_mode = 1;
                fontsize      = 22;
                markersize    = 3;
                grid_mode     = -3; 
                coast_wi      = 1.5; linewi = 2;
                map_proj      = 'ortho';
                map_lon_lower = -169;
                map_lon_upper = 190;
                map_lat_lower = -55;
                map_lat_upper = 89;
                map_center    = [84.8421, 20];
                title_pos = [0.12,1.01,1.0,0];
                str_regional  = '_Eurasia_spherical';

            elseif slct_region == 0        %/ global (centerd at 0E) 
                glb_data_mode = 1;
                fontsize      = 18;
                markersize    = 2; 
                grid_mode     = 0;
                coast_wi      = 1.5;
                map_proj      = 'robin';
                map_lon_lower = -169;
                map_lon_upper = 190;
                map_lat_lower = -89;
                map_lat_upper = 89;
                title_pos = [0.12,1.01,1.0,0];
                str_regional  = '_global';
                
            elseif slct_region == 1        %/ global (centerd at dateline) 
                glb_data_mode = 1;
                fontsize      = 18;
                markersize    = 2; 
                grid_mode     = 0;
                coast_wi      = 1.5;
                map_proj      = 'robin';
                map_lon_lower = 0;
                map_lon_upper = 359;
                map_lat_lower = -89;
                map_lat_upper = 89;
                title_pos = [0.12,1.01,1.0,0];
                str_regional  = '_global';

            elseif slct_region == 2        %/ North Pacific Rim
                glb_data_mode = 1;
                fontsize      = 18;
                markersize    = 2; 
                grid_mode     = 0;
                coast_wi      = 1.5;
                map_proj      = 'Miller Cylindrical';
                map_lon_lower = 60;
                map_lon_upper = 300;
                map_lat_lower = -10;
                map_lat_upper = 75;
                title_pos = [0.12,0.8,1.0,0];
                str_regional  = '_global';

            elseif slct_region == 3   
                glb_data_mode = 0;
                fontsize      = 18;
                markersize    = 2;
                grid_mode     = 0;
                coast_wi      = 1.5;
                map_proj      = 'robin';
                map_lon_lower = 0;
                map_lon_upper = 359;
                map_lat_lower = -89;
                map_lat_upper = 89;
                str_regional = '_global2';

            elseif slct_region == 4   %/ TP (zoom-in), no coastline
                linewi        = 3.5; %3;
                glb_data_mode = 0;
                fontsize      = 24;
                markersize    = 10;
                grid_mode     = -1;
                coast_wi      = 2;
%                 map_proj      = 'robin';
%                 map_proj      = 'Miller Cylindrical';
                map_proj      = 'equidistant';
                map_lon_lower = 63;
%                 map_lon_lower = 63.5;
                map_lon_upper = 108;
                map_lat_lower = 24.5;
                map_lat_upper = 44.5; 
%                 map_lat_upper = 43.5; 
                title_pos = [0.12,1.01,1.0,0];
                str_regional  = '_TP';
                coast_col     = 'none'; 
                
            elseif slct_region == 5   %/ Asia + IO (for PakCase)
                glb_data_mode = 0;
                fontsize      = 14;
                title_fontsize = 9;
                markersize    = 4;
                grid_mode     = 3;
                coast_wi      = 1.5;
                map_proj      = 'Miller Cylindrical';
                map_lat_lower = -33;
                map_lat_upper = 48; 
                map_lon_lower = 15;
                map_lon_upper = 130;
                title_pos = [0.12,1.01,1.0,0];
                str_regional  = '_IO';
                coast_col     = 'k';

            elseif slct_region == 5.5   %/ Asia + IO (for PakCase, Zoom-in)
                glb_data_mode = 0;
                fontsize      = 14;
                title_fontsize = 9;
                markersize    = 4;
                grid_mode     = 3;
                coast_wi      = 1.5;
                map_proj      = 'Miller Cylindrical';
                map_lat_lower = 18;
                map_lat_upper = 35; 
                map_lon_lower = 60;
                map_lon_upper = 78;
                title_pos = [0.12,1.01,1.0,0];
                str_regional  = '_IO_ZoomIn';
                coast_col     = 'k';

            elseif slct_region == 6   %/ Australia 
                glb_data_mode = 0;
                fontsize       = 18;
                title_fontsize = 12;
                markersize    = 4;
                grid_mode     = 3;
                coast_wi      = 1.5;
                map_proj      = 'Miller Cylindrical';
                map_lat_lower = -60;
                map_lat_upper = 20; 
                map_lon_lower = 110;
                map_lon_upper = 260;
                title_pos = [0.12,0.9,1.0,0];
                str_regional  = '_Aus';
                coast_col     = 'k';

            elseif slct_region == 6.5   %/ Australia (Zoom-in)
                glb_data_mode = 0;
                fontsize       = 18;
                title_fontsize = 12;
                markersize    = 4;
                grid_mode     = 3;
                coast_wi      = 1.5;
                map_proj      = 'Miller Cylindrical';
                map_lat_lower = -40;
                map_lat_upper = -15; 
                map_lon_lower = 135;
                map_lon_upper = 172;
                title_pos = [0.12,0.99,1.0,0];
                str_regional  = '_Aus_ZoomIn';
                coast_col     = 'k';

            elseif slct_region == 7   %/ Scotland 
                glb_data_mode = 0;  
                fontsize       = 18;
                title_fontsize = 12;
                markersize    = 4;
                grid_mode     = 1;
                coast_wi      = 1.5;
                map_proj      = 'lambert'; 
                % map_proj      = 'Miller Cylindrical';
                map_lat_lower = 0;
                map_lat_upper = 80; 
                map_lon_lower = -100;
                map_lon_upper = 30;
                title_pos = [0.12,0.99,1.0,0];
                str_regional  = '_Scot';
                coast_col     = 'k';

            elseif slct_region == 7.5   %/ Scotland (Zoom-In)
                glb_data_mode = 0;  
                fontsize       = 18;
                title_fontsize = 12;
                markersize    = 4;
                grid_mode     = 1;
                coast_wi      = 1.5;
                map_proj      = 'lambert'; 
                % map_proj      = 'Miller Cylindrical';
                map_lat_lower = 48;
                map_lat_upper = 63; 
                map_lon_lower = -20;
                map_lon_upper = 10;
                title_pos = [0.12,0.99,1.0,0];
                str_regional  = '_Scot_Zoomin';
                coast_col     = 'k';

            else
                error('slct_region not set!');
            end
            
            %/ Markersize should also consider grid res
            res = diff(contf_lon(1:2));
            markersize = markersize * 2.6^(res-1);
            create_fig    = 1;  
            FigName = strcat(str_all_vars, str_land_or_ocean, str_marker, str_trend, str_hatch, '_', str_dates, str_basin, str_drywet, str_alpha, str_regional, str_TP_only);
            titlename = strrep(strcat(str_basin, {' '}, str_all_vars, str_trend, str_drywet, str_scheme, '_', str_dates), '_', ' ');                                

            if ~isempty(SSIM)
                str_diff_metrics = sprintf(' (SSIM=%.2f r=%.2f RMSE=%.2f)', SSIM, r, RMSE); %/ This str will be used.
                titlename        = strcat(titlename, str_diff_metrics);
            end
            savepath = [];

            %/ Map
            if plot_map 
                plot_contfmap('contf_data', (contf_data + contf_plus)*contf_multiply, 'contf_lon', contf_lon, 'contf_lat', contf_lat, 'contf_levels', contf_levels,...
                              'contf_unit', contf_unit, 'colmap', colmap, 'cbar_interval', cbar_interval, 'pcolor_mode', pcolor_mode,...
                              'cont_data',  (cont_data + cont_plus)*cont_multiply,  'cont_lon', cont_lon, 'cont_lat', cont_lat, 'cont_levels', cont_levels, 'cont_colmap', cont_colmap, 'cont_linewi', cont_linewi, 'cont_labelsize', cont_labelsize, 'cont_label_col', cont_label_col, 'skip_zero_cont', skip_zero_cont,...
                              'point_data', point_data, 'marker', marker, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'linewi', linewi, 'color', color,...
                              'hatch_data', hatch_data, 'hatch_lon', hatch_lon, 'hatch_lat', hatch_lat, 'hatch_thres_pve', hatch_thres_pve, 'hatch_thres_nve', hatch_thres_nve,...
                              'color_hatch_pve', color_hatch_pve, 'color_hatch_nve', color_hatch_nve, 'hatch_mode', hatch_mode, 'hatch_linewi', hatch_linewi, 'hatch_intvl', hatch_intvl,...
                              'bndry_data', bndry_data, 'bndry_patch_mode', bndry_patch_mode, 'titlename', titlename, 'title_pos', title_pos, 'savepath', savepath, 'fig_fmt', fig_fmt,...
                              'Udata', (Udata + U_plus)*U_multiply, 'Vdata', (Vdata + V_plus)*V_multiply, 'uv_lon', uv_lon, 'uv_lat', uv_lat, 'vec_step_lon', vec_step_lon, 'vec_step_lat', vec_step_lat,...
                              'vector_color', vector_color, 'vector_edgecolor', vector_edgecolor, 'vecscale', vecscale, 'vecscale2', vecscale2, 'shaftwidth', shaftwidth, 'headlength', headlength, 'vec_lbs', vec_lbs, 'vec_mag_ref', vec_mag_ref, 'vec_ref_fontsize', vec_ref_fontsize,...
                              'U2data', (U2data + U2_plus)*U2_multiply, 'V2data', (V2data + V2_plus)*V2_multiply, 'uv2_lon', uv2_lon, 'uv2_lat', uv2_lat, 'vec2_step_lon', vec2_step_lon, 'vec2_step_lat', vec2_step_lat,...
                              'vector2_color', vector2_color, 'vector2_edgecolor', vector2_edgecolor, 'vec2scale', vec2scale, 'vec2scale2', vec2scale2, 'shaft2width', shaft2width, 'head2length', head2length, 'vec2_lbs', vec2_lbs, 'vec2_mag_ref', vec2_mag_ref, 'vec2_ref_fontsize', vec2_ref_fontsize,...
                              'glb_data_mode', glb_data_mode, 'glb_plateau_mode', glb_plateau_mode, 'plateau_hgt', plateau_hgt, 'plateau_col', plateau_col,...
                              'map_proj', map_proj, 'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper, 'map_center', map_center,'coast_col', coast_col, 'coast_wi', coast_wi, 'backcolor', backcolor,...
                              'fontsize', fontsize,  'title_fontsize', title_fontsize, 'create_fig', create_fig, 'grid_mode', grid_mode, 'cbar_mode', cbar_mode, 'cbar_location', cbar_location, 'cbar_position', cbar_position, 'cbar_YTick', cbar_YTick, 'cbar_YTickLabel', cbar_YTickLabel, 'cbar_fontsize', cbar_fontsize,...
                              'text_data', text_data, 'text_backgroundcolor', text_backgroundcolor, 'text_edgecolor', text_edgecolor, 'text_fontsize', text_fontsize,...
                              'trans_bg', trans_bg, 'coast_patch_col', coast_patch_col, 'patch_color', patch_color, 'patch_alpha', patch_alpha, 'draw_country', draw_country, 'draw_cbar_only', draw_cbar_only, 'draw_refvec_only', draw_refvec_only, 'curvature', curvature)
                fprintf('done \n')
            end

            %/ Compute area-weighted mean (AWM) or sum of contf_data
            if compute_contf_AWM
                [contf_AWM, ~, ~, area_unit, str_extent] = compute_area_wgted_meansum('data', (contf_data + contf_plus)*contf_multiply, 'lon', contf_lon, 'lat', contf_lat,...
                                                                                      'lon_extent', AWM_lon_extent, 'lat_extent', AWM_lat_extent, 'mean_or_sum', mean_or_sum);
                
                str_AWM = sprintf('AWM%s: %.5G %s %s', str_extent, contf_AWM, contf_unit, area_unit); 

                AWM_pos = title_pos; 
                AWM_pos(2) = AWM_pos(2) + 0.04;
                annotation('textbox', 'String', str_AWM, 'Color', 'k', 'FontSize', title_fontsize,...
                           'Units', 'normalized', 'EdgeColor', 'none',  'Position', AWM_pos)
            end
            
            %/ Save figue here to include AWM
            if savefig  
                final_FigName = char(fullfile(plotting_folder, strrep(FigName, ' ', '_')));
                export_fig(final_FigName,['-r',num2str(png_dpi)],'-png', '-opengl', '-nocrop', '-transparent');
            end

            %/ Some statistics
            if ismember(slct_contf_str, {'optimal_trajtime'}) && trend_mode == 0
                %/ unweighted
                tracking_time = (contf_data + contf_plus)*contf_multiply;
                tracking_time_gavg = mean(mean(tracking_time, 'omitnan'), 'omitnan'); %/ not considering grid area weighting (for now)
                fprintf('*** Global mean (unweighted)     tracking time = %.2f days ***\n', tracking_time_gavg);

                if ismember(ldirect, {'bwd'})
                    %/ areal-weighted
                    wgt = cosd(WSV_lat);                                %/ Create area weights first because the grid area decreases with latitude
                    wgt = wgt/sum(wgt);                             %/ Normalize wgt so that its sum is 1 
                    wgt = reshape(wgt, 1, [], 1);                   %/ Reshape to perform 3D element-wise mult.
                    tracking_time_gavg_wgt = squeeze(sum(mean(wgt.*tracking_time, 1, 'omitnan'), 2, 'omitnan'));
                    fprintf('*** Global mean (areal-weighted) trakcing time = %.2f days ***\n', tracking_time_gavg_wgt);
                end
            end
            
            %/ Save the difference map for later use (e.g., design RHc map)
            if save_diffmap && contains(slct_contf_str, '_diff')
                diffmap = contf_data;
                matfilename_diffmap = strcat(data_folder, sprintf('scheme%d_frombasin%d_%s%s%s.mat',...
                                             RHc_dqc_scheme, from_basin, slct_contf_str{:}, str_TP_only, str_dates));
                save(matfilename_diffmap, 'diffmap', '-v7.3')
                fprintf('*** diffmap has been saved to %s ***\n', matfilename_diffmap)
            end
            contf_data_daily = [];
        end
    end
end

%% [Step 5] Plot trajectories (fwd/bwd)

savefig       = 0;              
savemat       = 1;
recompute     = 0;
cbar_location = 'eastoutside'; 
fig_fmt       = 'png';
trans_bg      = 1;
NumWorkers    = [];
year_list_bc  = year_list;
WSV_name      = 'traj';

%==========================================================================
plot_traj           = 1;
traj_var            = 'q';               %/ Show 'Z', 'q', or 'T' as the colors of the trajs 
show_max_ntraj      = 200;               %/ Select a max (random) number of trajs to plot
show_max_trajtime   = 24*10;             %/ Select the backtrack time (in hours) of the trajs
% show_max_trajtime   = 24*maxtraj_day;  %/ Select the backtrack time (in hours) of the trajs

if contains(expmnt, 'PakCase') 
    st_month = 8; st_day = 10; ed_month = 8; ed_day = 24;   %/ Pakistan case (whole period mean)
    slct_region = 5; 

elseif contains(expmnt, 'AusCase')
    % st_month = 2; st_day = 22; ed_month = 2; ed_day = 28;   %/ Australian case (whole period mean)
    st_month = 2; st_day = 25; ed_month = 2; ed_day = 25;   %/ Australian case (the most extreme flood day)
    slct_region = 6;  
    
elseif contains(expmnt, 'ScotCase')
    if isequal(expmnt, 'domfill_EA_0.5deg_MPI_202208_ScotCase_test')
        st_month = 8; st_day = 10; ed_month = 8;  ed_day = 12;   %/ testing
    else
        st_month = 10; st_day = 6; ed_month = 10; ed_day = 8;    
    end
    slct_region = 7;
else
    %/ [By default] Use the tracking period 
    if numel(num2str(stdate)) == 8
        year_list = floor(stdate./1e4):floor(eddate./1e4); 
        
        st_month = floor(mod(stdate, 1e4)./1e2); st_day = mod(stdate, 1e2); 
        ed_month = floor(mod(eddate, 1e4)./1e2); ed_day = mod(eddate, 1e2); 
    end
    if contains(expmnt, 'ALA')
        slct_region = 2;
    else
        slct_region = 0;
    end
end

%==========================================================================
% basin_list = 1;                            %/ Choose the source region where traj is tracked from. See basin_catalog.
basin_list = 1:length(basin_catalog);  %/ Minus 1 is to avoid IPCC-PAK-rest (out of memory)
%==========================================================================

for top = basin_list
    basin_name = basin_catalog(top).name{:};
    bndry_data = basin_catalog(top).bndry;   %/ Outline the target region
    color      = [1, 0, 0];
    str_basin  = sprintf('_%s', basin_name);
    
    %/ Dates
    if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
        str_dates = sprintf('%d%02d%02d-%d%02d%02d', year_list_bc(1), st_month, st_day, year_list_bc(end), ed_month, ed_day);
    else
        if mth == 0        str_dates = strcat(str_years_bc);
        else               str_dates = strcat(str_years_bc,'_', str_mth{mth});    end
    end
    
    WSV_data_filename = string(strcat(data_folder, 'WSV_', WSV_name, str_basin,...
                           strrep(str_expmntinfo, ' ', '_'), str_dates, '.mat'));
    
    if isfile(WSV_data_filename) && recompute == 0
        fprintf('!!! The queried %s is found. Loading... !!! \n', WSV_name);
        traj_catalog = par_load(WSV_data_filename, 'traj_catalog');
    else
        fprintf('*** %s Not Found ***\n', WSV_name)
        fprintf('*** Now implementing retrieve_WSV... ***\n')
        [~, ~, ~, traj_catalog] = retrieve_WSV('WSV_name', WSV_name,  'WSV_dir',  WSV_dir, 'year_list', year_list_bc, 'mth', mth,...
                                                'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day, 'str_RHc_dqc', str_RHc_dqc,...
                                                'ldirect', ldirect, 'from_basin', from_basin, 'maxtraj_day', maxtraj_day, 'str_optimal', str_optimal,...
                                                'dt_slct', dt_slct, 'str_traj_rm_jump', str_traj_rm_jump, 'str_BLH_factor', str_BLH_factor, 'str_remark', str_remark,...
                                                'str_src_z', str_src_z, 'str_domain', str_domain, 'str_domain_trajfile', str_domain_trajfile, 'str_sharpcut', str_sharpcut, 'forcing', forcing,...
                                                'NumWorkers', NumWorkers, 'basin_name', basin_name, 'basin_catalog', basin_catalog);
        %/ Save only the whole year or seasonal WSV.
        if savemat
            fprintf('*** Saving traj_catalog: %s *** \n', WSV_data_filename{:})
            save(WSV_data_filename, 'traj_catalog', '-v7.3');
        end
    end
    
    %===== Plot trajectories =====%
    if plot_traj
        fprintf('*** Plotting trajs... ***\n');
        close all
    
        %/ Randomly select XXX trajs
        rng("default")
        ntraj         = size(traj_catalog.traj,1);
        ind_slct_traj = randperm(ntraj, show_max_ntraj); %/ create a list of non-repeating indices with number of slct_maxntraj ranging from 1 to ntraj
        
        traj_data = permute(traj_catalog.traj(ind_slct_traj,:,:), [1,3,2]);  %/ permute to [ntraj, TrajTime, XYZ]
        if isequal(traj_var, 'Z')
            traj_levels = -500:500:7000;  %/ Altitude (m) 
            traj_unit   = 'm'; 
            traj_colmap = viridis(length(traj_levels)-1);
    
        elseif isequal(traj_var, 'q')  
            traj_data(:,:,3) = permute(traj_catalog.(traj_var)(ind_slct_traj,:,:), [1,3,2]); %/ Replace Z with q
            traj_levels = 0:2:22;  %/ Specific humidty (g/kg)
            traj_unit   = 'g/kg'; 
            % traj_colmap = nclCM('precip4_11lev', length(traj_levels)-1);
            traj_colmap = nclCM('cmocean_ice', length(traj_levels)-1);
    
        elseif isequal(traj_var, 'T')
            traj_data(:,:,3) = permute(traj_catalog.(traj_var)(ind_slct_traj,:,:), [1,3,2]); %/ Replace Z with T
            traj_levels = 265:5:315;  %/ Temperature (K)
            traj_unit   = 'K'; 
            traj_colmap = nclCM('cmocean_thermal', length(traj_levels)-1); %/ MPL_RdYlBu
        else
            error('Invalid input of ''traj_var''!');
        end
    
        if mod(length(traj_levels), 2) == 1  %/ odd
            cbar_interval = 2;
        else
            cbar_interval = 3;
        end
    
        %/ Specify the max length of trajs
        traj_time             = 1:size(traj_data, 2); 
        ind_slct_max_trajtime = 1:floor(show_max_trajtime/dt);
    
        traj_time = traj_time(ind_slct_max_trajtime);
        traj_data = traj_data(:,ind_slct_max_trajtime,:);
    
        %/ Map setting
        fontsize        = 14;
        cbar_fontsize   = fontsize;
        traj_linewi     = 1;
        coast_wi        = 2;
        coast_col       = 'none';
        coast_patch_col = [240 239 221]./255;
        backcolor       = 'none';
        % backcolor       = [151 182 226]./255;

        %/ Make sure the region we plot covers all trajectories to avoid edge issue
        % linewi          = 1.5;
        % markersize      = 1;
        % sz = size(traj_data);
        % margin_lon = 5;
        % margin_lat = 3;
        % map_lon_lower = min(reshape(traj_data(:,:,1), sz(1)*sz(2), []))-margin_lon; 
        % map_lon_upper = max(reshape(traj_data(:,:,1), sz(1)*sz(2), []))+margin_lon; 
        % map_lat_lower = min(reshape(traj_data(:,:,2), sz(1)*sz(2), []))-margin_lat; 
        % map_lat_upper = max(reshape(traj_data(:,:,2), sz(1)*sz(2), []))+margin_lat; 

        if contains(expmnt, 'PakCase') 
            glb_data_mode = 0;
            grid_mode     = 3;
            linewi        = 2.5;
            markersize    = 1;
            map_lat_lower = -33;
            map_lat_upper = 48; 
            map_lon_lower = 15;
            map_lon_upper = 130;

        elseif contains(expmnt, 'AusCase') 
            glb_data_mode = 0;
            grid_mode     = 3;
            linewi        = 2.5;
            markersize    = 1;
            % map_lat_lower = 10;
            % map_lat_upper = 50; 
            % map_lon_lower = 190;
            % map_lon_upper = 300;
            map_lat_lower = -60;
            map_lat_upper = 20; 
            map_lon_lower = 110;
            map_lon_upper = 260;

        elseif contains(expmnt, 'ScotCase') 
            glb_data_mode = 0;
            grid_mode     = 3;
            linewi        = 2.5;
            markersize    = 1;
            map_lat_lower = 0;
            map_lat_upper = 80; 
            map_lon_lower = -100;
            map_lon_upper = 30;
        end

        titlename = sprintf('traj_%s_%s_%d_%dh_%s', traj_var, basin_name, show_max_ntraj, show_max_trajtime, str_dates);
        savepath = [];
        if savefig
            savepath = strcat(plotting_folder, titlename);
        end
        
        draw_cbar_only = 0;
        plot_contfmap('map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper,...
                      'bndry_data', bndry_data, 'color', color, 'bndry_patch_mode', bndry_patch_mode, 'coast_wi', coast_wi, 'fontsize', fontsize, 'create_fig', 1, 'grid_mode', grid_mode,...
                      'traj_data', traj_data, 'traj_time', traj_time, 'traj_levels', traj_levels, 'traj_colmap', traj_colmap, 'traj_linewi', traj_linewi, 'traj_unit', traj_unit,...
                      'glb_data_mode', glb_data_mode, 'markersize', markersize, 'linewi', linewi, 'coast_col', coast_col, 'coast_patch_col', coast_patch_col, 'backcolor', backcolor,...
                      'titlename', titlename, 'savepath', savepath, 'cbar_location', cbar_location, 'draw_cbar_only', draw_cbar_only, 'cbar_fontsize', cbar_fontsize, 'cbar_interval', cbar_interval, 'fig_fmt', fig_fmt, 'trans_bg', trans_bg);
    
    end
end

%% [Step 6] SR network (Pie Chart / Bar plots / Line plots / Table) (Figs. 4, 5)
close all; plot_SR_pie = 0; plot_SR_chord = 0; plot_SR_bar = 0; plot_lines = 0; plot_legend_only = 0; plot_table = 0; SR = []; anom_base_period = []; slct_xdata = [];
select_field = 'daily'; x_intvl = []; model_list = []; exp_list = []; ens_list = []; DivByAreaOf = [];  %/default
map_xticks = []; map_xticklabels = []; show_y1_as_bar   = 0;  N_running_avg = []; SR = []; yshift = []; line_data_xlim = []; show_diagonal = 0;
st_month = []; st_day = []; ed_month = []; ed_day = []; bar_pve_col = []; bar_nve_col = []; regime_ver = []; NumWorkers = []; derive_SR = 1; nc_remark = [];

savefig               = 1;
savemat               = 1;
savenc                = 1;     %/ whether to save SR time series into nc files
recompute_SR          = 0;     %<-- mind this!    
thres_agglomerate     = 3;     %/ group those contributions < x% into 'Others'
trend_mode            = 1;
trend_alpha           = 0.1; 
trend_test            = 'MK';   %/ [], 'MK'
asterisk              = {'(*)', '(**)', '(***)'}; 
alpha_level           = [  0.1,  0.05,   0.01]; 
show_regress_line     = 1;
box_mode              = 'on';
plot_or_not           = 1; 
fig_fmt = 'png'; png_dpi = 200;   %/ Output pdf may cause PPT to crash!

%==========================================================================
Njob = 3; job_id = 3;   %/ r147, r146, r145, 
%==========================================================================

% plot_SR_bar = 1;   %/ bar plots in Fig. 2
% plot_SR_pie = 1; 
plot_lines = 1; plot_legend_only = 0; show_legend = 1; lgd_position = 'northwest'; %'eastoutside'; southoutside
% plot_lines = 1; plot_legend_only = 1; show_legend = 1; lgd_position = 'eastoutside'; %'eastoutside'; southoutside
% plot_table = 1;  %/ get the numbers for Fig. 1 

if plot_SR_pie || plot_SR_chord || plot_SR_bar
    %======================================================================
    line_setname = 'VarSet_Cf_map_frac_real_Pakistan'; 
    %======================================================================
    if isequal(line_setname, 'VarSet_Pm_frac_real_Pakistan')
        recompute_SR          = 1;
        slct_data_list        = {'Pm_frac_real'};   %/ will compute ratio contribution
        area_mean_or_sum_list = {'sum'};
        ins_or_acc_list       = {'acc'};
        grouping_mode_list    = {'IPCC-PAK'};    
        DivByAreaOf           = 'Pakistan_box';   
        st_month = 8; st_day = 10; ed_month = 8; ed_day = 24; season = []; mth_list = 0; %/ Pakistan case (whole period mean)
        basin_list            = 1; %:length(basin_catalog); 
    end    
    if isequal(line_setname, 'VarSet_Cf_map_frac_real_Pakistan')
        recompute_SR          = 0;
        slct_data_list        = {'Cf_map_frac_real'};   %/ will compute ratio contribution
        area_mean_or_sum_list = {'sum'};
        ins_or_acc_list       = {'acc'};
        grouping_mode_list    = {'Pakistan_box'};   
        DivByAreaOf           = 'Pakistan_box';   
        st_month = 8; st_day = 10; ed_month = 8; ed_day = 24; season = []; mth_list = 0; %/ Pakistan case (whole period mean)
        basin_list            = 1:length(basin_catalog); 
    end    
    if isequal(line_setname, 'VarSet_land-ocean')   %/ Compute land, ocean contribution and traceability
        slct_data_list        = {'Pm_frac_real'};   %/ will compute ratio contribution
        area_mean_or_sum_list = {'sum'};
        ins_or_acc_list       = {'acc'};
        grouping_mode_list    = {'land-ocean'};    %/ 'basin-wise', 'domain-wise', 'land-ocean', 'local-nonlocal'
        mth_list              = [0,13:16];   
        basin_list            = 0; %:length(basin_catalog); 
    end
end
if plot_lines 
    fontsize         = 15.5;
    linewidth        = 2.5;
    markersize       = 10;
    markeredgecolor  = 'none';
    yyaxis1_col      = 'k'; %/ default
    yyaxis2_col      = 'k'; %/ default
    if show_legend  fig_width  = 900; fig_height = 600; 
    else            fig_width  = 650; fig_height = 600;     end
    %======================================================================
    line_setname = 'VarSet_ERA5_P_P_LA_Pm_attri';    %/ VarSet_ERA5_P_P_LA_Pm_attri   VarSet_ERA5_P_P_LA
    %======================================================================
    if isequal(line_setname, 'VarSet_ERA5_P_P_LA_Pm_attri')    %/ P & total moisture attribution
        trend_mode            = 0;      %/ do not compute trend for daily time series
        slct_data_list        = {'ERA5_P_05',   'P_LA',     'Pm'};   %/ Use the same resolution (0.5deg pr) to compute moisture attribution ratios
        area_mean_or_sum_list = {      'sum',    'sum',    'sum'};  
        ins_or_acc_list       = {      'acc',    'acc',    'acc'};  
        grouping_mode_list    = {         '',       '', 'global'}; 
        slct_src_list         = {         '',       '', 'global'};  
        unit_conv_list        = ones(1, length(slct_src_list));
        DivByAreaOf           = 'basin_itself';     
        marker                = [repmat({'^'},    1, length(slct_src_list))];
        year_range_bc         = [repmat([year_list(1), year_list(end)], length(slct_data_list), 1);];
        shared_year_list      = year_range_bc;
        select_field_list     = repmat({'daily'}, 1, length(slct_src_list));
        use_yyaxis            = 1;
        show_y1_as_bar        = 0;              %/ display y1 as bar
        bar_pve_col           = [102 204 255]./255;
        show_regress_line     = 0;
        anom_mode             = 0;   %/ whether to show absolute value or anomaly (i.e. abs - clim) 
        anom_base_period      = []; 
        
        season = []; mth_list = 0;
        if contains(expmnt, 'PakCase') 
            st_month = 8; st_day = 10; ed_month = 8; ed_day = 24;  %/ Pakistan case (whole period mean)
            line_data_ylim        = [0, 50];
            fontsize = 14; fig_width  = 1500; fig_height = 600;

        elseif contains(expmnt, 'AusCase')
            st_month = 2; st_day = 22; ed_month = 2; ed_day = 28;   %/ Australian case (whole period mean)
            line_data_ylim        = [0, 80];
            fontsize = 14; fig_width  = 1500; fig_height = 600;

        elseif contains(expmnt, 'ScotCase')
            if isequal(expmnt, 'domfill_EA_0.5deg_MPI_202208_ScotCase_test')
                line_data_ylim        = [0, 6];
                st_month = 8; st_day = 10; ed_month = 8;  ed_day = 12;    
            else
                line_data_ylim        = [0, 30];
                st_month = 10; st_day = 6; ed_month = 10; ed_day = 8;    
            end
            
            fontsize = 14; fig_width  = 1500; fig_height = 600;
            % fontsize = 24; fig_width  = 700; fig_height = 600;
        else
            error('code not set!');
        end
        basin_list  = 1;  
        yyaxis1_col = 'k';
        yyaxis2_col = 'k';
        line_colmap = [  0   0   0;
                       102 204 255;
                         0 255 151;]./255;
    end  
    if isequal(line_setname, 'VarSet_ERA5_P_P_LA')        %/ P & P_LA
        trend_mode            = 0;      %/ do not compute trend for daily time series
        slct_data_list        = {'ERA5_P_05',    'P_LA'};   %/ Use the same resolution (0.5deg pr) to compute moisture attribution ratios
        area_mean_or_sum_list = {      'sum',     'sum'};  
        ins_or_acc_list       = {      'acc',     'acc'};  
        grouping_mode_list    = {         '',  'global'}; 
        slct_src_list         = {         '',  'global'};  
        unit_conv_list        = ones(1, length(slct_src_list));
        DivByAreaOf           = 'basin_itself';     
        marker                = [repmat({'^'},    1, length(slct_src_list))];
        year_range_bc         = [repmat([year_list(1), year_list(end)], length(slct_data_list), 1);];
        shared_year_list      = year_range_bc;
        select_field_list     = repmat({'daily'}, 1, length(slct_src_list));
        use_yyaxis            = 1;
        show_y1_as_bar        = 0;              %/ display y1 as bar
        bar_pve_col           = [102 204 255]./255;
        show_regress_line     = 0;
        anom_mode             = 0;   %/ whether to show absolute value or anomaly (i.e. abs - clim) 
        anom_base_period      = []; 
        
        season = []; mth_list = 0;
        if contains(expmnt, 'PakCase') 
            st_month = 8; st_day = 10; ed_month = 8; ed_day = 24;  %/ Pakistan case (whole period mean)
            lline_data_ylim        = [0, 50];
            fontsize = 14; fig_width  = 1500; fig_height = 600;

        elseif contains(expmnt, 'AusCase')
            st_month = 2; st_day = 22; ed_month = 2; ed_day = 28;   %/ Australian case (whole period mean)
            line_data_ylim        = [0, 80];
            fontsize = 24; fig_width  = 700; fig_height = 600;

        elseif contains(expmnt, 'ScotCase')
            if isequal(expmnt, 'domfill_EA_0.5deg_MPI_202208_ScotCase_test')
                st_month = 8; st_day = 10; ed_month = 8;  ed_day = 12;    
            else
                st_month = 10; st_day = 6; ed_month = 10; ed_day = 8;    
            end
            line_data_ylim        = [0, 30];
            fontsize = 24; fig_width  = 700; fig_height = 600;
        else
            error('code not set!');
        end
        basin_list  = 1;  
        yyaxis1_col = 'k';
        yyaxis2_col = 'k';
        line_colmap = [  0   0   0;
                       102 204 255]./255;
    end  
    if isequal(line_setname, 'VarSet_Pm_frac_real_Pakistan_box')
        trend_mode            = 0;      %/ do not compute trend for daily time series
        slct_src_list         = {'global'};
        slct_data_list        = [repmat({'Pm_frac_real'}, 1, length(slct_src_list))];   %/ compute total water budget
        unit_conv_list        = ones(1, length(slct_data_list));
        area_mean_or_sum_list = [repmat({'sum'}, 1, length(slct_src_list))];
        ins_or_acc_list       = [repmat({'acc'}, 1, length(slct_src_list))];
        grouping_mode_list    = slct_src_list;
        DivByAreaOf           = 'basin_itself';   
        marker                = [repmat({'o'}, 1, length(slct_src_list))];
        year_range_bc         = [repmat([year_list(1), year_list(end)], length(slct_data_list), 1);];
        shared_year_list      = year_range_bc;
        select_field_list     = repmat({'daily'}, 1, length(slct_data_list));
        use_yyaxis            = 1;
        anom_mode             = 0;   %/ whether to show absolute value or anomaly (i.e. abs - clim) 
        anom_base_period      = []; 
        st_month = 8; st_day = 10; ed_month = 8; ed_day = 24; season = []; mth_list = 0; %/ Pakistan case (whole period mean)
        basin_list            = 1; 
        yyaxis1_col = 'k';
        yyaxis2_col = 'k';
        line_data_ylim = [0, 200];
        % line_colmap           = nclCM('Cat12', length(slct_data_list));
        % line_colmap           = nclCM('NCV_banded', length(slct_data_list));
        line_colmap = nclCM('amwg_blueyellowred', length(slct_data_list)); line_colmap(6,:) = [0 204 153]./255;

        fontsize = 14; fig_width  = 1500; fig_height = 600;
    end 
    if isequal(line_setname, 'VarSet_Pm_frac_real_Pakistan_IPCC')
        trend_mode            = 0;      %/ do not compute trend for daily time series
        slct_src_list         = labels_IPCC_PAK;
        slct_data_list        = [repmat({'Pm_frac_real'}, 1, length(slct_src_list))];   %/ compute total water budget
        unit_conv_list        = ones(1, length(slct_data_list));
        area_mean_or_sum_list = [repmat({'sum'}, 1, length(slct_src_list))];
        ins_or_acc_list       = [repmat({'acc'}, 1, length(slct_src_list))];
        grouping_mode_list    = [repmat({'IPCC-PAK'}, 1, length(slct_src_list))];
        DivByAreaOf           = 'basin_itself';   
        marker                = [repmat({'o'}, 1, length(slct_src_list))];
        year_range_bc         = [repmat([year_list(1), year_list(end)], length(slct_data_list), 1);];
        shared_year_list      = year_range_bc;
        select_field_list     = repmat({'daily'}, 1, length(slct_data_list));
        use_yyaxis            = 1;
        anom_mode             = 0;   %/ whether to show absolute value or anomaly (i.e. abs - clim) 
        anom_base_period      = []; 
        st_month = 8; st_day = 10; ed_month = 8; ed_day = 24; season = []; mth_list = 0; %/ Pakistan case (whole period mean)
        basin_list            = 1; 
        yyaxis1_col = 'k';
        yyaxis2_col = 'k';
        line_data_ylim = [0, 200];
        % line_colmap           = nclCM('Cat12', length(slct_data_list));
        % line_colmap           = nclCM('NCV_banded', length(slct_data_list));
        line_colmap = nclCM('amwg_blueyellowred', length(slct_data_list)); line_colmap(6,:) = [0 204 153]./255;

        fontsize = 14; fig_width  = 1500; fig_height = 600;
    end 
    if isequal(line_setname, 'VarSet_Cf_map_frac_real_Pakistan')    %/ P time series and total moisture attribution
        recompute_SR          = 1;
        savenc                = 0;      %/ whether to save SR time series into nc files
        trend_mode            = 0;      %/ do not compute trend for daily time series
        slct_src_list         = {'Pakistan_box'};
        slct_data_list        = {'Cf_map_frac_real'};   %/ compute total water budget
        unit_conv_list        = ones(1, length(slct_src_list));
        area_mean_or_sum_list = {'sum'};  
        ins_or_acc_list       = {'acc'};  
        grouping_mode_list    = {'Pakistan_box'}; 
        DivByAreaOf           = 'Pakistan_box';   
        marker                = [repmat({'^'},    1, length(slct_src_list))];
        year_range_bc         = [repmat([year_list(1), year_list(end)], length(slct_data_list), 1);];
        shared_year_list      = year_range_bc;
        select_field_list     = repmat({'daily'}, 1, length(slct_src_list));
        use_yyaxis            = 1;
        show_y1_as_bar        = 0;              %/ display y1 as bar
        bar_pve_col           = [102 204 255]./255;
        show_regress_line     = 0;
        anom_mode             = 0;   %/ whether to show absolute value or anomaly (i.e. abs - clim) 
        anom_base_period      = []; 
        line_data_ylim        = [0, 50];
        st_month = 8; st_day = 10; ed_month = 8; ed_day = 24; season = []; mth_list = 0; %/ Pakistan case (whole period mean)
        basin_list            = 7;  
        yyaxis1_col = 'k';
        yyaxis2_col = 'k';
        line_colmap = [  0   0   0;
                       102 204 255]./255;
        if show_legend  fontsize = 12; fig_width  = 1500; fig_height = 600;
        else            fontsize = 16; fig_width  = 600;  fig_height = 520;     end
    end  
   % if isequal(line_setname, 'VarSet_ERA5_E_Cf_attri_Pakistan')    %/ Because of the way we computed Cf_map, we can't measure its attribution ratio.
        % recompute_SR          = 1;
        % show_legend           = 1;
        % slct_data_list        = {'ERA5_E',  'Cf_map'};   %/ compute total water budget
        % area_mean_or_sum_list = {   'sum',    'sum'};  
        % ins_or_acc_list       = {   'acc',    'acc'};  
        % grouping_mode_list    = {      '',  'global'}; 
        % slct_src_list         = {      '',  'global'};  
        % unit_conv_list        = ones(1, length(slct_src_list));
        % DivByAreaOf           = 'basin_itself';    
        % marker                = [repmat({'^'},    1, length(slct_src_list))];
        % year_range_bc         = [repmat([year_list(1), year_list(end)], length(slct_data_list), 1);];
        % shared_year_list      = year_range_bc;
        % select_field_list     = repmat({'daily'}, 1, length(slct_src_list));
        % use_yyaxis            = 1;
        % show_y1_as_bar        = 0;              %/ display y1 as bar
        % bar_pve_col           = [102 204 255]./255;
        % show_regress_line     = 0;
        % % anom_mode             = 1;   %/ whether to show absolute value or anomaly (i.e. abs - clim) 
        % % anom_base_period      = 1971:2010; 
        % % line_data_ylim        = [-300, 300];
        % anom_mode             = 0;   %/ whether to show absolute value or anomaly (i.e. abs - clim) 
        % anom_base_period      = []; 
        % line_data_ylim        = [[0, 4];[0, 4];];
        % st_month = 8; st_day = 10; ed_month = 8; ed_day = 24; season = []; mth_list = 0; %/ Pakistan case (whole period mean)
        % basin_list            = 11;  
        % yyaxis1_col = 'k';
        % yyaxis2_col = 'k';
        % line_colmap = [  0   0   0;
        %                102 204 255]./255;
        % if show_legend  fontsize = 12; fig_width  = 1500; fig_height = 600;
        % else            fontsize = 16; fig_width  = 600;  fig_height = 520;     end
    % end  
end
if plot_table
    %======================================================================
    line_setname = 'VarSet_Pm'; %<- mind this!
    %======================================================================
    if isequal(line_setname, 'VarSet_Pm')
        slct_data_list        = [    'P',       'E',   repmat({'Pm'}, 1, 10)];   %/ compute total water budget
        unit_conv_list        = ones(1, length(slct_data_list));
        area_mean_or_sum_list = repmat({'sum'}, 1, length(slct_data_list));
        ins_or_acc_list       = repmat({'acc'}, 1, length(slct_data_list));       
        grouping_mode_list    = [   {''},      {''},   repmat({'domain-wise-exact'}, 1, length(slct_data_list)-2)];
        slct_src_list         = {    '',        '',    'local',  'nonlocal TP',  'NH_MidLat_Westerly_noTP', 'NH_Polar_Westerly_noTP',  'IM_noTP', 'EAM_noTP', 'Overturning_elsewhere_noTP', 'Easterly_noTP', 'SH_Westerly', 'Elsewhere_noTP'};
        year_range_bc         = repmat([year_list(1), year_list(end)], length(slct_data_list), 1);
        shared_year_list      = year_list;
        select_field_list     = repmat({'daily'}, 1, length(slct_data_list));
        anom_mode             = 0;   %/ whether to show absolute value or anomaly (i.e. abs - clim) 
        anom_base_period      = []; 
        mth_list              = 0; %/ load seasonal values, we will to recalcuate the annual value.
%         basin_list            = 0;  %0:length(basin_catalog);
        basin_list            = 1:length(basin_catalog);
        recompute_SR          = 1;
    end
end

%----------------------------------------
if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
    shared_dates = date_array_gen('year_list', shared_year_list, 'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day);
end

%/ Display the RHc scheme if plotting WSV variables
if any(ismember(slct_data_list, valid_WSV_dataname))
    str_RHc_dqc_scheme  = sprintf('_scheme%d',RHc_dqc_scheme);
else
    str_RHc_dqc_scheme = '';
end

dataset = []; dataset.placeholder = []; WSV = []; %/ clear out data to avoid bugs

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
            if top == -1  %/ global 
                basin_name = 'global';
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
            if contains(slct_data, 'EP_ratio')
                raw_slct_data_list = {'E', 'P'};   %/ Mind the ordering!  
            elseif contains(slct_data, 'evspsbl_pr_ratio')
                raw_slct_data_list = {'evspsbl', 'pr'};    %/ Mind the ordering!  
            else
                raw_slct_data_list = {slct_data};
            end

            %/ Due to my dumb coding, no need to recompute SR each time. Once suffices.
            if ~(isequal(slct_data, slct_data_old) && isequal(model, model_old) && isequal(exp, exp_old) && isequal(ens, ens_old) && ...
                     isequal(grouping_mode, grouping_mode_old) && isequal(area_mean_or_sum, area_mean_or_sum_old) && isequal(ins_or_acc, ins_or_acc_old))
                for jj = 1:length(raw_slct_data_list)
                    raw_slct_data = raw_slct_data_list{jj};
                    %/ [IMPORTANT] Get SR Matrix
                    [SR, SR_filename_suffix, dataset] = get_SR('SR', SR, 'project_name', project_name, 'param', param, 'dataset', dataset, 'dataname', dataname, 'data_folder', data_folder, 'valid_WSV_dataname', valid_WSV_dataname,...
                                            'optimal_rr', optimal_rr, 'expmnt', expmnt, 'ldirect', ldirect, 'output_res', output_res, 'dt', dt, 'RHc_dqc_scheme', RHc_dqc_scheme, 'from_basin', from_basin, 'masterfolder', masterfolder,...
                                            'slct_data', raw_slct_data, 'model', model, 'exp', exp, 'ens', ens, 'area_mean_or_sum', area_mean_or_sum, 'ins_or_acc', ins_or_acc, 'grouping_mode', grouping_mode, 'grouping_mode_labels', grouping_mode_labels,...
                                            'regime_ver', regime_ver, 'basin_name', basin_name, 'basin_catalog', basin_catalog, 'thres_agglomerate', thres_agglomerate,...
                                            'select_field', select_field, 'DivByAreaOf', DivByAreaOf, 'mth', mth, 'str_mth', str_mth,...
                                            'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day, 'year_list', year_list_bc, 'shared_year_list', shared_year_list, 'anom_base_period', anom_base_period,...
                                            'recompute_SR', recompute_SR, 'savemat', savemat, 'savemat_prefix', savemat_prefix, 'savenc', savenc, 'nc_remark', nc_remark, 'NumWorkers', NumWorkers, 'derive_SR', derive_SR);
                end
                
                %/ Further post-processing 
                if contains(slct_data, 'EP_ratio') || contains(slct_data, 'evspsbl_pr_ratio')
                    if contains(slct_data, 'evspsbl_pr_ratio')
                        str_cmip6 = strcat('_',model,'_',exp);
                    else
                        str_cmip6 = '';
                    end

                    SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale))          = SR.(basin_name_strrep).(strcat(raw_slct_data_list{1}, str_cmip6, str_timescale))./SR.(basin_name_strrep).(strcat(raw_slct_data_list{2}, str_cmip6, str_timescale))*100;
                    SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale, '_clim')) = mean(SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale)), 'omitnan');
                    SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale, '_sd'))   = std(SR.(basin_name_strrep).(strcat(raw_slct_data_list{1}, str_cmip6, str_timescale)), 0, 2, 'omitnan'); 
                    SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, '_date', str_timescale)) = SR.(basin_name_strrep).(strcat(raw_slct_data_list{1}, str_cmip6, '_date', str_timescale)); 
                    SR.(basin_name_strrep).(strcat(slct_data, str_cmip6, str_timescale, '_unit')) = '%';
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
        
        %/ Plote Pie/Chord Chart (after all basins loaded to get consistent colorings!)
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

            elseif isequal(grouping_mode, 'IPCC-PAK') 
            else
                sort_mode             = [];       %/ no sorting
                sig_src_ordered       = sig_src;  %/ no reordering is needed
                % sig_src_colorings = nclCM('NCV_banded', length(sig_src));
                sig_src_colorings = nclCM('amwg_blueyellowred', length(sig_src)); sig_src_colorings(6,:) = [0 204 153]./255;
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

                 %/ For forward contr to Pakistan_box, append contr from all basins before plotting.
                if contains(slct_data, 'Cf_map') && isequal(grouping_mode, 'Pakistan_box')
                    if cnt == 1
                        data4plot = []; errdata4plot = []; rowName = [];
                    end
                    if cnt == length(basin_list)
                        flag_to_plot = 1;
                    else
                        flag_to_plot = 0;
                    end

                    dataname4plot   = SR.(basin_name_strrep).([slct_data, str_grouping_mode, '_src_name']);
                    ind             = findismember_loop(dataname4plot, sig_src_ordered);
                    ind_local       = findismember_loop(dataname4plot, {basin_name, 'local'});
                    data4plot       = [data4plot; SR.(basin_name_strrep).([slct_data, str_grouping_mode, str_timescale, '_clim'])(ind) ]; 
                    errdata4plot    = [errdata4plot; SR.(basin_name_strrep).([slct_data, str_grouping_mode, str_timescale, '_sd'])(ind)];
                    rowName         = [rowName; strcat(pad(num2roman(top), 'right')',{': '}, basin_name)];  %/ rowName refers to the pie labels
                    colName         = {basin_name};
                    colmap          = nclCM('amwg_blueyellowred', length(basin_list)); colmap(6,:) = [0 204 153]./255;
                    % disp(data4plot)
                else
                    flag_to_plot = 1;
                    dataname4plot   = SR.(basin_name_strrep).([slct_data, str_grouping_mode, '_src_name']);
                    ind             = findismember_loop(dataname4plot, sig_src_ordered);
                    ind_local       = findismember_loop(dataname4plot, {basin_name, 'local'});
                    data4plot       = SR.(basin_name_strrep).([slct_data, str_grouping_mode, str_timescale, '_clim']); 
                    errdata4plot    = SR.(basin_name_strrep).([slct_data, str_grouping_mode, str_timescale, '_sd']);
                    data4plot       = data4plot(ind);       %<- rmb to update the ordering!
                    errdata4plot    = errdata4plot(ind);    %<- rmb to update the ordering!
                    rowName         = sig_src_ordered_label;
                    colName         = {basin_name};
                    colmap          = sig_src_colorings;
                end

                if flag_to_plot
                    if plot_SR_chord        %/ Not recommended
                        chord_data = data4plot;
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
                            [~, I]          = sort(data4plot(ind_others), 'descend');
                            ind_ordered     = [I; ind_unattr]; 
                            T_rowname_bc    = rowName(ind_ordered); 
                            pie_data        = data4plot(ind_ordered);
                            pie_colmap      = colmap(ind_ordered,:);
                        else
                            sort_mode       = 'descend';   
                            ind_ordered     = ind;
                            T_rowname_bc    = rowName;
                            pie_data        = data4plot;
                            pie_colmap      = colmap;
                        end
                        
                        T_varname_bc  = colName;
                        str_numbering = num2roman(ind_ordered);
                        pie_label     = strcat(str_numbering, {': '}, string(round(pie_data, 1)), '%');
    
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
                            pie_data_demo   = (1:length(pie_data))/sum(1:length(pie_data))*100;
                            pie_label_demo  = strcat(string(round(pie_data_demo, 1)), '%');
                            pie_legend_demo = strcat({' '}, strrep(rowName, '_', ' '));  %/ add a small placeholder
                            pie_colmap_demo = colmap;
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
                    if plot_SR_bar          %/ Preferred (climatological mean)
                        bar_mode     = 'grouped';
                        bar_data     = data4plot;
                        bar_error    = errdata4plot;
                        bar_labels   = num2roman(1:length(rowName));
                        stack_labels = rowName;
                        bar_colmap   = colmap;
                        
                        %/ remove 'Unattributed'
                        if isequal(grouping_mode, 'domain-wise-exact')
                            ind_unattr = find(contains(stack_labels, {'Unattributed'}));
                            bar_data(ind_unattr) = [];
                            bar_labels(ind_unattr) = [];
                            stack_labels(ind_unattr) = [];
                            bar_colmap(ind_unattr,:) = [];
                        end
                        
                        if contains(grouping_mode, 'domain-wise')
                            bar_range = [0, 55];
                        else
                            bar_range = [0, 100];
                        end
                        bar_ticks = [];
                        sort_mode = 'descend';
                        show_value = 1; 
                        edgecolor = 'k';
                        y_labels  = [];
                        orientation = [];
                        box_on    = 0;
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
                        plot_stackedbar('bar_mode', bar_mode, 'bar_data', bar_data, 'bar_error', bar_error, 'bar_range', bar_range, 'bar_ticks', bar_ticks,...
                                        'sort_mode', sort_mode, 'show_value', show_value, 'bar_labels', bar_labels, 'stack_labels', stack_labels,...
                                        'bar_colmap', bar_colmap, 'edgecolor', edgecolor, 'y_labels', y_labels, 'orientation', orientation, 'box_on', box_on,...
                                        'panel_pos', panel_pos, 'fontsize', fontsize, 'linewidth', linewidth,...
                                        'titlename', titlename, 'FigName_underscore', FigName_underscore, 'plotting_folder', plotting_folder, 'savefig', savefig);
                    end
                end
                if isequal(line_setname, 'VarSet_land-ocean')
                    fprintf('%s: traceability = %.1f %s %.1f%%\n', strrep(str_timescale, '_', ''), sum(data4plot), char(177), sqrt(sum(errdata4plot.^2)));
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
                basin_name  = 'global';
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
            %===== line_ydata (ntime, nvar) =====%
            if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
                line_ydata            = nan(length(shared_dates), length(slct_data_list));
            else
                line_ydata            = nan(length(shared_year_list), length(slct_data_list));
            end
            line_lgd_str          = cell(length(slct_data_list), 1);
            line_axislabel        = cell(length(slct_data_list), 1);
            ind_ydata_firstNonNaN = nan(length(slct_data_list),  1);
            str_period            = cell(length(line_lgd_str), 1);
            for j = 1:length(slct_data_list)
                slct_data = slct_data_list{j};
                if ~isempty(model_list)  model = model_list{j};  else  model = '';  end
                if ~isempty(exp_list)    exp   = exp_list{j};    else  exp   = '';  end
                if ~isempty(ens_list)    ens   = ens_list{j};    else  ens   = '';  end
                grouping_mode = grouping_mode_list{j};

                %/ Set slct_data_callfile, slct_data_new to handle CMIP6 data    
                if ~isempty(model) && ~isempty(exp) && ~isempty(ens)
                    slct_data_new      = strcat(slct_data,'_',model,'_',exp);
                else
                    slct_data_new      = slct_data;
                end

                if ~isempty(grouping_mode)
                    str_grouping_mode = strcat('_', strrep(grouping_mode, '-', '_'));
                else
                    str_grouping_mode = '';
                end
                if contains(slct_data_new, {'Pm', 'Cf_map'})
                    ind       = findismember_loop(SR.(basin_name_strrep).([slct_data_new, str_grouping_mode, '_src_name']), slct_src_list{j});
                    line_lgd_str{j} = sprintf('%s %s', slct_src_list{j}, slct_data_new);
                else
                    ind       = 1;
                    line_lgd_str{j} = slct_data_new;
                end
                if isequal(basin_name, 'TP') && isequal(slct_src_list{j}, 'nonlocal TP')
                    %/ No need to show 'nonlocal TP' contribution for 'TP'. 
                    continue;
                else
                    if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
                        %/ Plot the daily/monthly time series
                        ind_date = findismember_loop(shared_dates, SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, '_', 'date_yyyymmdd_AllYr')));
                        line_ydata(ind_date,j) = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, '_', select_field))(ind, :);
                        line_axislabel{j} = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, '_', select_field, '_unit'));
                    else
                        %/ Plot the annual/seasonal/monthly time series
                        ind_date = findismember_loop(shared_year_list, SR.(basin_name_strrep).([slct_data_new, str_grouping_mode, '_date', str_timescale]));
                        if anom_mode
                            line_ydata(ind_date,j) = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale, str_anom_base))(ind, :);
                        else
                            line_ydata(ind_date,j) = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale))(ind, :);
                        end
                        line_axislabel{j} = SR.(basin_name_strrep).(strcat(slct_data_new, str_grouping_mode, str_timescale, '_unit'));
                    end
                end
                
                %/ Unit Conversion (if any)
                line_ydata(:,j) = line_ydata(:,j)*unit_conv_list(j);
                
                ind_ydata_firstNonNaN(j) = find(~isnan(line_ydata(:,j)), 1, 'first');
                
                str_period{j} = convertStringsToChars( strjoin(string(year_range_bc(j,:)), '-') );
            end

            if anom_mode || ~isempty(yshift)
                show_zero_line = 1;
            else
                show_zero_line = 0;
            end
            %/ Finally, do moving averaging if required.
            if ~isempty(N_running_avg)
                % line_ydata = movmean(line_ydata, N_running_avg, 1);
                % line_ydata = movmean(line_ydata, N_running_avg, 1, 'omitnan');
                line_ydata = my_movmean('data', line_ydata, 'n', N_running_avg, 'dim', 1);

                str_N_moving_avg = sprintf('_%drunningavg', N_running_avg);
            else
                str_N_moving_avg = '';
            end
            line_lgd_str = strrep(line_lgd_str, '_', ' ');
      
            %===== line_xdata (ntime, 1) ======%
            if ~isempty(slct_xdata)
                %/ E.g. we can set xdata as 'T2m_global' to compute scaling
                ind = findismember_loop(slct_data_list, slct_xdata); 
                line_xdata = line_ydata(:,ind);
                line_ydata(:,ind)           = [];         %/ remove 
                line_lgd_str(ind)           = [];         %/ remove 
                line_axislabel(ind)         = [];         %/ remove 
                ind_ydata_firstNonNaN(ind)  = [];         %/ remove 
                str_period(ind)             = [];         %/ remove 
                linestyle = 'none';
                if isequal(slct_xdata, 'T2m_global')
                    xgap_left  = 0.05;            %/ Avoid clipping or a gap at the 1st bar. 
                    xgap_right = 0.05;            %/ Avoid clipping or a gap at the last bar. 
                else
                    xgap_left  = 0.1;            %/ Avoid clipping or a gap at the 1st bar. 
                    xgap_right = 0.1;            %/ Avoid clipping or a gap at the last bar. 
                end
                map_xtickangle = 0;
            else
                if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
                    shared_dates_mmdd = strcat(string(floor(mod(shared_dates, 1e4)/1e2)), '-', string(mod(shared_dates, 1e2)));  %/ Show mmdd is enough
                    line_xdata        = 1:length(shared_dates_mmdd);
                    x_intvl           = 2;   %/ interval: 2 days
                    map_xticks        = 1:x_intvl:length(line_xdata);
                    map_xticklabels   = string(shared_dates_mmdd(map_xticks)); %/ Convert numeric dates into strings
                    map_xtickangle    = 0;    
                    xgap_left         = 1;            %/ Avoid clipping or a gap at the 1st bar. 
                    xgap_right        = 1;            %/ Avoid clipping or a gap at the last bar. 
                    
                else
                    line_xdata = 1:length(shared_year_list);
                    if contains(line_setname, 'CMIP6_MME')
                        x_intvl = 20;  %/ interval: 20 years
                    else
                        x_intvl = 10;  %/ interval: 10 years
                    end
                    map_xticks       = 1:x_intvl:length(line_xdata);
                    map_xticklabels  = shared_year_list(map_xticks); 
                    map_xtickangle   = 0;    
                    xgap_left        = 0.5;            %/ Avoid clipping or a gap at the 1st bar. 
                    xgap_right       = 0.5;            %/ Avoid clipping or a gap at the last bar.  
                end
                linestyle        = '-';
            end
            
            %/ Check and reverse the dim if line_xdata is not strictly increasing
            if diff(line_xdata(1:2)) < 0
                line_xdata = flip(line_xdata, 1);
                line_ydata = flip(line_ydata, 1);
            end

            %/ Error shading
            line_yerror_upper = []; line_yerror_lower = []; lineerror_colmap = []; lineerror_alpha = [];

            %/ Grid and label settings
            title_pos = [0.15, 1.02, 1, 0];
            map_xlabel          = '';
            map_xticks_minor    = [];
            vert_lines          = []; 
            vert_lines_col      = 'r';
            str_dates_fig       = sprintf('_%d-%d%s', shared_year_list(1), shared_year_list(end), str_timescale);

            titlename           = strcat(basin_name, str_RHc_dqc_scheme, str_dates_fig, str_trend_test, str_N_moving_avg, {' '}, line_setname);
            FigName_underscore  = strrep(strcat('multilines_', line_setname, str_regime_ver, str_anom_base, str_trend_test, str_N_moving_avg, '_', basin_name, str_dates_fig), ' ', '_');

            [trend_data, pval_data, intercept_data, CI95_data, regress_x, regress_y] = plot_multi_ylines('line_ydata', line_ydata, 'line_xdata', line_xdata, 'line_lgd_str', line_lgd_str, 'line_axislabel', line_axislabel, 'line_colmap', line_colmap,...
                      'line_yerror_upper', line_yerror_upper, 'line_yerror_lower', line_yerror_lower, 'lineerror_colmap', lineerror_colmap, 'lineerror_alpha', lineerror_alpha, 'show_y1_as_bar', show_y1_as_bar, 'pve_nve_bar_mode', pve_nve_bar_mode, 'bar_pve_col', bar_pve_col, 'bar_nve_col', bar_nve_col, 'scatter_ydata', [], 'vert_lines', vert_lines, 'vert_lines_col', vert_lines_col,...
                      'yshift', yshift, 'show_y_max', show_y_max, 'show_zero_line', show_zero_line, 'show_axislabel', show_axislabel, 'map_xlabel', map_xlabel, 'map_xticks', map_xticks, 'map_xticklabels', map_xticklabels, 'map_xtickangle', map_xtickangle,...
                      'use_yyaxis', use_yyaxis, 'line_data_ylim', line_data_ylim, 'line_data_xlim', line_data_xlim, 'marker', marker, 'markersize', markersize, 'markeredgecolor', markeredgecolor, 'fontsize', fontsize, 'linewidth', linewidth, 'linestyle', linestyle,...
                      'titlename', titlename, 'FigName_underscore', FigName_underscore, 'show_legend', show_legend, 'lgd_position', lgd_position, 'plot_legend_only', plot_legend_only,...
                      'trend_mode', trend_mode, 'trend_alpha', trend_alpha, 'trend_test', trend_test, 'asterisk', asterisk, 'alpha_level', alpha_level, 'show_regress_line', show_regress_line, 'xgap_left', xgap_left, 'xgap_right', xgap_right,...
                      'fig_width', fig_width, 'fig_height', fig_height, 'plotting_folder', plotting_folder,...
                      'title_pos', title_pos, 'yyaxis1_col', yyaxis1_col, 'yyaxis2_col', yyaxis2_col, 'show_diagonal', show_diagonal, 'box_mode', box_mode, 'fig_fmt', fig_fmt, 'png_dpi', png_dpi, 'savefig', savefig);           
            
            if contains(line_setname, 'Pm_attri')
                attri = line_ydata(:,2)./line_ydata(:,1)*100;
                VariableNames = {'Attribution (%)'};
                RowNames      = string(shared_dates);
                T = array2table(round(attri, 2), 'RowNames', RowNames, 'VariableNames', VariableNames);
                disp(T)
                fprintf('*** Total attribution: %.2f%% ***\n', sum(line_ydata(:,2))./sum(line_ydata(:,1))*100);
            end
            %/ Print the table for trend information
            if trend_mode
                %/ Append the EM data 
                if isequal(line_setname, 'MultiProduct_P_E')
                    trend_data     = [trend_data;   EM_P_trend_data(end);   EM_E_trend_data(end)];
                    pval_data      = [pval_data;    EM_P_pval_data(end);    EM_E_pval_data(end)];
                    CI95_data      = [CI95_data;    EM_P_CI95_data(end);    EM_E_CI95_data(end)];
                    line_ydata     = [line_ydata,   EM_P_data(:,end),       EM_E_data(:,end)];
                    line_lgd_str   = [line_lgd_str; {'Ensemble P'; 'Ensemble E'}];
                    line_axislabel = [line_axislabel; line_axislabel(end-1:end)];
                    str_period     = [str_period;   str_period(end-1:end,:)];
                end
                
                if isempty(slct_xdata)              %/ that means we computed the trend of the time series
                    per_x_unit = {' decade^{-1}'};
                    trend_multiplier = 10;          %/ per yr -> per dec
                elseif contains(slct_xdata, 'T2m')  %/ that means we computed regress y onto T2m
                    per_x_unit = {' K^{-1}'};
                    trend_multiplier = 1;  
                elseif contains(slct_xdata, 'E')  %/ that means we computed regress y onto E
                    per_x_unit = {' mm yr^{-1}'};
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
%                 trend_per_dec = strcat(line_lgd_str, {': '}, num2str(round(trend_data*trend_multiplier, 2)),......
%                                        char(177), num2str(round(CI95_data*trend_multiplier, 2)),...
%                                        line_axislabel, per_x_unit, {' (p = '}, num2str(round(pval_data, 3)), ')', str_sig, {' - '}, str_period);
                trend_per_dec = strcat(line_lgd_str, {': '}, num2str(round(trend_data*trend_multiplier, 2, 'significant')),......
                                       char(177), num2str(round(CI95_data*trend_multiplier, 2, 'significant')),...
                                       line_axislabel, per_x_unit, {' (p = '}, num2str(round(pval_data, 3)), ')', str_sig, {' - '}, str_period);
                fprintf('*** %s, %s, trend test = %s ***\n', basin_name, strrep(str_timescale, '_', ' '), trend_test)
                disp(trend_per_dec)
                
                %/ Get the start point of the regression line at which the true value is not a NaN
                line_ydata_ref = nan(size(line_ydata,2), 1);
                for j = 1:size(line_ydata,2)
                    cond_nonNaN_y = ~isnan(line_ydata(:,j));
                    line_xdata_nonNaN_y = line_xdata;
                    line_xdata_nonNaN_y(~cond_nonNaN_y) = nan;
                    [B,I]=sort(line_xdata_nonNaN_y, 'ascend');
                    line_ydata_ref(j) = line_ydata(I(1),j);
                end
%                 trend_per_dec_inpercentage = strcat(line_lgd_str, {': '}, num2str(round(trend_data*trend_multiplier./line_ydata_ref*100, 2)),...
%                                                     char(177), num2str(round(CI95_data*trend_multiplier./line_ydata_ref*100, 2)), {'%'}, per_x_unit,...
%                                                     {' (p = '}, num2str(round(pval_data, 3)), ')', str_sig, {' - '}, str_period);
                trend_per_dec_inpercentage = strcat(line_lgd_str, {': '}, num2str(round(trend_data*trend_multiplier./line_ydata_ref*100, 2, 'significant')),...
                                                    char(177), num2str(round(CI95_data*trend_multiplier./line_ydata_ref*100, 2, 'significant')), {'%'}, per_x_unit,...
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
                A = strcat(string(round(R,2)), str_sig);
                T = array2table(A, 'RowNames', RowNames, 'VariableNames', VariableNames);
                disp(T)
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

                if contains(slct_data_new, {'Pm','Cf_map'})
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


%% Quick check of FLEXPART particle number (should be invariate with time)

% flexpart_output_dir = '/disk/r149/tfchengac/flexpart_output_test/domfill_EA_MPI_202207-08_2022_test/';
flexpart_output_dir = '/disk/r149/tfchengac/flexpart_output_test2/domfill_EA_0.5deg_MPI_202207-08_2022_test';
% numpart_esti = 24000000;

date_flag_list = {'20220701030000';
                  '20220701060000';
                  '20220701090000'};

nest     = 0;
readp    = 1; %/ read release points (0/1)
calcarea = 1;

[header, fail]   = flex_header(flexpart_output_dir, nest, readp, calcarea); %/ Postprocessing tool from flexpart
        
for t = 1:length(date_flag_list)
    date_flag = date_flag_list{t};

%     isfile(fullfile(flexpart_output_dir, 'partposit_', date_flag))
    numpart_esti = 'unknown';
    if ismember(numpart_esti, 'unknown')
        numpart_esti = readpart10(flexpart_output_dir, date_flag, numpart_esti); %/ If numpart is unknown, input 'unknown' to estimate it.
    end

    partoutput = readpart10(flexpart_output_dir, date_flag, numpart_esti);    
    disp(numpart_esti)
    
end
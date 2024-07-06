function S = quickplot_weather(varargin)

    pnames = {'project_name', 'compo_date', 'alpha', 'bndry_data', 'contf_dataname', 'hatch_dataname', 'cont_dataname', 'U_dataname', 'V_dataname',...
              'contf_data', 'point_data', 'contf_lon', 'contf_lat', 'contf_levels', 'contf_unit', 'colmap', 'cont_data', 'cont_data_raw', 'cont_lon', 'cont_lat', 'Udata', 'Vdata', 'uv_lon', 'uv_lat', 'U2data', 'V2data',...
              'Rossby_westerly', 'Kelvin_easterly', 'RK_ratio', 'Rossby_westerly_sig', 'Kelvin_easterly_sig', 'RK_ratio_sig', 'Rossby_westerly_error', 'Kelvin_easterly_error', 'RK_ratio_error', 'RK_diff_or_not',...
              'value_in_diff', 'sig_contf_as_stipp', 'pcolor_mode', 'auto_vec_step', 'draw_refvec_only', 'draw_cbar_only', 'cbar_mode', 'cbar_location', 'cbar_fontsize', 'dataset', 'select_field', 'sig_mode_contf', 'sig_mode_cont', 'sig_mode_vec',...
              'glb_data_mode', 'glb_plateau_mode', 'plateau_hgt', 'mask_mountains', 'coast_wi', 'coast_col', 'coast_patch_col', 'backcolor', 'show_topo', 'draw_country', 'province_on', 'draw_ZZ_county', 'map_proj', 'map_lon_lower', 'map_lon_upper', 'map_lat_lower', 'map_lat_upper', 'grid_mode',...
              'compute_RK_ratio', 'show_RK_box', 'show_MJO_centroid', 'MJO_centr_filename', 'recompute_MJO_centr', 'compute_contf_AWM', 'compute_AWM_lon_extent', 'compute_AWM_lat_extent', 'show_AMW_reg', 'AMW_reg_color', 'contf_AWM', 'str_contf_AWM', 'contf_AWM_error',...
              'all_in_one_format', 'fontsize', 'title_fontsize', 'color', 'markersize', 'markerfacecolor', 'linewi', 'cont_linewi', 'grid_linewi', 'plot_or_not', 'ax_panel', 'data_folder', 'plotting_folder', 'titlename', 'title_pos', 'savefig'};
    
    dflts  = cell(length(pnames), 1);
    
    [          project_name,   compo_date,   alpha,  bndry_data,   contf_dataname,   hatch_dataname,  cont_dataname,    U_dataname,   V_dataname,...  
               contf_data,   point_data,    contf_lon,   contf_lat,   contf_levels,  contf_unit, colmap,   cont_data,   cont_data_raw,   cont_lon,   cont_lat,   Udata,   Vdata,   uv_lon,   uv_lat,   U2data,   V2data,...
               Rossby_westerly,   Kelvin_easterly,   RK_ratio,   Rossby_westerly_sig, Kelvin_easterly_sig, RK_ratio_sig, Rossby_westerly_error, Kelvin_easterly_error, RK_ratio_error, RK_diff_or_not,...
               value_in_diff, sig_contf_as_stipp,   pcolor_mode,   auto_vec_step, draw_refvec_only, draw_cbar_only,    cbar_mode,  cbar_location,    cbar_fontsize,    dataset,   select_field,   sig_mode_contf,   sig_mode_cont,   sig_mode_vec,...
               glb_data_mode,   glb_plateau_mode,   plateau_hgt,   mask_mountains,   coast_wi,   coast_col,  coast_patch_col, backcolor, show_topo,    draw_country, province_on,   draw_ZZ_county,   map_proj,  map_lon_lower,   map_lon_upper,   map_lat_lower,   map_lat_upper, grid_mode,...
               compute_RK_ratio,   show_RK_box,   show_MJO_centroid, MJO_centr_filename, recompute_MJO_centr, compute_contf_AWM, compute_AWM_lon_extent, compute_AWM_lat_extent, show_AMW_reg, AMW_reg_color, contf_AWM, str_contf_AWM, contf_AWM_error,...
               all_in_one_format,   fontsize,   title_fontsize,  color, markersize,   markerfacecolor, linewi,   cont_linewi, grid_linewi,   plot_or_not,   ax_panel,   data_folder, plotting_folder,   titlename,   title_pos,   savefig] = ...
                                            internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Apr 25, 2024
    %/
    %/ Description: This function aims to quickly plot weather map based on 
    %/              dataset (only when contf_data, cont_data, etc. are not given)
    %/=====================================================================
    
    % ax_panel = []; all_in_one_format = []; 
    % pcolor_mode = 0; province_on = 0; draw_ZZ_county = 0;  
    % contf_data = []; contf_lon = []; contf_lat = []; contf_levels = []; colmap = []; contf_unit = '';
    % cont_data = []; cont_data_raw = []; cont_lon = []; cont_lat =[]; 
    % Udata=[]; Vdata=[]; uv_lon=[]; uv_lat=[]; U2data = []; V2data =[];
    % contf_AWM = []; contf_AWM_error = []; Rossby_westerly = [];
    % Kelvin_easterly = []; RK_ratio = [];

    %/ Initialization
    contf_data_events = []; cont_data_events = []; Udata_events = []; Vdata_events = [];
    cont_levels = []; cont_colmap = [];
    cont_unit = ''; cbar_interval = []; cbar_YTick = []; cbar_YTickLabel = [];
    hatch_data = []; hatch_lon = []; hatch_lat = []; hatch_thres_pve  = []; hatch_thres_nve  = []; color_hatch_pve = []; color_hatch_nve = [];   hatch_mode = [];    hatch_linewi  = []; hatch_intvl= [];
    cont_labelsize=[]; cont_label_col=[]; skip_zero_cont = 1; 
    vec_step_lon=[]; vec_step_lat=[]; vector_levels = []; vector_color=[]; vector_edgecolor=[]; vecscale=[]; vecscale2=[]; shaftwidth=[];
    vector2_color = []; vector2_edgecolor = [];
    headlength=[]; vec_lbs=[]; vec_mag_ref=[]; vec_ref_fontsize=[]; vec_ref_lat_shift=[]; str_show_topo = ''; 
    contf_multiply = 1;  cont_multiply = 1; contf_add = 0;  cont_add = 0;  fig_fmt = 'pdf';
    
    if isempty(select_field)
        select_field = 'subdaily';
    end

    %/ Make sure select_field is a string and only one field is quiried
    if iscell(select_field)
        if length(select_field) ~= 1
            error('Multiple ''select_field'' is not allowed!');
        else
            select_field = select_field{:};
        end
    end
    
    %/ Formatting of dates / timescales
    if contains(select_field, 'subdaily')
        time_unit     = 'subdaily';
        fld_date      = 'date_UTC_AllYr';
        fld_date_clim = [];
        str_flds_date = [];
        Lr            = 12;
        date_fmt      = 'yyyymmddHHMM';

    elseif contains(select_field, 'daily')
        time_unit     = 'daily';
        fld_date      = 'date_yyyymmdd_AllYr';
        fld_date_clim = 'date_yyyymmdd_AllYr_clim';
        str_flds_date = 'str_daily_dates';
        Lr            = 8;
        date_fmt      = 'yyyymmdd';

    elseif contains(select_field, 'pentad')
        time_unit     = 'pentad';
        fld_date      = 'date_yyyyptd';
        fld_date_clim = 'date_yyyyptd_clim';
        str_flds_date = 'str_ptd_dates';
        Lr            = 6;
        date_fmt      = 'index';  %/ Assume the input compo_date denote time indices, cos no date format recognized

    elseif contains(select_field, 'monthly')
        time_unit     = 'monthly';
        fld_date      = 'date_yyyymm';
        fld_date_clim = 'date_yyyymm_clim';
        str_flds_date = 'str_monthly_dates';
        Lr            = 6;
        date_fmt      = 'yyyymm';
    else
        error('code not ready!');
    end

    %/ Defaults
    plateau_col = [255 51 204]./255;
    
    marker = 'o';
    if isempty(color)           
        if ~isempty(AMW_reg_color)
            color = AMW_reg_color;
        else
            color = [0 0 0];             
        end
    else
        if ~isempty(AMW_reg_color)
            color = [color; AMW_reg_color];
        end
    end
    if isempty(contf_unit)      contf_unit = '';                      end
    if isempty(glb_data_mode)   glb_data_mode = 0;                    end
    if isempty(fontsize)        fontsize      = 20;                   end
    if isempty(title_fontsize)  title_fontsize = fontsize*0.8;        end
    if isempty(markersize)      markersize  = 3;                      end
    if isempty(linewi)          linewi      = 2.5;                    end
    if isempty(pcolor_mode)     pcolor_mode = 0;                      end
    if isempty(grid_mode)       grid_mode   = 2;                      end
    if isempty(plot_or_not)     plot_or_not = 1;                      end
    if isempty(cbar_mode)       cbar_mode     = 1;                    end 
    if isempty(cbar_location)   cbar_location = 'eastoutside';        end
    if isempty(alpha)           alpha = 0.05;                         end  %/ sig. level for one-sample t-test 
    if isempty(coast_wi)        coast_wi        = 2.5;                end
    if isempty(coast_col)       coast_col       = [.4 .4 .4];         end
    if isempty(backcolor)       backcolor       = 'none';             end
    if isempty(coast_patch_col) coast_patch_col = [];                 end
    if isempty(map_lon_lower)   map_lon_lower = 0;                    end
    if isempty(map_lon_upper)   map_lon_upper = 359;                  end
    if isempty(map_lat_lower)   map_lat_lower = -90;                  end
    if isempty(map_lat_upper)   map_lat_upper = 90;                   end
    if isempty(draw_cbar_only)  draw_cbar_only = 0;                   end
    if isempty(draw_refvec_only)draw_refvec_only = 0;                 end
    if isempty(contf_dataname)  contf_dataname = '';                  end
    if isempty(hatch_dataname)  hatch_dataname = '';                  end
    if isempty(cont_dataname)   cont_dataname = '';                   end
    if isempty(U_dataname)      U_dataname = '';                      end
    if isempty(V_dataname)      V_dataname = '';                      end
    if draw_cbar_only           str_cbar = '_cbar';      else str_cbar = '';    end
    if draw_refvec_only         str_refvec = '_refvec';  else str_refvec = '';  end

    %====== contf ======%
    if ~isempty(contf_dataname) || ~isempty(contf_data)
        if isempty(contf_data)
            flag_input = 0;
            if contains(select_field, '_clim')
                date = dataset.(contf_dataname).(fld_date_clim);
            else
                date = dataset.(contf_dataname).(fld_date);
            end

            if numel(num2str(compo_date(1))) == 8 && isequal(time_unit, 'monthly')
                ind  = findismember_loop(date, floor(compo_date./1e2));  %/ Then we will subset the month of the date in the monthly data
            else
                ind  = findismember_loop(date, compo_date);
            end

            if isempty(ind) 
                error('[quickplot_weather] Empty ind! Check compo_date!');
            elseif length(ind) ~= length(compo_date)
                warning('[quickplot_weather] Not all compo_date are available from the data. Filling missing data with NaN (useful when showing near-real time data).');
            end
    
            contf_lon = dataset.(contf_dataname).lon;
            contf_lat = dataset.(contf_dataname).lat;

            if isequal(sig_mode_contf, '_ttest')
                contf_data_events = dataset.(contf_dataname).(select_field)(:,:,ind); %/ To be output for other use (e.g., ttest2_sig_fn)
                [contf_data_sig, contf_data, ~] = ttest_sig_fn(contf_data_events, alpha, 3, 0);
                if sig_contf_as_stipp == 1
                    point_data = datasig2point(contf_data_sig, contf_lon, contf_lat); %/ My function to convert contf_data_sig to point_data
                elseif sig_contf_as_stipp == 2
                    contf_data = contf_data_sig;
                    point_data = [];
                else
                    point_data = [];
                end
                % size(point_data)
            elseif ~isempty(sig_mode_contf) 
                contf_data        = mean(dataset.(contf_dataname).(select_field)(:,:,ind), 3, 'omitnan');
                contf_data_sig    = mean(dataset.(contf_dataname).(strcat(select_field, sig_mode_contf))(:,:,ind), 3, 'omitnan');

                if sig_contf_as_stipp == 1
                    point_data  = contf_data_sig;
                elseif sig_contf_as_stipp == 2
                    contf_data  = contf_data_sig;
                    point_data  = [];
                else
                    point_data  = [];
                end
            else
                contf_data_events = dataset.(contf_dataname).(select_field)(:,:,ind);
                [~, contf_data, ~] = ttest_sig_fn(contf_data_events, alpha, 3, 0);
                point_data = []; 
            end
        else
            if sig_contf_as_stipp == 0
                point_data = [];   %/ make sure stippling is turned off as inquried
            end
            flag_input = 1;
        end

        if isempty(markerfacecolor)
            markerfacecolor = 'k';
        end

        cbar_interval   = 2;
        cbar_YTick      = [];
        cbar_YTickLabel = []; 
        if isempty(contf_levels) && isempty(colmap) 
            if ismember(contf_dataname, {'ERA5_P'}) || isequal(contf_dataname, 'pr')
                if contains(select_field, 'monthly')
                    contf_unit = 'mm/month';
                    CONV = 30;  
                else
                    contf_unit = 'mm/day';
                    CONV = 1;   
                end 
                if contains(select_field, {'_MJO'})
                    if value_in_diff
                        contf_levels  = CONV*(-3.6:0.6:3.6)/2;
                    else
                        contf_levels  = CONV*(-3.6:0.6:3.6);
                    end
                    colmap = nclCM('precip_diff_12lev', length(contf_levels)-1);
                    % colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                else
                    if value_in_diff
                        contf_levels  = CONV*(-3.6:0.6:3.6);
                        colmap = nclCM('precip_diff_12lev', length(contf_levels)-1);
                    else
                        contf_levels  = CONV*(0:2:22);
                        colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                    end
                    cbar_interval = 1;
                end
    
            elseif ismember(contf_dataname, {'ERA5_LcP'}) || isequal(contf_dataname, 'Lcpr')
                
                % Lc = 2.26e6;  %/ Latent heat for condensation (J kg-1)  
                % if contains(select_field, 'daily')
                %     CONV = 1/24/3600;   
                % elseif contains(select_field, 'monthly')
                %     CONV = 1/30/24/3600;   
                % else
                %     error('code not set!');
                % end
                % contf_multiply = Lc*CONV;
                contf_unit = 'W m^{-2}';  %/ NOTE: W m-2 == J s-1 m-2
    
                if contains(select_field, {'_MJO'})
                    if value_in_diff
                        contf_levels  = -45:5:45;
                    else
                        contf_levels  = -140:20:140;
                    end
                else
                    if value_in_diff
                        contf_levels  = -160:20:160;
                    else
                        contf_levels  = 0:50:850;
                    end
                    cbar_interval = 1;
                end
                colmap = nclCM('BlueWhiteOrangeRed', length(contf_levels)-1);
    
            elseif ismember(contf_dataname, {'ERA5_E'}) || isequal(contf_dataname, 'evspsbl')
                if contains(select_field, 'monthly')
                    contf_unit = 'mm/month';
                    CONV = 30;  
                else
                    contf_unit = 'mm/day';
                    CONV = 1;   
                end 
                if contains(select_field, {'_MJO'})
                    if value_in_diff
                        contf_levels  = CONV*(-3.6:0.6:3.6)/2;
                    else
                        contf_levels  = CONV*(-3.6:0.6:3.6);
                    end
                    colmap = nclCM('precip_diff_12lev', length(contf_levels)-1);
                    % colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                else
                    if value_in_diff
                        contf_levels  = CONV*(-3.6:0.6:3.6)/5; %/ Seems like change in evaporation is small (?)
                        colmap = nclCM('precip_diff_12lev', length(contf_levels)-1);
                        cbar_interval = 2;
                    else
                        contf_levels  = CONV*(0:2:22);
                        colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                        cbar_interval = 1;
                    end
                end
    
            elseif contains(contf_dataname, {'S500'}) || contains(contf_dataname, {'S5'})  %/ dry static stability (assuming in K/Pa)
                contf_multiply = 10^4;
                contf_unit     = 'K (100 hPa)^{-1}';
                if value_in_diff
                    % contf_levels   = (-1.6:0.2:1.8);
                    % contf_levels   = (-1.3:0.1:1.3);
                    % contf_levels   = (-1.05:0.05:1.05);
                    contf_levels   = (0.25:0.05:1.35);
                    NoOfColors     = length(contf_levels)-1; 
                    colmap         = nclCM('temp_19lev', NoOfColors);
                    % contf_levels(1:11) = [];  %/ Needless to show the full scale
                    % colmap(1:11,:)     = [];  %/ Needless to show the full scale
                    
                    cbar_interval = 4;
                    % disp(contf_levels(2:cbar_interval:end-1))
                else
                    contf_levels   = (2.5:0.5:8.5); 
                    NoOfColors    = length(contf_levels)-1; 
                    colmap         = nclCM('temp_19lev', NoOfColors);
                end
    
            elseif ismember(contf_dataname, {'pw'})     %/ pw == tcwv
                contf_unit    = 'mm';
                contf_levels  = 0:5:50;               %/ in mm (instantaneous)
                colmap        = brewermap(length(contf_levels)-1, 'BuGn');
    
            elseif contains(contf_dataname, {'omega500', 'W500', 'wap', 'W5'}) 
                contf_multiply  = 24*3600/100;  %/ Pa/s -> hPa/day
                contf_unit      = 'hPa/day';  
                if contains(select_field, 'anom')
                    % if value_in_diff
                    contf_levels = (-6:1:6)*2;
                    % else
                        % contf_levels = (-6:1:6)*2;
                    % end
                else
                    if value_in_diff
                        contf_levels = -28:4:28;
                    else
                        contf_levels = -70:10:70;
                    end
                end
                % contf_levels    = [flip([40:40:280])*-1, [40:40:280]];   %/ in hPa/day
                colmap = flip(nclCM('sunshine_diff_12lev', length(contf_levels)-1), 1);
    
            elseif ismember(contf_dataname, {'MC1000to850'})
                contf_unit      = 'mm/day';  
                if contains(select_field, 'anom')
                    contf_levels =  -0.9:0.15:0.9;

                elseif contains(select_field, 'daily_MJO')
                    if value_in_diff
                        % contf_levels = -0.6:0.1:0.6;
                        contf_levels = -0.9:0.15:0.9;
                        % contf_levels = -1.2:0.2:1.2;
                    else
                        contf_levels = -2.4:0.4:2.4;
                        % contf_levels = -3.6:0.6:3.6;
                        % contf_levels = -4.2:0.7:4.2;
                    end
                else
                    if value_in_diff
                        contf_levels = -24:3:24;
                    else
                        contf_levels = -70:10:70;
                    end
                end
                % colmap = flip(nclCM('WhiteBlueGreenYellowRed', length(contf_levels)-1), 1); %/ Recommended?
                % colmap = flip(nclCM('cmp_flux', length(contf_levels)-1), 1);
                % colmap = flip(nclCM('hotcolr_19lev', length(contf_levels)-1), 1);
                colmap = nclCM('CBR_drywet', length(contf_levels)-1);

            elseif ismember(contf_dataname, {'IVT'})    
                contf_unit   = 'kg/m/s';
                contf_levels = 0:50:500;               %/ in kg/m/s 
                colmap       = jet(length(contf_levels)-1);
                
            elseif isequal(contf_dataname, 'Z850')
                contf_unit    = 'm';
                if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                    contf_levels = -11:1:11;
                    colmap       = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
    
                elseif ismember(select_field, {'pentad_clim_AC'})
                    contf_levels = -70:10:70;
                    colmap       = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
    
                else
                    contf_levels = 1420:10:1540;
                    colmap       = brewermap(length(contf_levels), '*RdYlBu');
                end
    
            elseif ismember(contf_dataname, {'uv850mag'})   
                contf_unit    = 'm/s';
                contf_levels  = 0:2:22;               %/ in m/s 
                colmap        = brewermap(length(contf_levels)-1, '*Spectral');
    
            elseif contains(contf_dataname, 'OLR') || isequal(contf_dataname, 'rlut')
                contf_unit    = 'W m^{-2}';
                if ismember(select_field, {'pentad_clim_CISO_AC_abs_ratio'})
                    contf_levels  = 0:.2:2.2;               
                    NoOfColors     = length(contf_levels)-1;
                    colmap        = brewermap(NoOfColors, '*PiYG');
    %                 colmap(1:(NoOfColors-1)/2,:) = flip(colmap(1:(NoOfColors-1)/2,:), 1);
    %                 colmap((NoOfColors-1)/2+1:end,:) = flip(colmap((NoOfColors-1)/2+1:end,:), 1);
                    cbar_interval = 2;
                    backcolor     = [0.5 0.5 0.5]; %/ to indicate NaN.
                    pcolor_mode   = 1;
                    
                elseif contains(select_field, {'_MS'})
                    contf_levels  = 190:10:310;             %/ in W m^-2 
                    colmap = my_colormap(length(contf_levels)-1, 'radar_12lev');
                    cbar_interval = 2;
                    
                elseif contains(select_field, {'_AC'})
                    contf_levels   = -60:10:60;                %/ in W m^-2 
                    NoOfColors     = length(contf_levels)-1;
                    colmap         = nclCM('hotcold_18lev', NoOfColors);
    
                elseif contains(select_field, {'_MJO'})
                    if value_in_diff
                        % contf_levels   = [-9:1:9];                %/ in W m^-2 
                        contf_levels   = -14:2:14;                %/ in W m^-2 
                    else
                        contf_levels   = -35:5:35;                %/ in W m^-2 
                    end
                    NoOfColors     = length(contf_levels)-1;
                    colmap         = nclCM('BlueYellowRed', NoOfColors);
                    % colmap         = nclCM('hotcold_18lev', NoOfColors);
                    
                elseif contains(select_field, {'_CISO', '_TE', '_IAID'})
                    contf_levels  = -18:3:18;               %/ in W m^-2 
                    NoOfColors     = length(contf_levels)-1;
                    colmap = nclCM('hotcold_18lev', NoOfColors);
    
                else
                    contf_levels  = 190:10:310;             %/ in W m^-2 
                    colmap = my_colormap(length(contf_levels)-1, 'radar_12lev');
                    cbar_interval = 2;
                end
                
            elseif contains(contf_dataname, {'I_Wang1988'})   
                contf_unit    = ''; 
                if contains(select_field, {'_anom'})
                    contf_levels = (-2.4:0.4:2.4)./100;
                else
                    if value_in_diff
                        contf_levels = -0.18:0.03:0.18;
                    else
                        contf_levels = 0:0.1:1;
                    end
                end
                NoOfColors = length(contf_levels)-1;
                colmap = nclCM('BlueDarkRed18', NoOfColors);
    
            elseif contains(contf_dataname, {'alpha_Wang1988'})   
                contf_multiply = 1/1e-3;
                contf_unit    = 'x10^{-3}';
                if contains(select_field, {'_anom'})
                    contf_levels = -0.24:0.04:0.24;
                else
                    if value_in_diff
                        contf_levels = (0.8:0.1:2);
                    else
                        contf_levels = 7.5:0.5:14;
                    end
                end
                NoOfColors = length(contf_levels)-1;
                colmap = nclCM('temp_19lev', NoOfColors);
    
            elseif contains(contf_dataname, {'C0_Wang1988'})   
                contf_unit    = 'm s^{-1}';
                if contains(select_field, {'_anom'})
                    contf_levels = -0.6:0.1:0.96;
                else
                    if value_in_diff
                        contf_levels = 2.2:0.2:4.8;
                    else
                        contf_levels = 36:2:60;
                    end
                end
                NoOfColors = length(contf_levels)-1;
                colmap = nclCM('GreenMagenta16', NoOfColors);
    
            elseif contains(contf_dataname, {'C1_Wang1988'})   %/ == C0*sqrt(1-I)
                % if ~isempty(find(~isreal(contf_data), 1))
                %     contf_data(~isreal(contf_data)) = nan; %/ set complex number to nan 
                % end
                contf_unit    = 'm s^{-1}';
                if contains(select_field, {'_anom'})
                    contf_levels = -0.6:0.1:0.96;
                else
                    if value_in_diff
                        contf_levels = -3.6:0.4:3.6;
                    else
                        contf_levels = 24:2:48;
                    end
                end
                NoOfColors = length(contf_levels)-1;
                colmap = nclCM('GreenMagenta16', NoOfColors);
    
            elseif contains(contf_dataname, {'dq_bar_Wang1988'})   
                contf_multiply = 1/1e-3;
                contf_unit    = 'g/kg';
                if contains(select_field, {'_anom'})
                    contf_levels = (-0.3:0.05:0.3);
                else
                    if value_in_diff
                        contf_levels  = -0.2:0.2:2.2;
                    else
                        contf_levels  = 1.5:0.5:8.5;
                    end
                end
                NoOfColors = length(contf_levels)-1;
                colmap = nclCM('precip_diff_12lev', NoOfColors);
    
            elseif ismember(contf_dataname, {'TSR', 'SSR'})    
                contf_unit    = 'W m^{-2}';
                contf_levels  = 160:20:400;               %/ in W m^-2 
                colmap = my_colormap(length(contf_levels)-1, 'radar_12lev');
                
            elseif ismember(contf_dataname, {'TISR'})    
                contf_unit    = 'W m^{-2}';
                contf_levels  = 100:50:600;               %/ in W m^-2 
                colmap        = cmocean('solar', length(contf_levels)-1);
    
            elseif ismember(contf_dataname, {'dE_TOA'}) 
                contf_unit    = 'W m^{-2}';
    	        contf_levels  = (-220:20:220)/2;               %/ in W m^-2 
                colmap        = brewermap(length(contf_levels)-1, '*PuOr');
                NoOfColors    = length(contf_levels)-1;
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [1 1 1; 1 1 1];
    
            elseif ismember(contf_dataname, {'D_vm850to1000', 'D_vm850to1000_sm'})
                contf_unit = '10^{-6} s^{-1}';
                contf_multiply = 1e6;
                contf_levels  = -8:1:8;  %/ in 1e-6 s^-1 
                colmap = my_colormap(length(contf_levels)-1, 'amwg_blueyellowred_16lev');
                    
            elseif ismember({contf_dataname}, {'MSE_vi'})
                contf_multiply = 1e-9;
                contf_unit     = '10^{9} J m^{-2}';
                contf_levels   = linspace(2.925, 3.2, 12);             %/ in J m^-2 
                colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
    
            elseif ismember({contf_dataname}, {'dMSE_dt_vi', 'MSE_Res'})    
                contf_unit    = 'W m^{-2}';
                if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                    contf_levels  = -55:5:55;               %/ in W m^-2, -ve -> convection
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdYlBu');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
                    
                elseif ismember(select_field, {'pentad_clim', 'pentad_clim_AC'})
                    contf_levels  = -55:5:55;             %/ in W m^-2 
                    NoOfColors    = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdYlBu');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
                end
                
            elseif ismember({contf_dataname}, {'U_dMSE_dx_vi',   'V_dMSE_dy_vi',  'W_dMSE_dP_vi'})
                contf_multiply = -1;   %/ times a minus sign -> advection term.
                contf_unit     = 'W m^{-2}';
                if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                    contf_levels  = [-55:5:55];               %/ in W m^-2, -ve -> convection
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdYlBu');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
                    
                elseif ismember(select_field, {'pentad_clim', 'pentad_clim_AC'})
                    contf_levels  = [-110:10:110];             %/ in W m^-2 
                    NoOfColors    = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdYlBu');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
                end
    
            elseif ismember({contf_dataname}, {'netSW'})            
                contf_unit    = 'W m^{-2}';
                contf_levels  = 70:5:125;             %/ in W m^-2 
                colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
    
            elseif ismember({contf_dataname}, {'netLW'})
                contf_unit    = 'W m^{-2}';
                contf_levels  = -230:10:-120;             %/ in W m^-2 
                colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
    
            elseif ismember({contf_dataname}, {'netSW_LW', 'QR_vi'})
                contf_unit    = 'W m^{-2}';
                contf_levels  = -220:20:220;              %/ in W m^-2 
                NoOfColors    = length(contf_levels)-1;
                colmap = brewermap(NoOfColors, '*RdYlBu');
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
    
            elseif ismember({contf_dataname}, {'Q1_vi', 'Q2_vi'})
                contf_unit    = 'W m^{-2}';
                contf_levels  = -440:40:440;              %/ in W m^-2 
                NoOfColors    = length(contf_levels)-1;
                colmap = brewermap(NoOfColors, '*RdYlBu');
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
    
            elseif ismember(contf_dataname, {'SSHF', 'SLHF', 'SSHFLand', 'SLHFLand', 'SSHFOcean', 'SLHFOcean'})    
                contf_unit    = 'W m^{-2}';
                contf_levels  = 0:20:220;               %/ in W m^-2 
                NoOfColors    = length(contf_levels)-1;
    %             colmap        = brewermap(NoOfColors, '*PiYG');
                colmap        = my_colormap(NoOfColors, 'precip3_11lev');
    
            elseif contains(contf_dataname, {'T2m'}) || isequal(contf_dataname, 'tas')
                contf_unit    = sprintf('%sC', char(176));
                contf_add     = -273.15;   %/ K to deg C
                if value_in_diff
                    contf_levels  = (-6:1:6);
                else
                    contf_levels  = (10:2:36);
                end
                colmap        = nclCM('GMT_panoply', length(contf_levels)-1);
                % colmap        = brewermap(length(contf_levels)-1, '*spectral');
    
            elseif contains(contf_dataname, {'SST'}) || ismember(contf_dataname, {'tos'})
                % coast_patch_col = 'k';  %/ IMPORTANT: Set continent color to overlay the SST data (some CMIP6 models may unrealistically output non-NaN SST values over continent)
                coast_patch_col = [.5 .5 .5];  %/ IMPORTANT: Set continent color to overlay the SST data (some CMIP6 models may unrealistically output non-NaN SST values over continent)
                % coast_patch_col = [94, 132, 89]./255;
                % coast_patch_col = [141, 131, 97]./255;
    
                contf_unit    = sprintf('%sC', char(176));
                if contains(select_field, '_MJO')
                    pcolor_mode = 1;
                    if value_in_diff
                        contf_levels  = (-0.12:0.02:0.12);
                    else
                        contf_levels  = -0.21:0.03:0.21;
                    end
                elseif contains(select_field, {'monthly', 'daily'})
                    if contains(select_field, 'anom')
                        if value_in_diff
                            contf_levels  = (-0.64:0.08:0.64)/2;
                            cbar_interval = 2;
                        else
                            contf_levels  = -0.64:0.08:0.64;
                        end
                    else
                        contf_add     = -273.15;   %/ K to deg C (only for raw data)
                        if value_in_diff
                            contf_levels  = -0.25:0.25:4.25;
                            cbar_interval = 4;
                        else
                            contf_levels  = (21:1:33);
                        end
                    end
                end
                colmap = nclCM('temp_diff_18lev', length(contf_levels)-1);
                % if contains(select_field, 'anom') || value_in_diff
                %     n = length(contf_levels)-1;
                %     colmap(n/2:n/2+1,:) = [1 1 1; 1 1 1];
                % end
                % disp(contf_levels(2:2:end-1))
    
            elseif contains(contf_dataname, {'q1000', 'q850', 'hus', 'qs'}) 
                % coast_patch_col = 'k';  %/ IMPORTANT: Set continent color to overlay the SST data (some CMIP6 models may unrealistically output non-NaN SST values over continent)
                % coast_patch_col = [.15 .15 .15];  %/ IMPORTANT: Set continent color to overlay the SST data (some CMIP6 models may unrealistically output non-NaN SST values over continent)
    
                contf_unit        = 'g/kg';
                contf_multiply    = 1000;   %/ kg/kg to g/kg
    
                if contains(select_field, {'_MJO'})
                    if value_in_diff
                        % contf_levels  = (-0.36:0.06:0.36);
                        contf_levels  = (-0.30:0.05:0.30);
                    else
                        % contf_levels  = (-2.4:0.4:2.4)./2;
                        contf_levels  = (-0.6:0.1:0.6);
                    end
                else
                    if value_in_diff
                        contf_levels  = -4.5:0.5:4.5;
                        % contf_levels(1:8) = [];
                        % colmap(1:8, :) = [];
                    else
                        if contains(contf_dataname, {'1000'})
                            contf_levels  = 1.25:1.25:21.25;
                        else
                            contf_levels  = 0:1:16;
                        end
                    end
                end
                % colmap        = nclCM('CBR_drywet', length(contf_levels)-1);
                colmap        = nclCM('precip_diff_12lev', length(contf_levels)-1);
                % disp(contf_levels(2:2:end-1))
    
            elseif ismember(contf_dataname, {'psi200'})
                %/ NOTE:
                %/ In geostrophic balance (See Holton or notes in Notion):
                %/   High (low) stream function -> high (low) geopotential height in N.H.
                %/                              -> low (high) geopotential height in S.H.
                %/ https://www.cpc.ncep.noaa.gov/products/CDB/Tropics/figt22.shtml
                
                contf_multiply = 10^-6;
                contf_unit     = '10^{6} m^{2} s^{-1}';
                contf_levels   = -20:2:20; %/ This makes sure level interval to be const. 
                NoOfColors      = length(contf_levels)-1; 
                colmap         = brewermap(NoOfColors,'*RdBu');
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [1 1 1; 1 1 1;];
                
            elseif ismember(contf_dataname, {'S','S1','S2'}) %/ Static stability (computed by read_CMIP.m)
                            
                contf_multiply = 10^11; %/ scaling
                contf_unit     = '10^{-11} s^{-2}';
                
                if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                    contf_levels   = -30:3:30;  %/ +ve: source, -ve: sink.
                elseif ismember(select_field, {'pentad_clim_AC'})
                    contf_levels   = -50:5:50;  %/ +ve: source, -ve: sink.
                else
                    contf_levels   = -50:5:50;  %/ +ve: source, -ve: sink.
                end
                NoOfColors      = length(contf_levels)-1; 
                colmap         = brewermap(NoOfColors, '*RdGy');
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [1 1 1; 1 1 1;];
            else
                warning('The levels of the request data are not specified! plot_contfmap will Auto-generate contf_levels with the default colmap.');
            end
        end
        
        %/ Avoid performing unit conversion when contf_data is directly input
        if flag_input == 0
            contf_data = contf_data*contf_multiply + contf_add;
            contf_data_events = contf_data_events*contf_multiply + contf_add; %/ for consistency
        end
    end

    % %====== hatching =====%
    % if ~isempty(hatch_dataname) || ~isempty(hatch_data)
    %     hatch_mode = 1;
    % 
    %     if ~isempty(hatch_data)   
    %         hatch_lon       = dataset.(contf_dataname).lon;
    %         hatch_lat       = dataset.(contf_dataname).lat;
    %         hatch_thres_pve = 0;
    %         hatch_thres_nve = 0;
    %         color_hatch_pve = [0 0 0]; 
    %         color_hatch_nve = [0 0 0]; 
    %         hatch_linewi    = 3;
    %         hatch_intvl     = 10;
    % 
    %     else
    %         if isequal(date_fmt, 'yyyymmddHHMM')        
    %             X_dates     = datetime2int('X_datetime', dataset.(hatch_dataname).date_UTC_AllYr, 'format', 'yyyymmddHHMM');
    %             ind_date    = findismember_loop(X_dates, compo_date);
    % 
    %         elseif isequal(date_fmt, 'yyyymmdd')        
    %             X_dates     = dataset.(hatch_dataname).date_yyyymmdd_AllYr; 
    %             ind_date    = findismember_loop(X_dates, compo_date);
    % 
    %         elseif isequal(date_fmt, 'index')          
    %             ind_date    = compo_date;  %/ assume compo_date store pentad.
    %         end
    % 
    %         hatch_data      = mean(dataset.(hatch_dataname).(select_field)(:,:,ind_date), 3, 'omitnan');
    %         hatch_lon       = dataset.(hatch_dataname).lon;
    %         hatch_lat       = dataset.(hatch_dataname).lat;
    % 
    %         if isequal(hatch_dataname, 'OLR') && isempty(sig_mode_contf)
    %             hatch_thres_pve = 300;                 %/ surrogates for heat sinks               (Webster et al. 1998)
    %             hatch_thres_nve = 220;                 %/ surrogates for tropical deep convection (Webster et al. 1998)
    %             color_hatch_pve = [204 153 0]./255;    %? [171 165  84]./255;  %/ dark brown
    %             color_hatch_nve = [35 218 254]./255;  %/ cyan
    %             hatch_linewi    = 3;
    %             hatch_intvl     = 14;
    %         else
    %             hatch_thres_pve = 0;
    %             hatch_thres_nve = 0;
    %             color_hatch_pve = [0 0 0]; 
    %             color_hatch_nve = [0 0 0]; 
    %             hatch_linewi    = 1.5;
    %             hatch_intvl     = 14;
    %         end
    %     end
    % end

    %====== show_topo =====% (Caveat: will overlap the earlier contf.) 
    if show_topo
        if ~isempty(contf_dataname)
            error('''show_topo'' is on while contf_dataname is given! Turn off either one to avoid overlapping.');
        end
        fig_fmt = 'png';  %/ to minimize fig size
        str_show_topo = '_topo';
        %/ NOTE: m_elev  outputs 1 degree res. 
        %/       m_tbase outputs much higher res. (Default)
        
%         [contf_data, lon_2D_topo, lat_2D_topo]=m_elev([map_lon_lower, map_lon_upper, map_lat_lower, map_lat_upper]); %REGION =[west east south north];
%         pcolor_mode     = 0;
%         contf_lon       = lon_2D_topo(1,:);
%         contf_lat       = lat_2D_topo(:,1);
        if isempty(map_proj) map_proj = 'Miller Cylindrical';  end
        
        if isequal(map_proj, 'Miller Cylindrical')
            m_proj(map_proj,'longitudes',[map_lon_lower  map_lon_upper], 'latitudes', [map_lat_lower  map_lat_upper]);
        elseif isequal(map_proj, 'ortho')
            m_proj(map_proj,'lat', mean([map_lat_lower  map_lat_upper]),'long',mean([map_lon_lower, map_lon_upper]));
        else
            m_proj(map_proj,'longitudes',[map_lon_lower  map_lon_upper], 'latitudes', [map_lat_lower  map_lat_upper]);
        end
%         [contf_data, lon_2D_topo, lat_2D_topo]=m_tbase([map_lon_lower, map_lon_upper, map_lat_lower, map_lat_upper]); %REGION =[west east south north];
        [contf_data, lon_2D_topo, lat_2D_topo]=m_tbase([0, 360, -90, 90]); %REGION =[west east south north];
        
%         %/ 2D interpolation (nearest) to xx deg res 
%         new_res = 0.25;
%         lon_topo_new = [map_lon_lower:new_res:map_lon_upper];
%         lat_topo_new = [map_lat_lower:new_res:map_lat_upper];
%         [lon_2D_topo, lat_2D_topo] = meshgrid(lon_topo_new, lat_topo_new);
%         contf_data = interp2(lon_2D_topo, lat_2D_topo, contf_data, lon_2D_topo_new, lat_2D_topo_new, 'nearest');
        
        pcolor_mode     = 1;
        contf_data(contf_data < 0) = -9999;   %/ replace ocean basin with one value.
        contf_lon       = lon_2D_topo(1,:);
        contf_lat       = lat_2D_topo(:,1);
%         contf_levels    = [-250:250:6250];
%         colmap          = [[0 45 100]./255; m_colmap('gland',length(find(contf_levels > 0)))];
%         cbar_YTick      = [-4000:1000:7000];
%         cbar_YTickLabel = cbar_YTick./1000; %/ show in km.
        topo_intvl      = 225;
        contf_levels    = [-topo_intvl:topo_intvl:topo_intvl*26];
        contf_unit       = 'm';
        colmap          = my_colormap(length(contf_levels)-1,'topo_land_27lev');
        colmap(1,:)     = [0 45 100]./255;
        cbar_YTick      = contf_levels(2:4:end-1);
        cbar_YTickLabel = cbar_YTick;
        cbar_mode       = 1;
        cbar_interval   = 2;
        glb_plateau_mode = 0;
    end
    
    %====== cont ======%
    if ~isempty(cont_dataname) || ~isempty(cont_data)
        if isempty(cont_data) && isempty(cont_data_raw)
            flag_input = 0;
            if contains(select_field, '_clim')
                date = dataset.(cont_dataname).(fld_date_clim);
            else
                date = dataset.(cont_dataname).(fld_date);
            end

            if numel(num2str(compo_date(1))) == 8 && isequal(time_unit, 'monthly')
                ind  = findismember_loop(date, floor(compo_date./1e2));  %/ Then we will subset the month of the date in the monthly data
            else
                ind  = findismember_loop(date, compo_date);
            end

            if isempty(ind) 
                error('[quickplot_weather] Empty ind! Check compo_date!');
            elseif length(ind) ~= length(compo_date)
                warning('[quickplot_weather] Not all compo_date are available from the data. Filling missing data with NaN (useful when showing near-real time data).');
            end
    
            if isequal(sig_mode_cont, '_ttest')
                cont_data_events = dataset.(cont_dataname).(select_field)(:,:,ind);
                if isequal(contf_dataname, cont_dataname)  %/ Then just to outline contouf.
                    cont_data = mean(cont_data_events, 3, 'omitnan');
                    cont_data_raw = [];
                else
                    [cont_data, cont_data_raw, ~] = ttest_sig_fn(cont_data_events, alpha, 3, 0);
                end
            elseif ~isempty(sig_mode_cont)
                cont_data_events = dataset.(cont_dataname).(strcat(select_field, sig_mode_cont))(:,:,ind);
                cont_data     = mean(cont_data_events, 3, 'omitnan');
                cont_data_raw = mean(dataset.(cont_dataname).(select_field)(:,:,ind), 3, 'omitnan');
            else 
                cont_data_events = dataset.(cont_dataname).(select_field)(:,:,ind);
                [~, cont_data, ~] = ttest_sig_fn(cont_data_events, alpha, 3, 0);
                cont_data_raw = [];
            end
            
            fprintf('*** Mean cont data = %.2f ***\n', mean(cont_data, 'all', 'omitnan'));
            cont_lon       = dataset.(cont_dataname).lon;
            cont_lat       = dataset.(cont_dataname).lat;
        else
            flag_input = 1;
        end

        if isempty(cont_linewi) 
            if isequal(project_name, 'ITCC')   
                cont_linewi = 2;  
            else  
                cont_linewi = 4; 
            end
        end
        % cont_labelsize = fontsize*0.8;
        cont_labelsize = [];  %/ do not show contour label
        cont_label_col = 'k';
        skip_zero_cont = 1;
        cont_unit      = '';

        if contains(cont_dataname, {'prcp'}) || ismember(cont_dataname, {'ERA5_P'}) || isequal(cont_dataname, 'pr')
            if contains(select_field, 'monthly')
                CONV = 30;  %/ Convert into mm/month
            else
                CONV = 1;   %/ Default magnitude: mm/day
            end 

            if contains(select_field, {'_MJO'})
                if value_in_diff
                    cont_levels  = CONV*(-3.6:0.6:3.6)/2;
                else
                    cont_levels  = CONV*(-3.6:0.6:3.6);
                end
            else
                if value_in_diff
                    cont_levels  = CONV*(-3.6:0.6:3.6);
                else
                    cont_levels  = CONV*(0:2:24);
                end
            end
            if isequal(contf_dataname, cont_dataname)
                cont_colmap   = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);
            else
                cont_colmap   = create_div_colmap([42, 203, 238]./255, [238, 77, 42]./255, cont_levels);
            end

        elseif ismember(cont_dataname, {'MC1000to850'})
            if contains(select_field, 'anom')
                cont_levels = (-6:1:6)*2;

            elseif contains(select_field, 'daily_MJO')
                if value_in_diff
                    cont_levels = (-2.4:0.3:2.4)/3;
                else
                    % cont_levels = -3.6:0.6:3.6;
                    cont_levels = -4.2:0.7:4.2;
                end
            else
                if value_in_diff
                    cont_levels = -24:3:24;
                else
                    cont_levels = -70:10:70;
                end
            end
            cont_colmap   = create_div_colmap([173, 124, 0]./255, [0, 145, 92]./255, cont_levels);
            % cont_colmap   = create_div_colmap([204, 147, 14]./255, [14, 204, 131]./255, cont_levels);
            % colmap = nclCM('precip_diff_12lev', length(contf_levels)-1);

        elseif ismember(cont_dataname, {'Z850'})
            if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                cont_levels = -11:1:11;
                cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
%                 cont_colmap = brewermap(length(cont_levels), '*PuOr');              
%                 cont_colmap = create_div_colmap([33 121 255]./255, [255 102 51]./255, cont_levels);
                
            elseif ismember(select_field, {'pentad_clim_AC'})
                cont_levels    = -70:10:70;
                cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);

            else
                cont_levels    = 1420:10:1540;
%                 cont_colmap    = brewermap(length(cont_levels), '*RdYlBu');
                cont_colmap    = repmat([.1 .1 .1], length(cont_levels), 1);
            end
            cont_unit = 'm';
            
        elseif ismember(cont_dataname, {'Z500'})
            if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                cont_levels = -11:1:11;
                cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
%                 cont_colmap = brewermap(length(cont_levels), '*PuOr');              
%                 cont_colmap = create_div_colmap([33 121 255]./255, [255 102 51]./255, cont_levels);
                
            elseif ismember(select_field, {'pentad_clim_AC'})
                cont_levels    = -70:10:70;
                cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                
            else
                cont_levels    = 5400:50:6000;
                cont_colmap    = brewermap(length(cont_levels), '*RdYlBu');
            end
            cont_unit = 'm';

        elseif ismember({cont_dataname}, {'dMSE_dt_vi', 'MSE_Res'})
            cont_levels  = -55:5:55;             %/ in W m^-2 
            cont_levels(end) = [];
            skip_zero_cont   = 0;
            cont_colmap = brewermap(length(cont_levels), '*RdYlBu');
            cont_unit = 'W m^{-2}';
            
        elseif ismember(cont_dataname, {'psi200'})
            %/ NOTE:
            %/ In geostrophic balance (See Holton or notes in Notion):
            %/   High (low) stream function -> high (low) geopotential height in N.H.
            %/                              -> low (high) geopotential height in S.H.
            %/ https://www.cpc.ncep.noaa.gov/products/CDB/Tropics/figt22.shtml
            
            if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                cont_levels   = -10:1:10;   %/ x 10^6
                cont_colmap   = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                
            elseif ismember(select_field, {'pentad_clim_AC'})
                cont_levels   = -100:10:100;   %/ x 10^6
                cont_colmap   = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
            else
                cont_levels   = -100:10:100;   %/ x 10^6
                cont_colmap   = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
            end
            cont_multiply = 10^-6;
            cont_unit = '10^{-6} m^{2} s^{-1}';
        
        elseif ismember(cont_dataname, {'Z200'})
            if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                cont_levels   = -50:5:50;  
                cont_colmap   = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                
            elseif ismember(select_field, {'pentad_clim_AC'})
                cont_levels   = -200:20:200;   
                cont_colmap   = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
            else
                cont_levels   = 11550:50:13050;
                cont_colmap   = repmat([20 20 20]./255, length(cont_levels), 1);
            end
            cont_unit = 'm';

        elseif ismember(cont_dataname, {'slp'})       
            cont_levels    = 1000:2:1020;       %/ in hPa
            cont_colmap    = jet(length(cont_levels));
            cont_unit      = 'hPa';

        elseif contains(cont_dataname, {'SST'}) || ismember(cont_dataname, {'tos'})
            cont_unit    = sprintf('%sC', char(176));
            if contains(select_field, '_MJO')
                if value_in_diff
                    cont_levels  = (-0.3:0.05:0.3)/2;
                else
                    cont_levels  = -0.3:0.05:0.3;
                end
            elseif contains(select_field, {'monthly', 'daily'})
                if contains(select_field, 'anom')
                    if value_in_diff
                        cont_levels  = (-0.64:0.08:0.64)/2;
                    else
                        cont_levels  = -0.64:0.08:0.64;
                    end
                else
                    cont_add     = -273.15;   %/ K to deg C (only for raw data)
                    if value_in_diff
                        cont_levels  = -0.25:0.25:4.25;
                    else
                        cont_levels  = (22:1:33);
                    end
                end
            end
            % cont_colmap = create_div_colmap([255 255 255]./255, [255 255 255]./255, cont_levels);
            cont_colmap = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);

        elseif contains(cont_dataname, {'T2m'}) || isequal(cont_dataname, 'tas')
            cont_add       = -273.15;   %/ K to deg C
            cont_unit      = sprintf('%sC', char(176));
            if value_in_diff
                cont_levels  = (-6:1:6);
            else
                cont_levels    = (10:2:36);
            end
            if isequal(contf_dataname, cont_dataname)
                cont_colmap   = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);
            else
                cont_colmap    = brewermap(length(cont_levels), '*spectral');
            end
            
        elseif any(ismember(cont_dataname, {'uWS_LL','uWS_UL','vWS_LL','vWS_UL'}))
            cont_unit      = 'm/s';
            cont_levels    = -50:5:50;   %/ in m/s
            cont_colmap    = jet(length(cont_levels));
            
        elseif contains(cont_dataname, {'omega500', 'W500', 'wap', 'W5'}) 
            cont_multiply  = 24*3600/100;  %/ Pa/s -> hPa/day
            cont_unit      = 'hPa/day';
            if contains(select_field, 'anom')
                % if value_in_diff
                %     cont_levels = (-6:1:6)*2;
                % else
                cont_levels = (-6:1:6)*2;
                % end
            else
                if value_in_diff
                    cont_levels = -28:4:28;
                else
                    cont_levels = -70:10:70;
                end
            end
            if isequal(contf_dataname, cont_dataname)
                cont_colmap   = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);
            else
                cont_colmap   = nclCM('sunshine_diff_12lev', length(cont_levels));
                % cont_colmap   = create_div_colmap([58 187 0]./255, [187 0  187]./255, cont_levels);
            end

        elseif ismember(cont_dataname, {'pw'})     %/ pw == tcwv
            cont_unit      = 'mm';
            cont_levels    = 10:10:60;               %/ in mm (instantaneous)
            cont_colmap    = cmocean('algae', length(cont_levels));

        elseif ismember(cont_dataname, {'TSR', 'SSR'})
            cont_unit      = 'W m^{-2}';
            if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                cont_levels  = -40:4:40;               %/ in W m^-2
                cont_colmap   = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                
            elseif ismember(select_field, {'pentad_clim_AC'})
                cont_levels  = -100:20:100;               %/ in W m^-2 
                cont_colmap   = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                
            else
                cont_levels    = 190:10:330;               %/ in W m^-2 
%                 contf_levels   = [190:10:360];             %/ in W m^-2 
                cont_colmap    = brewermap(length(cont_levels), '*RdYlBu');
            end
            
        elseif contains(cont_dataname, {'OLR'}) || ismember(cont_dataname, {'rlut'})
            cont_unit      = 'W m^{-2}';
            if contains(select_field, {'_MS'})
                cont_levels  = 190:10:310;             %/ in W m^-2 
                cont_colmap  = my_colormap(length(cont_levels), 'radar_12lev');

            elseif contains(select_field, {'_AC'})
                cont_levels   = -60:10:60;                %/ in W m^-2 
                cont_colmap   = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);

            elseif contains(select_field, {'_MJO'})
                if value_in_diff
                    % cont_levels   = [-9:1:9];                %/ in W m^-2 
                    cont_levels   = -14:2:14;                %/ in W m^-2 
                else
                    % cont_levels   = [-35:5:35];                %/ in W m^-2 
                    cont_levels   = -40:4:40;                %/ in W m^-2 
                end
                cont_colmap   = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);
                % if isequal(contf_dataname, cont_dataname)
                %     cont_colmap   = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);
                % else
                %     cont_colmap   = create_div_colmap([42, 203, 238]./255, [238, 77, 42]./255, cont_levels);
                % end
            elseif contains(select_field, {'_CISO', '_TE', '_IAID'})
                cont_levels   = -18:3:18;               %/ in W m^-2 
                cont_colmap   = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);

            elseif contains(select_field, {'_anom'})
                if value_in_diff
                    cont_levels   = -6:1:6;                %/ in W m^-2 
                else
                    cont_levels   = -12:2:12;                %/ in W m^-2 
                end
                cont_colmap   = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);
            else
                cont_levels   = 190:10:310;             %/ in W m^-2 
                cont_colmap   = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);
            end

        elseif contains(cont_dataname, {'q1000', 'q850', 'hus'}) 
            cont_unit        = 'g/kg';
            cont_multiply    = 1000;   %/ kg/kg to g/kg

            if contains(select_field, {'_MJO'})
                if value_in_diff
                    cont_levels  = (-0.36:0.06:0.36);
                else
                    cont_levels  = (-2.4:0.4:2.4)./2;
                end
            else
                if value_in_diff
                    cont_levels  = -4.5:0.5:4.5;
                else
                   if contains(cont_dataname, {'1000'})
                        cont_levels  = 1.25:1.25:21.25;
                   else
                        cont_levels  = 0:1:16;
                   end
                end
            end
            cont_colmap = repmat([0 0 0], length(cont_levels), 1);

        elseif contains(cont_dataname, {'S500'}) || contains(cont_dataname, {'S5'})  %/ dry static stability (assuming in K/Pa)
            cont_multiply = 10^4;
            cont_unit     = 'K (100 hPa)^{-1}';
            if value_in_diff
                cont_levels   = (0.25:0.05:1.35);
            else
                cont_levels   = (2.5:0.5:8.5); 
            end
            cont_colmap = repmat([0 0 0], length(cont_levels), 1);

        elseif ismember(cont_dataname, {'dE_TOA'})    
            cont_unit      = 'W m^{-2}';
            if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                cont_levels   = -40:4:40;               %/ in W m^-2
                cont_colmap   = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                
            elseif ismember(select_field, {'pentad_clim_AC'})
                cont_levels   = -100:20:100;               %/ in W m^-2 
                cont_colmap   = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);

            else
                cont_levels    = -200:20:200;               %/ in W m^-2 
                cont_colmap    = brewermap(length(cont_levels), '*RdGy');
            end
            
        elseif ismember(cont_dataname, {'D_vm850to1000', 'D_vm850to1000_sm'})
            cont_multiply  = 1e6;
            cont_unit      = '10^{-6} s^{-1}';
            cont_levels  = -10:2:10;             %/ in 1e-6 s^-1 
            cont_colmap  = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);

        elseif ismember({cont_dataname}, {'MSE_vi'})
            cont_multiply  = 1e-9;
            cont_unit      = '10^{9} J m^{-2}';
            cont_levels = linspace(2.925, 3.2, 12);             %/ in 10^9 J m^-2 
            cont_colmap = jet(length(cont_levels));
        else
            error('The cont_levels of the request data has not been defined!')
        end

        %/ Avoid performing unit conversion when cont_data is directly input
        if flag_input == 0
            cont_data     = cont_data * cont_multiply + cont_add;
            cont_data_raw = cont_data_raw * cont_multiply + cont_add;
            cont_data_events = cont_data_events * cont_multiply + cont_add;
        end
    end
    
    %====== vector ======%
    if (~isempty(U_dataname) && ~isempty(V_dataname)) || (~isempty(Udata) && ~isempty(Vdata)) 
        if isempty(Udata) && isempty(Vdata) 
            if contains(select_field, '_clim')
                date = dataset.(U_dataname).(fld_date_clim);
            else
                date = dataset.(U_dataname).(fld_date);
            end

            if numel(num2str(compo_date(1))) == 8 && isequal(time_unit, 'monthly')
                ind  = findismember_loop(date, floor(compo_date./1e2));  %/ Then we will subset the month of the date in the monthly data
            else
                ind  = findismember_loop(date, compo_date);
            end

            if isempty(ind) 
                error('[quickplot_weather] Empty ind! Check compo_date!');
            elseif length(ind) ~= length(compo_date)
                warning('[quickplot_weather] Not all compo_date are available from the data. Filling missing data with NaN (useful when showing near-real time data).');
            end
    
            % if isequal(sig_mode_cont, '_ttest')
            %     cont_data_events = dataset.(cont_dataname).(select_field)(:,:,ind);
            %     if isequal(contf_dataname, cont_dataname)  %/ Then just to outline contouf.
            %         cont_data = mean(cont_data_events, 3, 'omitnan');
            %         cont_data_raw = [];
            %     else
            %         [cont_data, cont_data_raw, ~] = ttest_sig_fn(cont_data_events, alpha, 3, 0);
            %     end
            % elseif ~isempty(sig_mode_cont)
            %     cont_data_events = dataset.(cont_dataname).(strcat(select_field, sig_mode_cont))(:,:,ind);
            %     cont_data     = mean(cont_data_events, 3, 'omitnan');
            %     cont_data_raw = mean(dataset.(cont_dataname).(select_field)(:,:,ind), 3, 'omitnan');
            % else 
            %     cont_data_events = dataset.(cont_dataname).(select_field)(:,:,ind);
            %     [~, cont_data, ~] = ttest_sig_fn(cont_data_events, alpha, 3, 0);
            %     cont_data_raw = [];
            % end

            if isequal(sig_mode_vec, '_ttest')
                Udata_events = dataset.(U_dataname).(select_field)(:,:,ind);
                Vdata_events = dataset.(V_dataname).(select_field)(:,:,ind);
                [U2data, Udata, ~] = ttest_sig_fn(Udata_events, alpha, 3, 0);
                [V2data, Vdata, ~] = ttest_sig_fn(Vdata_events, alpha, 3, 0);
            elseif ~isempty(sig_mode_vec)
                Udata_events = dataset.(U_dataname).(strcat(select_field, sig_mode_vec))(:,:,ind);
                Vdata_events = dataset.(V_dataname).(strcat(select_field, sig_mode_vec))(:,:,ind);
                Udata = mean(Udata_events, 3, 'omitnan');
                Vdata = mean(Vdata_events, 3, 'omitnan');
                U2data = [];
                V2data = [];
            else
                Udata_events = dataset.(U_dataname).(select_field)(:,:,ind);
                Vdata_events = dataset.(V_dataname).(select_field)(:,:,ind);
                [~, Udata, ~] = ttest_sig_fn(Udata_events, alpha, 3, 0);
                [~, Vdata, ~] = ttest_sig_fn(Vdata_events, alpha, 3, 0);
                U2data = [];
                V2data = [];
            end
            uv_lon = dataset.(U_dataname).lon;
            uv_lat = dataset.(U_dataname).lat;
        end

        vector_color        = 'k'; 
        vector_edgecolor    = 'none';
        vector2_color       = 'r';    %/ default
        vector2_edgecolor   = 'none';

        if isequal(project_name, 'ITCC')
            if contains(U_dataname, {'IVT'})
                vector_color = 'k'; %blue 
                vector_edgecolor = 'w';
                if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                    vec_step_lon = 6;
                    vec_step_lat = 3;
                    vecscale     = 50;     % the smaller value the bigger vector. for winds
                    vecscale2    = 1;       % control shaft length.
                    shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                    headlength   = shaftwidth*4; % control length of the head, the larger value the longer  
                    vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                    vec_mag_ref  = 50;

                elseif ismember(select_field, {'pentad_clim_AC'})
                    vec_step_lon = 6;
                    vec_step_lat = 3;
                    vecscale     = 500;     % the smaller value the bigger vector. for winds
                    vecscale2    = 1;       % control shaft length.
                    shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                    headlength   = shaftwidth*4; % control length of the head, the larger value the longer  
                    vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                    vec_mag_ref = 500;
                else           
                    vector_levels = [0:50:500];    %/ in kg/m/s
                    vecscale2     = 1;           % control shaft length.
                    vector_color      = jet(length(vector_levels));
                    vector_edgecolor  = 'k';
                    vec_step_lon = 5;
                    vec_step_lat = 3;
                    vecscale     = 500;    % the smaller value the bigger vector. for winds
                    shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                    headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                    vec_mag_ref = 500;
                end

                %/ vector reference
                vec_lbs = strcat(num2str(vec_mag_ref), {' kg/m/s'});

            elseif contains(U_dataname, {'850'})
                vector_color = 'k';
                vector_edgecolor = 'w';
                if contains(select_field, {'_TE', '_CISO', '_AC'})
                    if map_lon_upper - map_lon_lower > 350  %/ assume global
                        vec_step_lon = 10;
                        vec_step_lat = 5;
                        vecscale     = 25;     % the smaller value the bigger vector. for winds
                        vecscale2    = 3;       % control shaft length.
                        shaftwidth   = 3;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    else
                        vec_step_lon = 4;
                        vec_step_lat = 2;
                        vecscale     = 5;     % the smaller value the bigger vector. for winds
                        vecscale2    = 1;       % control shaft length.
                        shaftwidth   = 3;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    end
                    vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                    vec_mag_ref  = 10;

                else
                    if map_lon_upper - map_lon_lower > 350  %/ assume global
                        vec_step_lon = 10;
                        vec_step_lat = 5;
                    else
                        vec_step_lon      = 6;           %/ then, no need to have high density vectors.
                        vec_step_lat      = 3;
                    end
                    vector_color      = 'k';
                    vector_edgecolor  = 'w';         
                    vecscale2         = 1;  
                    vecscale          = 20;             % the smaller value the bigger vector. for winds
                    shaftwidth        = 2.8;            % control width of the shaft, the larger value the thicker
                    headlength        = shaftwidth*3.5;   % control length of the head, the larger value the longer  
                    vec_ref_lat_shift = 0.08;           % the smaller, the less the ref vector shift to the south
                    vec_mag_ref       = 20;
                end

                %/ vector reference
                vec_lbs = strcat(num2str(vec_mag_ref), {' m/s'});
                
            elseif contains(U_dataname, {'10m'})
                vector_color = 'k';
                vector_edgecolor = 'w';
                if contains(select_field, {'_TE', '_CISO', '_AC'})

                    if map_lon_upper - map_lon_lower > 350  %/ assume global
                        vec_step_lon = 10;
                        vec_step_lat = 5;
                        vecscale     = 25;     % the smaller value the bigger vector. for winds
                        vecscale2    = 3;       % control shaft length.
                        shaftwidth   = 3;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    else
                        vec_step_lon = 4;
                        vec_step_lat = 2;
                        vecscale     = 5;     % the smaller value the bigger vector. for winds
                        vecscale2    = 1;       % control shaft length.
                        shaftwidth   = 3;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    end
                    vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                    vec_mag_ref  = 10;

                else
                    if map_lon_upper - map_lon_lower > 350  %/ assume global
                        vec_step_lon = 10;
                        vec_step_lat = 5;
                    else
                        vec_step_lon  = 5;           %/ then, no need to have high density vectors.
                        vec_step_lat  = 2;
                    end
                    vector_color      = 'k';
                    vector_edgecolor  = 'w';         
                    vecscale2         = 1;  
                    vecscale          = 20;             % the smaller value the bigger vector. for winds
                    shaftwidth        = 2.8;            % control width of the shaft, the larger value the thicker
                    headlength        = shaftwidth*3.5;   % control length of the head, the larger value the longer  
                    vec_ref_lat_shift = 0.08;           % the smaller, the less the ref vector shift to the south
                    vec_mag_ref       = 20;
                end

                %/ vector reference
                vec_lbs = strcat(num2str(vec_mag_ref), {' m/s'});
                
            elseif contains(U_dataname, {'200', '200chi'})
                vector_color = [147 255 50]./255; 
                vector_edgecolor = 'k';

                if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})
                    if map_lon_upper - map_lon_lower > 350  %/ assume global
                        vec_step_lon = 10;
                        vec_step_lat = 5;
                        vecscale     = 10;     % the smaller value the bigger vector. for winds
                        vecscale2    = 5;       % control shaft length.
                        shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    else
                        vec_step_lon = 6;
                        vec_step_lat = 3;
                        vecscale     = 10;     % the smaller value the bigger vector. for winds
                        vecscale2    = 5;       % control shaft length.
                        shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    end
                    vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                    vec_mag_ref  = 2;

                elseif ismember(select_field, {'pentad_clim_AC'})
                    if map_lon_upper - map_lon_lower > 350  %/ assume global
                        vec_step_lon = 10;
                        vec_step_lat = 5;
                        vecscale     = 25;     % the smaller value the bigger vector. for winds
                        vecscale2    = 4;       % control shaft length.
                        shaftwidth   = 3;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    else
                        vec_step_lon = 6;
                        vec_step_lat = 3;
                        vecscale     = 25;     % the smaller value the bigger vector. for winds
                        vecscale2    = 4;       % control shaft length.
                        shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    end
                    vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                    vec_mag_ref  = 10;
                else
                    if map_lon_upper - map_lon_lower > 350  %/ assume global
                        vec_step_lon = 10;
                        vec_step_lat = 5;
                    else
                        vec_step_lon      = 6;           %/ then, no need to have high density vectors.
                        vec_step_lat      = 4;
                    end
%                     vector_color      = jet(length(vector_levels));
                    vector_color      = [102 0 204]./255;
                    vector_edgecolor  = 'w';         
                    vecscale2         = 1;  
                    vecscale          = 30;             % the smaller value the bigger vector. for winds
                    shaftwidth        = 2.8;            % control width of the shaft, the larger value the thicker
                    headlength        = shaftwidth*3.5;   % control length of the head, the larger value the longer  
                    vec_ref_lat_shift = 0.08;           % the smaller, the less the ref vector shift to the south
                    vec_mag_ref       = 10;
                end

                %/ vector reference
                vec_lbs = strcat(num2str(vec_mag_ref), {' m/s'});

            elseif contains(U_dataname, {'WAF200'})  %/ order of magnitude 10^2 m^2 s^-1
                vector_color = [147 255 50]./255; 
                vector_edgecolor = 'k';
                if ismember(select_field, {'pentad_clim_TE', 'pentad_clim_CISO'})

                    if map_lon_upper - map_lon_lower > 350  %/ assume global
                        vec_step_lon = 10;
                        vec_step_lat = 5;
                        vecscale     = 50;     % the smaller value the bigger vector. for winds
                        vecscale2    = 2;       % control shaft length.
                        shaftwidth   = 3;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    else
                        vec_step_lon = 6;
                        vec_step_lat = 3;
                        vecscale     = 50;     % the smaller value the bigger vector. for winds
                        vecscale2    = 2;       % control shaft length.
                        shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    end
                    vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                    vec_mag_ref  = 2;

                elseif ismember(select_field, {'pentad_clim_AC'})

                    if map_lon_upper - map_lon_lower > 350  %/ assume global
                        vec_step_lon = 10;
                        vec_step_lat = 5;
                        vecscale     = 600;     % the smaller value the bigger vector. for winds
                        vecscale2    = 2;       % control shaft length.
                        shaftwidth   = 3;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    else
                        vec_step_lon = 6;
                        vec_step_lat = 3;
                        vecscale     = 600;     % the smaller value the bigger vector. for winds
                        vecscale2    = 2;       % control shaft length.
                        shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                        headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                    end
                    vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                    vec_mag_ref  = 10;
                else
                    vec_step_lon = 10;
                    vec_step_lat = 6;
                    vecscale     = 100;     % the smaller value the bigger vector. for winds
                    vecscale2    = 3;       % control shaft length.
                    shaftwidth   = 3;    % control width of the shaft, the larger value the thicker
                    headlength   = shaftwidth*4; % control length of the head, the larger value the longer  
                    vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                    vec_mag_ref  = 10;
                end

                %/ vector reference
                vec_lbs = strcat(num2str(vec_mag_ref), {' m**2 s**-1'});

            elseif contains(U_dataname, {'WS_LL', 'WS_UL'})
    %             vector_color = [.3 .3 .3]; % grey
                vector_edgecolor = 'w';
                vec_step_lon = 2;
                vec_step_lat = 2;
                vecscale     = 30;     % the smaller value the bigger vector. for winds
                vecscale2    = 1;       % control shaft length.
                shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                headlength   = shaftwidth*4; % control length of the head, the larger value the longer  
                vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south

                %/ vector reference
                vec_mag_ref = 20;
                vec_lbs = strcat(num2str(vec_mag_ref), {' m/s'});
            end    
        else
            if contains(U_dataname, {'IVT'})
                vector_levels = 0:50:500;    %/ in kg/m/s
                vecscale2     = 1;           % control shaft length.
                vector_color      = jet(length(vector_levels));
                vector_edgecolor  = 'k';
                vec_step_lon = 5;
                vec_step_lat = 3;
                vecscale     = 500;    % the smaller value the bigger vector. for winds
                shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                headlength   = shaftwidth*4; % control length of the head, the larger value the longer
                vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                vec_mag_ref = 500;

                %/ vector reference
                vec_lbs = strcat(num2str(vec_mag_ref), {' kg/m/s'});

            elseif contains(U_dataname, {'850'})
                vector_levels = [];    
                vector_color = [.25 .25 .25];
                vector_edgecolor = 'none';
                % vector_levels = [0:1:10];    %/ in m/s
                % vector_color      = jet(length(vector_levels));
                % vector_edgecolor  = 'k';

                vecscale2     = 1;           
                vec_step_lon  = 5;           %/ then, no need to have high density vectors.
                vec_step_lat  = 4;

                if contains(select_field, '_MJO') 
                    if value_in_diff
                        vecscale  = 1;             % the smaller value the bigger vector. for winds
                    else
                        vecscale  = 3;             % the smaller value the bigger vector. for winds
                    end
                elseif contains(select_field, 'monthly') && contains(select_field, '_anom')
                    if value_in_diff
                        vecscale  = 1;             % the smaller value the bigger vector. for winds
                    else
                        vecscale  = 2;             % the smaller value the bigger vector. for winds
                    end
                else
                    vecscale  = 10;
                end
                vecscale
                
                %/ Vector shape depends on the map size
                if map_lon_upper-map_lon_lower >= 250
                    shaftwidth        = 1.5;    % control width of the shaft, the larger value the thicker
                    headlength        = shaftwidth*4;   % control length of the head, the larger value the longer  
                else
                    shaftwidth        = 1.5;    % control width of the shaft, the larger value the thicker
                    headlength        = shaftwidth*6;   % control length of the head, the larger value the longer  
                end
                % shaftwidth        = 2;            % control width of the shaft, the larger value the thicker
                % headlength        = shaftwidth*4;   % control length of the head, the larger value the longer  
                vec_ref_lat_shift = 0.08;           % the smaller, the less the ref vector shift to the south
                vec_mag_ref       = 10;

                %/ vector reference
                vec_lbs = strcat(num2str(vec_mag_ref), {' m/s'});
            
            elseif contains(U_dataname, {'200', '200chi'})
%                 vector_color = [147 255 50]./255; 
%                 vector_edgecolor = 'k';
                vector_levels = [0:1:10];    %/ in m/s
                vecscale2     = 1;           % control shaft length.
                vector_color      = jet(length(vector_levels));
                vector_edgecolor  = 'k';
                vec_step_lon      = 5;
                vec_step_lat      = 3;
                vecscale          = 10;             % the smaller value the bigger vector. for winds
                shaftwidth        = 3.5;            % control width of the shaft, the larger value the thicker
                headlength        = shaftwidth*4;   % control length of the head, the larger value the longer  
                vec_ref_lat_shift = 0.08;           % the smaller, the less the ref vector shift to the south
                vec_mag_ref       = 10;

                %/ vector reference
                vec_lbs = strcat(num2str(vec_mag_ref), {' m/s'});

            elseif contains(U_dataname, {'WAF200'})  %/ order of magnitude 10^2 m^2 s^-1
                vector_color = [147 255 50]./255; 
                vector_edgecolor = 'k';
                vec_step_lon = 10;
                vec_step_lat = 6;
                vecscale     = 100;     % the smaller value the bigger vector. for winds
                vecscale2    = 3;       % control shaft length.
                shaftwidth   = 3;    % control width of the shaft, the larger value the thicker
                headlength   = shaftwidth*4; % control length of the head, the larger value the longer  
                vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                vec_mag_ref  = 10;

                %/ vector reference
                vec_lbs = strcat(num2str(vec_mag_ref), {' m**2 s**-1'});

            elseif contains(U_dataname, {'WS_LL', 'WS_UL'})
                vector_color = [.3 .3 .3]; % grey
                vector_edgecolor = 'w';
                vec_step_lon = 2;
                vec_step_lat = 2;
                vecscale     = 30;     % the smaller value the bigger vector. for winds
                vecscale2    = 1;       % control shaft length.
                shaftwidth   = 3.5;    % control width of the shaft, the larger value the thicker
                headlength   = shaftwidth*4; % control length of the head, the larger value the longer  
                vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south

                %/ vector reference
                vec_mag_ref = 20;
                vec_lbs = strcat(num2str(vec_mag_ref), {' m/s'});

            else
                vector_levels = [];    
                % vector_color = 'k';
                vector_color = [.25 .25 .25];
                % vector_edgecolor = 'w';
                vector_edgecolor = 'none';
                % vector_levels = [0:1:10];    %/ in m/s
                % vector_color      = jet(length(vector_levels));
                % vector_edgecolor  = 'k';

                vecscale2     = 1;           
                vec_step_lon  = 5;           %/ then, no need to have high density vectors.
                vec_step_lat  = 4;

                if contains(select_field, '_MJO')
                    if value_in_diff && contains(select_field, '_anom')
                        vecscale  = 1;             % the smaller value the bigger vector. for winds
                    else
                        vecscale  = 3;             % the smaller value the bigger vector. for winds
                    end
                else
                    vecscale  = 10;
                end

                shaftwidth        = 1.5;    % control width of the shaft, the larger value the thicker
                headlength        = shaftwidth*6;   % control length of the head, the larger value the longer  
                % shaftwidth        = 2;            % control width of the shaft, the larger value the thicker
                % headlength        = shaftwidth*4;   % control length of the head, the larger value the longer  
                vec_ref_lat_shift = 0.08;           % the smaller, the less the ref vector shift to the south
                vec_mag_ref       = 10;

                %/ vector reference
                vec_lbs = strcat(num2str(vec_mag_ref), {' m/s'});
            end
        end
        
        %/ Auto adjust vec_step_lon and vec_step_lat based on hori resolution
        if auto_vec_step
            res_lon = abs(diff(uv_lon(1:2)));
            res_lat = abs(diff(uv_lat(1:2)));
            vec_step_lon = round(1/res_lon*4); %/ round up to integer
            vec_step_lat = round(1/res_lat*4); %/ round up to integer
        end
    end

    %====== Additional data processing (project-oriented) ======%    
    %/ Show the centroid of the MJO
    if show_MJO_centroid
        if isfile(MJO_centr_filename) && recompute_MJO_centr == 0
            %/ Load centroid
            fprintf('*** Loading MJO_centr from %s ***\n', MJO_centr_filename)
            load(MJO_centr_filename, 'MJO_centr');

        elseif (contains(contf_dataname, 'OLR') || isequal(contf_dataname, 'rlut')) && contains(select_field, 'daily_MJO')
            %/ 1. Subset data on tropical Indo-Pacific Ocean
            dom_name = 'TIPO';       
            [lon_range, lat_range, ~, ~] = hovmoller_box_reg('project_name', project_name, 'dom_name', dom_name, 'zm_or_mm', 2, 'select_field', select_field); 
            
            ind_lon = contf_lon >= lon_range(1) & contf_lon <= lon_range(2);
            ind_lat = contf_lat >= lat_range(1) & contf_lat <= lat_range(2);
            contf_lon_subset = contf_lon(ind_lon);
            contf_lat_subset = contf_lat(ind_lat);
        
            [contf_lon_subset_2D, contf_lat_subset_2D] = meshgrid(contf_lon_subset, contf_lat_subset);
            contf_lon_subset_2D      = contf_lon_subset_2D'; 
            contf_lat_subset_2D      = contf_lat_subset_2D'; 
            contf_lon_subset_2Dto1D  = reshape(contf_lon_subset_2D, [], 1);
            contf_lat_subset_2Dto1D  = reshape(contf_lat_subset_2D, [], 1);
            contf_data_2Dto1D = reshape(contf_data(ind_lon,ind_lat), [], 1);

            %/ 2. Get the centroid by the weighted average (Recommended)
            wgt = contf_data_2Dto1D;
            wgt(wgt > -15 | isnan(wgt)) = 0;  %/ No weighting for OLR > -15 W m^-2
            wgt = -wgt;
            wgt = (wgt-min(wgt))./(max(wgt)-min(wgt)); %/normalization

            MJO_centr = [sum(contf_lon_subset_2Dto1D.*wgt)./sum(wgt), sum(contf_lat_subset_2Dto1D.*wgt)./sum(wgt)];
            
            %/ Save centroid
            fprintf('*** Saving MJO_centr to %s ***\n', MJO_centr_filename)
            save(MJO_centr_filename, 'MJO_centr');
        else
            error('Set contf_dataname to OLR or rlut and select_field to daily_MJO (or the like) in order to compute MJO_centr!');
        end
        fprintf('*** MJO_centr: %.2f, %.2f ***\n', MJO_centr(1), MJO_centr(2));
    end

    %/ Rossby-Kelvin ratio (based on 850-hPa winds)
    if compute_RK_ratio
        %/===================================================================================
        %/ Dissect two tropical lon bands following longitudinal position of the MJO centroid
        %/ NOTE: The planetary wave (Rossby) to the west of the convection is
        %/       only one-third of the size of the easterly wave to the east,
        %/       because of the lower phase speed (Gill 1980)
        %/===================================================================================
        if isempty(U_dataname)
            compute_RK_ratio = 0;
            warning('The Udata is not given! Auto setting compute_RK_ratio == 0...')
        else
            if show_MJO_centroid == 0
               error('Set show_MJO_centroid == 1 (or load from ''MJO_centr_filename'' if contf_dataname is not OLR or rlut'); 
            end
    
            lon_size_R   = 40;
            lon_size_K   = lon_size_R*3;  %/ See the note above
            lon_range_R  = [MJO_centr(1)-lon_size_R, MJO_centr(1)           ];
            lon_range_K  = [MJO_centr(1),            MJO_centr(1)+lon_size_K];
            lat_range_RK = [-15 15]; %/ See Fig. 4 in Wang et al (2018)
            ind_lat_RK   = uv_lat >= lat_range_RK(1) & uv_lat <= lat_range_RK(2);
            ind_lon_R    = uv_lon >= lon_range_R(1)  & uv_lon <= lon_range_R(2);
            ind_lon_K    = uv_lon >= lon_range_K(1)  & uv_lon <= lon_range_K(2);
    
            %/ If to show the box region of where Rossby westerler/Kelvin easterly are averaged
            if show_RK_box 
                if show_MJO_centroid == 0
                   error('Set show_MJO_centroid == 1 (or load from ''MJO_centr_filename'' if contf_dataname is not OLR or rlut');  
                end
                bndry_data_R = get_box_bndry(lon_range_R, lat_range_RK);
                bndry_data_K = get_box_bndry(lon_range_K, lat_range_RK);
                bndry_data   = [bndry_data, {bndry_data_R}, {bndry_data_K}];
            end
            
            if isempty(Rossby_westerly) && isempty(Kelvin_easterly) && isempty(RK_ratio)
                %/ [IMPORTANT] Separate westerly from easterly before taking an average (Not recommended)
                % westerly = Udata; 
                % easterly = Udata; 
                % westerly(westerly <= 0) = nan;  
                % easterly(easterly >= 0) = nan;
                % Rossby_westerly = mean(westerly(ind_lon_R, ind_lat_RK), 'all', 'omitnan');
                % Kelvin_easterly = mean(easterly(ind_lon_K, ind_lat_RK), 'all', 'omitnan');
    
                Rossby_westerly = mean(Udata(ind_lon_R, ind_lat_RK), 'all', 'omitnan');
                Kelvin_easterly = mean(Udata(ind_lon_K, ind_lat_RK), 'all', 'omitnan') * -1;  %/ Convert the sign; do not take absolute value
                if Rossby_westerly < 0
                    warning('The regional mean Rossby_westerly is -ve (i.e., easterly)! Check if your wind data is really the one coupled with MJO!')
                end
                if Kelvin_easterly < 0
                    warning('The regional mean Kelvin_easterly is -ve (i.e., westerly)! Check if your wind data is really the one coupled with MJO!')
                end
                RK_ratio        = Rossby_westerly./Kelvin_easterly;  
                RK_diff_or_not  = 0; 
    
            elseif isempty(RK_diff_or_not)
                error('Set ''RK_diff_or_not'' to indicate whether ''Rossby_westerly'', ''Kelvin_easterly'', and ''RK_ratio'' are differences or not.');
            end
    
            if RK_diff_or_not
                str_RK_diff_or_not = '\Delta';
            else
                str_RK_diff_or_not = '';
            end

            if ~isempty(Rossby_westerly_sig) && ~isnan(Rossby_westerly_sig)
                str_RW_sig = {' (*)'};
            else
                str_RW_sig = {''};
            end
            if ~isempty(Kelvin_easterly_sig) && ~isnan(Kelvin_easterly_sig)
                str_KE_sig = {' (*)'};
            else
                str_KE_sig = {''};
            end
            if ~isempty(RK_ratio_sig) && ~isnan(RK_ratio_sig)
                str_RK_sig = {' (*)'};
            else
                str_RK_sig = {''};
            end
    
            if ~isempty(Rossby_westerly_error) && ~isempty(Kelvin_easterly_error) && ~isempty(RK_ratio_error)
                %/ Check if it is a two-element error (e.g., 25th to 75th percentile), otherwise assume SD
                if numel(Rossby_westerly_error) == 2 
                    str_RW_error = sprintf('[%.2f, %.2f]', Rossby_westerly_error(1), Rossby_westerly_error(2));
                else  
                    str_RW_error = sprintf('%s%.2f', char(177), Rossby_westerly_error);
                end
                if numel(Kelvin_easterly_error) == 2  
                    str_KE_error = sprintf('[%.2f, %.2f]', Kelvin_easterly_error(1), Kelvin_easterly_error(2));
                else 
                    str_KE_error = sprintf('%s%.2f', char(177), Kelvin_easterly_error);
                end
                if numel(RK_ratio_error) == 2  
                    str_ratio_error = sprintf('[%.2f, %.2f]', RK_ratio_error(1), RK_ratio_error(2));
                else 
                    str_ratio_error = sprintf('%s%.2f', char(177), RK_ratio_error);
                end
            else
                str_RW_error = ''; str_KE_error = ''; str_ratio_error = '';
            end
            str_RK_ratio = sprintf('%sRossby westerly: %.2f%s%s m/s, %sKelvin easterly: %.2f%s%s m/s, %sRK ratio: %.2f%s%s',...
                                    str_RK_diff_or_not, Rossby_westerly, str_RW_sig{:}, str_RW_error, str_RK_diff_or_not, Kelvin_easterly, str_KE_sig{:}, str_KE_error, str_RK_diff_or_not, RK_ratio, str_RK_sig{:}, str_ratio_error);
            % str_RK_ratio = sprintf('%sRossby westerly: %.2f m/s, %sKelvin easterly: %.2f m/s, %sRK ratio: %.2f', str_RK_diff_or_not, Rossby_westerly, str_RK_diff_or_not, Kelvin_easterly, str_RK_diff_or_not, RK_ratio);
            % disp(str_RK_ratio);
        end
    end

    %/ Area-weighted mean of contf_data shown in the domain
    if compute_contf_AWM   
        if ~isempty(compute_AWM_lon_extent) && ~isempty(compute_AWM_lat_extent)  
            AWM_lon_extent = compute_AWM_lon_extent;
            AWM_lat_extent = compute_AWM_lat_extent;
            AWM_bndry = get_box_bndry(AWM_lon_extent, AWM_lat_extent);

            if show_AMW_reg
                bndry_data = {bndry_data, AWM_bndry};  %/ Outline the region where the AWM is computed.

                if isempty(AMW_reg_color)
                    AMW_reg_color = [0 0 0];
                end
                color = [color; AMW_reg_color];
            end
        else
            AWM_lon_extent = [map_lon_lower, map_lon_upper];
            AWM_lat_extent = [map_lat_lower, map_lat_upper];
        end

        %/ Compute it directly if the AWM of contf_data is not given
        if isempty(contf_AWM)
            mean_or_sum = 'mean'; %/ default
            [contf_AWM, ~, ~, ~, str_extent] = compute_area_wgted_meansum('data', contf_data, 'lon', contf_lon, 'lat', contf_lat,...
                                                           'lon_extent', AWM_lon_extent, 'lat_extent', AWM_lat_extent, 'mean_or_sum', mean_or_sum);
        end
        % str_extent = sprintf(' (%d-%dE, %d-%dN)', AWM_lon_extent(1), AWM_lon_extent(2), AWM_lat_extent(1), AWM_lat_extent(2));

        if ~isempty(contf_AWM_error)
            if numel(contf_AWM_error) == 2  
                str_AWM_error = sprintf('[%.2f, %.2f]', contf_AWM_error(1), contf_AWM_error(2));
            else 
                str_AWM_error = sprintf('%s%.2f', char(177), contf_AWM_error);
            end
        else
            str_AWM_error = ''; 
        end
        
        if ~isempty(str_contf_AWM)
            str_AWM = sprintf('AWM%s: %s %s %s', str_extent, str_contf_AWM, contf_unit, str_AWM_error); 
        else
            str_AWM = sprintf('AWM%s: %.5G %s %s', str_extent, contf_AWM, contf_unit, str_AWM_error); 
        end
        % str_AWM    = sprintf('Area-weighted mean: %.3G %s, Mean %.3G %s', contf_AWM, contf_unit, contf_M, contf_unit);  %/ almost the same 
    end
    
    %/ Output a struct of the plotted data for post-processing 
    S = [];
    S.contf_data_events = contf_data_events; %/ for other use, e.g., ttest2_sig_fn
    S.cont_data_events  = cont_data_events;  %/ for other use, e.g., ttest2_sig_fn
    S.Udata_events      = Udata_events;      %/ for other use, e.g., ttest2_sig_fn
    S.Vdata_events      = Vdata_events;      %/ for other use, e.g., ttest2_sig_fn
    S.contf_data        = contf_data;
    S.contf_lon         = contf_lon;
    S.contf_lat         = contf_lat;
    S.contf_unit        = contf_unit;
    S.point_data        = point_data;
    S.cont_data         = cont_data;
    S.cont_lon          = cont_lon;
    S.cont_lat          = cont_lat;
    S.cont_unit         = cont_unit;
    S.Udata             = Udata;
    S.Vdata             = Vdata;
    S.U2data            = U2data;
    S.V2data            = V2data;
    S.uv_lon            = uv_lon;
    S.uv_lat            = uv_lat;
    S.Rossby_westerly   = Rossby_westerly;
    S.Kelvin_easterly   = Kelvin_easterly;
    S.RK_ratio          = RK_ratio;
    S.contf_AWM         = contf_AWM;
    S.AWM_lon_extent = compute_AWM_lon_extent;
    S.AWM_lat_extent = compute_AWM_lat_extent;
  
    %====== to plot or not to plot ======%
    if plot_or_not
        %/ axis/f/figname
        titlename = strrep(titlename, '_', ' ');
        
        if all_in_one_format == 1  %/ tentative formatting for a plot with 4x2 panels
            fontsize       = fontsize*1.5;
            markersize     = markersize*1.5;
            linewi         = linewi*1.5;
            coast_wi       = coast_wi*1.5;
            cont_linewi    = cont_linewi*1.5;
            cont_labelsize = cont_labelsize*1.5;
            shaftwidth   = shaftwidth*1.5;    % control width of the shaft, the larger value the thicker
            headlength   = headlength*1.5; % control length of the head, the larger value the longer  
        end

        plot_contfmap('contf_data', contf_data, 'contf_lon', contf_lon, 'contf_lat', contf_lat, 'contf_levels', contf_levels,...
                      'contf_unit', contf_unit, 'colmap', colmap, 'cbar_interval', cbar_interval, 'pcolor_mode', pcolor_mode,...
                      'cont_data',  cont_data,  'cont_data_raw', cont_data_raw, 'cont_lon', cont_lon, 'cont_lat', cont_lat, 'cont_levels', cont_levels, 'cont_colmap', cont_colmap, 'cont_linewi', cont_linewi, 'cont_labelsize', cont_labelsize, 'cont_label_col', cont_label_col, 'skip_zero_cont', skip_zero_cont,...
                      'Udata', Udata, 'Vdata', Vdata, 'uv_lon', uv_lon, 'uv_lat', uv_lat, 'vec_step_lon', vec_step_lon, 'vec_step_lat', vec_step_lat,...
                      'vector_levels', vector_levels, 'vector_color', vector_color, 'vector_edgecolor', vector_edgecolor, 'vecscale', vecscale, 'vecscale2', vecscale2, 'shaftwidth', shaftwidth, 'headlength', headlength, 'vec_lbs', vec_lbs, 'vec_mag_ref', vec_mag_ref, 'vec_ref_fontsize', vec_ref_fontsize, 'vec_ref_lat_shift', vec_ref_lat_shift,...
                      'U2data', U2data, 'V2data', V2data, 'vector2_color', vector2_color, 'vector2_edgecolor', vector2_edgecolor,...
                      'hatch_data', hatch_data, 'hatch_lon', hatch_lon, 'hatch_lat', hatch_lat, 'hatch_thres_pve', hatch_thres_pve, 'hatch_thres_nve', hatch_thres_nve,...
                      'color_hatch_pve', color_hatch_pve, 'color_hatch_nve', color_hatch_nve, 'hatch_mode', hatch_mode, 'hatch_linewi', hatch_linewi, 'hatch_intvl', hatch_intvl,...
                      'point_data', point_data, 'bndry_data', bndry_data, 'marker', marker, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'linewi', linewi, 'grid_linewi', grid_linewi, 'color', color, 'draw_country', draw_country,...
                      'titlename', titlename, 'title_fontsize', title_fontsize, 'title_pos', title_pos, 'savepath', [], ...
                      'glb_data_mode', glb_data_mode, 'mask_mountains', mask_mountains, 'glb_plateau_mode', glb_plateau_mode, 'plateau_hgt', plateau_hgt, 'plateau_col', plateau_col,...
                      'map_proj', map_proj, 'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper, 'coast_col', coast_col, 'coast_wi', coast_wi, 'coast_patch_col', coast_patch_col, 'backcolor', backcolor,...
                      'fontsize', fontsize,  'create_fig', 1, 'grid_mode', grid_mode, 'draw_refvec_only', draw_refvec_only, 'draw_cbar_only', draw_cbar_only, 'cbar_mode', cbar_mode, 'cbar_location', cbar_location, 'cbar_fontsize', cbar_fontsize, 'cbar_YTick', cbar_YTick, 'cbar_YTickLabel', cbar_YTickLabel,...
                      'ax_panel', ax_panel);
        fprintf('done \n')
        
        if draw_cbar_only == 0
            %/ Show the centroid of the MJO
            if show_MJO_centroid
                markerfacecolor_MJO = [0 0 255]./255;
                m_line(MJO_centr(1), MJO_centr(2), 'marker', 'o', 'markersize', 16, 'markerfacecolor', markerfacecolor_MJO,...
                            'linest', 'none', 'color', 'none', 'linewi', linewi, 'clip', 'point');
                hold on;
            end
            if compute_RK_ratio
                % RK_pos = [0.7,0,1.0,0];
                RK_pos = title_pos; RK_pos(2) = RK_pos(2) + 0.04;
    
                annotation('textbox', 'String', str_RK_ratio, 'Color', 'k', ...
                        'FontSize', title_fontsize, 'Units', 'normalized', 'EdgeColor', 'none', ...
                        'Position', RK_pos)
            end
            if compute_contf_AWM   %/ Area-weighted mean of the domain shown
                % RK_pos = [0.7,0,1.0,0];
                AWM_pos = title_pos; 
                if compute_RK_ratio
                    AWM_pos(2) = AWM_pos(2) + 0.08;
                else
                    AWM_pos(2) = AWM_pos(2) + 0.04;
                end
                
                annotation('textbox', 'String', str_AWM, 'Color', 'k', ...
                        'FontSize', title_fontsize, 'Units', 'normalized', 'EdgeColor', 'none', ...
                        'Position', AWM_pos)
            end
            if province_on
                %/ load Chinese province borders
                %/ https://www.mathworks.com/matlabcentral/fileexchange/29462-china-province-mat
                filename =  strcat('/disk/r059/tfchengac/henan_floods_2021/prcssd_data_4plotting/china.province.mat');
                L = load(filename);
                m_line(L.long, L.lat, 'marker', 'none', 'markersize', 0.001, 'markerfacecolor', 'none',...
                                  'linest', '-', 'color', 'g', 'linewi', 1);
                drawnow; 
            end
            if draw_ZZ_county
                %/ load ZZ county borders
                %/ http://gaohr.win/site/blogs/2017/2017-04-18-GIS-basic-data-of-China.html
                filename =  strcat('/disk/r059/tfchengac/henan_floods_2021/prcssd_data_4plotting/Henan/Henan_county.shp');
                [shp,A] = shaperead(filename,'UseGeoCoords',true);
    
                ind_Zhengzhou = find([A.CENTROID_Y] > 34.4 & [A.CENTROID_Y] < 34.9 & ...
                                     [A.CENTROID_X] > 113 & [A.CENTROID_X] < 114.1);
    
                henan_county_bndry_data = [cat(2, shp(ind_Zhengzhou).Lon)', cat(2, shp(ind_Zhengzhou).Lat)'];
    
                m_line(henan_county_bndry_data(:,1), henan_county_bndry_data(:,2), 'marker', 'none', 'markersize', 0.001, 'markerfacecolor', 'none',...
                                      'linest', '-', 'color', 'r', 'linewi', 2);
                drawnow; 
            end
        end
        if savefig
            figname = char(strcat(plotting_folder, titlename, str_show_topo, str_cbar, str_refvec, '.', fig_fmt));
            figname = strrep(figname, ' ', '_');
            
            if isequal(fig_fmt, 'pdf')
                export_fig(figname,'-pdf','-painters', '-nocrop', '-transparent'); %, '-nocrop'); % '-c[inf, inf, inf, inf]', '-transparent');
            elseif isequal(fig_fmt, 'png')
                export_fig(figname,'-r300','-png','-opengl', '-nocrop', '-transparent'); %, '-nocrop');
            else
                error('code not set for the given fig_fmt = %s', fig_fmt);
            end
            fprintf('*** Figure saved into %s ***\n', figname)
        end
    end
end
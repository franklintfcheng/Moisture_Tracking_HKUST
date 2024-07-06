function S = quickplot_hovmoller(varargin)
    
    pnames = {'project_name', 'slct_date',  'zm_or_mm', 'lat_range', 'contf_dataname', 'cont_dataname', 'vector_dataname',...
              'contf_data', 'contf_hori', 'cont_data', 'cont_hori', 'stipp_data', 'stipp_hori', 'Udata', 'Vdata', 'U2data', 'V2data', 'uv_hori', 'contf_levels', 'colmap', 'cont_levels', 'cont_colmap',...
              'dataset', 'select_field', 'select_field_cont', 'sig_mode_contf', 'sig_mode_cont', 'sig_mode_vec', 'alpha', 'sig_contf_as_stipp', 'skip_zero_cont', 'not_show_cont_label',...
              'compo_date', 'timelag', 'hori_refpt', 'max_leadlag', 'show_hori_range', 'noleap', 'border_lines', 'time_lines', 'border_color', 'border_linestyle', 'border_linewidth', 'stipp_markersize', 'stipp_marker', 'stipp_color',...
              'value_in_diff', 'scatter_data', 'scatter_size', 'scatter_edgecolor', 'scatter_facecolor', 'scatter_alpha',...
              'fit_SLR_on', 'SLR_rmval', 'SLR_min_or_max', 'SLR_largest_obj', 'SLR_largest_obj_time', 'SLR_largest_obj_lon', 'SLR_largest_obj_lat', 'SLR_fitting_range_lon', 'SLR_recur',  'recur_max_LX', 'recur_min_RX', 'recur_Rsquared', 'SLR_test', 'SLR_alpha', 'cp_error', 'SLR_linewi', 'SLR_color',...
              'line_data', 'line_color', 'line_linewi', 'cont_linewi', 'always_time_as_y', 'draw_cbar_only', 'cbar_mode', 'cbar_location', 'pcolor_mode', 'draw_refvec_only', 'fontsize', 'title_fontsize', 'linewidth', 'map_xticks',...
              'gcf_position', 'plotting_folder', 'titlename', 'title_ypos', 'FigName_underscore', 'savefig', 'fig_fmt',  'png_dpi', 'plot_or_not'};
    
    dflts  = cell(length(pnames), 1);
    
    [          project_name,   slct_date,    zm_or_mm,   lat_range, contf_dataname,  cont_dataname,   vector_dataname,...
               contf_data, contf_hori, cont_data, cont_hori, stipp_data, stipp_hori, Udata, Vdata, U2data, V2data, uv_hori, contf_levels, colmap, cont_levels, cont_colmap,...
               dataset,   select_field,   select_field_cont,  sig_mode_contf, sig_mode_cont, sig_mode_vec, alpha, sig_contf_as_stipp, skip_zero_cont, not_show_cont_label,...
               compo_date,  timelag  ,  hori_refpt,   max_leadlag,  show_hori_range, noleap,  border_lines,  time_lines,    border_color,   border_linestyle,  border_linewidth, stipp_markersize, stipp_marker, stipp_color,...
               value_in_diff, scatter_data, scatter_size,   scatter_edgecolor, scatter_facecolor,  scatter_alpha,...
               fit_SLR_on,  SLR_rmval,  SLR_min_or_max,   SLR_largest_obj,     SLR_largest_obj_time,  SLR_largest_obj_lon,   SLR_largest_obj_lat,  SLR_fitting_range_lon,  SLR_recur,  recur_max_LX,   recur_min_RX, recur_Rsquared,  SLR_test,   SLR_alpha,  cp_error, SLR_linewi,   SLR_color,...
               line_data,  line_color,   line_linewi, cont_linewi, always_time_as_y,   draw_cbar_only,   cbar_mode,   cbar_location,   pcolor_mode,   draw_refvec_only,   fontsize, title_fontsize, linewidth, map_xticks,...
               gcf_position, plotting_folder,   titlename,   title_ypos, FigName_underscore, savefig, fig_fmt, png_dpi, plot_or_not] = ...
                                            internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Mar 27, 2024
    %/
    %/ Description: This function aims to quickly plot hovmoller diagram 
    %/              based on dataset (only when contf_data, cont_data, etc. 
    %/              are not given)
    %/
    %/              
    %/   'fit_SLR_on == 'contf': fit a simple linear regression (SLR) on 
    %                          contf_data based on the relevant settings.
    %/
    %/                   'cp': phase speed
    %/             'cp_error': phase speed error, can be an input (useful when the
    %/                         inter-model spread is computed elsewhere)
    %/=====================================================================
    % 
    % contf_data = []; contf_hori = []; contf_levels = []; colmap = [];
    % cont_data = []; cont_hori = []; cont_levels = []; cont_colmap = [];
    % stipp_data = []; stipp_hori= []; Udata = []; Vdata = []; U2data = [];
    % V2data = []; uv_hori = []; cp_error = []; cont_linewi = [];

    contf_multiply = 1;  cont_multiply = 1; contf_add = 0;  cont_add = 0;
    stat_contf_range = [-inf inf];   %/ by default
    stat_contf_mean_or_sum = 'mean'; %/ by default
    if isempty(zm_or_mm)
        zm_or_mm = 2;
    end

    if zm_or_mm == 1
        str_zm_or_mm = '_zm'; 
    elseif zm_or_mm == 2 
        str_zm_or_mm = '_mm';  
    end

    if isempty(select_field)
        select_field = 'daily';
    end

    if contains(select_field, 'daily')
        time_unit     = 'daily';
        fld_date      = 'date_yyyymmdd_AllYr';
        fld_date_clim = 'date_yyyymmdd_AllYr_clim';
        str_flds_date = 'str_daily_dates';
        Lr            = 8;
    elseif contains(select_field, 'pentad')
        time_unit     = 'pentad';
        fld_date      = 'date_yyyyptd';
        fld_date_clim = 'date_yyyyptd_clim';
        str_flds_date = 'str_ptd_dates';
        Lr            = 6;
    elseif contains(select_field, 'monthly')
        time_unit     = 'monthly';
        fld_date      = 'date_yyyymm';
        fld_date_clim = 'date_yyyymm_clim';
        str_flds_date = 'str_monthly_dates';
        Lr            = 6;
    else
        error('code not ready!');
    end

    if isempty(plot_or_not)
        plot_or_not = 1;
    end

    if isempty(linewidth)
        linewidth = 2;
    end

    if ~isempty(max_leadlag)
        max_leadlag = 25;
    end

    grid_mode = 0;
    %===== contf setting =====%
    if isempty(pcolor_mode)    pcolor_mode = 0;     end
    if isempty(cbar_mode)      cbar_mode   = 0;     end  %/ not recommond to use this, unless it doesn't distort fig.
    cbar_interval     = 2;
    contf_data_events = [];
    contf_unit        = '';
    if ~isempty(contf_dataname)
        %/ Horizontal axis
        if isempty(contf_hori)
            contf_hori  = dataset.(contf_dataname).hori_array;
        end
        if isempty(stipp_hori)
            stipp_hori  = contf_hori;
        end

        if isempty(contf_data)  %/ If contf_data is not given, get it from dataset and compute also stipp_data if needed
            flag_input = 0;
            date        = dataset.(contf_dataname).(fld_date);
            date_clim   = dataset.(contf_dataname).(fld_date_clim);
            sz          = size(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm)));
            
            if contains(select_field, '_clim')   %/ If to plot climatology, simply show it.
                if ~isempty(sig_mode_contf)
                    if isequal(sig_mode_contf, '_localcorr')
                        if sig_contf_as_stipp == 1
                            contf_data = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, '_corr'));
                            stipp_data = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, '_corr_sig'));
                        elseif sig_contf_as_stipp == 2
                            contf_data = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, '_corr_sig'));
                        else
                            contf_data = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, '_corr_sig'));
                        end
                    elseif isequal(sig_mode_contf, '_ASsigv2')
                        if sig_contf_as_stipp == 1
                            contf_data = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm));
                            stipp_data = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, sig_mode_contf));
                        elseif sig_contf_as_stipp == 2
                            contf_data = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, sig_mode_contf));
                        else
                            contf_data = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, sig_mode_contf));
                        end
                    elseif isequal(sig_mode_contf, '_ttest')
                        if sig_contf_as_stipp == 1
                            contf_data = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm));
                            stipp_data = dataset.(contf_dataname).(strcat(select_field, '_sig', str_zm_or_mm));
                        elseif sig_contf_as_stipp == 2
                            contf_data = dataset.(contf_dataname).(strcat(select_field, '_sig', str_zm_or_mm));
                        else
                            contf_data = dataset.(contf_dataname).(strcat(select_field, '_sig', str_zm_or_mm));
                        end
                    else
                        error('Code not set for the given sig_mode_contf ''%s''!', sig_mode_contf);
                    end
                else 
                    contf_data = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm)); %/ Assume daily_clim/pentad_clim etc.
                end

            elseif ~isempty(compo_date) && ~isempty(timelag)  %/ then check if to plot lag-lon/lat diagram
                NoOfevents = length(compo_date);
                timelag_full = timelag;  
                
                NoOflags     = length(timelag_full);
                year_list_bc = unique(floor(date./(10^(Lr-4)))); 
                all_dates    = date_array_gen('year_list', (year_list_bc(1)-1):(year_list_bc(end)+1), 'noleap', noleap); %/ +- 1 year as a buffer
                contf_data_events   = zeros([sz(1), NoOflags, NoOfevents]);  %/ make zeros for reduction
                for i = 1:NoOfevents
                    I       = findismember_loop(all_dates, compo_date(i));   
                    datelag = all_dates((I+timelag_full(1)):(I+timelag_full(end)));    %/ 1. Get the correct date lags
                    ind     = findismember_loop(date, datelag);              %/ 2. Get the indices of date lags based on the input data's date
                    if length(ind) < length(datelag)
                        warning('No enough dates for the quiried timelag for MJO Event #%d on %d. Skip it.', i, compo_date(i));
                    else
                        contf_data_events(:,:,i) = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm))(:,ind);  %/ 3. Store into contf_data
                    end
                end
                dim_event = 3; %/ Always assume to be 3. Do not write length(size(contf_data_events)), as it does not work for 1 event
                
                %/ 3. Further processing
                if NoOfevents > 1
                    if isequal(sig_mode_contf, '_ttest')
                        if sig_contf_as_stipp == 1
                            [stipp_data, contf_data, ~] = ttest_sig_fn(contf_data_events, alpha, dim_event, 0);
                        elseif sig_contf_as_stipp == 2
                            [contf_data, ~, ~] = ttest_sig_fn(contf_data_events, alpha, dim_event, 0);
                        else
                            [contf_data, ~, ~] = ttest_sig_fn(contf_data_events, alpha, dim_event, 0);
                        end
    
                    elseif ~isempty(sig_mode_contf)
                        error('Code not set for sig_mode_contf = %s!', sig_mode_contf);
    
                    else
                        contf_data = mean(contf_data_events, dim_event, 'omitnan');   %/ Simply take average over events 
                    end
                else
                    contf_data = mean(contf_data_events, dim_event, 'omitnan');   %/ Simply take average over events 
                end
            elseif ~isempty(max_leadlag)   %/ This is for lead-lag correlation
                ind = findismember_loop(date, slct_date);
                if isempty(ind) 
                    error('[quickplot_hovmoller] Empty ind! Check slct_date!');
                elseif length(ind) ~= length(slct_date)
                    warning('[quickplot_hovmoller] Not all slct_date are available from the data. Filling missing data with NaN (useful when showing near-real time data).');
                end
                contf_data = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm))(:, ind); %/ Assume daily/pentad
                if find(isnan(contf_data), 1)
                    error('contf_data contains nan!');
                end

                ind_refpt = find((contf_hori >= hori_refpt(1) & contf_hori <= hori_refpt(end))); %/ In case hori_refpt is a range
                if isempty(ind_refpt)
                    error('ind_refpt is empty! Check hori_refpt!');
                end
                refpt = mean(contf_data(ind_refpt,:), 1, 'omitnan');  %/ average over the reference point (if in a range)
                
                [nhori, ~] = size(contf_data);
                xcorr_M = nan(nhori, 2*max_leadlag+1);
                for i = 1:nhori
                    %/ NOTE: The first vector is actually to be shifted
                    %/ (try the example below to prove yourself!)
                    % B = randn(40,1);
                    % A = [zeros(4,1) ; B]; %/ A is a delayed version of B, delayed by 4 samples
                    % [xc,lags] = xcorr(A,B,20);
                    % stem(lags,xc)  %/ delay is positive
                    % [xc,lags] = xcorr(B,A,20); %/ now reverse order of inputs
                    % figure;
                    % stem(lags,xc) %delay is negative
                    xcorr_M(i,:) = xcorr(contf_data(i,:), refpt, max_leadlag, 'coeff');  %/ RMB to input 'coeff', otherwise the output is not correlation!
                end
                contf_data = xcorr_M;  %/ Update

            elseif ~isempty(slct_date)   %/ Then, check if quiring a certain period only
                ind = findismember_loop(date, slct_date);
                if isempty(ind) 
                    error('[quickplot_hovmoller] Empty ind! Check slct_date!');
                elseif length(ind) ~= length(slct_date)
                    warning('[quickplot_hovmoller] Not all slct_date are available from the data. Filling missing data with NaN (useful when showing near-real time data).');
                end
                %/ Inversely, find the indices of the available date in slct_date
                ind_inverse = findismember_loop(date_clim, mod(date(ind), 1e4));
                contf_data = nan(length(contf_hori), length(date_clim));
                contf_data(:,ind_inverse) = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm))(:, ind); %/ Assume daily/pentad
            end
        else
            flag_input = 1;
        end

        %/ Set levels & colormaps
        if isempty(contf_levels) || isempty(colmap)
            if contains(contf_dataname, 'OLR') || isequal(contf_dataname, 'rlut')
                contf_unit = 'W/m^{2}';
                if isequal(sig_mode_contf, '_localcorr') || ~isempty(max_leadlag)   %/ This is for lead-lag correlation
                    if value_in_diff
                        contf_levels = -0.16:0.02:0.16;
                    else
                        contf_levels  = (-1:0.1:1);  
                    end
                    NoOfColors = length(contf_levels)-1;
                    colmap = nclCM('CBR_drywet', NoOfColors);
                    % colmap = nclCM('hotcold_18lev', NoOfColors);
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];
                else
                    if contains(select_field, {'_MS'})
                        contf_levels  = 190:10:310;               %/ in W m^-2 
                        colmap        = my_colormap(length(contf_levels)-1, 'radar_12lev');
                        cbar_interval = 2;
    
                    elseif contains(select_field, {'_AC'})
                        contf_levels = -36:6:36;                  %/ in W m^-2 
                        NoOfColors   = length(contf_levels)-1;
                        colmap = nclCM('hotcold_18lev', NoOfColors);
    
                    elseif contains(select_field, {'_MJO'})
                        stat_contf_range = [-inf 0];  %/ Compute average OLR of convective MJO
                        if value_in_diff
                            contf_levels = -5.6:0.8:5.6;
                        else
                            if ~isempty(slct_date)
                                contf_levels  = [-60:10:60];            %/ in W m^-2, -ve -> convection
                            else
                                if contains(select_field, '_clim')
                                    contf_levels = [-6:1:6];                  %/ in W m^-2 
                                else
                                    contf_levels = [-35:5:35];               %/ in W m^-2 
                                end
                            end
                        end
                        NoOfColors   = length(contf_levels)-1;
                        colmap = nclCM('BlueYellowRed', NoOfColors);
    
                    elseif contains(select_field, {'_CISO', '_TE', '_IAID'})
                        contf_levels = [-6:1:6];                  %/ in W m^-2 
                        NoOfColors   = length(contf_levels)-1;
                        colmap = nclCM('BlueYellowRed', NoOfColors);
                        % colmap = nclCM('hotcold_18lev', NoOfColors);
                    else
                        contf_levels  = [190:10:310];             %/ in W m^-2 
                        colmap        = my_colormap(length(contf_levels)-1, 'radar_12lev');
                        cbar_interval = 2;
                    end
                end
    
            elseif contains(contf_dataname, '_P') || isequal(contf_dataname, 'pr') 
                contf_unit = 'mm/day';
                if contains(select_field, {'_AC'})
                    cbar_interval = 2;
                    contf_levels  = -6:1:6;             %/ in mm/day
                    NoOfColors    = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, 'BrBG');
                    % colmap = nclCM('precip_diff_12lev', NoOfColors);
                    precip_diff_12lev
    
                elseif contains(select_field, {'_MJO'})
                    % stat_contf_range = [0 inf];  %/ Compute average prcp of convective MJO

                    stat_contf_range = [0 inf];  %/ Compute average prcp of convective MJO
                    stat_contf_mean_or_sum = 'sum';
                    cbar_interval = 2;
                    if value_in_diff
                        contf_levels = -1.2:.2:1.2; 
                        % contf_levels = -0.9:.15:0.9; 
                    else
                        if contains(select_field, '_clim')
                            contf_levels  = -1.2:.2:1.2;        %/ in mm/day
                        else
                            contf_levels  = -3.6:0.6:3.6;       %/ in mm/day
                        end
                    end
                    NoOfColors    = length(contf_levels)-1;
                    % colmap = brewermap(NoOfColors, 'BrBG'); 
                    colmap = nclCM('precip_diff_12lev', NoOfColors);
                    
                elseif contains(select_field, {'_CISO', '_TE', 'IAID'})
                    cbar_interval = 2;
                    contf_levels  = [-1.2:.2:1.2];             %/ in mm/day
                    % contf_levels  = [-3:.5:3];             %/ in mm/day
                    NoOfColors    = length(contf_levels)-1;
                    % colmap = brewermap(NoOfColors, 'BrBG'); 
                    colmap = nclCM('precip_diff_12lev', NoOfColors);
                else
                    if isequal(project_name, 'trimonsoon')
                        cbar_interval = 3;
                        contf_levels  = [0:1:11];             %/ in mm/day
                        NoOfColors    = length(contf_levels)-1;
                        colmap = brewermap(NoOfColors, 'Blues');
                    else
                        cbar_interval = 1;
                        contf_levels  = [0:1:11];             %/ in mm/day
                        NoOfColors    = length(contf_levels)-1;
                        colmap = my_colormap(NoOfColors, 'precip3_11lev');  
                    end
                end
    
            elseif contains(contf_dataname, {'q1000', 'q850', 'hus'}) 
                contf_unit        = 'g/kg';
                contf_multiply    = 1000;   %/ kg/kg to g/kg
                if contains(select_field, {'_MJO'})
                    if value_in_diff
                        contf_levels  = (-0.36:0.06:0.36);
                    else
                        contf_levels  = (-2.4:0.4:2.4)./4;
                    end
                else
                    if value_in_diff
                        contf_levels  = -4.5:0.5:4.5;
                    else
                        if contains(contf_dataname, {'1000'})
                            contf_levels  = 1.25:1.25:21.25;
                        else
                            contf_levels  = 0:1:16;
                        end
                    end
                end
                colmap        = nclCM('CBR_drywet', length(contf_levels)-1);

            elseif contains(contf_dataname, {'SST'}) || ismember(contf_dataname, {'tos'})
                contf_unit    = sprintf('%sC', char(176));
                if contains(select_field, '_MJO')
                    if value_in_diff
                        contf_levels  = (-0.12:0.02:0.12)/2;
                    else
                        % contf_levels  = (-0.3:0.05:0.3);
                        contf_levels  = (-0.24:0.04:0.24);
                    end
                elseif ismember(select_field, {'monthly', 'daily'})
                    
                    if value_in_diff
                        contf_levels  = -0.25:0.25:5.25;
                        cbar_interval = 4;
                    else
                        % contf_add     = -273.15;   %/ K to deg C (only for raw data)
                        % contf_levels  = (22:1:33);
                        contf_levels  = (22:1:33)+273.15;
                    end
                end
                colmap        = nclCM('temp_diff_18lev', length(contf_levels)-1);

            elseif contains({contf_dataname}, {'u10m', 'v10m', 'u850', 'v850'})
                if contains(select_field, {'_MJO'})
                    contf_levels = (-3.6:0.6:3.6);             %/ in m/s
                else
                    contf_levels = (-12:2:12);             %/ in m/s
                end
                NoOfColors   = length(contf_levels)-1;
                colmap       = nclCM('GreenMagenta16', NoOfColors);
    
            elseif ismember({contf_dataname}, {'dE_TOA'})
                if contains(select_field, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    contf_levels  = [-11:1:11];                   %/ in W m^-2, -ve -> convection
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*PiYG');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
    
                elseif contains(select_field, {'pentad_clim_AC'})
                    contf_levels  = [-44:4:44];               %/ in W m^-2, -ve -> convection
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*PiYG');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
    
                elseif ismember(select_field, {'pentad_clim'})
                    contf_levels  = [-110:10:110];             %/ in W m^-2 
                    colmap = brewermap(length(contf_levels)-1, '*PiYG');
    
                elseif ismember(select_field, {'daily', 'daily_clim'})
                    contf_levels  = [-140:20:140];
                    colmap = brewermap(length(contf_levels)-1, '*PiYG');
                end
    
            elseif ismember({contf_dataname}, {'MSE_vi'})
                contf_data   = contf_data/1e9;  %/ scale by 1e-9.
                contf_levels = [2.9:0.02:3.12];             %/ in J m^-2 
                colmap = brewermap(length(contf_levels)-1, '*spectral');
    %             contf_levels = [2.9:0.02:3.12];             %/ in J m^-2 
    %             colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
    
            elseif ismember({contf_dataname}, {'MSE_vi650to1000'})
                contf_levels  = linspace(2.925, 3.2, 12)*1e9*0.39;              %/ in J m^-2 
                colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
    
            elseif ismember({contf_dataname}, {'MSE_vi100to650'})
                contf_levels  = linspace(2.925, 3.2, 12)*1e9*0.61;              %/ in J m^-2 
                colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
    
            elseif contains({contf_dataname}, {'dMSE_dt_vi'})
                contf_levels  = -40:4:40;             %/ in W m^-2 
                NoOfColors    = length(contf_levels)-1;
                colmap = brewermap(NoOfColors, '*RdYlBu');
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
    
            elseif contains({contf_dataname}, {'U_dMSE_dx_vi',   'V_dMSE_dy_vi',  'W_dMSE_dP_vi'})
                contf_data = contf_data * -1; %/ times a minus sign -> advection term.
                contf_levels  = -110:10:110;             %/ in W m^-2 
                NoOfColors    = length(contf_levels)-1;
                colmap = brewermap(NoOfColors, '*RdYlBu');
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
    
            elseif ismember({contf_dataname}, {'E'})
                contf_levels  = 0:.5:5.5;             %/ in W m^-2 
                colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                cbar_interval = 1;
                
            elseif ismember({contf_dataname}, {'SLHF'})
                cbar_interval = 1;
                if contains(select_field, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    contf_levels  = -5:.5:5;               %/ in W m^-2, -ve -> convection
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdYlBu');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
    
                elseif contains(select_field, {'pentad_clim_AC'})
                    contf_levels  = -30:3:30;               %/ in W m^-2, -ve -> convection
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdYlBu');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
                else
                    contf_levels  = 50:10:160;             %/ in W m^-2 
                    colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                end
                
            elseif ismember({contf_dataname}, {'SSHF'})
                cbar_interval = 1;
                if contains(select_field, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    contf_levels  = -10:1:10;               %/ in W m^-2, -ve -> convection
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdYlBu');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
    
                elseif contains(select_field, {'pentad_clim_AC'})
                    contf_levels  = -20:2:20;               %/ in W m^-2, -ve -> convection
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdYlBu');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
                else
                    contf_levels  = [0:5:55];             %/ in W m^-2 
                    colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                end
                
            elseif ismember({contf_dataname}, {'netSW'})
                contf_levels  = [70:5:125];             %/ in W m^-2 
                colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
    
            elseif ismember({contf_dataname}, {'netLW'})
                contf_levels  = [-230:10:-120];             %/ in W m^-2 
                colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
    
            elseif contains({contf_dataname}, {'QR_vi', 'netSWLW'}) 
                contf_levels  = [-220:20:220];              %/ in W m^-2 
                NoOfColors    = length(contf_levels)-1;
                colmap = brewermap(NoOfColors, '*RdYlBu');
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
    
            elseif contains({contf_dataname}, {'Q1_vi', 'Q2_vi'})
                contf_levels  = [-440:40:440];              %/ in W m^-2 
                NoOfColors    = length(contf_levels)-1;
                colmap = brewermap(NoOfColors, '*RdYlBu');
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1]; 
    
            elseif contains({contf_dataname}, {'MSE_Res'})
                contf_levels  = [-40:4:40];             %/ in W m^-2 
                NoOfColors    = length(contf_levels)-1;
                colmap = brewermap(NoOfColors, '*RdYlBu');
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
    
            elseif ismember({contf_dataname}, {'IVTmag'})
                if contains(select_field, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    contf_levels  = [-11:1:11];               %/ in W m^-2, -ve -> convection
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*PiYG');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
    
                elseif contains(select_field, {'pentad_clim_AC'})
                    contf_levels  = [-55:5:55];               %/ in W m^-2, -ve -> convection
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*PiYG');
                    colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
                else
                    contf_levels  = [0:30:330];             %/ in W m^-2 
                    colmap = jet(length(contf_levels)-1);
                end
    
            elseif ismember({contf_dataname}, {'dEPTdy'})
                contf_levels  = [-4:0.5:4];             %/ in K/100km 
                NoOfColors    = length(contf_levels)-1;
                colmap = brewermap(NoOfColors, '*PiYG');
    
            elseif ismember({contf_dataname}, {'Z850'})
                contf_levels  = [1340:20:1560];             %/ in gpm
                NoOfColors    = length(contf_levels)-1;
                colmap = brewermap(NoOfColors, '*spectral');
    
            elseif ismember({contf_dataname}, {'omega500'})
                contf_levels  = [-220:20:220];             %/ in hPa/day
                NoOfColors    = length(contf_levels)-1;
                colmap = brewermap(NoOfColors, '*RdYlBu');
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1]; 
                
            elseif ismember({contf_dataname}, {'T2m', 'T2mLand', 'T2mOcean'})
                if contains(select_field, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    contf_levels  = [-0.33:0.06:0.33];               %/ in W m^-2, -ve -> convection
                    cbar_interval = 1;
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdBu');
                    colmap(length(contf_levels)/2,:) = [ 1 1 1;];  
    
                elseif contains(select_field, {'pentad_clim_AC'})
                    contf_levels  = [-4.4:.8:4.4];               %/ in W m^-2, -ve -> convection
                    cbar_interval = 1;
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdBu');
    %                     colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
                    colmap(length(contf_levels)/2,:) = [ 1 1 1;];
                else
                    contf_data    = contf_data - 273.15; 
                    cbar_interval = 1;
                    contf_levels  = [20:1:31];   
                    NoOfColors    = length(contf_levels)-1; 
                    colmap = brewermap(NoOfColors, '*spectral');
                end
    
            elseif ismember({contf_dataname}, {'SST'})
                if contains(select_field, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    contf_levels  = [-0.33:0.06:0.33];               %/ in W m^-2, -ve -> convection
                    cbar_interval = 1;
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdBu');
                    colmap(length(contf_levels)/2,:) = [ 1 1 1;];
    
                elseif contains(select_field, {'pentad_clim_AC'})
                    contf_levels  = [-4.4:.8:4.4];                %/ in W m^-2, -ve -> convection
                    cbar_interval = 1;
                    NoOfColors = length(contf_levels)-1;
                    colmap = brewermap(NoOfColors, '*RdBu');
    %                     colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  
                    colmap(length(contf_levels)/2,:) = [ 1 1 1;];
                else
                    contf_data    = contf_data - 273.15; 
                    cbar_interval = 1;
                    contf_levels  = [25:0.5:30.5];
                    NoOfColors    = length(contf_levels)-1;
    %                     colmap        = cmocean('thermal', NoOfColors);
                    colmap        = brewermap(NoOfColors, '*spectral');
                end
            else
                error('contf_levels not pre-set for the input contf data!')
            end
        end

        %/ Avoid performing unit conversion when contf_data is directly input
        if flag_input == 0
            contf_data = contf_data*contf_multiply + contf_add;
            contf_data_events = contf_data_events*contf_multiply + contf_add; %/ for consistency
        end
        
        % contf_levels
        fprintf('Max contf_data = %.2f\n',  max(contf_data, [], 'all', 'omitnan'))
        fprintf('Mean contf_data = %.2f\n', mean(contf_data,    'all', 'omitnan'))
        fprintf('Min contf_data = %.2f\n',  min(contf_data, [], 'all', 'omitnan'))
    end
    
    %===== cont setting =====%
    cont_data_raw       = [];
    cont_data_events    = [];
    cont_labelsize      = [];
    cont_unit           = '';
    if isempty(cont_linewi)
        if ~isempty(linewidth)
            cont_linewi = linewidth*0.8;
        else
            cont_linewi = 2;
        end
    end
    cont_linest         = 'auto';  %/ Set it empty to allow it show dashed contour for -ve values
    if ~isempty(cont_dataname)
        if isempty(select_field_cont)
            select_field_cont = select_field;  %/ if empty, then assume the contour var uses the same common field.
        end
        if isempty(cont_hori)
            cont_hori = dataset.(cont_dataname).hori_array;
        end

        if isempty(cont_data)
            flag_input = 0;

            date      = dataset.(cont_dataname).(fld_date);
            date_clim = dataset.(cont_dataname).(fld_date_clim);
            sz        = size(dataset.(cont_dataname).(strcat(select_field, str_zm_or_mm)));

            if contains(select_field, '_clim') %/ If to plot climatology, simply show it.
                if ~isempty(sig_mode_cont)
                    if isequal(sig_mode_cont, '_localcorr')
                        cont_data     = dataset.(cont_dataname).(strcat(select_field_cont, str_zm_or_mm, '_corr_sig'));
                        cont_data_raw = dataset.(cont_dataname).(strcat(select_field_cont, str_zm_or_mm, '_corr'));
        
                    elseif isequal(sig_mode_cont, '_ASsigv2')
                        cont_data     = dataset.(cont_dataname).(strcat(select_field_cont, str_zm_or_mm, sig_mode_cont));
                        cont_data_raw = dataset.(cont_dataname).(strcat(select_field_cont, str_zm_or_mm));
        
                    elseif isequal(sig_mode_cont, '_ttest')
                        cont_data     = dataset.(cont_dataname).(strcat(select_field_cont, str_zm_or_mm));
                        cont_data_raw = dataset.(cont_dataname).(strcat(select_field_cont, '_sig', str_zm_or_mm));
                    else
                        error('Code not set for the given sig_mode_cont ''%s''!', sig_mode_cont);
                    end
                else
                    cont_data     = dataset.(cont_dataname).(strcat(select_field_cont, str_zm_or_mm));
                end

            elseif ~isempty(compo_date) && ~isempty(timelag)  %/ then check if to plot lag-lon/lat diagram
                NoOfevents   = length(compo_date);
                timelag_full = timelag;  
                NoOflags     = length(timelag_full);
                year_list_bc = unique(floor(date./(10^(Lr-4)))); 
                all_dates    = date_array_gen('year_list', (year_list_bc(1)-1):(year_list_bc(end)+1), 'noleap', noleap);
                cont_data_events   = zeros([sz(1), NoOflags, NoOfevents]);  %/ make zeros for reduction
                for i = 1:NoOfevents
                    I       = findismember_loop(all_dates, compo_date(i));   
                    datelag = all_dates((I+timelag_full(1)):(I+timelag_full(end)));    %/ 1. Get the correct date lags
                    ind     = findismember_loop(date, datelag);              %/ 2. Get the indices of date lags based on the input data's date
                    if length(ind) < length(datelag)
                        warning('No enough dates for the quiried timelag for MJO Event #%d on %d. Skip it.', i, compo_date(i));
                    else
                        cont_data_events(:,:,i) = dataset.(cont_dataname).(strcat(select_field, str_zm_or_mm))(:,ind);  %/ 3. Store into cont_data
                    end
                end
                dim_event = 3; %/ Always assume to be 3. Do not write length(size(cont_data_events)), as it does not work for 1 event

                %/ 3. Further processing
                %     %/ Average over events first, then convert into leadlagcorr data
                %     hovdata    = mean(cont_data_events, dim_event, 'omitnan');   
                %     cont_data  = hovmoller2leadlagcorr('hovdata', hovdata, 'hori', cont_hori, 'hori_refpt', hori_refpt,...
                %                           'along_timelag_or_hori', along_timelag_or_hori, 'sliding_range', sliding_range,...
                %                           'zm_or_mm', zm_or_mm, 'timelag', timelag, 'timelag_full', timelag_full);
                if NoOfevents > 1
                    if isequal(sig_mode_cont, '_ttest')
                        [cont_data, ~, ~] = ttest_sig_fn(cont_data_events, alpha, dim_event, 0);
    
                    elseif ~isempty(sig_mode_cont)
                        error('Code not set for sig_mode_contf = %s!', sig_mode_cont);
                    else
                        cont_data = mean(cont_data_events, dim_event, 'omitnan');   %/ 4. Take average
                    end
                else
                    cont_data = mean(cont_data_events, dim_event, 'omitnan');   %/ 4. Take average
                end

            elseif ~isempty(max_leadlag)
                ind = findismember_loop(date, slct_date);
                if isempty(ind) 
                    error('[quickplot_hovmoller] Empty ind! Check slct_date!');
                elseif length(ind) ~= length(slct_date)
                    warning('[quickplot_hovmoller] Not all slct_date are available from the data. Filling missing data with NaN (useful when showing near-real time data).');
                end
                cont_data = dataset.(cont_dataname).(strcat(select_field, str_zm_or_mm))(:, ind); %/ Assume daily/pentad
                
                ind_refpt = (cont_hori >= hori_refpt(1) & cont_hori <= hori_refpt(end));
                refpt = mean(cont_data(ind_refpt,:), 1, 'omitnan');  %/ average over the reference point (if in a range)
                
                [nhori, ~] = size(cont_data);
                xcorr_M = nan(nhori, 2*max_leadlag+1);
                for i = 1:nhori
                    %/ NOTE: The first vector is actually to be shifted
                    %/ (try the example below to prove yourself!)
                    % B = randn(40,1);
                    % A = [zeros(4,1) ; B]; %/ A is a delayed version of B, delayed by 4 samples
                    % [xc,lags] = xcorr(A,B,20);
                    % stem(lags,xc)  %/ delay is positive
                    % [xc,lags] = xcorr(B,A,20); %/ now reverse order of inputs
                    % figure;
                    % stem(lags,xc) %delay is negative
                    xcorr_M(i,:) = xcorr(cont_data(i,:), refpt, max_leadlag, 'coeff');  %/ RMB to input 'coeff', otherwise the output is not correlation!
                end
                cont_data = xcorr_M;  %/ Update
            
            elseif ~isempty(slct_date)   %/ Then, check if quiring a certain period only
                ind = findismember_loop(date, slct_date);
                if isempty(ind) 
                    error('[quickplot_hovmoller] Empty ind! Check slct_date!');
                elseif length(ind) ~= length(slct_date)
                    warning('[quickplot_hovmoller] Not all slct_date are available from the data. Filling missing data with NaN (useful when showing near-real time data).');
                end
                ind_inverse              = findismember_loop(date_clim, mod(date(ind), 1e4));
                cont_data                = nan(length(cont_hori), length(date_clim));
                cont_data(:,ind_inverse) = dataset.(cont_dataname).(strcat(select_field_cont, str_zm_or_mm))(:, ind); %/ Assume daily/pentad
            end
        else
            flag_input = 1;
        end

        if isempty(cont_levels) || isempty(cont_colmap)
            if contains(cont_dataname, 'OLR') || contains(cont_dataname, 'rlut')
                cont_unit = 'mm/day';
                if isequal(sig_mode_cont, '_localcorr') || ~isempty(max_leadlag)
                    if value_in_diff
                        cont_levels = -0.16:0.02:0.16;
                    else
                        cont_levels = [-1:0.1:-0.2, 0.2:0.1:1];
                    end
                    cont_colmap = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels); %/ Set it black if contf and cont are identical.
                    % cont_colmap = create_div_colmap([255 255 255]./255, [255 255 255]./255, cont_levels); %/ Set it black if contf and cont are identical.
                else
                    if contains(select_field_cont, {'_Res', '_CISO'})
                        cont_levels  = -12:2:12;               %/ in W m^-2, -ve -> convection
                        
                    elseif contains(select_field_cont, {'_MJO'})
                        if value_in_diff
                            cont_levels = -5.6:0.8:5.6;
                        else
                            if ~isempty(slct_date)
                                cont_levels  = -60:10:60;               %/ in W m^-2, -ve -> convection
                            else
                                if contains(select_field_cont, '_clim')
                                    cont_levels = -6:1:6;                  %/ in W m^-2 
                                else
                                    cont_levels  = -35:5:35;               %/ in W m^-2, -ve -> convection
                                end
                            end
                        end

                        %/ Thicken the -5 W m^-2 contour
                        cont_linewi = repmat(cont_linewi, length(cont_levels), 1);
                        ind = find(cont_levels == -5);
                        cont_linewi(ind) = cont_linewi(ind)*1.5;

                    elseif contains(select_field_cont, {'_AC'})
                        cont_levels  = -48:8:48;               %/ in W m^-2, -ve -> convection
                        
                    else
                        cont_levels  = 190:10:360;             %/ in W m^-2 
                        if ismember(select_field_cont, {'daily_clim'})
                            not_show_cont_label = 1;  %/ not to draw labels for clarify.
                        end
                    end  
                    % if isequal(contf_dataname, cont_dataname) 
                    cont_colmap = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels); %/ Set it black if contf and cont are identical.
                    % else
                        % cont_colmap = create_div_colmap([29 68 131]./255, [172 72 31]./255, cont_levels);
                    % end
                end
    
            elseif contains(cont_dataname, '_P') || isequal(cont_dataname, 'pr')
                cont_unit = 'mm/day';
                if contains(select_field, {'_AC'})
                    cont_levels  = -6:1:6;                 %/ in mm/day
    
                elseif contains(select_field, {'_MJO'})
                    if contains(select_field, '_clim')
                        cont_levels  = -1.2:.2:1.2;        %/ in mm/day
                    else
                        cont_levels  = -6:1:6;             %/ in mm/day
                    end
    
                elseif contains(select_field, {'_CISO', '_TE', 'IAID'})
                    cont_levels  = -1.2:.2:1.2;            %/ in mm/day
                else
                    cont_levels  = 0:1:11;                 %/ in mm/day
                end
                if isequal(contf_dataname, cont_dataname) 
                    cont_colmap = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels); %/ Set it black if contf and cont are identical.
                else
                    cont_colmap = brewermap(length(cont_levels), '*BrBG');
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
                cont_colmap = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);

            elseif contains(cont_dataname, {'SST'}) || ismember(cont_dataname, {'tos'})
                cont_unit    = sprintf('%sC', char(176));
                if contains(select_field, '_MJO')
                    if value_in_diff
                        cont_levels  = (-0.3:0.05:0.3)/2;
                    else
                        cont_levels  = -0.3:0.05:0.3;
                    end
                elseif ismember(select_field, {'monthly', 'daily'})
                    cont_add     = -273.15;   %/ K to deg C (only for raw data)
                    if value_in_diff
                        cont_levels  = -0.25:0.25:5.25;
                    else
                        cont_levels  = (22:1:33);
                    end
                end
                cont_colmap = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels);

            elseif contains({cont_dataname}, {'u10m', 'v10m', 'u850', 'v850'})
                if contains(select_field_cont, {'_MJO'})
                    cont_levels = (-3.6:0.6:3.6);             %/ in m/s
                else
                    cont_levels = (-12:2:12);             %/ in m/s
                end
                if isequal(contf_dataname, cont_dataname) 
                    cont_colmap = create_div_colmap([0 0 0]./255, [0 0 0]./255, cont_levels); %/ Set it black if contf and cont are identical.
                else
                    cont_colmap = create_div_colmap([0 102 0]./255, [255 51 204]./255, cont_levels); %/ Set it black if contf and cont are identical.
                    % cont_colmap = nclCM('GreenMagenta16', length(cont_levels));
                end
    
            elseif ismember({cont_dataname}, {'MSE_vi'})
                cont_data   = cont_data/1e9;                 %/ scale by 1e-9.
                cont_levels = [2.9:0.02:3.1];      %/ in J m^-2 
                cont_colmap = brewermap(length(cont_levels), 'Purples');
    %             cont_colmap = repmat([204 51 255]./255, length(cont_levels), 1);
                
            elseif ismember({cont_dataname}, {'dE_TOA'})
                if contains(select_field_cont, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    cont_levels  = [-11:1:11];               %/ in W m^-2, -ve -> convection
                    cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
    
                elseif contains(select_field_cont, {'pentad_clim_AC'})
                    cont_levels  = [-80:10:80];               %/ in W m^-2, -ve -> convection
                    cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
    
                elseif ismember(select_field_cont, {'pentad_clim'})
                    cont_levels  = [-120:5:-80, -60:20:60, 80:5:120];             %/ in W m^-2 
        %             cont_levels  = [-110:20:-10, 10:20:110];             %/ in W m^-2 
        %             cont_colmap = repmat([0 0 0], length(cont_levels), 1);
        %             cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                    cont_colmap = brewermap(length(cont_levels), '*PuOr');
                                
                elseif ismember(select_field_cont, {'daily', 'daily_clim'})
                    cont_levels  = [-140:10:140];
    %                 cont_colmap = brewermap(length(cont_levels), '*PuOr');
                    cont_colmap = brewermap(length(cont_levels), '*PiYG');
                    
                end
            
            elseif ismember({cont_dataname}, {'U_dMSE_dx_vi',   'V_dMSE_dy_vi',  'W_dMSE_dP_vi'})
                cont_data = cont_data * -1; %/ times a minus sign -> advection term.
                
                cont_levels = [-110:10:110];             %/ in W m^-2 
                cont_colmap = brewermap(length(cont_levels), '*RdYlBu');
                
            elseif ismember({cont_dataname}, {'SSHF', 'SSHFLand', 'SSHFOcean'})
                if contains(select_field_cont, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    cont_levels  = [-5:0.5:5];               %/ in W m^-2, -ve -> convection
                    cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
    
                elseif contains(select_field_cont, {'pentad_clim_AC'})
                    cont_levels  = [-20:2:20];               %/ in W m^-2, -ve -> convection
                    cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                    
                else
                    cont_levels  = [5:5:60];             %/ in W m^-2 
        %             cont_levels  = [-110:20:-10, 10:20:110];             %/ in W m^-2 
        %             cont_colmap = repmat([0 0 0], length(cont_levels), 1);
        %             cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                    cont_colmap = brewermap(length(cont_levels), 'PuRd');
        %             cont_colmap = cmocean('thermal', length(cont_levels));
                end
                
            elseif ismember({cont_dataname}, {'SLHF', 'SLHFLand', 'SLHFOcean'})
                if contains(select_field_cont, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    cont_levels  = [-10:1:10];               %/ in W m^-2, -ve -> convection
                    cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
    
                elseif contains(select_field_cont, {'pentad_clim_AC'})
                    cont_levels  = [-20:2:20];               %/ in W m^-2, -ve -> convection
                    cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
    
                else
                    cont_levels = [50:10:150];             %/ in W m^-2 
                    cont_colmap = my_colormap(length(cont_levels), 'precip3_11lev');
        %             cont_levels  = [-110:20:-10, 10:20:110];             %/ in W m^-2 
        %             cont_colmap = repmat([0 0 0], length(cont_levels), 1);
        %             cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
    %                 cont_colmap = brewermap(length(cont_levels), 'PuRd');
        %             cont_colmap = cmocean('thermal', length(cont_levels));
                end
                
            elseif ismember({cont_dataname}, {'EPT'})
    %             cont_levels  = [270:5:345];             %/ in K
                  cont_levels  = [295:5:345];             %/ in K
    %             cont_levels  = [-110:20:-10, 10:20:110];             %/ in W m^-2 
    %             cont_colmap = repmat([0 0 0], length(cont_levels), 1);
    %             cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                cont_colmap = brewermap(length(cont_levels), 'BrBG');
    %             cont_colmap = brewermap(length(cont_levels), '*RdYlBu');
    %             cont_colmap = brewermap(length(cont_levels), '*PuOr');
    %             cont_colmap = cmocean('thermal', length(cont_levels));
    %             cont_colmap = my_colormap(length(cont_levels), 'precip3_11lev');
    %               cont_colmap = repmat([0 0 0], length(cont_levels), 1);
                
            elseif ismember({cont_dataname}, {'dEPTdy'})
                cont_levels  = [-4:0.5:4];             %/ in K/100km 
    %             cont_levels  = [-110:20:-10, 10:20:110];             %/ in W m^-2 
    %             cont_colmap = repmat([0 0 0], length(cont_levels), 1);
    %             cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
    %             cont_colmap = brewermap(length(cont_levels), 'PuRd');
    %             cont_colmap = cmocean('thermal', length(cont_levels));
                cont_colmap = brewermap(length(cont_levels), '*PuOr');
                
            elseif ismember({cont_dataname}, {'omega500'})
                cont_levels  = [-110:10:110];             %/ in hPa/day
                cont_colmap  = brewermap(length(cont_levels), '*RdYlBu');
                
            elseif ismember({cont_dataname}, {'Z850'})
                cont_levels  = [1440:10:1560];             %/ in gpm
                cont_colmap  = jet(length(cont_levels));
    %             cont_colmap = brewermap(length(cont_levels), '*spectral');
    
            elseif contains({cont_dataname}, {'QR_vi', 'netSWLW'})
                cont_levels  = [-220:20:220];              %/ in W m^-2 
                cont_colmap  = brewermap(length(cont_levels), '*RdYlBu');     
                
            elseif contains({cont_dataname}, {'MSE_Res'})
                cont_levels  = [-40:4:40];               %/ in W m^-2 
    %             cont_colmap  = brewermap(length(cont_levels), '*RdYlBu'); 
                cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                
            elseif ismember({cont_dataname}, {'T2m', 'T2mLand'})
                if contains(select_field, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    cont_levels  = [-0.33:0.06:0.33];               %/ in W m^-2, -ve -> convection
    %                 cont_colmap = brewermap(length(cont_levels), '*RdBu');
                    cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                    
                elseif contains(select_field, {'pentad_clim_AC'})
                    cont_levels  = [-4.4:.8:4.4];               %/ in W m^-2, -ve -> convection
    %                 cont_colmap = brewermap(length(cont_levels), '*RdBu'); 
                    cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                else
                    cont_data    = cont_data - 273.15; 
                    cont_levels  = [20:1:31];
    %                 cont_levels  = [25:0.5:30.5];
    %                 cont_colmap = brewermap(length(cont_levels), '*spectral');
    %                 cont_colmap = brewermap(length(cont_levels), 'PuRd');
                    cont_colmap = repmat([204 51 255]./255, length(cont_levels), 1);
                end
    
            elseif ismember({cont_dataname}, {'SST'})
                if contains(select_field, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    cont_levels  = [-0.33:0.06:0.33];               %/ in W m^-2, -ve -> convection
    %                 cont_colmap = brewermap(length(cont_levels), '*RdBu');
                    cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
    
                elseif contains(select_field, {'pentad_clim_AC'})
                    cont_levels  = [-4.4:.8:4.4];                %/ in W m^-2, -ve -> convection
    %                 cont_colmap = brewermap(length(cont_levels), '*RdBu');
                    cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
                else
                    cont_data    = cont_data - 273.15; 
                    cont_levels  = [25:0.5:30.5];
    %                 cont_colmap  = brewermap(length(cont_levels), '*spectral');
    %                 cont_colmap = brewermap(length(cont_levels), 'PuRd');
                    cont_colmap = repmat([204 51 255]./255, length(cont_levels), 1);
                end
            else
                error('cont_levels not pre-set for the input cont data!')
            end
        end

        %/ Avoid performing unit conversion when contf_data is directly input
        if flag_input == 0
            cont_data = cont_data*cont_multiply + cont_add;
            cont_data_events = cont_data_events*cont_multiply + cont_add; %/ for consistency
        end
        
        fprintf('Mean cont_data = %.2f\n', mean(cont_data, 'all', 'omitnan'))
    end
    
    %===== vector setting =====%
    qscale = []; vector_levels = []; vector_color = []; vector_linewidth = [];
    vector_step_time = []; vector_step_hori = []; vecref_mag = [];
    if isempty(draw_refvec_only)  draw_refvec_only = 0;   end
    if ~isempty(vector_dataname)
        if isempty(uv_hori)
            uv_hori = dataset.(strcat('u', vector_dataname)).hori_array;
        end
        
        if isempty(Udata) && isempty(Vdata)
            date      = dataset.(strcat('u', vector_dataname)).(fld_date);
            date_clim = dataset.(strcat('u', vector_dataname)).(fld_date_clim);
            sz        = size(dataset.(strcat('u', vector_dataname)).(strcat(select_field, str_zm_or_mm)));
    
            if contains(select_field, '_clim')   %/ If to plot climatology, simply show it.
                if ~isempty(sig_mode_vec)
                    if isequal(sig_mode_vec, '_localcorr')
                        Udata  = dataset.(strcat('u', vector_dataname)).(strcat(select_field, str_zm_or_mm));
                        Vdata  = dataset.(strcat('v', vector_dataname)).(strcat(select_field, str_zm_or_mm));
                        U2data = dataset.(strcat('u', vector_dataname)).(strcat(select_field, str_zm_or_mm, '_corr_sig'));
                        V2data = dataset.(strcat('v', vector_dataname)).(strcat(select_field, str_zm_or_mm, '_corr_sig'));
        
                    elseif isequal(sig_mode_vec, '_ASsigv2')
                        Udata  = dataset.(strcat('u', vector_dataname)).(strcat(select_field, str_zm_or_mm));
                        Vdata  = dataset.(strcat('v', vector_dataname)).(strcat(select_field, str_zm_or_mm));
                        U2data = dataset.(strcat('u', vector_dataname)).(strcat(select_field, str_zm_or_mm, sig_mode_vec));
                        V2data = dataset.(strcat('v', vector_dataname)).(strcat(select_field, str_zm_or_mm, sig_mode_vec));
        
                    elseif isequal(sig_mode_vec, '_ttest')
                        Udata  = dataset.(strcat('u', vector_dataname)).(strcat(select_field, str_zm_or_mm));
                        Vdata  = dataset.(strcat('v', vector_dataname)).(strcat(select_field, str_zm_or_mm));
                        U2data = dataset.(strcat('u', vector_dataname)).(strcat(select_field, '_sig', str_zm_or_mm));
                        V2data = dataset.(strcat('v', vector_dataname)).(strcat(select_field, '_sig', str_zm_or_mm));
                    else
                        error('Code not set for the given sig_mode_vec ''%s''!', sig_mode_vec);
                    end
        
                    %/ [CAVEAT]: Show sig vector only, since two sets of quivers
                    %/           NEVER match with each other.
                    cond_U_sig = (~isnan(U2data));
                    cond_V_sig = (~isnan(V2data));
                    Udata(~cond_U_sig & ~cond_V_sig) = nan;
                    Vdata(~cond_U_sig & ~cond_V_sig) = nan;
                    U2data = []; V2data = [];
                else
                    Udata = dataset.(strcat('u', vector_dataname)).(strcat(select_field, str_zm_or_mm));
                    Vdata = dataset.(strcat('v', vector_dataname)).(strcat(select_field, str_zm_or_mm));
                end

            elseif ~isempty(compo_date) && ~isempty(timelag) %/ then check if to plot lag-lon/lat diagram
                NoOfevents   = length(compo_date);
                NoOflags     = length(timelag);
                year_list_bc = unique(floor(date./(10^(Lr-4)))); 
                all_dates    = date_array_gen('year_list', year_list_bc, 'noleap', noleap);
                Udata        = zeros([sz(1), NoOflags, NoOfevents]);  %/ make zeros for reduction
                Vdata        = Udata;
                for i = 1:NoOfevents
                    I       = findismember_loop(all_dates, compo_date(i));   
                    datelag = all_dates((I+timelag(1)):(I+timelag(end)));    %/ 1. Get the correct date lags
                    ind     = findismember_loop(date, datelag);              %/ 2. Get the indices of date lags based on the input data's date
                    if length(ind) < length(datelag)
                        warning('No enough dates for the quiried timelag for MJO Event #%d on %d. Skip it.', i, compo_date(i));
                    else
                        Udata(:,:,i) = dataset.(strcat('u', vector_dataname)).(strcat(select_field, str_zm_or_mm))(:,ind);  %/ 3. Store into data
                        Vdata(:,:,i) = dataset.(strcat('v', vector_dataname)).(strcat(select_field, str_zm_or_mm))(:,ind);  %/ 3. Store into data
                    end
                end
                dim_event = 3; %/ Always assume to be 3. Do not write length(size(Udata)), as it does not work for 1 event

                Udata = mean(Udata, dim_event, 'omitnan');   %/ 4. Take average
                Vdata = mean(Vdata, dim_event, 'omitnan');   %/ 4. Take average

                if dim_event == 2
                    warning('Seems like only one event for Udata! No ttest is performed!');
                else
                    if isequal(sig_mode_cont, '_ttest')
                        [U2data, ~, ~] = ttest_sig_fn(Udata, alpha, dim_event, 0);
                        [V2data, ~, ~] = ttest_sig_fn(Vdata, alpha, dim_event, 0);
    
                        %/ [CAVEAT]: Show sig vector only, since two sets of quivers
                        %/           will NEVER be matched with each other.
                        cond_U_sig = (~isnan(U2data));
                        cond_V_sig = (~isnan(V2data));
                        Udata(~cond_U_sig & ~cond_V_sig) = nan;
                        Vdata(~cond_U_sig & ~cond_V_sig) = nan;
                    end
                end

            elseif ~isempty(slct_date)   %/ Then, check if quiring a certain period only
                ind = findismember_loop(date, slct_date);
                if isempty(ind) 
                    error('[quickplot_hovmoller] Empty ind! Check slct_date!');
                elseif length(ind) ~= length(slct_date)
                    warning('[quickplot_hovmoller] Not all slct_date are available from the data. Filling missing data with NaN (useful when showing near-real time data).');
                end
                ind_inverse          = findismember_loop(date_clim, mod(date(ind), 1e4));
                Udata                = nan(length(uv_hori), length(date_clim));
                Vdata                = nan(length(uv_hori), length(date_clim));
                Udata(:,ind_inverse) = dataset.(strcat('u', vector_dataname)).(strcat(select_field, str_zm_or_mm))(:, ind); %/ Assume daily/pentad
                Vdata(:,ind_inverse) = dataset.(strcat('v', vector_dataname)).(strcat(select_field, str_zm_or_mm))(:, ind); %/ Assume daily/pentad
            end
        end

        if isequal(project_name, 'ITCC')
            %/ determine uv step in time-axis based on pentad or daily
            if isequal(time_unit, 'pentad')           vector_step_time = 2;
            else                                      vector_step_time = 4;  end

%             %/ determine uv step in time-axis based on hori dim.
            if size(Udata,1) < 15                     vector_step_hori = 2;
            else                                      vector_step_hori = 3;  end

            if contains({vector_dataname}, {'IVT'})
                vector_levels    = [0:30:300];   %/ in kg/m/s
                vector_color     = jet(length(vector_levels));
                vector_linewidth = linewidth;
                vecref_mag       = 300;   %/ in kg/m/s

                if contains(select_field, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    qscale = 1;
                elseif contains(select_field, {'pentad_clim_AC'})
                    qscale = 0.02;
                else
                    qscale = 0.02;
                end

            elseif contains({vector_dataname}, {'850'})
                vector_levels    = [];   %/ in m/s
                vector_color     = [0 0 0]./255;
                vector_linewidth = linewidth;
                if contains(select_field, '_CISO')
                    qscale     = 2;
                    vecref_mag = 2.5;   %/ in m/s
                else
                    qscale     = 0.5;
                    vecref_mag = 10;   %/ in m/s
                end
                
            elseif contains({vector_dataname}, {'10m'})
                vector_levels    = [];   %/ in m/s
                vector_color     = [0 0 0]./255;
                vector_linewidth = linewidth;
                if contains(select_field, '_CISO')
                    qscale     = 2;
                    vecref_mag = 2.5;   %/ in m/s
                else
                    qscale     = 0.5;
                    vecref_mag = 10;   %/ in m/s
                end
            else
                error('cont_levels not pre-set for the input contf data!')
            end
            
        elseif isequal(project_name, 'trimonsoon')
            %/ determine uv step in time-axis based on pentad or daily
            if isequal(time_unit, 'pentad')       vector_step_time = 3;
            else                                      vector_step_time = 6;  end

            %/ determine uv step in time-axis based on hori dim.
            if size(Udata,1) < 15                     vector_step_hori = 2;
            else                                      vector_step_hori = 3;  end
            
            if ismember({vector_dataname}, {'IVT'})
                vector_levels    = [0:30:300];   %/ in kg/m/s
                vector_color     = jet(length(vector_levels));
                vector_linewidth = linewidth;
                vecref_mag       = 300;   %/ in kg/m/s

                if contains(select_field, {'pentad_clim_Res', 'pentad_clim_CISO'})
                    qscale = 1;
                elseif contains(select_field, {'pentad_clim_AC'})
                    qscale = 0.02;
                else
                    qscale = 0.02;
                end

            elseif ismember({vector_dataname}, {'850'})
                vector_levels    = [0:1:10];   %/ in m/s
                vector_color     = jet(length(vector_levels));
                vector_linewidth = linewidth;
                vecref_mag       = vector_levels;   %/ in m/s
%                 vecref_mag       = repmat(14, 1, length(vector_levels));   %/ in m/s
%                 vecref_mag       = [0:1:7, 14 14, 14];   %/ in m/s
%                 vecref_mag       = [300, 1:9, 14];   %/ in m/s
                
                if zm_or_mm == 1        qscale = 1;
                elseif zm_or_mm == 2    qscale = 0.5; end
            else
                error('vector_levels not pre-set for the input vector data!')
            end
        else
            error('vector parameters not set for project_name = %s!', project_name);
        end
        fprintf('Mean Udata = %.2f\n', mean(Udata, 'all'))
        fprintf('Mean Vdata = %.2f\n', mean(Vdata, 'all'))
    end
    
    %===== grids/labels setting =====%    
    %/ Time axis
    if ~isempty(timelag) || ~isempty(max_leadlag)
        if ~isempty(max_leadlag)
            contf_time = (-max_leadlag:max_leadlag)';
        else
            contf_time = timelag;
        end
        I                = find(contf_time == 0);  %/ locate day 0
        intvl            = 10;
        ntime            = length(contf_time);
        ind              = [flip(I:-intvl:1), (I+intvl:intvl:ntime)];
        map_yticks       = contf_time(ind);
        map_yticklabels  = map_yticks;               %/ label pentad axis by month. string() is better than num2str as it eliminates spaces.
        map_ytickangle   = 0; 

        if isempty(fontsize) && zm_or_mm == 1   fontsize = 30;          end 
        if isempty(fontsize) && zm_or_mm == 2   fontsize = 30;          end 
    else
        contf_time       = 1:size(dataset.(contf_dataname).(str_flds_date), 1); 
        contf_year       = floor(slct_date./1e4);

        str_time_label   = dataset.(contf_dataname).(str_flds_date)(:,2);  %/ Assume the last column is about the month 
        [~, ind, ~]      = unique(str_time_label, 'stable');  %/ get the first occurrence index of each month.
        map_yticks       = ind;
        map_yticklabels  = str_time_label(ind);               %/ label pentad axis by month. string() is better than num2str as it eliminates spaces.
        % map_yticklabels  = strcat(str_time_label(ind), {' '}, string(contf_year(ind)));               %/ label pentad axis by month. string() is better than num2str as it eliminates spaces.
        map_ytickangle   = 0; 

        if isequal(time_unit, 'pentad')
            % str_time_label   = dataset.(contf_dataname).(str_flds_date)(:,2);  %/ Assume the last column is about the month 
            % [~, ind, ~]      = unique(str_time_label, 'stable');  %/ get the first occurrence index of each month.
            % map_yticks       = ind;
            % map_yticklabels  = str_time_label(ind);               %/ label pentad axis by month. string() is better than num2str as it eliminates spaces.
            % map_ytickangle   = 0; 

            if isempty(fontsize) && zm_or_mm == 1   fontsize = 22;          end 
            if isempty(fontsize) && zm_or_mm == 2   fontsize = 20;          end 
            
        elseif isequal(time_unit, 'daily')
            % intvl = 30;
            % ntime = length(contf_time);
            % map_yticks       = (1:intvl:ntime);
            % map_yticklabels  = contf_time(map_yticks);               %/ label pentad axis by month. string() is better than num2str as it eliminates spaces.
            % map_ytickangle   = 0; 
    
            if isempty(fontsize) && zm_or_mm == 1   fontsize = 30;        end 
            if isempty(fontsize) && zm_or_mm == 2   fontsize = 26;        end 
        end
    end

    if max(contf_hori) - min(contf_hori) > 150
        intvl_hori       = 60;  
    else
        intvl_hori       = 30;  
    end
    stipp_time = contf_time;
    cont_time  = contf_time;
    uv_time    = contf_time;
    line_time  = contf_time;

    %/ Finally, subset the range of data to show in the hovmoller diagram (if inquried)
    if ~isempty(show_hori_range)
        if ~isempty(contf_data)
            ind = find(contf_hori >= show_hori_range(1) & contf_hori <= show_hori_range(end));
            contf_hori = contf_hori(ind);
            contf_data = contf_data(ind,:);
            if ~isempty(contf_data_events)
                contf_data_events = contf_data_events(ind,:,:);
            end
        end
        if ~isempty(stipp_data)
            ind = find(stipp_hori >= show_hori_range(1) & stipp_hori <= show_hori_range(end));
            stipp_hori = stipp_hori(ind);
            stipp_data = stipp_data(ind,:);
        end
        if ~isempty(cont_data)
            ind = find(cont_hori >= show_hori_range(1) & cont_hori <= show_hori_range(end));
            cont_hori = cont_hori(ind);
            cont_data = cont_data(ind,:);
            if ~isempty(cont_data_raw)
                cont_data_raw = cont_data_raw(ind,:);
            end
            if ~isempty(cont_data_events)
                cont_data_events = cont_data_events(ind,:,:);
            end
        end
        if ~isempty(Udata)
            ind = find(uv_hori >= show_hori_range(1) & uv_hori <= show_hori_range(end));
            uv_hori = uv_hori(ind);
            Udata = Udata(ind,:);
            Vdata = Vdata(ind,:);
            if ~isempty(U2data)
                U2data = U2data(ind,:);
                V2data = V2data(ind,:);
            end
        end
    end

    %/ Compute average contf_data (for reference/writing)
    stat_contf = contf_data;
    cond = stat_contf > min(stat_contf_range) & stat_contf < max(stat_contf_range);
    stat_contf(~cond) = nan;

    stat_contf_unit = contf_unit;
    if isequal(stat_contf_mean_or_sum, 'mean')
        stat_contf = mean(stat_contf, 'all', 'omitnan');
    elseif isequal(stat_contf_mean_or_sum, 'sum')
        stat_contf = sum(stat_contf, 'all', 'omitnan');
        if contains(contf_dataname, '_P') || isequal(contf_dataname, 'pr') 
            stat_contf_unit = 'mm'; %/ Since we accumulate over time as well, so the unit [mm/day] becomes [mm].
        end
    end
    
    titlename = strcat(titlename, sprintf(', %s contf: %.2f %s,', stat_contf_mean_or_sum, stat_contf, stat_contf_unit));
    % stat_contf_range
    % stat_contf
    % titlename

    % close all;
    cp = [];
    if plot_or_not
        cp = plot_hovmoller('contf_data', contf_data, 'contf_hori', contf_hori, 'contf_time', contf_time, 'colmap', colmap, 'contf_levels', contf_levels,...
                       'pcolor_mode', pcolor_mode, 'cbar_mode', cbar_mode, 'cbar_location', cbar_location, 'cbar_interval', cbar_interval, 'draw_cbar_only', draw_cbar_only,...
                       'stipp_data', stipp_data, 'stipp_hori', stipp_hori, 'stipp_time', stipp_time, 'stipp_markersize', stipp_markersize, 'stipp_marker', stipp_marker, 'stipp_color', stipp_color,...
                       'cont_data',  cont_data, 'cont_data_raw', cont_data_raw, 'cont_hori', cont_hori,  'cont_time', cont_time, 'cont_colmap', cont_colmap, 'cont_levels', cont_levels, 'skip_zero_cont', skip_zero_cont, 'cont_linewi', cont_linewi, 'cont_linest', cont_linest, 'cont_labelsize', cont_labelsize, 'not_show_cont_label', not_show_cont_label,...
                       'Udata', Udata, 'Vdata', Vdata, 'U2data', U2data, 'V2data', V2data, 'uv_hori', uv_hori, 'uv_time', uv_time, 'qscale', qscale, 'vector_levels', vector_levels, 'vector_color', vector_color, 'vector_linewidth', vector_linewidth,...
                       'vector_step_time', vector_step_time, 'vector_step_hori', vector_step_hori, 'draw_refvec_only', draw_refvec_only, 'vecref_mag', vecref_mag,...
                       'border_lines', border_lines, 'time_lines', time_lines, 'border_color', border_color, 'border_linestyle', border_linestyle, 'border_linewidth', border_linewidth, 'fontsize', fontsize,  'title_fontsize', title_fontsize, 'linewidth', linewidth, 'grid_mode', grid_mode, 'map_yticks', map_yticks, 'map_yticklabels', map_yticklabels, 'map_ytickangle', map_ytickangle, 'intvl_hori', intvl_hori, 'map_xticks', map_xticks,...
                       'scatter_data', scatter_data, 'scatter_size', scatter_size, 'scatter_edgecolor', scatter_edgecolor, 'scatter_facecolor', scatter_facecolor, 'scatter_alpha', scatter_alpha,...
                       'fit_SLR_on', fit_SLR_on, 'SLR_rmval', SLR_rmval, 'SLR_min_or_max', SLR_min_or_max, 'SLR_largest_obj', SLR_largest_obj, 'SLR_largest_obj_time', SLR_largest_obj_time, 'SLR_largest_obj_lon', SLR_largest_obj_lon, 'SLR_largest_obj_lat', SLR_largest_obj_lat, 'SLR_fitting_range_lon', SLR_fitting_range_lon, 'SLR_recur', SLR_recur,  'recur_max_LX', recur_max_LX, 'recur_min_RX', recur_min_RX, 'recur_Rsquared', recur_Rsquared,...
                       'SLR_test', SLR_test, 'SLR_alpha', SLR_alpha, 'SLR_linewi', SLR_linewi, 'SLR_color', SLR_color, 'cp_error', cp_error, 'time_unit', time_unit, 'lat_range', lat_range,...
                       'line_data', line_data, 'line_time', line_time, 'line_color', line_color, 'line_linewi', line_linewi,...
                       'zm_or_mm', zm_or_mm, 'always_time_as_y', always_time_as_y, 'gcf_position', gcf_position, 'plotting_folder', plotting_folder, 'titlename', titlename, 'title_ypos', title_ypos, 'FigName_underscore',  FigName_underscore, 'savefig', savefig, 'fig_fmt', fig_fmt, 'png_dpi', png_dpi);
    end

    %/ Output a struct of the plotted data for post-processing 
    S = [];
    S.contf_data        = contf_data;
    S.contf_hori        = contf_hori;
    S.contf_time        = contf_time;
    S.contf_unit        = contf_unit;
    S.cont_data         = cont_data;
    S.cont_hori         = cont_hori;
    S.cont_time         = cont_time;
    S.cont_unit         = cont_unit;
    S.contf_data_events = contf_data_events;
    S.cont_data_events  = cont_data_events;
    S.cp                = cp;
end
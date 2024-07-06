
function S = quickplot_hgtlonlat(varargin)
    
    pnames = {'project_name', 'compo_date', 'timelag', 'zm_or_mm', 'lat_range', 'contf_dataname', 'cont_dataname', 'vector_dataname_hori', 'vector_dataname_vert',...
              'contf_data', 'contf_hori', 'contf_hgt', 'cont_data', 'cont_hori', 'cont_hgt', 'stipp_data', 'stipp_hori', 'vec_hori_data', 'vec_hgt_data', 'vec_hori', 'vec_hgt', 'contf_levels', 'colmap', 'cont_levels', 'cont_colmap',...
              'dataset', 'select_field', 'select_field_cont', 'sig_mode_contf', 'sig_mode_cont', 'sig_mode_vec', 'alpha', 'sig_contf_as_stipp', 'skip_zero_cont', 'not_show_cont_label',...
              'hori_refpt', 'sliding_range', 'show_hori_range', 'noleap', 'border_lines', 'time_lines', 'border_color', 'border_linestyle', 'border_linewidth', 'stipp_markersize', 'stipp_marker', 'stipp_color',...
              'value_in_diff', 'scatter_data', 'scatter_size', 'scatter_edgecolor', 'scatter_facecolor', 'scatter_alpha',...
              'line_data', 'line_color', 'line_linewi', 'cont_linewi', 'always_time_as_y', 'draw_cbar_only', 'cbar_mode', 'cbar_location', 'pcolor_mode', 'draw_refvec_only', 'fontsize', 'title_fontsize', 'linewidth', 'map_xticks',...
              'gcf_position', 'plotting_folder', 'titlename', 'title_ypos', 'FigName_underscore', 'savefig', 'fig_fmt',  'png_dpi', 'plot_or_not'};
    
    dflts  = cell(length(pnames), 1);
    
    [          project_name,   compo_date,  timelag,  zm_or_mm,   lat_range, contf_dataname,  cont_dataname,   vector_dataname_hori,  vector_dataname_vert,...
               contf_data, contf_hori, contf_hgt, cont_data, cont_hori, cont_hgt, stipp_data, stipp_hori, vec_hori_data, vec_hgt_data, vec_hori, vec_hgt, contf_levels, colmap, cont_levels, cont_colmap,...
               dataset,   select_field,   select_field_cont,  sig_mode_contf, sig_mode_cont, sig_mode_vec, alpha, sig_contf_as_stipp, skip_zero_cont, not_show_cont_label,...
               hori_refpt,   sliding_range,  show_hori_range, noleap,  border_lines,  time_lines,    border_color,   border_linestyle,  border_linewidth, stipp_markersize, stipp_marker, stipp_color,...
               value_in_diff, scatter_data, scatter_size,   scatter_edgecolor, scatter_facecolor,  scatter_alpha,...
               line_data,  line_color,   line_linewi, cont_linewi, always_time_as_y,   draw_cbar_only,   cbar_mode,   cbar_location,   pcolor_mode,   draw_refvec_only,   fontsize, title_fontsize, linewidth, map_xticks,...
               gcf_position, plotting_folder,   titlename,   title_ypos, FigName_underscore, savefig, fig_fmt, png_dpi, plot_or_not] = ...
                                            internal.stats.parseArgs(pnames, dflts, varargin{:});
%%

    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Apr 10, 2024
    %/
    %/ Description: This function aims to quickly plot hovmoller diagram 
    %/              based on dataset (only when contf_data, cont_data, etc. 
    %/              are not given)
    %/=====================================================================

    % contf_data = []; contf_hori = []; cont_data = []; cont_hori = []; cont_hgt = []; cont_colmap = []; cont_levels = []; 
    % slct_year = []; line_data = [];  line_hori = []; map_xticks = []; map_xticklabels = [];

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

    if ~isempty(compo_date)
        if isempty(timelag)
            timelag = 0;  %/ Lag 0 by default
        elseif numel(timelag)
            error('timelag must be a constant!')
        end
    end

    %===== contf setting =====%
    if isempty(pcolor_mode)
        pcolor_mode     = 0;
    end
    cbar_interval   = 2;
    cbar_YTick      = [];
    cbar_YTickLabel = [];
    if ~isempty(contf_dataname)
        %/ Horizontal axis
        if isempty(contf_hori)
            contf_hori  = dataset.(contf_dataname).hori_array;
        end
        if isempty(contf_hgt)
            P           = dataset.(contf_dataname).level;
            nlevel      = length(P);
            contf_hgt   = p2z_hydrostatic('P', P, 'unit', P_unit)./1000;  %/ to km
        end
        stipp_hori  = contf_hori;
        stipp_hgt   = contf_hgt;

        if isempty(contf_data)
            date        = dataset.(contf_dataname).(fld_date);
            date_clim   = dataset.(contf_dataname).(fld_date_clim);
            sz          = size(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm))); %/ (hori, plevel, time)
            if ~isequal(nlevel, sz(2))
                error('Check the dimension of your zm/mm data! plevel dimension should be the 2nd one!')
            end

            if contains(select_field, '_clim')   %/ If to plot climatology, simply show it.
                if ~isempty(compo_date)   %/ Then, check if quiring a certain period only
                    ind = findismember_loop(date_clim, compo_date);
                    if isempty(ind) 
                        error('[quickplot_hgtlonlat] Empty ind! Check slct_date!');
                    elseif length(ind) ~= length(compo_date)
                        warning('[quickplot_hgtlonlat] Not all compo_date are available from the data. Filling missing data with NaN (useful when showing near-real time data).');
                    end
                else
                    ind = 1:length(date_clim); 
                end

                if ~isempty(sig_mode_contf)
                    if isequal(sig_mode_contf, '_localcorr')
                        if sig_contf_as_stipp == 1
                            contf_data = mean(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, '_corr'))(:,:,ind), 3, 'omitnan');
                            stipp_data = mean(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, '_corr_sig'))(:,:,ind), 3, 'omitnan');
                        elseif sig_contf_as_stipp == 2
                            contf_data = mean(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, '_corr_sig'))(:,:,ind), 3, 'omitnan');
                        else
                            contf_data = mean(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, '_corr_sig'))(:,:,ind), 3, 'omitnan');
                        end
                    elseif isequal(sig_mode_contf, '_ASsigv2')
                        if sig_contf_as_stipp == 1
                            contf_data = mean(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm))(:,:,ind), 3, 'omitnan');
                            stipp_data = mean(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, sig_mode_contf))(:,:,ind), 3, 'omitnan');
                        elseif sig_contf_as_stipp == 2
                            contf_data = mean(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, sig_mode_contf))(:,:,ind), 3, 'omitnan');
                        else
                            contf_data = mean(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm, sig_mode_contf))(:,:,ind), 3, 'omitnan');
                        end
                    elseif isequal(sig_mode_contf, '_ttest')
                        if sig_contf_as_stipp == 1
                            contf_data = mean(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm))(:,:,ind), 3, 'omitnan');
                            stipp_data = mean(dataset.(contf_dataname).(strcat(select_field, '_sig', str_zm_or_mm))(:,:,ind), 3, 'omitnan');
                        elseif sig_contf_as_stipp == 2
                            contf_data = mean(dataset.(contf_dataname).(strcat(select_field, '_sig', str_zm_or_mm))(:,:,ind), 3, 'omitnan');
                        else
                            contf_data = mean(dataset.(contf_dataname).(strcat(select_field, '_sig', str_zm_or_mm))(:,:,ind), 3, 'omitnan');
                        end
                    else
                        error('Code not set for the given sig_mode_contf ''%s''!', sig_mode_contf);
                    end
                else 
                    contf_data = mean(dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm))(:,:,ind), 3, 'omitnan'); %/ Assume daily_clim/pentad_clim etc.
                end
            else
                NoOfevents   = length(compo_date); %/ then check if to plot lag-lon/lat diagram
                year_list_bc = unique(floor(date./(10^(Lr-4)))); 
                all_dates    = date_array_gen('year_list', (year_list_bc(1)-1):(year_list_bc(end)+1), 'noleap', noleap); %/ +- 1 year as a buffer
                contf_data_events = zeros([sz(1), sz(2), NoOfevents]);  %/ make zeros for reduction
                for i = 1:NoOfevents
                    I       = findismember_loop(all_dates, compo_date(i));   
                    datelag = all_dates(I+timelag);    %/ 1. Get the correct date lags
                    ind     = findismember_loop(date, datelag);              %/ 2. Get the indices of date lags based on the input data's date
                    if length(ind) < length(datelag)
                        warning('No enough dates for the quiried timelag for MJO Event #%d on %d. Skip it.', i, compo_date(i));
                    else
                        contf_data_events(:,:,i) = dataset.(contf_dataname).(strcat(select_field, str_zm_or_mm))(:,:,ind);  %/ 3. Store into contf_data
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
            end
        end

        local_grad_contf   = [];   %/ Though not use, keep it for now
        unitconv_local_grad = 1;   %/ unit conversion for local gradient
        if ~isempty(local_grad_contf)
            if ismember(local_grad_contf, [1, 2]) 
               unitconv_local_grad = 1e5;   %/ per m  -> per 100km. (hori. local grad)
            elseif ismember(local_grad_contf, 3) 
               unitconv_local_grad = 100;   %/ per Pa -> per hPa.     (vert. local grad)
            else
               error('Wrong input of ''local_grad_contf''!') 
            end
        end
        
        %/ Set levels & colormaps
        if isempty(contf_levels) || isempty(colmap)
            if ismember(contf_dataname, {'MSE'})
                contf_data = contf_data * 1e-5 * unitconv_local_grad;    
                
                if contains(select_field, {'_AC', '_CISO', '_HF'})
                    if unitconv_local_grad ~= 1 
                        contf_levels  = (-0.16:0.02:0.16)/10;          %/ [10^5 J kg-1 100km-1]    
                    else
                        contf_levels  = -0.16:0.02:0.16;             %/ [10^5 J kg-1]
                    end
                    colmap        = my_colormap(length(contf_levels)-1, 'amwg_blueyellowred_16lev_white');
                else
                    if unitconv_local_grad ~= 1 
                        contf_levels  = (3.3:0.02:3.52)/10;            %/ [10^5 J kg-1 100km-1]    
                    else
                        contf_levels  = 3.3:0.02:3.52;               %/ [10^5 J kg-1]
                    end
                    
                    colmap        = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                end
                cbar_interval = 2;
                
            elseif ismember({contf_dataname}, {'Q1', 'Q2'})     
                cp         = 1004;                                   %/ [J kg-1 K-1]
                contf_data = contf_data./cp*24*3600 * unitconv_local_grad; %/ [W kg-1] / [J kg-1 K-1] = [K s-1] -> [K day-1]
                                                                     %/ If local grad is computed, then [K day-1 m-1] -> [K day-1 100km-1]
                if unitconv_local_grad ~= 1 
                    contf_levels = (-8:1:8)/10;           %/ [K day-1 100km-1]
                else
                    contf_levels = -8:1:8;              %/ [K day-1]
                end
                if ismember({contf_dataname}, {'Q1'})    
                    colmap        = my_colormap(length(contf_levels)-1, 'amwg_blueyellowred_16lev_white');
                elseif ismember({contf_dataname}, {'Q2'})    
                    colmap        = brewermap(length(contf_levels)-1, 'BrBG');
                end
                cbar_interval = 2;
                
            elseif ismember({contf_dataname}, {  'dMSE_dt',   'MSE_Res', 'Q1_minus_Q2',...
                                               'W_dMSE_dP', 'U_dMSE_dx',   'V_dMSE_dy', 'hori_MSE_adv'})    
                cp         = 1004;                                   %/ [J kg-1 K-1]
                contf_data = contf_data./cp*24*3600 * unitconv_local_grad; %/ [W kg-1] / [J kg-1 K-1] = [K s-1] -> [K day-1]
                                                                     %/ If local grad is computed, then [K day-1 m-1] -> [K day-1 100km-1]
                
                if ismember({contf_dataname}, {'W_dMSE_dP', 'U_dMSE_dx', 'V_dMSE_dy'})
                    contf_data = contf_data * -1; %/ advection term
                end
                
                if contains(select_field, {'_AC', '_CISO', '_HF'})
                    %/ For the large terms (RHS)
                    if ismember({contf_dataname}, {'Q1_minus_Q2', 'W_dMSE_dP', 'U_dMSE_dx', 'V_dMSE_dy', 'hori_MSE_adv', 'MSE_Res'})
                        if unitconv_local_grad ~= 1
                            contf_levels = (-3.2:0.4:3.2)/2/10;     %/ [K day-1 100km-1]
                        else
                            contf_levels = (-3.2:0.4:3.2)/2;        %/ [K day-1]
                        end
                    %/ For the small term (LHS)
                    elseif ismember({contf_dataname}, {'dMSE_dt'})
                        if unitconv_local_grad ~= 1
                            contf_levels = (-0.8:0.1:0.8)/2/5;     %/ [K day-1 100km-1]
                        else
                            contf_levels = (-0.8:0.1:0.8)/2;        %/ [K day-1]
                        end
                    end
                else
                    %/ For the large terms (RHS)
                    if ismember({contf_dataname}, {'Q1_minus_Q2', 'W_dMSE_dP', 'U_dMSE_dx', 'V_dMSE_dy', 'hori_MSE_adv', 'MSE_Res'})
                        if unitconv_local_grad ~= 1
                            contf_levels = (-3.2:0.4:3.2)/10;     %/ [K day-1 100km-1]
                        else
                            contf_levels = -3.2:0.4:3.2;        %/ [K day-1]
                        end
                    %/ For the small term (LHS)
                    elseif ismember({contf_dataname}, {'dMSE_dt'})
                        if unitconv_local_grad ~= 1
                            contf_levels = (-0.8:0.1:0.8)/5;     %/ [K day-1 100km-1]
                        else
                            contf_levels = -0.8:0.1:0.8;        %/ [K day-1]
                        end
                    end
                end
                colmap          = my_colormap(length(contf_levels)-1, 'amwg_blueyellowred_16lev');
%                 colmap          = my_colormap(length(contf_levels)-1, 'amwg_blueyellowred_16lev_white');
                cbar_interval   = 2;
                cbar_YTick      = cbar_YTick(2:end-1);
                cbar_YTickLabel = cbar_YTick;
            else
                error('contf_levels not pre-set for the input contf data!')
            end
        end
        fprintf('Mean contf_data = %.2f\n', mean(contf_data, 'all', 'omitnan'));
    end
    
    %===== cont setting =====%
    if isempty(cont_linewi)
        cont_linewi    = 3;
    end
    if isempty(skip_zero_cont)
        skip_zero_cont = 1; %/ by default, skip zero-contour
    end
    if ~isempty(cont_dataname)
        
        %/ Horizontal axis
        if isempty(cont_hori)
            cont_hori  = dataset.(cont_dataname).hori_array;
        end
        if isempty(cont_hgt)
            P           = dataset.(cont_dataname).level;
            nlevel      = length(P);
            cont_hgt   = p2z_hydrostatic('P', P, 'unit', P_unit)./1000;  %/ to km
        end
        if isempty(cont_data)
            date        = dataset.(cont_dataname).(fld_date);
            date_clim   = dataset.(cont_dataname).(fld_date_clim);
            sz          = size(dataset.(cont_dataname).(strcat(select_field, str_zm_or_mm))); %/ (hori, plevel, time)
            if ~isequal(nlevel, sz(2))
                error('Check the dimension of your zm/mm data! plevel dimension should be the 2nd one!')
            end

            if contains(select_field, '_clim')   %/ If to plot climatology, simply show it.
                if ~isempty(compo_date)   %/ Then, check if quiring a certain period only
                    ind = findismember_loop(date_clim, compo_date);
                    if isempty(ind) 
                        error('[quickplot_hgtlonlat] Empty ind! Check slct_date!');
                    elseif length(ind) ~= length(compo_date)
                        warning('[quickplot_hgtlonlat] Not all compo_date are available from the data. Filling missing data with NaN (useful when showing near-real time data).');
                    end
                else
                    ind = 1:length(date_clim); 
                end

                if ~isempty(sig_mode_cont)
                    if isequal(sig_mode_cont, '_localcorr')
                        cont_data = mean(dataset.(cont_dataname).(strcat(select_field, str_zm_or_mm, '_corr_sig'))(:,:,ind), 3, 'omitnan');
                    elseif isequal(sig_mode_cont, '_ASsigv2')
                        cont_data = mean(dataset.(cont_dataname).(strcat(select_field, str_zm_or_mm, sig_mode_cont))(:,:,ind), 3, 'omitnan');
                    elseif isequal(sig_mode_cont, '_ttest')
                        cont_data = mean(dataset.(cont_dataname).(strcat(select_field, '_sig', str_zm_or_mm))(:,:,ind), 3, 'omitnan');
                    else
                        error('Code not set for the given sig_mode_cont ''%s''!', sig_mode_cont);
                    end
                else 
                    cont_data = mean(dataset.(cont_dataname).(strcat(select_field, str_zm_or_mm))(:,:,ind), 3, 'omitnan'); %/ Assume daily_clim/pentad_clim etc.
                end
            else
                NoOfevents   = length(compo_date); %/ then check if to plot lag-lon/lat diagram
                year_list_bc = unique(floor(date./(10^(Lr-4)))); 
                all_dates    = date_array_gen('year_list', (year_list_bc(1)-1):(year_list_bc(end)+1), 'noleap', noleap); %/ +- 1 year as a buffer
                cont_data_events = zeros([sz(1), sz(2), NoOfevents]);  %/ make zeros for reduction
                for i = 1:NoOfevents
                    I       = findismember_loop(all_dates, compo_date(i));   
                    datelag = all_dates(I+timelag);    %/ 1. Get the correct date lags
                    ind     = findismember_loop(date, datelag);              %/ 2. Get the indices of date lags based on the input data's date
                    if length(ind) < length(datelag)
                        warning('No enough dates for the quiried timelag for MJO Event #%d on %d. Skip it.', i, compo_date(i));
                    else
                        cont_data_events(:,:,i) = dataset.(cont_dataname).(strcat(select_field, str_zm_or_mm))(:,:,ind);  %/ 3. Store into contf_data
                    end
                end
                dim_event = 3; %/ Always assume to be 3. Do not write length(size(contf_data_events)), as it does not work for 1 event
                
                %/ 3. Further processing
                if NoOfevents > 1
                    if isequal(sig_mode_cont, '_ttest')
                        [cont_data, ~, ~] = ttest_sig_fn(cont_data_events, alpha, dim_event, 0);

                    elseif ~isempty(sig_mode_cont)
                        error('Code not set for sig_mode_cont = %s!', sig_mode_cont);
                    else
                        cont_data = mean(cont_data_events, dim_event, 'omitnan');   %/ Simply take average over events 
                    end
                else
                    cont_data = mean(cont_data_events, dim_event, 'omitnan');   %/ Simply take average over events 
                end
            end
        end

        local_grad_cont     = [];
        unitconv_local_grad = 1;   %/ uc == unit conversion.
        if ~isempty(local_grad_cont)
            if ismember(local_grad_cont, [1, 2]) 
               unitconv_local_grad = 1e5;   %/ per m  -> per 100km.   (hori. local grad)
            elseif ismember(local_grad_cont, 3) 
               unitconv_local_grad = 100;   %/ per Pa -> per hPa.     (vert. local grad)
            else
               error('Wrong input of ''local_grad_contf''!') 
            end
        end
        
        if ismember(cont_dataname, {'MSE'}) 
            cont_data   = cont_data * 1e-5 * unitconv_local_grad;
            if contains(select_field, {'_AC', '_CISO', '_HF'})
                if unitconv_local_grad ~= 1 
                    cont_levels    = (-0.16:0.02:0.16)/10;         %/ [10^5 J kg-1 100km-1]    
                else
                    cont_levels  = [-0.16:0.02:0.16];            %/ [10^5 J kg-1]
                end
                cont_colmap  = create_div_colmap([0 43 204]./255, [204 0 47]./255, cont_levels);
            else
                if unitconv_local_grad ~= 1 
                    cont_levels  = (3.3:0.02:3.52)/10;           %/ [10^5 J kg-1 100km-1]    
                else
                    cont_levels  = 3.3:0.02:3.5;               %/ [10^5 J kg-1]
                end
                cont_colmap  = my_colormap(length(cont_levels), 'radar_11lev');
            end
        
        %/ For the large terms (RHS) of the MSE balance
        elseif ismember({cont_dataname}, {'Q1_minus_Q2', 'W_dMSE_dP', 'U_dMSE_dx', 'V_dMSE_dy'})  
            cp         = 1004;                                  %/ [J kg-1 K-1]
            cont_data = cont_data./cp*24*3600*unitconv_local_grad;    %/ [W kg-1] / [J kg-1 K-1] = [K s-1] -> [K day-1]
                                                                %/ If local grad is computed, then [K day-1 m-1] -> [K day-1 100km-1]            
            if ismember({cont_dataname}, {'W_dMSE_dP', 'U_dMSE_dx', 'V_dMSE_dy'})
                cont_data = cont_data * -1; %/ advection term
            end
            
            if contains(select_field, {'_AC', '_CISO', '_HF'})
                if unitconv_local_grad ~= 1 
                    cont_levels = (-3.2:0.4:3.2)/2/10;       %/ [K day-1 100km-1]
                else
                    cont_levels = (-3.2:0.4:3.2)/2;          %/ [K day-1]
                end
            else
                if unitconv_local_grad ~= 1 
                    cont_levels = (-3.2:0.4:3.2)/10;         %/ [K day-1 100km-1]
                else
                    cont_levels = -3.2:0.4:3.2;            %/ [K day-1]
                end
            end
            cont_colmap  = create_div_colmap([0 43 204]./255, [204 0 47]./255, cont_levels);

         %/ For the small term (LHS) of the MSE balance
        elseif ismember({cont_dataname}, {'MSE_Res', 'dMSE_dt'})   
            cp         = 1004;                                  %/ [J kg-1 K-1]
            cont_data = cont_data./cp*24*3600*unitconv_local_grad;    %/ [W kg-1] / [J kg-1 K-1] = [K s-1] -> [K day-1]
                                                                %/ If local grad is computed, then [K day-1 m-1] -> [K day-1 100km-1]            
%                 skip_zero_cont = 0;
            if contains(select_field, {'_AC', '_CISO', '_HF'})
                if unitconv_local_grad ~= 1
                    cont_levels = (-0.8:0.1:0.8)/2/5;       %/ [K day-1 100km-1]
                else
                    cont_levels = (-0.8:0.1:0.8)/2;          %/ [K day-1]
                end
            else
                if unitconv_local_grad ~= 1
                    cont_levels = (-0.8:0.1:0.8)/5;         %/ [K day-1 100km-1]
                else
                    cont_levels = -0.8:0.1:0.8;            %/ [K day-1]
                end
            end
            cont_colmap  = create_div_colmap([0 43 204]./255, [204 0 47]./255, cont_levels);
            
        elseif ismember({cont_dataname}, {'Q1', 'Q2'})       %/ [J kg-1 s-1]
            cp        = 1004;                                %/ [J kg-1 K-1]
            cont_data = cont_data./cp*24*3600*unitconv_local_grad; %/ [K day-1] or [K day-1 100km-1]
            
            if unitconv_local_grad ~= 1
                cont_levels      = (-8:1:8)/10;   %/ [K day-1 100km-1]
            else
                cont_levels      = -8:1:8;      %/ [K day-1]
            end
            cont_levels(end) = [];
            skip_zero_cont   = 0;
            cont_colmap = my_colormap(length(cont_levels), 'amwg_blueyellowred_16lev_white');

        elseif ismember({cont_dataname}, {'U', 'V'})     
            cont_data = cont_data*unitconv_local_grad;  %/ [m s-1] or [m s-1 100km-1]
            
            if unitconv_local_grad ~= 1
                cont_levels      = (-20:2:20)/10; %/ [m s-1 100km-1]
            else
                cont_levels      = -20:2:20;    %/ [m s-1]
            end
            skip_zero_cont   = 0;
            cont_colmap  = create_div_colmap([0 43 204]./255, [204 0 47]./255, cont_levels);
        else
            error('cont_levels not pre-set for the input contf data!')
        end
        fprintf('Mean cont_data = %.2f\n', mean(cont_data, 'all', 'omitnan'))
    end
    
    %===== vector setting =====%
    qscale = [];  vector_color = []; vector_linewidth = [];
    vector_step_hori = []; vector_step_hgt = []; refvec_hori_mag = []; refvec_hgt_mag = [];
    if ~isempty(vector_dataname)
        % vec_hori_data = dataset.(vector_dataname_hori).(strcat(select_field, str_zm_or_mm, str_local_grad{local_grad_vec}))(:, :, t);
        % vec_hgt_data  = dataset.(vector_dataname_vert).(strcat(select_field, str_zm_or_mm, str_local_grad{local_grad_vec}))(:, :, t);

        %/ Horizontal axis
        if isempty(vec_hori)
            vec_hori  = dataset.(vector_dataname_hori).hori_array;
        end
        if isempty(cont_hgt)
            P           = dataset.(vector_dataname_hori).level;
            nlevel      = length(P);
            vec_hgt   = p2z_hydrostatic('P', P, 'unit', P_unit)./1000;  %/ to km
        end

        if isempty(vec_hori_data) && isempty(vec_hgt_data)
            date        = dataset.(vector_dataname_hori).(fld_date);
            date_clim   = dataset.(vector_dataname_hori).(fld_date_clim);
            sz          = size(dataset.(vector_dataname_hori).(strcat(select_field, str_zm_or_mm))); %/ (hori, plevel, time)
            if ~isequal(nlevel, sz(2))
                error('Check the dimension of your zm/mm data! plevel dimension should be the 2nd one!')
            end

            if contains(select_field, '_clim')   %/ If to plot climatology, simply show it.
                if ~isempty(compo_date)   %/ Then, check if quiring a certain period only
                    ind = findismember_loop(date_clim, compo_date);
                    if isempty(ind) 
                        error('[quickplot_hgtlonlat] Empty ind! Check slct_date!');
                    elseif length(ind) ~= length(compo_date)
                        warning('[quickplot_hgtlonlat] Not all compo_date are available from the data. Filling missing data with NaN (useful when showing near-real time data).');
                    end
                else
                    ind = 1:length(date_clim); 
                end

                if ~isempty(sig_mode_vec)
                    if isequal(sig_mode_vec, '_localcorr')
                        vec_hori_data = mean(dataset.(vector_dataname_hori).(strcat(select_field, str_zm_or_mm, '_corr_sig'))(:,:,ind), 3, 'omitnan');
                        vec_hgt_data  = mean(dataset.(vector_dataname_vert).(strcat(select_field, str_zm_or_mm, '_corr_sig'))(:,:,ind), 3, 'omitnan');
                    elseif isequal(sig_mode_vec, '_ASsigv2')
                        vec_hori_data = mean(dataset.(vector_dataname_hori).(strcat(select_field, str_zm_or_mm, sig_mode_vec))(:,:,ind), 3, 'omitnan');
                        vec_hgt_data  = mean(dataset.(vector_dataname_vert).(strcat(select_field, str_zm_or_mm, sig_mode_vec))(:,:,ind), 3, 'omitnan');
                    elseif isequal(sig_mode_vec, '_ttest')
                        vec_hori_data = mean(dataset.(vector_dataname_hori).(strcat(select_field, '_sig', str_zm_or_mm))(:,:,ind), 3, 'omitnan');
                        vec_hgt_data  = mean(dataset.(vector_dataname_vert).(strcat(select_field, '_sig', str_zm_or_mm))(:,:,ind), 3, 'omitnan');
                    else
                        error('Code not set for the given sig_mode_vec ''%s''!', sig_mode_vec);
                    end
                else 
                    vec_hori_data = mean(dataset.(vector_dataname_hori).(strcat(select_field, str_zm_or_mm))(:,:,ind), 3, 'omitnan'); %/ Assume daily_clim/pentad_clim etc.
                end
            else
                NoOfevents   = length(compo_date); %/ then check if to plot lag-lon/lat diagram
                year_list_bc = unique(floor(date./(10^(Lr-4)))); 
                all_dates    = date_array_gen('year_list', (year_list_bc(1)-1):(year_list_bc(end)+1), 'noleap', noleap); %/ +- 1 year as a buffer
                cont_data_events = zeros([sz(1), sz(2), NoOfevents]);  %/ make zeros for reduction
                for i = 1:NoOfevents
                    I       = findismember_loop(all_dates, compo_date(i));   
                    datelag = all_dates(I+timelag);    %/ 1. Get the correct date lags
                    ind     = findismember_loop(date, datelag);              %/ 2. Get the indices of date lags based on the input data's date
                    if length(ind) < length(datelag)
                        warning('No enough dates for the quiried timelag for MJO Event #%d on %d. Skip it.', i, compo_date(i));
                    else
                        cont_data_events(:,:,i) = dataset.(vector_dataname_hori).(strcat(select_field, str_zm_or_mm))(:,:,ind);  %/ 3. Store into contf_data
                    end
                end
                dim_event = 3; %/ Always assume to be 3. Do not write length(size(contf_data_events)), as it does not work for 1 event
                
                %/ 3. Further processing
                if NoOfevents > 1
                    if isequal(sig_mode_vec, '_ttest')
                        [vec_hori_data, ~, ~] = ttest_sig_fn(cont_data_events, alpha, dim_event, 0);

                    elseif ~isempty(sig_mode_vec)
                        error('Code not set for sig_mode_vec = %s!', sig_mode_vec);
                    else
                        vec_hori_data = mean(cont_data_events, dim_event, 'omitnan');   %/ Simply take average over events 
                    end
                else
                    vec_hori_data = mean(cont_data_events, dim_event, 'omitnan');   %/ Simply take average over events 
                end
            end
        end

        vector_linewidth = 2.2;
        vector_step_hori    = 2;    %/ Reminder: here 'UV' is not zonal / merid, but horiz and vert winds.
        vector_step_hgt     = 2;
        
        vec_hori_factor = 1; vec_hgt_factor = 1;
        if isequal(vector_dataname, 'UVW')
            vec_hori_factor  = abs(diff(vec_hori(1:2)))/2;      %/If assume horizonal traj speed ~20 m/s. x-axis interval = 1 deg. (ignore unit conversion)
            vec_hgt_factor   = abs(diff(vec_hgt(1:2)))*-5;      %/NOTE: the negative sign is to get the direction of vertical velocity right. if assume vetical traj speed ~10 cm/s. z-axis interval divided by w magnitude = diff(vec_hgt(1:2))/w = 1/50; (ignore unit conversion)

            if zm_or_mm == 1
                qscale           = 2;
                refvec_hori_mag  = 3;         %/ [U] or [V]: m/s 
                refvec_hgt_mag   = -0.05;     %/ [W]: Pa/s  (NOTE: 1 hPa/h = 0.0278 Pa/s)
            elseif zm_or_mm == 2              %/ since zonal winds are much stronger.
                qscale           = 1;
                refvec_hori_mag  = 6;         %/ [U] or [V]: m/s
                refvec_hgt_mag   = -0.12;     %/ [W]: Pa/s  (NOTE: 1 hPa/h = 0.0278 Pa/s)
            end
            vector_color     = [0 0 0]./255; 
        else
            error('code not set');
        end
        fprintf('Mean vec_hori_data = %.2f, Mean vec_hori_data * vec_hori_factor = %.2f\n', mean(vec_hori_data, 'all', 'omitnan'), mean(vec_hori_data, 'all', 'omitnan')*vec_hori_factor)
        fprintf('Mean vec_hgt_data = %.2f,  Mean vec_hgt_data * vec_hgt_factor = %.2f\n',   mean(vec_hgt_data,  'all', 'omitnan'), mean(vec_hgt_data,   'all', 'omitnan')*vec_hgt_factor)
        fprintf('Max(abs(vec_hori_data)) = %.2f, Max(abs(vec_hori_data)) * vec_hori_factor = %.2f\n', max(abs(vec_hori_data), [], 'all'), max(abs(vec_hori_data), [], 'all')*vec_hori_factor)
        fprintf('Max(abs(vec_hgt_data)) = %.2f,  Max(abs(vec_hgt_data)) * vec_hgt_factor = %.2f\n',   max(abs(vec_hgt_data), [], 'all'),  max(abs(vec_hgt_data), [], 'all')*vec_hgt_factor)
    end
    
    %/ grids/labels
    grid_mode       = 0;  %/ use this default setting
    create_fig      = 1;

    if isempty(gcf_position)
        gcf_position    = [600 100 900 750];
    end
    if isemtpy(fontsize)
        fontsize        = 22; 
    end
    if isempty(linewidth)
        linewidth       = 2;
    end
    if isequal(P_unit, 'hPa')
        ind              = find(mod(P, 100) == 0);
        map_yticks       = contf_hgt(ind);  
        map_yticklabels  = P(ind);  %/ show at 1000, 900, 800 hPa, etc.
        map_ylabel       = '';      %/ Not to plot ylabel to avoid distortion of the fig. 'p (hPa)';
    else
        error('define P_unit!');
    end

    %/ plotting
    if plot_or_not
        plot_hgtlonlat('contf_data', contf_data, 'contf_hori', contf_hori, 'contf_hgt', contf_hgt, 'colmap', colmap, 'contf_levels', contf_levels, 'pcolor_mode', pcolor_mode,...
                      'cbar_mode', cbar_mode, 'cbar_location', cbar_location, 'cbar_interval', cbar_interval, 'cbar_YTick', cbar_YTick, 'cbar_YTickLabel', cbar_YTickLabel, 'draw_cbar_only', draw_cbar_only,...
                      'cont_data', cont_data, 'cont_hori', cont_hori, 'cont_hgt', cont_hgt, 'cont_colmap', cont_colmap, 'cont_levels', cont_levels, 'cont_linewi', cont_linewi, 'skip_zero_cont', skip_zero_cont,...
                      'vec_hori_data', vec_hori_data, 'vec_hgt_data', vec_hgt_data, 'vec_hori', vec_hori, 'vec_hgt', vec_hgt, 'vec_hori_factor', vec_hori_factor, 'vec_hgt_factor', vec_hgt_factor, 'qscale', qscale, 'vector_color', vector_color,...
                      'vector_linewidth', vector_linewidth, 'vector_step_hori', vector_step_hori, 'vector_step_hgt', vector_step_hgt, 'draw_refvec_only', draw_refvec_only, 'refvec_hori_mag', refvec_hori_mag, 'refvec_hgt_mag', refvec_hgt_mag,...
                      'line_data', line_data, 'line_hori', line_hori, 'scatter_data', scatter_data, 'scatter_hgt', scatter_hgt, 'border_lines', border_lines, 'border_color', border_color, 'border_linestyle', border_linestyle, 'border_linewidth', border_linewidth,...
                      'zm_or_mm', zm_or_mm, 'terrain_mode', terrain_mode, 'terrain_bndry', terrain_bndry,...
                      'grid_mode', grid_mode, 'map_xticks', map_xticks, 'map_xticklabels', map_xticklabels, 'map_yticks', map_yticks, 'map_yticklabels', map_yticklabels, 'map_ylabel', map_ylabel,...
                      'create_fig', create_fig, 'gcf_position', gcf_position, 'fontsize', fontsize, 'linewidth', linewidth, 'plotting_folder', plotting_folder, 'titlename', titlename, 'FigName_underscore', FigName_underscore, 'savefig', savefig)
    end

    %/ Output a struct of the plotted data for post-processing 
    S = [];
    S.contf_data        = contf_data;
    S.contf_hori        = contf_hori;
    S.contf_hgt         = contf_hgt;
    S.cont_data         = cont_data;
    S.cont_hori         = cont_hori;
    S.cont_hgt          = cont_hgt;
    S.vec_hori_data     = vec_hori_data;
    S.vec_hgt_data      = vec_hgt_data;
    S.vec_hori          = vec_hori;
    S.vec_hgt           = vec_hgt;
    S.vec_hori_factor   = vec_hori_factor;
    S.vec_hgt_factor    = vec_hgt_factor;
    S.contf_data_events = contf_data_events;
    S.cont_data_events  = cont_data_events;
    S.cp                = cp;
        

end
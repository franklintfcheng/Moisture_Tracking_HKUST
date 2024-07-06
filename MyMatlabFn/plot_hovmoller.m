function cp = plot_hovmoller(varargin)
      
    pnames       = {'contf_data',    'contf_hori',        'contf_time',           'colmap',                'contf_levels',     'contf_unit', 'pcolor_mode',  'cbar_mode',    'cbar_location', 'cbar_interval',  'draw_cbar_only',...
                    'stipp_data',    'stipp_hori',        'stipp_time',           'stipp_markersize',      'stipp_marker',     'stipp_color',...
                    'cont_data',     'cont_data_raw',     'cont_hori',            'cont_time',             'cont_colmap',      'cont_levels',  'skip_zero_cont', 'cont_linewi', 'cont_linest', 'cont_labelsize', 'not_show_cont_label', ...
                    'Udata',         'Vdata',             'U2data',               'V2data',                'uv_hori',          'uv_time',      'qscale',  'vector_levels',       'vector_color', 'vector_linewidth', 'vector_step_time', 'vector_step_hori',  'draw_refvec_only', 'vecref_mag',...
                    'border_lines',  'time_lines',        'border_color',         'border_linestyle',      'border_linewidth',...
                    'line_data',     'line_time',         'line_color',           'line_linewi',           'fontsize',         'title_fontsize', 'linewidth',    'grid_mode',    'map_yticks',    'map_yticklabels',  'map_ytickangle', 'intvl_hori',  'map_xticks',...
                    'scatter_data',  'scatter_size',      'scatter_edgecolor',    'scatter_facecolor',     'scatter_alpha',...
                    'fit_SLR_on', 'SLR_rmval',   'SLR_min_or_max', 'SLR_largest_obj', 'SLR_largest_obj_time', 'SLR_largest_obj_lon', 'SLR_largest_obj_lat', 'SLR_fitting_range_lon', 'SLR_recur', 'recur_max_LX', 'recur_min_RX', 'recur_Rsquared',...
                    'SLR_test',      'SLR_alpha',         'SLR_linewi',           'SLR_color', 'cp_error', 'time_unit', 'lat_range',...
                    'zm_or_mm',      'always_time_as_y',  'gcf_position',         'plotting_folder',  'titlename',    'title_ypos', 'FigName_underscore', 'savefig', 'fig_fmt',  'png_dpi'};
    
    dflts        = cell(length(pnames), 1);
    
    [                contf_data,     contf_hori,         contf_time,           colmap,                contf_levels,        contf_unit,     pcolor_mode,    cbar_mode,      cbar_location,     cbar_interval,    draw_cbar_only,...
                     stipp_data,     stipp_hori,         stipp_time,           stipp_markersize,      stipp_marker,        stipp_color,...
                     cont_data,      cont_data_raw,      cont_hori,            cont_time,             cont_colmap,         cont_levels,    skip_zero_cont, cont_linewi, cont_linest,   cont_labelsize,    not_show_cont_label,...
                     Udata,          Vdata,              U2data,               V2data,                uv_hori,             uv_time,        qscale,  vector_levels,          vector_color,   vector_linewidth,  vector_step_time,   vector_step_hori,    draw_refvec_only,   vecref_mag,...
                     border_lines,   time_lines,         border_color,         border_linestyle,      border_linewidth,...
                     line_data,      line_time,          line_color,           line_linewi,           fontsize,            title_fontsize, linewidth,    grid_mode,      map_yticks,     map_yticklabels,  map_ytickangle,    intvl_hori,     map_xticks,...
                     scatter_data,   scatter_size,       scatter_edgecolor,    scatter_facecolor,     scatter_alpha,...
                     fit_SLR_on,  SLR_rmval,    SLR_min_or_max, SLR_largest_obj, SLR_largest_obj_time,  SLR_largest_obj_lon, SLR_largest_obj_lat, SLR_fitting_range_lon, SLR_recur, recur_max_LX,  recur_min_RX,   recur_Rsquared,...
                     SLR_test,       SLR_alpha,          SLR_linewi,           SLR_color,             cp_error,                    time_unit,                 lat_range,...
                     zm_or_mm,       always_time_as_y,   gcf_position,         plotting_folder,       titlename,                   title_ypos,                FigName_underscore,       savefig,  fig_fmt,   png_dpi] = ...
            internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
%%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Apr 2, 2024
    %/
    %/ This function is designed for plotting Homvoller diagrams by overlaying 
    %/ contourf, contour, vector and line plots using MATLAB built-in functions.
    %/=====================================================================
    cp = nan;

    %/ transpose
    if isempty(draw_cbar_only)       cbar_mode  = 1;                 end
    if isempty(grid_mode)            grid_mode  = 0;                 end
    if isempty(skip_zero_cont)       skip_zero_cont = 0;             end
    if isempty(not_show_cont_label)  not_show_cont_label = 0;        end
    if size(contf_hori, 1) == 1      contf_hori = contf_hori';       end
    if size(cont_hori, 1)  == 1      cont_hori  = cont_hori';        end
    if size(uv_hori, 1)    == 1      uv_hori    = uv_hori';          end
    if isempty(cbar_location)        cbar_location = 'eastoutside';  end
%     if size(line_hori, 1)  == 1      line_hori  = line_hori';        end
    
    %/ make contf_hori strictly increasing (e.g., lat)
    if any(unique(diff(contf_hori)) < 0)     contf_hori = flip(contf_hori);  contf_data = flip(contf_data, 1);                                         end
    if any(unique(diff(stipp_hori)) < 0)     stipp_hori = flip(stipp_hori);  stipp_data = flip(stipp_data, 1);                                         end
    if any(unique(diff(cont_hori))  < 0)     cont_hori  = flip(cont_hori);   cont_data  = flip(cont_data, 1);  cont_data_raw = flip(cont_data_raw, 1); end
    if any(unique(diff(uv_hori))    < 0)     uv_hori    = flip(uv_hori);     Udata      = flip(Udata, 1);      Vdata = flip(Vdata, 1);                 end

    if isempty(gcf_position)   
        if zm_or_mm == 1 && always_time_as_y == 0
            gcf_position = [200 50 1225 400];   % x,y,w,h; sets figure siz, landscape if always_time_as_y == 0 
        else
            gcf_position = [1400 100 525 1225]; % x,y,w,h; default
        end
    end   
    str_cbar = []; str_qscale = []; str_vecref = []; res_x = 0; res_y = 0; 
    contains_transparent = 0;

    %=============
    figure
    set(gcf, 'color','w');
    set(gcf,'position', gcf_position) % x,y,w,h; 
    ax_ori = gca;
    
    %======== contf  ========%    
    % if isequal(cbar_location, 'eastoutside')
    %     cb_fontsize     = fontsize*0.5;
    % else
    cb_fontsize     = fontsize*0.5;
    % end
                
    if draw_refvec_only    %/ if draw_refvec_only, we plot vectors only.
        contf_data   = nan(length(contf_hori), length(contf_time));
        cont_data    = [];
        scatter_data = [];
        line_data    = [];
        border_lines = [];
        time_lines   = [];
    end  
    
    if ~isempty(contf_data)
        if isempty(contf_unit)
            contf_unit      = '';
        end
        if draw_cbar_only
            str_cbar    = '_cbar';
            set(ax_ori, 'visible', 'off');
            
            clim([min(contf_levels) max(contf_levels)]);
            colormap(ax_ori, colmap)
            % drawnow; pause(0.05);
            
            cb              = colorbar;
            cbar_YTick      = contf_levels(2:cbar_interval:end-1);
            cbar_YTickLabel = round(cbar_YTick,2);

            set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', cb_fontsize,...
                'Location', cbar_location)   
            set(get(cb,'Title'),'String',contf_unit, 'Fontsize', cb_fontsize)
%             cb.Position = cb.Position + 1e-10;  
%             pos = get(cb, 'Position');
%             if isequal(cbar_location, 'southoutside')
%                 pos(2) = pos(2)*5;
%                 set(cb, 'Position', pos);
%             elseif isequal(cbar_location, 'eastoutside')
%                 pos(1) = pos(1)*0.8;
%                 set(cb, 'Position', pos);
%             end
        else
            res_x = unique(round(diff(contf_hori), 3));  %/ round up to 3 decimal points to avoid numerical errors
            res_y = unique(round(diff(contf_time), 3));
            if pcolor_mode
                h = pcolor(contf_hori-res_x/2, contf_time-res_y/2, contf_data');     %/ shift by half grid to correct the bug in pcolor.
                set(h, 'edgecolor','none');
                shading flat
            else
%                 warning('contourf may not show the discontinued grids.');
                contourf(contf_hori, contf_time, contf_data', [-99999999, contf_levels], 'edgecolor', 'none') %/ IMPORTANT: input contf_levels, otherwise the results are completely wrong!!
                shading interp
            end
            
            if cbar_mode
                cb              = colorbar;
                cbar_YTick      = contf_levels(2:cbar_interval:end-1);
                cbar_YTickLabel = round(cbar_YTick,2);

                set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', cb_fontsize, 'Location', cbar_location)   
                set(get(cb,'Title'),'String',contf_unit, 'Fontsize', cb_fontsize)
                
                cb.Position = cb.Position + 1e-10;  %/ IMPORTANT: First stop auto resizing by setting cb pos. 
                
                %/ Fine-tune cb pos
                if isequal(cbar_location, 'eastoutside')
                    %/ Shift it eastward by a distance 3 times its width.
                    cb.Position(1) = cb.Position(1) + cb.Position(3)*3;  %/ x, y, width, height
                    
                elseif isequal(cbar_location, 'southoutside')
                    %/ Shift it southward by a distance 3 times its height.
                    cb.Position(2) = cb.Position(2) - cb.Position(4)*4;  %/ x, y, width, height
                end
            end
            
            clim([min(contf_levels) max(contf_levels)]);
            colormap(ax_ori, colmap)
            % drawnow; pause(0.05);
        end
    end
    hold on;
    

    if draw_cbar_only == 0  %/ only cont, vector only if draw_cbar_only = 0.

        %======== cont ========%
        if ~isempty(cont_data)
            if isempty(cont_labelsize)      cont_labelsize = fontsize*0.5;     end

            for cc = 1:length(cont_levels)
                %/ skip drawing the zero level.
                if cont_levels(cc) == 0 
                    if skip_zero_cont
                        continue;  
                    else
                        a = 1.5;  %/ thicken the zero-contour line
                    end
                else
                    a = 1;
                end
                cont_color = cont_colmap(cc,:);
                if numel(cont_linewi) ~= 1
                    cont_wi = cont_linewi(cc);
                else
                    cont_wi = cont_linewi;
                end

                %/ show raw data thinner contours.
                if ~isempty(cont_data_raw)
                    [~, ~] = contour(cont_hori, cont_time, cont_data_raw', [cont_levels(cc), cont_levels(cc)],...
                                     'LineWidth', cont_wi*0.3*a, 'LineStyle','--', 'Color',cont_color);
                    hold on;
                    cont_linest_bc = '-'; 
                
                elseif isequal(cont_linest, 'auto') || isempty(cont_linest)
                    if cont_levels(cc) < 0
                        cont_linest_bc = '--'; %/ Normally, use dashed contours to indicate -ve values
                    else
                        cont_linest_bc = '-';
                    end
                else
                    cont_linest_bc = cont_linest;
                end
                
                [CS_cont, h_cont] = contour(cont_hori, cont_time, cont_data', [cont_levels(cc), cont_levels(cc)],...
                                            'LineWidth', cont_wi*a, 'LineStyle', cont_linest_bc, 'Color',cont_color);
                hold on;

                if not_show_cont_label == 0
                    clabel(CS_cont, h_cont, 'Color', 'k', 'FontSize',cont_labelsize,'FontWeight','bold', 'labelspacing', 1e4); %, 'BackgroundColor',[1 1 1]);
%                     clabel(CS_cont, 'Color', 'k', 'FontSize',cont_labelsize,'FontWeight','bold', 'labelspacing', 2000, 'Rotation', 0); %, 'BackgroundColor',[1 1 1]);
                    hold on; %/ if not hold on, then changes in clabel() will not be retained.
                end
            end
        end

        %======== Stippling (put it after cont to avoid being blocked =======%
        if ~isempty(stipp_data)
            %/ scatter (Simple, complete and error-free)
            [ind_hori, ind_time] = find(~isnan(stipp_data));

            stipp_x = stipp_hori(ind_hori);
            stipp_y = stipp_time(ind_time);

            if isempty(stipp_markersize) stipp_markersize = 3;                   end
            if isempty(stipp_color)      stipp_color      = [255 255 255]./255;  end
            if isempty(stipp_marker)     stipp_marker     = 'o';                 end
            
            % stipp_markersize = 3;
            % stipp_color      = [0 0 0]./255;
            % stipp_marker     = 'o';

            % h_sc = scatter(stipp_x, stipp_y, stipp_markersize, stipp_color, stipp_marker);
            h_sc = scatter(stipp_x, stipp_y, stipp_markersize, stipp_color, 'filled',  stipp_marker);
            set(h_sc, 'linewidth', 1);
            hold on;

            %/ hatching (WARNING: hatching will miss all isolated points, and is error-prone. Do NOT use it...)
%             hatch_data = stipp_data;
%             hatch_data(~isnan(hatch_data)) = 1;   %/ set whatever sig values (both +ve and -ve) to 1.
% %             hatch_data( isnan(hatch_data)) = 0;
%             
%             [~, h2] = contourf(stipp_hori, stipp_time, hatch_data', [1 1]); %/ IMPORTANT: input contf_levels, otherwise the results are completely wrong!!
%             set(h2,'linestyle','none','Tag','HatchingRegion');
%             hold off;                                 % if you want to have more control
%             ax1 = gca;
% 
%             hp = findobj(ax1,'Tag','HatchingRegion');
%             hh = hatchfill2(hp, 'cross','Fill','off', 'HatchDensity', 100, 'HatchLineWidth', 1, 'HatchColor', [0 0 0]);
        end

        %======== Simple Linear Regression based on contf_data or cont_data ========%
        ax_SLR = []; str_SLR = ''; 
        if ~isempty(fit_SLR_on) && ~isequal(fit_SLR_on, 'off')
            ax_SLR = axes;                             %/ IMPORTANT: axes is to create a new axes!
            set(ax_SLR, 'Visible', 'off')              %/ IMPORTANT: hide the 2nd axes 
            hold on;

            if isempty(SLR_rmval)      
                warning('Detected ''SLR_rmval'' is empty. No values will be removed before the fitting.'); 
            end
            if isempty(SLR_min_or_max)   
                error('Specify ''SLR_min_or_max'' to indicate your choice to fit the regression based on min or max value!'); 
            elseif ~isequal(SLR_min_or_max, 'min') && ~isequal(SLR_min_or_max, 'max')   
                error('Invalid ''SLR_min_or_max''! Input ''min'' or ''max''.');
            end
            if isempty(SLR_largest_obj)  
                SLR_largest_obj = 1;       %/ By default, fit on the longest continuous segment (if contains NaNs)
            end   
            
            %/==== Data processing based on the input criteria ====%
            if isequal(fit_SLR_on, 'contf')
                if zm_or_mm == 2 && ~isempty(SLR_fitting_range_lon)
                    ind = find(contf_hori >= SLR_fitting_range_lon(1) & contf_hori <= SLR_fitting_range_lon(end));
                    horidata4SLR = contf_hori(ind);
                    timedata4SLR = contf_time;
                    Zdata4SLR    = contf_data(ind,:);  %/ hori * time
                else
                    horidata4SLR = contf_hori;
                    timedata4SLR = contf_time;
                    Zdata4SLR    = contf_data;  %/ hori * time
                end 
            elseif isequal(fit_SLR_on, 'cont')
                if zm_or_mm == 2 && ~isempty(SLR_fitting_range_lon)
                    ind = find(cont_hori >= SLR_fitting_range_lon(1) & cont_hori <= SLR_fitting_range_lon(end));
                    horidata4SLR = cont_hori(ind);
                    timedata4SLR = cont_time;
                    Zdata4SLR    = cont_data(ind,:);  %/ hori * time
                else
                    horidata4SLR = cont_hori;
                    timedata4SLR = cont_time;
                    Zdata4SLR    = cont_data;  %/ hori * time
                end
            else
                error('Set fit_SLR_on = ''contf'', ''cont'' or []!');
            end

            %/ 1. Replace values greater than SLR_rmval with nan
            if ~isempty(SLR_rmval)
                if isequal(SLR_min_or_max, 'min')
                    Zdata4SLR(Zdata4SLR > SLR_rmval) = nan; %/ 2D, hori x time
                elseif isequal(SLR_min_or_max, 'max')
                    Zdata4SLR(Zdata4SLR < SLR_rmval) = nan; %/ 2D, hori x time
                end
            end

            %/ 2. Find the largest object by bwareafilt
            if SLR_largest_obj
                logi_data = Zdata4SLR;  %/ hori * time
                logi_data(~isnan(logi_data)) = 1;
                logi_data(isnan(logi_data)) = 0;
                logi_data = logical(logi_data);

                if isempty(SLR_largest_obj_lon) && isempty(SLR_largest_obj_lat)
                    BW2 = bwareafilt(logi_data, 1, 'largest', 8); %/ Get the largest object by default. Set connectivity = 8
                else
                    N   = 20;   %/ dummy
                    BW2 = false(size(logi_data)); %/ Initialization (useful when no object at all based on the criteria)
                    for n = 1:N
                        BW2_1ton = bwareafilt(logi_data, n, 'largest', 8); %/ Get n number of largest object by default. Set connectivity = 8
                        if n == 1
                            BW2_1tonm1 = false(size(logi_data)); %/ dummy
                        else
                            BW2_1tonm1 = bwareafilt(logi_data, n-1, 'largest', 8); %/ Get n-1 number of largest object by default. Set connectivity = 8
                        end
                        BW2      = (BW2_1ton & ~BW2_1tonm1);   %/ "substract" to get the real n-th largest object
                        % contf_time
                        ind_t    = contf_time >= SLR_largest_obj_time(1) & contf_time <= SLR_largest_obj_time(end);
                        BW2_t    = BW2(:,ind_t);  
                        obj_hori = horidata4SLR(logical(BW2_t));
                        if zm_or_mm == 1
                            b1 = length(find(obj_hori     >= SLR_largest_obj_lat(1) & obj_hori     <= SLR_largest_obj_lat(end)));
                            b2 = length(find(horidata4SLR >= SLR_largest_obj_lat(1) & horidata4SLR <= SLR_largest_obj_lat(end)));
                        elseif zm_or_mm == 2
                            b1 = length(find(obj_hori     >= SLR_largest_obj_lon(1) & obj_hori     <= SLR_largest_obj_lon(end)));
                            b2 = length(find(horidata4SLR >= SLR_largest_obj_lon(1) & horidata4SLR <= SLR_largest_obj_lon(end)));
                        end
                        occupy_ratio = b1/b2;
                        if occupy_ratio > 0.5
                            break;  %/ if occupying more than half of the hori range, break the loop and select the object.
                        end
                    end
                end
                Zdata4SLR(~BW2) = nan;  %/ Remove other all smaller objects
            end

            if size(horidata4SLR, 1) == 1  horidata4SLR = horidata4SLR'; end
            if size(timedata4SLR, 1) == 1  timedata4SLR = timedata4SLR'; end

            %/ 3. Find hori positions of min/max along the *hori* dim
            %/    min(A,[],2) returns a column vector containing the minimum value of each row.
            if isequal(SLR_min_or_max, 'min')
                [B, I] = min(Zdata4SLR, [], 2, 'omitnan');
            elseif isequal(SLR_min_or_max, 'max')
                [B, I] = max(Zdata4SLR, [], 2, 'omitnan');
            end
            timedata4SLR = timedata4SLR(I);
            timedata4SLR(isnan(B)) = nan;             %/ As min/max() will output a dummy index position of 1 even if the vector contains nan only. We must distinguish it by nan.
            horidata4SLR(isnan(timedata4SLR)) = nan;  %/ For consistency

            %/ 4. Perform Simple Linear Regression (with t-test or MK test)
            if isempty(SLR_test)   SLR_test = 'MK';   end
            if isempty(SLR_alpha)  SLR_alpha = 0.05;  end

            %/ Sort the xdata before computing trends
            [horidata4SLR, I] = sort(horidata4SLR, 'ascend'); %/ [nan nan 5 10 3 nan]  -> [3     5    10   NaN   NaN]
            timedata4SLR = timedata4SLR(I); %/ because xdata has been sorted, we need to sort ydata accordingly
            flag_small_obj = 0;
            if SLR_recur == 1
                ind_RX = find(horidata4SLR >= recur_min_RX); %/ recursively fitting SLR until X is shortened to a certain point
                ind_RX = flip(ind_RX);                       %/ Shift the rightmost boundary from right to left

                ind_LX = find(horidata4SLR <= recur_max_LX); %/ Shift the leftmost boundary from left to right
                if isempty(ind_RX)
                    flag_small_obj = 1;
                    warning('ind_till_X is empty! Either the target does not propagate to recur_min_RX or invalid recur_min_RX.');
                end
            else
                ind_RX = length(horidata4SLR); %/ just use all points
                ind_LX = 1;                    %/ just use all points
            end
            if isempty(recur_Rsquared)
                recur_Rsquared = 0.8;
            end
            % if isempty(recur_RMSE)
            %     recur_RMSE = inf;  %/ if not given, then regardless of RMSE associated with the fitting
            % end

            pval_data = nan;
            trend_data_temp     = nan(length(ind_LX), length(ind_RX));
            pval_data_temp      = nan(length(ind_LX), length(ind_RX));
            intercept_data_temp = nan(length(ind_LX), length(ind_RX));
            Rsquared_temp       = nan(length(ind_LX), length(ind_RX));
            RMSE_temp           = nan(length(ind_LX), length(ind_RX));
            for i = 1:length(ind_LX)
                % break_outer_loop = 0;
                for j = 1:length(ind_RX)   
                    horidata4SLR_recur = horidata4SLR(ind_LX(i):ind_RX(j));         
                    timedata4SLR_recur = timedata4SLR(ind_LX(i):ind_RX(j));         
    
                    %/ NOTE: 'data_yearly' -> Y 
                    %/         'year_list' -> X
                    [trend_data, ~, ~, pval_data, intercept_data, ~, Rsquared_data, Rsquared_adj_data, RMSE_data] = ...
                            compute_trend('data_yearly', timedata4SLR_recur', 'year_list', horidata4SLR_recur', 'trend_test', SLR_test, 'trend_alpha', SLR_alpha);

                    if SLR_recur == 1
                        fprintf('Recursive fitting (leftmost: %2.fE rightmost: %.2fE)...\n', horidata4SLR(ind_LX(i)), horidata4SLR(ind_RX(j)));
                        % 
                        % if (Rsquared_data >= recur_Rsquared) %&& (RMSE_data <= recur_RMSE)
                        %     horidata4SLR = horidata4SLR_recur; %/ Update
                        %     % timedata4SLR = timedata4SLR_recur; %/ not used later
                        % 
                        %     break_outer_loop = 1; %/ This is to break the outer loop
                        %     fprintf('Recursive fitting stopped: R-squared = %.2f (>= %.2f), RMSE = %.2f\n', Rsquared_data, recur_Rsquared, RMSE_data);
                        %     break;
                        % else
                        trend_data_temp(i,j)     = trend_data;
                        pval_data_temp(i,j)      = pval_data;
                        intercept_data_temp(i,j) = intercept_data;
                        Rsquared_temp(i,j)       = Rsquared_data;  
                        RMSE_temp(i,j)           = RMSE_data;  

                        if i == length(ind_LX) && j == length(ind_RX)  
                            %/ 1. find those with R-2 > recur_Rsquared
                            %/ 2. select the one with the minimum RMSE
                            
                            [ind_row, ind_col] = find(Rsquared_temp >= recur_Rsquared);
                            if isempty(ind_row)
                                % [ind_row, ind_col]  = find(max(Rsquared_temp, [], 'all'));
                                pval_data = nan; %/ Set pval_data to nan to not show the fitting
                                warning('No Rsquared meets the criteria (%.2f), no fitting will be shown', recur_Rsquared);
                            else

                                RMSE_temp_list = nan(length(ind_row), 1);
                                for ii = 1:length(ind_row)
                                    RMSE_temp_list(ii) = RMSE_temp(ind_row(ii), ind_col(ii));
                                end
                                [~, I] = min(RMSE_temp_list);
                                ind_row = ind_row(I);  %/ Update
                                ind_col = ind_col(I);  %/ Update
                                
                                trend_data          = trend_data_temp(ind_row, ind_col);  
                                pval_data           = pval_data_temp(ind_row, ind_col);  
                                intercept_data      = intercept_data_temp(ind_row, ind_col);  
                                Rsquared_data       = Rsquared_temp(ind_row, ind_col);  
                                RMSE_data           = RMSE_temp(ind_row, ind_col);  
                                horidata4SLR        = horidata4SLR(ind_LX(ind_row):ind_RX(ind_col)); %/ Update
                                
                                % if isempty(ind_row)
                                    % warning('No Rsquared meets the criteria (%.2f), then we choose the largest one (%.2f) (leftmost: %2.fE rightmost: %.2fE)', recur_Rsquared, Rsquared_data, horidata4SLR(ind_LX(ind_row)), horidata4SLR(ind_RX(ind_col)));
                                % else
                                fprintf('Recursive fitting stopped: R-squared = %.2f (>= %.2f), RMSE = %.2f\n', Rsquared_data, recur_Rsquared, RMSE_data);
                            end
                        end
                        % end
                    end
                end
                % if break_outer_loop
                    % break;
                % end
            end

            %/ [Optional] Use asterisks to indicate significance level 
            asterisk    = {'(*)', '(**)', '(***)'};    
            alpha_level = [0.05,    0.01,   0.001];    
            str_sig = repmat({''}, length(pval_data), 1);
            for i = 1:length(alpha_level)
                ind = find(pval_data <= alpha_level(i));
                str_sig(ind) = repmat(asterisk(i), length(ind), 1);
            end

            if (flag_small_obj == 0) && (pval_data <= SLR_alpha)
                draw_SLR_line = 1;
            else
                draw_SLR_line = 0;
            end

            %/ Fit a simple linear regression line (if the column is not entirely NaN, or if the trend is significant)
            %/ If not significant, regress_x may be unrealistic.
            if draw_SLR_line
                %/ Time conversion to daily
                if isequal(time_unit, 'daily')
                    conv2daily = 1;
                elseif isequal(time_unit, 'pentad')
                    conv2daily = 1/5;
                elseif isequal(time_unit, 'monthly')
                    conv2daily = 1/30;
                end
    
                %/ Distance conversion
                a = 6371e3;                      %/ Earth's radius (m)
                if zm_or_mm == 1
                    deg2m = a*pi/180;            %/ 1 deg of latitude is always ~111 km
                elseif zm_or_mm == 2
                    mean_lat = mean(lat_range);  %/ Mean lat of the meridional mean
                    deg2m    = a*cosd(mean_lat)*pi/180;
                end
                cp = 1./trend_data * deg2m * conv2daily /24/3600;      %/ Propagation speed (day/deg -> m/s)

                SLR_linest = '-';    %/ significant
                if isempty(SLR_linewi)       SLR_linewi = linewidth*2;    end
                if isempty(SLR_color)        SLR_color = [0 0 0];         end

                %/============================================================================================
                %/ Since we fit the regression with X = time, Y = lon,
                %/ while the diagram shows x = lon, y = time, and also because we take min/max along time dim,
                %/ we have a fuller spectrum of hori data (Y) to *retrospectively* deduce the time (X).
                %/ Finally, we swtich the two axis when calling plot().
                %/============================================================================================
                ind_nonnan = nan(2,1);
                ind_nonnan(1) = find(horidata4SLR == min(horidata4SLR, [], 'omitnan'), 1, 'first'); %/ Min values can be multiple. Choose the first one.
                ind_nonnan(2) = find(horidata4SLR == max(horidata4SLR, [], 'omitnan'), 1, 'first'); %/ Max values can be multiple. Choose the first one.
                regress_x = horidata4SLR(ind_nonnan([1,end]));
                regress_y = trend_data.*regress_x + intercept_data;
                plot(regress_x, regress_y, 'color', SLR_color, 'linewidth', SLR_linewi, 'linestyle', SLR_linest, 'marker', 'none');
                hold on;

                if ~isempty(cp_error)
                    str_SLR = sprintf(' R2: %.2f, RMSE: %.2f, cp: %.1f%s%.1f m/s', Rsquared_data, RMSE_data, cp, char(177), cp_error); 
                else
                    str_SLR = sprintf(' R2: %.2f, RMSE: %.2f, cp: %.1f m/s', Rsquared_data, RMSE_data, cp); 
                end
                % cp_CI95        = CI95_data * deg2m * conv2daily /24/3600;       %/ Propagation speed uncertainty (m/s)
                % str_SLR_contf = sprintf('_%.1f%s%.1f m/s', cp, char(177), cp_CI95); 
            else
                warning('Insignificant SLR slope. No drawing of it. cp will be nan.');
            end
            linkaxes([ax_ori, ax_SLR], 'xy');                  %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
        end

        %======== scatter plot ========%
        ax_scatter = [];
        if ~isempty(scatter_data)
            ax_scatter = axes;                             %/ IMPORTANT: axes is to create a new axes!
            ax_scatter.Clipping = 'off';                   %/ IMPORTANT: turn clipping off. (allow to draw a point outside the axes!
            set(ax_scatter, 'Visible', 'off')              %/ IMPORTANT: hide the 2nd axes 
            hold on;

            scatter_linewi     = linewidth*0.7;

%             scatter_edgecolor  = [255 51 204]./255;
%             scatter_facecolor  = scatter_edgecolor;

%             scatter_edgecolor  = [153 153 153]./255;
%             scatter_facecolor  = [153 153 153]./255;

            if isempty(scatter_edgecolor)  scatter_edgecolor  = [255 255 255]./255;   end
            if isempty(scatter_facecolor)  scatter_facecolor  = scatter_edgecolor;    end
            if isempty(scatter_alpha)      scatter_alpha      = 0.4;                  end

            if scatter_alpha == 1
                contains_transparent = 0;
            else
                contains_transparent = 1;
            end

            if ~isempty(scatter_size)
                %/ Normalize the size
                scatter_size = (scatter_size - min(scatter_size))/(max(scatter_size) - min(scatter_size));
                scatter_size = 50 + round(1000*scatter_size);
                scatter_size(isnan(scatter_size)) = 1.e-15;
                scatter_size(scatter_size == 0)   = 1.e-15;  %/ since size must not be zero or contain NaN.
            else
                scatter_size = 300;
            end
%                 plot(scatter_data, 1:length(scatter_data), 'marker', 'none', 'color', track_color,...
%                      'linestyle', '-', 'linewidth', linewidth*1.5);

            sc = scatter(scatter_data, 1:length(scatter_data), scatter_size, 'MarkerFaceColor', scatter_facecolor, 'MarkerEdgeColor', scatter_edgecolor); 
            sc.MarkerFaceAlpha = scatter_alpha;
            sc.LineWidth       = scatter_linewi;
            linkaxes([ax_ori, ax_scatter], 'xy');                  %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
        end
        
        %======== 1D line plot ========%
        ax_line = [];
        if ~isempty(line_data)
            ax_line = axes;                             %/ IMPORTANT: axes is to create a new axes!
            ax_line.Clipping = 'off';                   %/ IMPORTANT: turn clipping off. (allow to draw a point outside the axes!
            set(ax_line, 'Visible', 'off')              %/ IMPORTANT: hide the 2nd axes 
            hold on;

            if isempty(line_linewi)    line_linewi = repmat(linewidth, size(line_data,1), 1);  end
            
            for k = 1:size(line_data,1)
                plot(line_data(k,:), line_time, 'marker', 'none', 'color', line_color(k,:),...
                     'linestyle', '-', 'linewidth', line_linewi(k));
                hold on;
            end
            linkaxes([ax_ori, ax_line], 'xy');                  %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
        end
        
        %======== border / time lines =======%
        %/ Do NOT put this after set(gca), or will cause annoying margins of contf!
        ax_border = [];
        if ~isempty(border_lines) || ~isempty(time_lines)
            ax_border = axes;                                             %/ IMPORTANT: axes is to create a new axes!
            set(ax_border, 'Visible', 'off')                              %/ IMPORTANT: hide the 2nd axes 
            hold on;
%             bar(topo_hori, topo, 'Barwidth', 1, 'edgecolor', 'none', 'facecolor', [.6 .6 .6], 'FaceAlpha', 1); %'facecolor', [.3 .3 .3], 'FaceAlpha', 0.9);

            if isempty(border_color)       border_color = 'k';           end
            if isempty(border_linestyle)   border_linestyle = '--';      end
            if isempty(border_linewidth)   border_linewidth = linewidth; end
            for i = 1:length(border_lines)
                xline(border_lines(i), 'color', border_color, 'linestyle', border_linestyle, 'linewidth', border_linewidth);
            end
            hold on;
            for i = 1:length(time_lines)
                yline(time_lines(i),   'color', border_color, 'linestyle', border_linestyle, 'linewidth', border_linewidth);
            end
            hold on;
            linkaxes([ax_ori ax_border], 'xy');                                 %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
        end
        
        %======== vector ========%
        ax_vec = [];
        if ~isempty(Udata) && ~isempty(Vdata)
            ax_vec = axes;                             %/ IMPORTANT: axes is to create a new axes!
            ax_vec.Clipping = 'on';                   %/ IMPORTANT: turn clipping off. (allow to draw a point outside the axes!
            set(ax_vec, 'Visible', 'off')              %/ IMPORTANT: hide the 2nd axes 
            hold on;
            
            [uv_hori_2D, uv_time_2D] = meshgrid(uv_hori, uv_time);  %/ (time x hori) -> as required by quiver().
            if size(Udata, 1) == length(uv_hori)    Udata = Udata';   end
            if size(Vdata, 1) == length(uv_hori)    Vdata = Vdata';   end
%             
%             warning('testing the correctness of UV vector..')
%             Udata(:,:) = 7;  
%             Vdata(:,:) = 7; 
            
            %/ IMPORTANT: quiver() somehow does NOT follow view(90) ->
            %             rotate counterclockwise. We need to correct it manually.
            if zm_or_mm == 1  && always_time_as_y == 0
                Udata_bc = Vdata;  %/ swtich u and v.
                Vdata_bc = Udata;
            else
                Udata_bc = Udata;  %/ keep it the same
                Vdata_bc = Vdata;
            end

            if ~isempty(U2data) && ~isempty(V2data)
                vector_color_bc = [.4 .4 .4];  %/ If vector 2 exists, make vector 1 grey
            else
                vector_color_bc = vector_color;
            end
            
            %/ Check if plotting the vector scale only.
            if draw_refvec_only
                str_vecref = '_vecref';
                disp(vecref_mag)

                X = uv_hori_2D;
                Y = uv_time_2D;
                %/ [IMPORTANT]: Make a U / V matrix with the reference vector and zeros elsewhere to get the correct quiver scale!
                U = zeros(size(Udata_bc)); 
                V = zeros(size(Vdata_bc)); 
                
                if size(vector_color_bc, 1) > 1 %/ if for multiple colorbars.
                    error('Somehow I cannot use m_quiver to plot consistent scale for the reference vector. Consider using a colorbar to indicate the magnitude of the colored vectors.')
                    
                    if size(vector_color_bc, 1) ~= length(vecref_mag)
                        error('Inconsistent data length between vector_color_bc and vecref_mag!');
                    end
%                     ind_hori_mid = round(length(uv_hori)/2);
%                     ind_time_mid = round(linspace(1, length(uv_time), size(vector_color_bc, 1)*2 ));

                    ind_hori_mid = round(length(uv_hori)/2);
                    ind_time_mid = round(linspace(floor(length(uv_time))*0.1, floor(length(uv_time))*0.9, size(vector_color_bc, 1)));   
                else
                    ind_hori_mid = round(length(uv_hori)/2);
                    ind_time_mid = round(length(uv_time)/2);
                end
                
                %/ IMPORTANT: quiver() somehow does NOT follow view(90) ->
                %             rotate counterclockwise. We need to correct it manually.
                for i = 1:length(ind_hori_mid)
                    if zm_or_mm == 1  && always_time_as_y == 0
                        V(ind_time_mid, ind_hori_mid(i)) = vecref_mag;
                    else
                        U(ind_time_mid, ind_hori_mid(i)) = vecref_mag;
                    end

                    text(X(ind_time_mid+1,ind_hori_mid+1), Y(ind_time_mid+1,ind_hori_mid+1),...
                         string(vecref_mag), 'FontSize', fontsize, 'verticalalignment','bottom',...
                         'horizontalalignment','left');
                    hold on; 
                end
            else
                X = uv_hori_2D(2:vector_step_time:end-1, 2:vector_step_hori:end-1);
                Y = uv_time_2D(2:vector_step_time:end-1, 2:vector_step_hori:end-1);
                U = Udata_bc  (2:vector_step_time:end-1, 2:vector_step_hori:end-1);
                V = Vdata_bc  (2:vector_step_time:end-1, 2:vector_step_hori:end-1);
            end

            %/ If vector_levels is given, draw vectors with different colors to indicate the magnitude level.
            if ~isempty(vector_levels)  
                if length(vector_levels) ~= size(vector_color_bc, 1)
                    error('The length of vector_levels has to be consistent with the rows of vector_color_bc!');
                end
                
                UVmag = sqrt(U.^2 + V.^2);
                for k = 1:length(vector_levels)

                    %/ a faster way (should be)
                    if k == length(vector_levels)
                        cond = (vector_levels(k) <= UVmag); %/ UVmag >= upper bound.
                    else
                        cond = (vector_levels(k) <= UVmag & UVmag < vector_levels(k+1)); %/ eg. 0 <= UVmag < 20.
                    end
                    X_bc = X;  X_bc(~cond) = nan;  
                    Y_bc = Y;  Y_bc(~cond) = nan;  
                    U_bc = U;  U_bc(~cond) = nan;  
                    V_bc = V;  V_bc(~cond) = nan;  
                    h_quiver = quiver(X_bc, Y_bc, U_bc, V_bc, 0, 'color', vector_color_bc(k,:), 'linewidth', vector_linewidth); %, 'Clipping', 'on');
                    hU = get(h_quiver,'UData') ;
                    hV = get(h_quiver,'VData') ;
                    set(h_quiver,'UData',qscale*hU,'VData',qscale*hV) 
                    hold on; %drawnow;
                end
            else
                h_quiver = quiver(X, Y, U, V, 0, 'color', vector_color_bc, 'linewidth', vector_linewidth);
                hU = get(h_quiver,'UData') ;
                hV = get(h_quiver,'VData') ;
                set(h_quiver,'UData',qscale*hU,'VData',qscale*hV) 
                hold on; %drawnow;
            end
            
            fprintf('*** Spatial avg of time-mean Udata = %.2f ***\n', my_nanmean(Udata, 'all'));
            fprintf('*** Spatial avg of time-mean Vdata = %.2f ***\n', my_nanmean(Vdata, 'all'));
            fprintf('*** Max of time-mean Udata = %.2f ***\n', max(Udata, [], 'all'));
            fprintf('*** Max of time-mean Vdata = %.2f ***\n', max(Vdata, [], 'all'));
            linkaxes([ax_ori, ax_vec], 'xy');                  %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
        end

        %======== vector2 ========%
        ax_vec2 = []; 
        % U2data = []; V2data = []; %/ testing
        if ~isempty(U2data) && ~isempty(V2data)
            ax_vec2 = axes;                             %/ IMPORTANT: axes is to create a new axes!
            ax_vec2.Clipping = 'on';                   %/ IMPORTANT: turn clipping off. (allow to draw a point outside the axes!
            set(ax_vec2, 'Visible', 'off')              %/ IMPORTANT: hide the 2nd axes 
            hold on;
            
            [uv_hori_2D, uv_time_2D] = meshgrid(uv_hori, uv_time);  %/ (time x hori) -> as required by quiver().
            if size(U2data, 1) == length(uv_hori)    U2data = U2data';   end
            if size(V2data, 1) == length(uv_hori)    V2data = V2data';   end
            
            %/ IMPORTANT: quiver() somehow does NOT follow view(90) ->
            %             rotate counterclockwise. We need to correct it manually.
            if zm_or_mm == 1  && always_time_as_y == 0
                U2data_bc = V2data;  %/ swtich u and v.
                V2data_bc = U2data;
            else
                U2data_bc = U2data;  %/ keep it the same
                V2data_bc = V2data;
            end
            vector_color_bc = vector_color;
            
            %/ Check if plotting the vector scale only.
            if draw_refvec_only
                str_vecref = '_vecref';
                disp(vecref_mag)

                X = uv_hori_2D;
                Y = uv_time_2D;
                %/ [IMPORTANT]: Make a U / V matrix with the reference vector and zeros elsewhere to get the correct quiver scale!
                U = zeros(size(U2data_bc)); 
                V = zeros(size(V2data_bc)); 

                if size(vector_color_bc, 1) > 1 %/ if for multiple colorbars.
                    error('Somehow I cannot use m_quiver to plot consistent scale for the reference vector. Consider using a colorbar to indicate the magnitude of the colored vectors.')
                    
                    if size(vector_color_bc, 1) ~= length(vecref_mag)
                        error('Inconsistent data length between vector_color_bc and vecref_mag!');
                    end
%                     ind_hori_mid = round(length(uv_hori)/2);
%                     ind_time_mid = round(linspace(1, length(uv_time), size(vector_color_bc, 1)*2 ));

                    ind_hori_mid = round(length(uv_hori)/2);
                    ind_time_mid = round(linspace(floor(length(uv_time))*0.1, floor(length(uv_time))*0.9, size(vector_color_bc, 1)));   
                else
                    ind_hori_mid = round(length(uv_hori)/2);
                    ind_time_mid = round(length(uv_time)/2);
                end
                
                %/ IMPORTANT: quiver() somehow does NOT follow view(90) ->
                %             rotate counterclockwise. We need to correct it manually.
                for i = 1:length(ind_hori_mid)
                    if zm_or_mm == 1  && always_time_as_y == 0
                        V(ind_time_mid, ind_hori_mid(i)) = vecref_mag;
                    else
                        U(ind_time_mid, ind_hori_mid(i)) = vecref_mag;
                    end

                    text(X(ind_time_mid+1,ind_hori_mid+1), Y(ind_time_mid+1,ind_hori_mid+1),...
                         string(vecref_mag), 'FontSize', fontsize, 'verticalalignment','bottom',...
                         'horizontalalignment','left');
                    hold on; 
                end
            else
                X = uv_hori_2D(2:vector_step_time:end-1, 2:vector_step_hori:end-1);
                Y = uv_time_2D(2:vector_step_time:end-1, 2:vector_step_hori:end-1);
                U = U2data_bc  (2:vector_step_time:end-1, 2:vector_step_hori:end-1);
                V = V2data_bc  (2:vector_step_time:end-1, 2:vector_step_hori:end-1);
            end

            %/ If vector_levels is given, draw vectors with different colors to indicate the magnitude level.
            if ~isempty(vector_levels)  
                if length(vector_levels) ~= size(vector_color_bc, 1)
                    error('The length of vector_levels has to be consistent with the rows of vector_color_bc!');
                end
                
                UVmag = sqrt(U.^2 + V.^2);
                for k = 1:length(vector_levels)
                    %/ a faster way (should be)
                    if k == length(vector_levels)
                        cond = (vector_levels(k) <= UVmag); %/ UVmag >= upper bound.
                    else
                        cond = (vector_levels(k) <= UVmag & UVmag < vector_levels(k+1)); %/ eg. 0 <= UVmag < 20.
                    end
                    X_bc = X;  X_bc(~cond) = nan;  
                    Y_bc = Y;  Y_bc(~cond) = nan;  
                    U_bc = U;  U_bc(~cond) = nan;  
                    V_bc = V;  V_bc(~cond) = nan;  
                    h_quiver = quiver(X_bc, Y_bc, U_bc, V_bc, 0, 'color', vector_color_bc(k,:), 'linewidth', vector_linewidth); %, 'Clipping', 'on');
                    hU = get(h_quiver,'UData') ;
                    hV = get(h_quiver,'VData') ;
                    set(h_quiver,'UData',qscale*hU,'VData',qscale*hV) 
                    hold on; %drawnow;
                end
            else
                h_quiver = quiver(X, Y, U, V, 0, 'color', vector_color_bc, 'linewidth', vector_linewidth);
                hU = get(h_quiver,'UData') ;
                hV = get(h_quiver,'VData') ;
                set(h_quiver,'UData',qscale*hU,'VData',qscale*hV) 
                hold on; %drawnow;
            end
            
            fprintf('*** Spatial avg of time-mean U2data = %.2f ***\n', my_nanmean(U2data, 'all'));
            fprintf('*** Spatial avg of time-mean V2data = %.2f ***\n', my_nanmean(V2data, 'all'));
            fprintf('*** Max of time-mean U2data = %.2f ***\n', max(U2data, [], 'all'));
            fprintf('*** Max of time-mean V2data = %.2f ***\n', max(V2data, [], 'all'));
            linkaxes([ax_ori, ax_vec2], 'xy');                  %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
        end

        %======== axis/title/figname ========%
        if grid_mode == 0               
            map_xlabel = '';  map_ylabel = '';
            if zm_or_mm == 1                 %/ Time-lat   
                if isempty(intvl_hori) 
                    intvl_hori = 5;
                end

                %/ make use of the min and max contf_hori, create an array that start from the multiplier of the interval
                if isempty(map_xticks)
                    map_xticks = ((min(contf_hori) - mod(min(contf_hori), intvl_hori)):intvl_hori:max(contf_hori))'; %/ since contf_hori is not index array but real lon or lat array.
                elseif size(map_xticks, 1) == 1
                    map_xticks = map_xticks';
                end
                map_xticklabels                      = map_xticks;
                hori_unit                            = repmat({'N'}, length(map_xticklabels), 1);
                deg                                  = repmat(char(176), length(map_xticklabels), 1);
                hori_unit(map_xticklabels == 0)      = {''};
                hori_unit(map_xticklabels < 0)       = {'S'};
                map_xticklabels(map_xticklabels < 0) = abs(map_xticklabels(map_xticklabels < 0));
                map_xticklabels                      = strcat(string(map_xticklabels), deg, hori_unit);

            elseif zm_or_mm == 2              %/ Time-lon        
                if isempty(intvl_hori) 
                    if max(contf_hori) - min(contf_hori) <= 40   
                        intvl_hori = 10;
                    else                                         
                        intvl_hori = 20;  
                    end
                end

                %/ make use of the min and max contf_hori, create an array that start from the multiplier of the interval
                if isempty(map_xticks)
                    map_xticks = ((min(contf_hori) - mod(min(contf_hori), intvl_hori)):intvl_hori:max(contf_hori))'; %/ since contf_hori is not index array but real lon or lat array.
                elseif size(map_xticks, 1) == 1
                    map_xticks = map_xticks';
                end
                map_xticklabels                 = map_xticks;
                hori_unit                       = repmat({'E'}, length(map_xticklabels), 1);  %/ assume lon goes from 0 to 360.
                deg                             = repmat(char(176), length(map_xticklabels), 1);
                hori_unit(map_xticklabels == 0) = {''};
                map_xticklabels                 = strcat(string(map_xticklabels), deg, hori_unit);
            end
        end  
        
        xlabel(ax_ori,  map_xlabel);
        ylabel(ax_ori,  map_ylabel);     
        xticks(ax_ori,  map_xticks) 
        yticks(ax_ori,  map_yticks)
        xticklabels(ax_ori,  map_xticklabels)
        yticklabels(ax_ori,  map_yticklabels)
          
        if pcolor_mode 
            map_xlim        = [min(contf_hori)-res_x/2  max(contf_hori)-res_x/2];   %/ shift to get good visualization
            map_ylim        = [min(contf_time)-res_y/2  max(contf_time)-res_y/2];   %/ shift to get good visualization
            linewidth_multiplier = 3;  %/ somehow pcolor will overlay the panel edge. So make it 2 times thicker.
        else
%             map_xlim        = [min(contf_hori)         max(contf_hori)      ];  
%             map_ylim        = [min(contf_time)+res_y   max(contf_time)+res_y];  
            map_xlim        = [min(contf_hori)   max(contf_hori)];  
            map_ylim        = [min(contf_time)   max(contf_time)];   
            linewidth_multiplier = 1.5;
        end
        xlim(ax_ori,  map_xlim)  %<-- this causes scatter plot to disappear.
        ylim(ax_ori,  map_ylim)
        
        if zm_or_mm == 1        
            if always_time_as_y == 0
                %/ view(az, el)
                %/ az > 0   -> turn counterclockwise with z-axis.
                %/ el = 90  -> plan view.
                %/ el = -90 -> bottom view
                az = 90;  el = -90;
                view(ax_ori, [az el]);    %/ convenient way to swap x & y axis (but have to do it for each axis).
                if ~isempty(ax_border)     view(ax_border,  [az el]);    end
                if ~isempty(ax_scatter)    view(ax_scatter, [az el]);    end
                if ~isempty(ax_line)       view(ax_line,    [az el]);    end
                if ~isempty(ax_vec)        view(ax_vec,     [az el]);    end
                if ~isempty(ax_vec2)       view(ax_vec2,    [az el]);    end
                if ~isempty(ax_SLR)        view(ax_SLR,    [az el]);     end
                ytickangle(ax_ori, map_ytickangle);
            end
        elseif zm_or_mm == 2      
%             set(ax_ori, 'YDir','reverse');   %/ better not to        
        end
        
        set(ax_ori, 'linewidth', linewidth*linewidth_multiplier, 'fontsize', fontsize*0.6)
        hold on;
        box on
        grid off     %/ grid off only works without set() function after it!
        
        %/ Set title after set gca!!
        if isempty(title_fontsize)
            title_fontsize = fontsize*0.6; 
        end
        if isempty(title_ypos)  %/ Standard position for title is 1.
            title_ypos = 1.02;
        end
        titlename = strrep(strcat(titlename, str_SLR), '_', ' ');
        title_pos = [0.12,title_ypos,0.9,0];
        % title(ax_ori,  titlename, 'Fontsize', fontsize*0.8);
        
        %/ Annotation is indep. of axes -> not distorting the aspect ratio.
        annotation( 'textbox', 'String', titlename, 'Color', 'k', ...
                    'FontSize', title_fontsize, 'Units', 'normalized', 'EdgeColor', 'none', ...
                    'Position', title_pos)

%         titlePos = get( th , 'position');
%         titlePos(2) = titlePos(2)*1.4;
%         set(th , 'Fontsize', fontsize*0.7, 'position' , titlePos);
%         drawnow; pause(0.05);
    end
    hold on;
    drawnow;   %/ Show all elements immediately
    pause(5);  %/ 5-sec pause is needed to allow elements to update (after testing)

    FigName_underscore = strcat(FigName_underscore, str_qscale, str_cbar, str_vecref); %/ append vars to FigName_underscore.
    if savefig
        if isempty(fig_fmt)  fig_fmt = 'pdf';   end
        if contains_transparent || isequal(fig_fmt, 'png')
            if isempty(png_dpi)   png_dpi = 300;   end
            final_figname = char(strcat(plotting_folder, FigName_underscore,'.png'));
            fprintf('*** Saving Fig to: %s ***\n', final_figname);
            export_fig(final_figname,['-r',num2str(png_dpi)], '-r300','-opengl', '-nocrop', '-transparent'); %/ '-c[inf, inf, inf, inf]'
            
        elseif isequal(fig_fmt, 'pdf')
            final_figname = char(strcat(plotting_folder, FigName_underscore,'.pdf'));
            fprintf('*** Saving Fig to: %s ***\n', final_figname);
            export_fig(final_figname,'-pdf','-painters', '-nocrop', '-transparent'); %/ '-c[inf, inf, inf, inf]', ); %, 
        else
            error('Invalid input of ''fig_fmt''!');
        end
    end

end
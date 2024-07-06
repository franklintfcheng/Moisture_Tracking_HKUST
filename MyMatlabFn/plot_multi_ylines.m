function [trend_data, pval_data, intercept_data, CI95_data, regress_x, regress_y] = plot_multi_ylines(varargin)

    pnames = { 'line_ydata',       'line_ydata_ref',       'line_xdata',       'line_yerror_upper', 'line_yerror_lower', 'line_yerror_type',   'line_lgd_str',     'line_axislabel',   'line_colmap',...
               'marker',           'markersize',           'markeredgecolor',  'fontsize',          'title_fontsize',    'linewidth',          'linestyle',        'panel_linewidth',...
               'lineerror_colmap', 'lineerror_alpha',      'use_yyaxis',       'yyaxis1_col',       'yyaxis2_col',       'show_y1_as_bar',     'pve_nve_bar_mode', 'bar_pve_col',      'bar_nve_col',      'barwidth',         'scatter_ydata',...
               'yshift',           'show_zero_line',       'hori_lines',       'hori_lines_col',    'hori_lines_st',     'vert_lines',         'vert_lines_col',   'vert_lines_st',    'vert_patches',     'patch_color',      'show_axislabel',...
               'map_yticks',       'map_xlabel',           'map_xticks',       'map_xticklabels',   'map_xticks_minor',  'map_xtickangle',     'xaxis_fontsize',...
               'fig_width',        'fig_height',           'xgap_left',        'xgap_right',        'show_y_max',        'show_diagonal',      'box_mode',...
               'line_data_ylim',   'line_data_xlim',       'trend_mode',       'trend_abs_or_rel',  'trend_alpha',       'trend_test',         'asterisk',         'alpha_level',      'show_regress_line','regress_linewidth',...
               'boxplot_data',     'boxplot_show_whisker', 'lower_whisker',    'upper_whisker',     'vertprofile_mode', 'titlename',         'title_pos',         'FigName_underscore', 'show_legend',      'lgd_position',     'plot_legend_only', 'plotting_folder', 'fig_fmt', 'png_dpi',  'savefig'};
          
    dflts  = cell(1, length(pnames));
    
    [          line_ydata,          line_ydata_ref,        line_xdata,         line_yerror_upper,   line_yerror_lower,   line_yerror_type,    line_lgd_str,       line_axislabel,     line_colmap,...
               marker,              markersize,            markeredgecolor,    fontsize,            title_fontsize,      linewidth,           linestyle,          panel_linewidth,...
               lineerror_colmap,    lineerror_alpha,       use_yyaxis,         yyaxis1_col,         yyaxis2_col,         show_y1_as_bar,      pve_nve_bar_mode,   bar_pve_col,        bar_nve_col,        barwidth,           scatter_ydata,...
               yshift,              show_zero_line,        hori_lines,         hori_lines_col,      hori_lines_st,       vert_lines,          vert_lines_col,     vert_lines_st,      vert_patches,       patch_color,        show_axislabel,...
               map_yticks,          map_xlabel,            map_xticks,         map_xticklabels,     map_xticks_minor,    map_xtickangle,      xaxis_fontsize,...
               fig_width,           fig_height,            xgap_left,          xgap_right,          show_y_max,          show_diagonal,       box_mode,...
               line_data_ylim,      line_data_xlim,        trend_mode,         trend_abs_or_rel,    trend_alpha,         trend_test,          asterisk,           alpha_level,        show_regress_line,  regress_linewidth,...
               boxplot_data,        boxplot_show_whisker,  lower_whisker,      upper_whisker,       vertprofile_mode,   titlename,           title_pos,           FigName_underscore,  show_legend,        lgd_position,       plot_legend_only,   plotting_folder,    fig_fmt,    png_dpi,     savefig]...
                        = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: May 22, 2024
    %/
    %/ DESCRIPTION:
    %/       This function is designed for drawing lines with different
    %/       scales on the same plot.
    %/
    %/       NEW feature: allow to imbed a boxplot on the right
    %/                    allow to plot vertical profile
    %/
    %/ INPUT:
    %/      'line_ydata': A column vector / matrix      (ntime, nvar)
    %/      'line_xdata': Assumes to be a column vector (ntime, 1)
    %/=====================================================================
    
    %------------------------
    if isempty(fontsize)        fontsize = 14;                                  end
    if isempty(xaxis_fontsize)  xaxis_fontsize = fontsize;                      end
    if isempty(linewidth)       linewidth = 2;                                  end
    if isempty(marker) 
        marker = repmat({'o'}, size(line_ydata,2), 1);  
    elseif ~iscell(marker)
        marker = repmat({marker}, size(line_ydata,2), 1);  
    end
    if isempty(markersize)      markersize = 5;                                 end   %/ markersize for line plots
    if size(line_xdata, 1)== 1  line_xdata = line_xdata';                       end   %/ make it a column vector.
    if isempty(use_yyaxis)      use_yyaxis = 0;                                 end
    if isempty(linestyle)       linestyle = '-';                                end                 %/ not a good idea to set 'none' if for time series
    if isempty(box_mode)        box_mode = 'on';                                end
    if isempty(patch_color)     patch_color = repmat([.9 .9 .9], length(vert_patches), 1);     end
    if isempty(line_colmap)     line_colmap = jet(size(line_ydata,2));          end
    if isempty(fig_fmt)         fig_fmt = 'pdf';                                end
    if isempty(barwidth)        
        barwidth = 1;
    elseif barwidth <= 0 || barwidth > 1
        error('Set barwidth to be between 0 and 1!');
    end
    if plot_legend_only         
        show_regress_line = 0;
    elseif isempty(show_regress_line)
        show_regress_line = 1;
    end
    if isempty(show_axislabel)
        show_axislabel = 1;  %/ show yaxis labels by default
    end
    if isempty(line_data_ylim)  %/ Auto set a single ylim
        data_max = max(line_ydata, [], 'all');
        data_min = min(line_ydata, [], 'all');
        data_range = data_max - data_min;
        line_data_ylim = [data_min-data_range*0.1, data_max+data_range*0.1];

        % warning('please input ''line_data_ylim''!');  
    end
    if isempty(line_lgd_str) 
        line_lgd_str = strcat('Var', string((1:size(line_ydata,2))'))';
    else
        if ischar(line_lgd_str)
            line_lgd_str = {line_lgd_str};
        end
        if size(line_lgd_str, 2) ~= 1      %/ make line_lgd_str a column vector to avoid bug
            line_lgd_str = line_lgd_str';  
        end
    end
    
    if vertprofile_mode
        if use_yyaxis == 1
            use_yyaxis = 0;  
            warning('vertprofile_mode requires to set use_yyaxis == 0. Auto setting it to 0!')
        end
        if size(line_data_ylim, 1) ~= 1
            error('vertprofile_mode requires to input only a single row of line_data_ylim!')
        end
    end
    
    %/ To compute relative trend, we may simply take natural log of line_ydata
    if isempty(trend_abs_or_rel)  
        trend_abs_or_rel = 'abs';   %/ do nothing.
    end

    if ~iscell(line_axislabel)
        line_axislabel = {line_axislabel};
    end

    %/ If to show trend in relative change (%), divide it by the value at the smallest x
    if ismember(trend_abs_or_rel, {'rel'})
        if isempty(line_ydata_ref)
            error('Input line_ydata_ref for trend_abs_or_rel = ''rel''! ');
        end
        line_ydata_multiply = 1./line_ydata_ref*100;   
        line_ydata          = line_ydata        .* line_ydata_multiply;
        line_yerror_upper   = line_yerror_upper .* line_ydata_multiply;
        line_yerror_lower   = line_yerror_lower .* line_ydata_multiply;
        line_axislabel = repmat({'%'}, length(line_axislabel), 1);  %/ update unit
        
    elseif ~ismember(trend_abs_or_rel, {'abs'})
        error('''trend_abs_or_rel'' can only be ''abs'' or ''rel''!');
    end
    
    %/ reference line
    if ~isempty(yshift)  %/ then we substract yshift from the original line_ydata
        line_ydata     = line_ydata - yshift;
        line_data_ylim = line_data_ylim - yshift;
    else
        yshift = 0;
    end
    
    %/ Handle error bars / shadings
    if ~isempty(line_yerror_upper) && ~isempty(line_yerror_lower) 
        if isempty(line_yerror_type) || ~ismember(line_yerror_type, {'sd', 'percentile'})
            error('Please clarify ''line_yerror_type''!! E.g., ''sd'' or ''percentile''?');
        end
        if line_yerror_upper < line_yerror_lower
            error('Detected that line_yerror_upper < line_yerror_lower! Please check!');
        end
    end
    
    trend_data = []; pval_data = []; intercept_data = []; regress_x =[]; str_sig = []; %/ as they are output, initialize them first.
    str_legend_only = '';
    if plot_legend_only
        str_legend_only         = '_lgd';
        line_ydata_bc_trendmode = line_ydata;  %/ Keep the data for processing trends -> generate str_sig for legend
        line_ydata_bc           = nan(size(line_ydata));  %/ make dummy plot 
    else
        line_ydata_bc_trendmode = line_ydata;
        line_ydata_bc           = line_ydata;
    end
    Nvar = size(line_ydata_bc, 2);
    line_xdata_bc   = line_xdata;

    if isempty(fig_width)       fig_width  = 1200;          end
    if isempty(fig_height)      fig_height = 600;           end
    str_repeated    = [];

    figure
    set(gcf, 'color','w');
    set(gcf, 'position', [0 0 fig_width fig_height])
    set(gca, 'linewidth', 2.5, 'fontsize', fontsize) %/ set it first
    hold on; 

    %/ find the unique set of ylims if multiple rows of line_data_ylim are given
    if size(line_data_ylim, 1) ~= 1
        [ylim_sets, ~, ind_which] = unique(line_data_ylim, 'rows', 'stable');
        if size(ylim_sets, 1) == 1
            ylim_sets = [ylim_sets; ylim_sets];
            if use_yyaxis 
                %/ nothing needed to change.
            else
                ind_which = 1:length(line_data_ylim);
            end
        end
    else
        ylim_sets   = repmat(line_data_ylim, Nvar, 1);
        ind_which   = ones(Nvar, 1);
    end

    if ~isempty(map_yticks)
        if ~iscell(map_yticks) %/ For it can only be cell to indicate multiple y-axis ticks
            map_yticks = {map_yticks};                %/ Update
            map_yticks = repmat(map_yticks, Nvar, 1); %/ Update
        end
    end

    %/ Marker for the maxima (if show_y_max == 1)
    scatter_marker = 's'; 
    scatter_alpha  = 0.3;
    scatter_size   = 200;

    if isempty(bar_pve_col) 
        bar_pve_col = [255 176 127]./255; %/ orange
    end
    if isempty(bar_nve_col)
        bar_nve_col = [142 127 255]./255; %/ purple
    end

    trend_data = []; pval_data = []; intercept_data = []; CI95_data = [];
    h_line = nan(Nvar, 1);
    regress_y = nan(2, Nvar);  %/ ([start,end], nvar)
    for k = 1:Nvar
        if ismember(marker(k), {'d','s'})  %/ make it a hollow diamond
            markerfacecolor = 'none'; 
            markeredgecolor_bc = line_colmap(k,:);
        else
            markerfacecolor = line_colmap(k,:);
            if ~isempty(markeredgecolor) 
                markeredgecolor_bc = markeredgecolor;          
            else
                markeredgecolor_bc = 'none';
            end
        end
        
        if k == 1 && trend_mode
            [trend_data, ~, ~, pval_data, intercept_data, CI95_data] = ...
                    compute_trend('data_yearly', line_ydata_bc_trendmode', 'year_list', line_xdata_bc, 'trend_test', trend_test, 'trend_alpha', trend_alpha);
            
            if isempty(asterisk)    asterisk    = {'(*)', '(**)', '(***)'};    end
            if isempty(alpha_level) alpha_level = [0.05,    0.01,   0.001];    end
            
            %/ Use asterisks to indicate significance level
            str_sig = repmat({''}, length(pval_data), 1);
            for i = 1:length(alpha_level)
                ind = find(pval_data <= alpha_level(i));
                str_sig(ind) = repmat(asterisk(i), length(ind), 1);
            end
            for i = 1:length(pval_data)
                if isempty(str_sig{i})
                    str_sig{i} = sprintf(' (p=%.2f)',pval_data(i));
                end
            end

            %/ Sort the xdata, then get min and max of it (useful when xdata is a parameter other than time)
            [line_xdata_bc_sorted, I] = sort(line_xdata_bc, 'ascend'); %/ [nan nan 5 10 3 nan]  -> [3     5    10   NaN   NaN]
            line_ydata_bc_trendmode   = line_ydata_bc_trendmode(I,:); %/ because xdata has been sorted, we need to sort ydata accordingly
        end

        if use_yyaxis 
            if ind_which(k) == 1
                yyaxis left
                if k == 1 && show_y1_as_bar
                    %/ Whether to draw +ve and -ve bars with two different colors.
                    if pve_nve_bar_mode
                        bar_data_pve = line_ydata_bc(:,k); bar_data_pve(bar_data_pve <  0) = nan;
                        bar_data_nve = line_ydata_bc(:,k); bar_data_nve(bar_data_nve >= 0) = nan;

                        bar(line_xdata_bc, bar_data_pve, 'facecolor', bar_pve_col, 'edgecolor', 'k', 'linewidth', linewidth*0.5, 'BarWidth', barwidth); %/ make the first lineplot thicker.
                        h_line(k) = bar(line_xdata_bc, bar_data_nve, 'facecolor', bar_nve_col, 'edgecolor', 'k', 'linewidth', linewidth*0.5, 'BarWidth', barwidth); %/ make the first lineplot thicker.
                    else
                        h_line(k) = bar(line_xdata_bc, line_ydata_bc(:,k), 'facecolor', line_colmap(k,:), 'edgecolor', 'k', 'linewidth', linewidth*0.5, 'BarWidth', barwidth); %/ make the first lineplot thicker.
                    end
                else
                    %/ Show errorbar as *shading* 
                    if ~isempty(line_yerror_upper) && ~isempty(line_yerror_lower)
                        ind_nonnan = find(~isnan(line_ydata_bc(:,k)));
                        x_vector = [line_xdata_bc(ind_nonnan)', fliplr(line_xdata_bc(ind_nonnan)')];

                        if ismember(line_yerror_type, {'sd'}) %/ for sd, no need to substract it from yshift 
                            y_vector = [line_ydata_bc(ind_nonnan,k)'+line_yerror_upper(ind_nonnan,k)',fliplr(line_ydata_bc(ind_nonnan,k)'-line_yerror_lower(ind_nonnan,k)')];

                        elseif ismember(line_yerror_type, {'percentile'})
                            y_vector = [line_yerror_upper(ind_nonnan,k)',fliplr(line_yerror_lower(ind_nonnan,k)')] - yshift;
                        end
                        error_patch = fill(x_vector, y_vector, lineerror_colmap(k,:));
                        set(error_patch, 'edgecolor', 'none');
                        set(error_patch, 'FaceAlpha', lineerror_alpha);
                        hold on;
                    end
                    h_line(k) = plot(line_xdata_bc, line_ydata_bc(:,k), 'color', line_colmap(k,:), 'linewidth', linewidth*0.5, 'linestyle', linestyle, 'marker', marker{k}, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'markeredgecolor', markeredgecolor_bc);
                end
            
            elseif ind_which(k) == 2
                yyaxis right
                %/ Show errorbar as *shading* 
                if ~isempty(line_yerror_upper) && ~isempty(line_yerror_lower)
                    ind_nonnan = find(~isnan(line_ydata_bc(:,k)));
                    x_vector = [line_xdata_bc(ind_nonnan)', fliplr(line_xdata_bc(ind_nonnan)')];

                    if ismember(line_yerror_type, {'sd'}) %/ for sd, no need to substract it from yshift 
                        y_vector = [line_ydata_bc(ind_nonnan,k)'+line_yerror_upper(ind_nonnan,k)',fliplr(line_ydata_bc(ind_nonnan,k)'-line_yerror_lower(ind_nonnan,k)')];

                    elseif ismember(line_yerror_type, {'percentile'})
                        y_vector = [line_yerror_upper(ind_nonnan,k)',fliplr(line_yerror_lower(ind_nonnan,k)')] - yshift;
                    end

                    error_patch = fill(x_vector, y_vector, lineerror_colmap(k,:));
                    set(error_patch, 'edgecolor', 'none');
                    set(error_patch, 'FaceAlpha', lineerror_alpha);
                    hold on;
                end

                h_line(k) = plot(line_xdata_bc, line_ydata_bc(:,k), 'color', line_colmap(k,:), 'linewidth', linewidth*0.5, 'linestyle', linestyle, 'marker', marker{k}, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'markeredgecolor', markeredgecolor_bc);
            
                ylim(ylim_sets(ind_which(k), :));
                if show_axislabel
                    ylabel(line_axislabel{k});
                end
            end
            
            %/ Set yyaxis color here
            yyaxis left
            ax = gca;
            if isempty(yyaxis1_col)  
                yyaxis1_col = line_colmap(k,:);   
            end 
            ax.YAxis(1).Color = yyaxis1_col;
            yyaxis right
            ax = gca;
            if isempty(yyaxis2_col)  
                yyaxis2_col = line_colmap(k,:);   
            end 
            ax.YAxis(2).Color = yyaxis2_col;
            
            %/ Set 'ylim' and 'ylabel'
            %/ if each var is assigned 'ylim' and 'ylabel', we will follow that
            %/ if only two ylims and ylabels are provided, then we assign
            %/ each to left and right yyaxis, respectively.
            if ind_which(k) == 1
                yyaxis left
            else
                yyaxis right
            end
            if length(ylim_sets) == Nvar
                ylim(ylim_sets(k, :));
            else
                ylim(ylim_sets(ind_which(k), :));
            end

            if ~isempty(map_yticks)
                yticks(map_yticks{k});
            end

            if show_axislabel
                if length(line_axislabel) == Nvar
                    ylabel(line_axislabel{k});
                else
                    ylabel(line_axislabel{ind_which(k)});
                end
            end
            if ~isempty(scatter_ydata)
                sc = scatter(line_xdata_bc, scatter_ydata(:,k), scatter_size, 'Marker', scatter_marker,...
                            'MarkerFaceColor', 'none', 'MarkerEdgeColor', line_colmap(k,:), 'linewidth', linewidth*0.5); 
                sc.MarkerFaceAlpha = scatter_alpha;
            end
            
            %/ Add trend to each var if trend_mode = 1
            if trend_mode
                % if k == 1
                %     [trend_data, ~, ~, pval_data, intercept_data, CI95_data] = ...
                %             compute_trend('data_yearly', line_ydata_bc_trendmode', 'year_list', line_xdata_bc, 'trend_test', trend_test, 'trend_alpha', trend_alpha);
                % 
                %     if isempty(asterisk)    asterisk    = {'(*)', '(**)', '(***)'};    end
                %     if isempty(alpha_level) alpha_level = [0.05,    0.01,   0.001];    end
                % 
                %     %/ Use asterisks to indicate significance level
                %     str_sig = repmat({''}, length(pval_data), 1);
                %     for i = 1:length(alpha_level)
                %         ind = find(pval_data <= alpha_level(i));
                %         str_sig(ind) = repmat(asterisk(i), length(ind), 1);
                %     end
                %     for i = 1:length(pval_data)
                %         if isempty(str_sig{i})
                %             str_sig{i} = sprintf(' (p=%.2f)',pval_data(i));
                %         end
                %     end
                % 
                %     %/ Sort the xdata, then get min and max of it (useful when xdata is a parameter other than time)
                %     [line_xdata_bc_sorted, I] = sort(line_xdata_bc, 'ascend'); %/ [nan nan 5 10 3 nan]  -> [3     5    10   NaN   NaN]
                %     line_ydata_bc_trendmode   = line_ydata_bc_trendmode(I,:); %/ because xdata has been sorted, we need to sort ydata accordingly
                % end
                
                if ind_which(k) == 1
                    yyaxis left
                elseif ind_which(k) == 2
                    yyaxis right
                end
                
                %/ Fit a simple linear regression line (if the column is not entirely NaN)
                if ~all(isnan(line_ydata_bc_trendmode(:,k)))
                    ind_nonnan    = nan(2,1);
                    ind_nonnan(1) = find(~isnan(line_ydata_bc_trendmode(:,k)), 1, 'first');
                    ind_nonnan(2) = find(~isnan(line_ydata_bc_trendmode(:,k)), 1, 'last');

                    if pval_data(k) <= trend_alpha
                        trend_linestyle = '-';    %/ significant
                    else
                        trend_linestyle = '--';
                    end

                    %/ We only need two point to draw a regression line
                    regress_x      = line_xdata_bc_sorted(ind_nonnan([1,end]));
                    regress_y(:,k) = trend_data(k).*regress_x + intercept_data(k);

%                     regress_line_color = 'k';
                    regress_line_color = line_colmap(k,:);
                    if show_regress_line
                        if isempty(regress_linewidth)   
                            regress_linewidth = linewidth;   
                        end
                        plot(regress_x, regress_y(:,k), 'color', regress_line_color, 'linewidth', regress_linewidth,...
                             'linestyle', trend_linestyle, 'marker', 'none');
                    end
                end
            end
            
            %/ if all ylim_sets are identical, we copy the axis from left to right
            if k == Nvar && size(ylim_sets, 1) == 2
                if all(ind_which == 1)
                    ax.YAxis(2).Limits = ax.YAxis(1).Limits;
                end
            end

        elseif use_yyaxis == 0 && size(line_data_ylim, 1) == 1  %/ Just a normal line plot without multiple axis (useful for vertical profile!)
            if length(ylim_sets) == Nvar
                ylim(ylim_sets(k, :));
                if ~isempty(map_yticks)
                    yticks(map_yticks{k});
                end
            else
                ylim(ylim_sets(ind_which(k), :));
                if ~isempty(map_yticks)
                    yticks(map_yticks{ind_which(k)});
                end
            end
            if show_axislabel
                if length(line_axislabel) == Nvar
                    ylabel(line_axislabel{k});
                else
                    ylabel(line_axislabel{ind_which(k)});
                end
            end
            hold on;
            if show_y1_as_bar
                %/ Whether to draw +ve and -ve bars with two different colors.
                if pve_nve_bar_mode
                    bar_data_pve = line_ydata_bc(:,k); bar_data_pve(bar_data_pve <  0) = nan;
                    bar_data_nve = line_ydata_bc(:,k); bar_data_nve(bar_data_nve >= 0) = nan;
                    
                    bar(line_xdata_bc, bar_data_pve, 'facecolor', bar_pve_col, 'edgecolor', 'k', 'linewidth', linewidth, 'BarWidth', barwidth); %/ make the first lineplot thicker.
                    bar(line_xdata_bc, bar_data_nve, 'facecolor', bar_nve_col, 'edgecolor', 'k', 'linewidth', linewidth, 'BarWidth', barwidth); %/ make the first lineplot thicker.
                else
                    bar(line_xdata_bc, line_ydata_bc(:,k), 'facecolor', line_colmap(k,:), 'edgecolor', 'k', 'linewidth', linewidth, 'BarWidth', barwidth); %/ make the first lineplot thicker.
                end

                %/ Draw a dummy line (will remove later) to enable addaxis() to work!
                h_line_dummy = plot(line_xdata_bc, line_ydata_bc(:,k), 'color', line_colmap(k,:), 'linewidth', linewidth*1.5, 'marker', marker{k}, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'markeredgecolor', markeredgecolor_bc); %/ dummy. Just for addaxis.
            else
                %/ Show errorbar as *shading* 
                if ~isempty(line_yerror_upper) && ~isempty(line_yerror_lower)
                    ind_nonnan = find(~isnan(line_ydata_bc(:,k)));
                    x_vector = [line_xdata_bc(ind_nonnan)', fliplr(line_xdata_bc(ind_nonnan)')];
        
                    if ismember(line_yerror_type, {'sd'}) %/ for sd, no need to substract it from yshift 
                        y_vector = [line_ydata_bc(ind_nonnan,k)'+line_yerror_upper(ind_nonnan,k)',fliplr(line_ydata_bc(ind_nonnan,k)'-line_yerror_lower(ind_nonnan,k)')];
        
                    elseif ismember(line_yerror_type, {'percentile'})
                        y_vector = [line_yerror_upper(ind_nonnan,k)',fliplr(line_yerror_lower(ind_nonnan,k)')] - yshift;
                    end
                    error_patch = fill(x_vector, y_vector, lineerror_colmap(k,:));
                    set(error_patch, 'edgecolor', 'none');
                    set(error_patch, 'FaceAlpha', lineerror_alpha);
                    hold on;
                end
                h_line(k) = plot(line_xdata_bc, line_ydata_bc(:,k), 'color', line_colmap(k,:), 'linewidth', linewidth*0.5, 'linestyle', linestyle, 'marker', marker{k}, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'markeredgecolor', markeredgecolor_bc);
            end

            %/ Fit a simple linear regression line (if the column is not entirely NaN)
            if ~all(isnan(line_ydata_bc_trendmode(:,k)))
                ind_nonnan    = nan(2,1);
                ind_nonnan(1) = find(~isnan(line_ydata_bc_trendmode(:,k)), 1, 'first');
                ind_nonnan(2) = find(~isnan(line_ydata_bc_trendmode(:,k)), 1, 'last');

                if pval_data(k) <= trend_alpha
                    trend_linestyle = '-';    %/ significant
                else
                    trend_linestyle = '--';
                end

                %/ We only need two point to draw a regression line
                regress_x      = line_xdata_bc_sorted(ind_nonnan([1,end]));
                regress_y(:,k) = trend_data(k).*regress_x + intercept_data(k);

%                     regress_line_color = 'k';
                regress_line_color = line_colmap(k,:);
                if show_regress_line
                    if isempty(regress_linewidth)   
                        regress_linewidth = linewidth;   
                    end
                    plot(regress_x, regress_y(:,k), 'color', regress_line_color, 'linewidth', regress_linewidth,...
                         'linestyle', trend_linestyle, 'marker', 'none');
                end
            end

        else  %/ Otherwise, add an axis to each var
            if trend_mode   
                error('trend_mode has not yet worked with use_yyaxis!');   
            end
            if k == 1
                if show_y1_as_bar
                    %/ Whether to draw +ve and -ve bars with two different colors.
                    if pve_nve_bar_mode
                        bar_data_pve = line_ydata_bc(:,k); bar_data_pve(bar_data_pve <  0) = nan;
                        bar_data_nve = line_ydata_bc(:,k); bar_data_nve(bar_data_nve >= 0) = nan;
                        
                        bar(line_xdata_bc, bar_data_pve, 'facecolor', bar_pve_col, 'edgecolor', 'k', 'linewidth', linewidth, 'BarWidth', barwidth); %/ make the first lineplot thicker.
                        bar(line_xdata_bc, bar_data_nve, 'facecolor', bar_nve_col, 'edgecolor', 'k', 'linewidth', linewidth, 'BarWidth', barwidth); %/ make the first lineplot thicker.
                    else
                        bar(line_xdata_bc, line_ydata_bc(:,k), 'facecolor', line_colmap(k,:), 'edgecolor', 'k', 'linewidth', linewidth, 'BarWidth', barwidth); %/ make the first lineplot thicker.
                    end

                    %/ Draw a dummy line (will remove later) to enable addaxis() to work!
                    h_line_dummy = plot(line_xdata_bc, line_ydata_bc(:,k), 'color', line_colmap(k,:), 'linewidth', linewidth*1.5, 'marker', marker{k}, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'markeredgecolor', markeredgecolor_bc); %/ dummy. Just for addaxis.
                else
                    %/ Show errorbar as *shading* if 'line_yerror_bc' is given
                    if ~isempty(line_yerror_upper) && ~isempty(line_yerror_lower)
                        ind_nonnan = find(~isnan(line_ydata_bc(:,k)));
                        x_vector = [line_xdata_bc(ind_nonnan)', fliplr(line_xdata_bc(ind_nonnan)')];
                        
                        if ismember(line_yerror_type, {'sd'}) %/ for sd, no need to substract it from yshift 
                            y_vector = [line_ydata_bc(ind_nonnan,k)'+line_yerror_upper(ind_nonnan,k)',fliplr(line_ydata_bc(ind_nonnan,k)'-line_yerror_lower(ind_nonnan,k)')];
    
                        elseif ismember(line_yerror_type, {'percentile'})
                            y_vector = [line_yerror_upper(ind_nonnan,k)',fliplr(line_yerror_lower(ind_nonnan,k)')] - yshift;
                        end
                        error_patch = fill(x_vector, y_vector, lineerror_colmap(k,:));
                        set(error_patch, 'edgecolor', 'none');
                        set(error_patch, 'FaceAlpha', lineerror_alpha);
                        hold on;
                    end
                    
                    plot(line_xdata_bc, line_ydata_bc(:,k), 'color', line_colmap(k,:), 'linewidth', linewidth*1.5); %/ make the first lineplot thicker.
                end
                set(gca, 'ylim', ylim_sets(ind_which(k), :));
            else
                %/ Show errorbar as *shading* if 'line_yerror_bc' is given
                if ~isempty(line_yerror_upper) && ~isempty(line_yerror_lower)
                    ind_nonnan = find(~isnan(line_ydata_bc(:,k)));
                    x_vector = [line_xdata_bc(ind_nonnan)', fliplr(line_xdata_bc(ind_nonnan)')];

                    if ismember(line_yerror_type, {'sd'}) %/ for sd, no need to substract it from yshift 
                        y_vector = [line_ydata_bc(ind_nonnan,k)'+line_yerror_upper(ind_nonnan,k)',fliplr(line_ydata_bc(ind_nonnan,k)'-line_yerror_lower(ind_nonnan,k)')];

                    elseif ismember(line_yerror_type, {'percentile'})
                        y_vector = [line_yerror_upper(ind_nonnan,k)',fliplr(line_yerror_lower(ind_nonnan,k)')] - yshift;
                    end

                    error_patch = fill(x_vector, y_vector, lineerror_colmap(k,:));
                    set(error_patch, 'edgecolor', 'none');
                    set(error_patch, 'FaceAlpha', lineerror_alpha);
                    hold on;
                end
                
                %/ The 3rd argument is for ylim!
                addaxis(line_xdata_bc, line_ydata_bc(:,k), ylim_sets(ind_which(k), :), 'color', line_colmap(k,:), 'linewidth', linewidth*1.5, 'marker', marker{k}, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'markeredgecolor', markeredgecolor_bc); %/ NOTE: addaxis forces the original axes to show in black.
            end
            
            if ~isempty(map_yticks)
                warning('Adjusting yticks for each addaxis may not be possible! Please check.');
                yticks(map_yticks{k});
            end

            if show_axislabel
                addaxislabel(k, strrep(line_axislabel{k}, '_', ' '), 'linewidth', linewidth, 'fontsize', fontsize*0.5);
            else
                addaxislabel(k, '', 'linewidth', linewidth, 'fontsize', fontsize*0.5);
            end
            hold on;
        end
        
        %/ Show maximum of the lines.
        if show_y_max
            if k == 1 && show_y1_as_bar
                continue;
            else
                [MAX, ind_MAX] = max(line_ydata_bc(:,k));

                %/ IMPORTANT: addaxisplot(x,y,axis_number,...);
                %/            you must indicate axis_number to overlay the marker
                %/            on the correct axis of the line plot!!
                addaxisplot(line_xdata_bc(ind_MAX), MAX, k, 'linestyle', 'none', 'Marker', 'p', 'linewidth', linewidth*0.5,...
                            'MarkerFaceColor', line_colmap(k,:), 'MarkerEdgeColor', 'w',  'markersize', fontsize*3);
                hold on;
            end
        end
    end

    if use_yyaxis == 0 && show_y1_as_bar    
        h_line_dummy.Color = 'none';   
        h_line_dummy.Marker = 'none';  
    end  %/ Remove line plot of y1 when bar plot is shown.

    %/ Add vertical patches (it should be a cell with a pair of x values indicating the range)
    for i = 1:length(vert_patches)
        xl = vert_patches{i};
        yl = ylim; %/ get the current y limits
        patch_x = [xl(1) xl(2) xl(2) xl(1)];
        patch_y = [yl(1) yl(1) yl(2) yl(2)];
        p = patch(patch_x,patch_y, patch_color(i,:), 'edgecolor', 'none');
        uistack(p,'bottom')
    end
    %/ [IMPORTANT] Move the axis to the top of all object (otherwise patches will overlap it!)
    set(gca, 'Layer', 'top');   


    %/ Add horizontal lines if queried
    for i = 1:length(hori_lines)
        if isempty(hori_lines_col)  hori_lines_col = 'k';         end
        if isempty(hori_lines_st)   hori_lines_st  = '--';        end
        plot([-inf, inf], [hori_lines(i), hori_lines(i)], 'color', hori_lines_col, 'linestyle', hori_lines_st, 'linewidth', linewidth*0.75);
        % yline(hori_lines(i), 'color', hori_lines_col, 'linestyle', hori_lines_st, 'linewidth', linewidth*0.75);
    end


    %/ Add vertical lines if queried
    for i = 1:length(vert_lines)
        if isempty(vert_lines_col)  vert_lines_col = 'k';         end
        if isempty(vert_lines_st)   vert_lines_st  = '--';        end
        xline(vert_lines(i), 'color', vert_lines_col, 'linestyle', vert_lines_st, 'linewidth', linewidth*0.75);
    end
    
    %/ Add a zero line (if no bar plot)
    if show_zero_line
        yline(0, 'color', [0 0 0]./255, 'linestyle', '--', 'linewidth', linewidth*0.25);
        hold on;
    end

    %/ Overlay a boxplot on the right
    if ~isempty(boxplot_data)
        %/ Default
        boxplot_gap       = 8;
        boxplot_width     = 5;
        boxplot_shift     = 0;
        boxplot_linewidth = 2;
        boxplot_color     = line_colmap;
        boxplot_outlier_marker       = 'none';
        boxplot_mean_markersize      = 8;
        boxplot_mean_markeredgecolor = 'w';
        boxplot_mean_markerfacecolor = [204 0 0]./255;
        

        %/ Clear the column with all NaNs
        boxplot_data_bc = boxplot_data;
        ind = find(all(isnan(boxplot_data_bc)));
        boxplot_data_bc(:,ind) = [];
        boxplot_color(ind,:) = [];
        y_mean = mean(boxplot_data_bc, 1, 'omitnan');  %/ Get mean from the original data distribution

        %/ When showing whisker, show 10th-90th instead of min and max!
        if isempty(boxplot_show_whisker)
            boxplot_show_whisker = 0;   %/ 0]: show the 25th-75th range only; 1]: show also the non-outliner max and min
        end
        if isempty(lower_whisker)
            lower_whisker = 10;
        end
        if isempty(upper_whisker)
            upper_whisker = 90;
        end
        
        for i = 1:size(boxplot_data_bc,2)
            %/ Whether to show whiskers or not
            if boxplot_show_whisker   
                WhiskerLineColor = 'k';
                %/ 1. Remove outliers using the "quartiles" method, the detection threshold factor replaces the number of interquartile ranges, which is 1.5 by default.
                y = rmoutliers(boxplot_data_bc(:,i), "quartiles", 1); 

                %/ 2. Remove unwanted data (e.g., <10th or >90th) 
                LB = prctile(y, lower_whisker);  %/ Lower bound (by percentile)
                UB = prctile(y, upper_whisker);  %/ Upper bound (by percentile)
                y(y < LB | y > UB) = [];         %/ Remove those outside the bounds
            else
                y = boxplot_data_bc(:,i);
                WhiskerLineColor = 'none';
            end

            boxplot_pos = line_xdata(end)+boxplot_gap*(1:size(boxplot_data_bc,2))+boxplot_shift;
            xdata = repelem(boxplot_pos(i),size(y,1),1); % Specify x=coordinate of box

            %/ NOTE: In boxchart function, Outliers are values that are more than 1.5 · IQR away from the top or bottom of the box
            boxchart(xdata, y, 'BoxWidth', boxplot_width, 'WhiskerLineColor', WhiskerLineColor, 'MarkerStyle', boxplot_outlier_marker,...
                     'BoxFaceColor', boxplot_color(i,:), 'BoxEdgeColor', boxplot_color(i,:), 'BoxMedianLineColor', boxplot_color(i,:),...
                     'LineWidth',  boxplot_linewidth, 'BoxFaceAlpha', lineerror_alpha);
            hold on;

            plot(boxplot_pos(i), y_mean(i), 'marker', 'o', 'markerfacecolor', boxplot_mean_markerfacecolor, 'MarkerEdgeColor', boxplot_mean_markeredgecolor,...
            'linewidth', boxplot_linewidth, 'markersize', boxplot_mean_markersize);
            hold on;
        end
    end
    
    %===== grids/labels setting =====%
    if isempty(xgap_left)           
        xgap_left  = 0.5;                
    else
        xgap_left = abs(xgap_left);  %/ make sure it's always in positive number
    end
    if isempty(xgap_right)          
        xgap_right = 0.5;                
    else
        xgap_right = abs(xgap_right);  %/ make sure it's always in positive number
    end
    
    clear xlim;  %/ make sure xlim is a built-in function, not a variable that was mistakenly defined
    if isempty(line_data_xlim)
        if ~isempty(boxplot_data)
            line_data_xlim = [min(line_xdata_bc)-xgap_left max(boxplot_pos)+boxplot_width/2+xgap_right]; %/ add space in the rightmost to show boxplots (if present)
        else
            line_data_xlim = [min(line_xdata_bc)-xgap_left max(line_xdata_bc)+xgap_right];
        end
    end
    % line_xdata_bc
    % line_data_xlim
    xlim(line_data_xlim);
    if ~isempty(map_xticks)         
        [map_xticks,I] = sort(map_xticks, 'ascend'); %/ it must be in ascending order, otherwise will encounter a bug
        if ~isempty(map_xticklabels)
            map_xticklabels = map_xticklabels(I);
        end
        xticks(map_xticks);
    end
    if ~isempty(map_xticklabels)    xticklabels(map_xticklabels);   end
    if ~isempty(map_xtickangle)     xtickangle(map_xtickangle);     end
    if ~isempty(map_xlabel)         xlabel(map_xlabel);             end
    
    %/ Add minor xticks if queried
    if ~isempty(map_xticks_minor)
        hA = get(gca);
        hA.XAxis.MinorTickValues = map_xticks_minor;
        hA.XAxis.MinorTick='on';
    end
    if isempty(lgd_position)   lgd_position = 'southeast';    end
    if show_legend 
        legend(h_line(:), strcat(strrep(line_lgd_str, '_', ' '), str_sig), 'Location',  lgd_position,  'linewidth', linewidth, 'fontsize', fontsize);
        legend boxoff
    end
    if plot_legend_only
        legend(h_line(:), strcat(strrep(line_lgd_str, '_', ' '), str_sig), 'Location',  lgd_position,  'linewidth', linewidth, 'fontsize', fontsize);
        axis off %hide axis
        set(gca,'visible','off')
        legend boxoff
    end
    ax=gca;
    ax.XAxis.FontSize = xaxis_fontsize;  %/ enlarge xtick fontsize!
    
    %/ grid off only works without set() function after it!
    grid off;
    hold on;
    
    %/ add a diagonal line
    if show_diagonal
        plot([0 ax.XLim(end)],[0 ax.YLim(end)], 'linestyle', '-', 'color', 'k', 'linewidth', linewidth*0.5); 
        hold on;
    end

    %/ Control the box here!
    if isequal(box_mode, 'on') || isequal(box_mode, 1)
        box on;  
    else
        box off;
        if use_yyaxis
            ax.YAxis(2).Visible='off';
        end
    end

    %/ Write titlename using Annotation is indep. of axes -> not distorting the aspect ratio.
    if isempty(title_fontsize)   title_fontsize = fontsize;  end
    % if  isempty(shift_title_y)   shift_title_y = 0.02;     end

    if isempty(title_pos)  title_pos = [0, 1.01, 1, 0];   end

    annotation( 'textbox', 'String', strrep(titlename, '_', ' '), 'Color', 'k', ...
                'FontSize', title_fontsize, 'Units', 'normalized', 'EdgeColor', 'none', ...
                'Position', title_pos, 'FontWeight', 'Bold'); %'HorizontalAlignment', 'left', 

    if isempty(panel_linewidth)
        panel_linewidth = linewidth;
    end

    set(gca, 'linewidth', panel_linewidth);
    set(gca, 'Color', 'none');
    set(gcf, 'Color', 'none');

    if vertprofile_mode %/ Rotate the line plot into a vertical profile (doesn't work with multiple y-axis!!)
        view([90 90])  %/ This works better than camroll
    end
    drawnow; pause(3); 
    
    if savefig
        fig = get(0, 'CurrentFigure');
        hasTransparency = ~isempty(findall(fig,'-property','FaceAlpha','-and','-not','FaceAlpha',1)); %/ code script borrowed from 'export_fig'!
        if hasTransparency
            set(fig, 'InvertHardCopy', 'off'); %/ [KEEP THIS LINE! OR YOU WON"T GET ANY PLOTS!] transparent background 
            set(fig, 'Color', 'None');         %/ [KEEP THIS LINE! OR YOU WON"T GET ANY PLOTS!] transparent background 
            if isempty(png_dpi)
                png_dpi = 300;
            end
            final_FigName = char(strcat(plotting_folder, FigName_underscore, str_repeated, str_legend_only, '.png'));
            export_fig(final_FigName,['-r',num2str(png_dpi)],'-png', '-opengl', '-nocrop', '-transparent');

            % final_FigName = char(strcat(plotting_folder, FigName_underscore, str_repeated, str_legend_only, '.pdf'));
            % print(gcf, '-dpng',final_FigName); % '-fillpage'); 
        else
            % final_FigName = char(strcat(plotting_folder, FigName_underscore, str_repeated, str_legend_only, '.pdf'));
            % print(gcf, '-dpdf',final_FigName, '-bestfit'); % '-fillpage'); 
    
            set(fig, 'InvertHardCopy', 'off'); %/ [KEEP THIS LINE! OR YOU WON"T GET ANY PLOTS!] transparent background 
            set(fig, 'Color', 'None');         %/ [KEEP THIS LINE! OR YOU WON"T GET ANY PLOTS!] transparent background 
            if isequal(fig_fmt, 'pdf')
                final_FigName = char(strcat(plotting_folder, FigName_underscore, str_repeated, str_legend_only, '.pdf'));
                export_fig(final_FigName,'-pdf','-painters', '-nocrop', '-transparent'); %, '-c[inf, inf, inf, inf]'  ); %'-transparent');
            elseif isequal(fig_fmt, 'png')
                final_FigName = char(strcat(plotting_folder, FigName_underscore, str_repeated, str_legend_only, '.png'));
                export_fig(final_FigName,['-r',num2str(png_dpi)],'-png', '-opengl', '-nocrop', '-transparent');
            else
                error('invalid fig_fmt = %s!', fig_fmt);
            end
        end
        fprintf('*** Fig is saved: %s ***\n', final_FigName);
    end
    set(gca, 'Color', 'w'); %/ resume the white background for preview
    set(gcf, 'Color', 'w');

    %/ Finally, show the mean of each variable for checking
    fprintf('\n*** [plot_multi_ylines] Mean of each variable in line_ydata ***\n')
    RowNames = {'mean'};
    VariableNames = line_lgd_str;
    A = mean(line_ydata, 1, 'omitnan');
    A = string(round(A,2,'significant'));
    T = array2table(A, 'RowNames', RowNames, 'VariableNames', VariableNames);
    disp(T)
    
end
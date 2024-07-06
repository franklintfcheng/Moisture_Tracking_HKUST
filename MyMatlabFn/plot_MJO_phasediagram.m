function plot_MJO_phasediagram(varargin)
        
    pnames = {'RMM_source', 'RMM1', 'RMM2', 'dates', 'colmap', 'x_lim', 'y_lim', 'labelgap', 'textshift',...
              'fontsize', 'fontsize_phase', 'fontsize_lgd', 'fontsize_title', 'lgd_location',...
              'pos_title', 'savefig', 'trans_bg', 'plotting_folder'};

    dflts  = cell(1, length(pnames));

             [RMM_source, RMM1,  RMM2, dates, colmap, x_lim, y_lim, labelgap, textshift,...
              fontsize, fontsize_phase, fontsize_lgd, fontsize_title, lgd_location,...
              pos_title, savefig, trans_bg, plotting_folder]...
                           = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 20 Feb 2024
    %/=====================================================================


    %/ Defaults
    ticklength = [0.01, 0.01];
    color_quad = [0 0 0];

    %/ Makesure dates in the format of datetime
    if isempty(dates)
        error('dates is empty!');
    elseif ~isdatetime(dates) && isinteger(dates)
        dates_dt = int2datetime(dates, 'yyyyMMdd');
    else
        dates_dt = dates;
    end

    if isempty(RMM1) || isempty(RMM2)
        error('RMM1 or RMM2 is empty!');
    end

    uni_mths = unique(dates_dt.Month,'stable');
    

    if isempty(pos_title)
        pos_title = [0.2, 1.005, 1, 0];
    end


    figure
    set(gcf, 'Color', 'w');
    set(gca, 'Color', 'w');

    %/ Set the quadrants
    plot(x_lim, y_lim,      'linestyle', '--', 'color', color_quad); hold on;
    plot(x_lim, flip(y_lim),'linestyle', '--', 'color', color_quad); hold on;
    plot(x_lim, [0, 0],     'linestyle', '--', 'color', color_quad); hold on;  %/ Avoid using xline/yline, since it will always overlay of other layers
    plot([0, 0], y_lim,     'linestyle', '--', 'color', color_quad); hold on;  %/ Avoid using xline/yline, since it will always overlay of other layers
    poly = draw_polyshape('type', 'circle', 'center', [0,0], 'radius', 1);
    plot(poly, 'FaceColor','w','FaceAlpha', 1); hold on;

    %/ Label the day numbers
    text(RMM1-textshift, RMM2+textshift, string(dates_dt.Day), 'HorizontalAlignment', 'center');

    %/ Color the monthly RMM data points (NOTE: uni_mths can be from 11 to 4)
    hseg = [];
    for m = 1:length(uni_mths)
        if m == length(uni_mths)
            ind = find(dates_dt.Month == uni_mths(m));
        else
            ind = find(dates_dt.Month == uni_mths(m)); 
            ind = [ind; ind(end)+1];  %/ Append the next index so as to link the monthly data segments
        end
        hseg{m} = plot(RMM1(ind), RMM2(ind), 'linestyle', '-', 'marker', 'o', 'color', colmap(m,:), 'markerfacecolor', colmap(m,:));    
        % hseg{m} = plot(RMM1(ind), RMM2(ind), 'linestyle', '-', 'marker', 'o', 'color', colmap(uni_mths(m),:), 'markerfacecolor', colmap(uni_mths(m),:));
        drawnow;
    end
    box on;
    axis equal;
    xlim(x_lim);
    ylim(y_lim);
    xlabel('RMM1');
    ylabel('RMM2');
    set(gca,'fontsize', fontsize, 'TickLength',ticklength);  %/ Make the tick length consistent

    %/ Add xticks on the top and yticks on the right (no need now)
    % a2 = axes;
    % set(a2,'Color','none')
    % set(a2,'XAxisLocation','top', 'xlim', x_lim);  
    % set(a2,'YAxisLocation','right', 'ylim', y_lim);  
    % set(a2,'TickLength',ticklength/2);  %/ Somehow the added axes always show a doubly longer length of ticks. Fix it manually.
    % hold on;
    % xlabel('Phase 7 (Western Pacific) Phase 6');
    % ylabel('Phase 4 (Maritime) Phase 5');

    %/ Label the phase quadrants
    text(0, y_lim(1)+labelgap,       'Indian Ocean',       'fontsize', fontsize, 'Color', 'k', 'HorizontalAlignment', 'center');
    text(0, y_lim(2)-labelgap,       'Western Pacific',    'fontsize', fontsize, 'Color', 'k', 'HorizontalAlignment', 'center');
    hx1 = text(x_lim(1)+labelgap, 0, 'Western Hem.',       'fontsize', fontsize, 'Color', 'k', 'HorizontalAlignment', 'center');
    hx2 = text(x_lim(2)-labelgap, 0, 'Maritime Continent', 'fontsize', fontsize, 'Color', 'k', 'HorizontalAlignment', 'center');
    set(hx1,'Rotation',90);
    set(hx2,'Rotation',270);
    
    text(x_lim(1)/3, y_lim(1)+labelgap*2, '2', 'fontsize', fontsize_phase, 'Color', 'k', 'HorizontalAlignment', 'center');
    text(x_lim(2)/3, y_lim(1)+labelgap*2, '3', 'fontsize', fontsize_phase, 'Color', 'k', 'HorizontalAlignment', 'center');
    text(x_lim(1)/3, y_lim(2)-labelgap*2, '7', 'fontsize', fontsize_phase, 'Color', 'k', 'HorizontalAlignment', 'center');
    text(x_lim(2)/3, y_lim(2)-labelgap*2, '6', 'fontsize', fontsize_phase, 'Color', 'k', 'HorizontalAlignment', 'center');

    text(x_lim(2)-labelgap*2, y_lim(1)/3, '4', 'fontsize', fontsize_phase, 'Color', 'k', 'HorizontalAlignment', 'center');
    text(x_lim(2)-labelgap*2, y_lim(2)/3, '5', 'fontsize', fontsize_phase, 'Color', 'k', 'HorizontalAlignment', 'center');
    text(x_lim(1)+labelgap*2, y_lim(1)/3, '1', 'fontsize', fontsize_phase, 'Color', 'k', 'HorizontalAlignment', 'center');
    text(x_lim(1)+labelgap*2, y_lim(2)/3, '8', 'fontsize', fontsize_phase, 'Color', 'k', 'HorizontalAlignment', 'center');

    %/ Legend showing the correspoding month of the segment color
    str_mths = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
    legend([hseg{1:end}], str_mths(uni_mths), 'fontsize', fontsize_lgd, 'location', lgd_location,...
            'color', 'none', 'edgecolor', 'none');
    
    %/ Title
    titlename = sprintf('%s to %s', string(dates_dt(1), 'yyyy-MMM-dd'), string(dates_dt(end), 'yyyy-MMM-dd'));
    figname   = sprintf('phase_diag_%s_%s_to_%s', RMM_source, string(dates_dt(1), 'yyyy-MMM-dd'), string(dates_dt(end), 'yyyy-MMM-dd'));

    annotation( 'textbox', 'String', titlename, 'Color', 'k', ...
                    'FontSize', fontsize_title, 'Units', 'normalized', 'EdgeColor', 'none', ...
                    'Position', pos_title)
    hold on; drawnow; pause(3);

    if savefig
        final_FigName = fullfile(plotting_folder, strcat(figname, '.pdf'));
        fprintf('*** Saving figure into: %s ***\n', final_FigName);

        %/ Cropping often causes unnecessary inconsistency in fontsize. Disable it.
        if trans_bg
            export_fig(final_FigName,'-pdf','-painters', '-nocrop', '-transparent'); %/ transparent bg (do not auto-crop; it would cut out the annotation object.
        else
            export_fig(final_FigName,'-pdf','-painters', '-c[inf, inf, inf, inf]'); %, '-nocrop'); 
        end
    end


end
function plot_hist_cdf(varargin)
    
    % create a set of valid parameters and their default value
    pnames = {'y1', 'y2', 'nbins', 'y1_binlimit', 'y2_binlimit', 'var_label', 'y1name', 'y2name', 'y_lim',...
              'perform_KS_test', 'set_prcntl', 'basic_stat', 'zoomin_mode', 'y1_col', 'y2_col', 'facealpha', 'draw_pdf', 'draw_cdf',...
              'fontsize', 'title_fontsize', 'linewidth', 'titlename', 'FigName_underscore', 'plotting_folder', 'savefig'};
    
    dflts  = cell(1, length(pnames));
    
    %/ parse function arguments
    [         y1,    y2,   nbins,   y1_binlimit,   y2_binlimit,   var_label,   y1name,   y2name,   y_lim,...
              perform_KS_test,  set_prcntl,    basic_stat,    zoomin_mode, y1_col,   y2_col,    facealpha,   draw_pdf,   draw_cdf,...
              fontsize,    title_fontsize,   linewidth,   titlename,   FigName_underscore,   plotting_folder,   savefig] ...
                 = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%

    %/ NOTE:   y1 = target var
    if isempty(y1_col)
        % y1_col = [0 0 0]./255;
        y1_col = [0, 32, 63]./255;
    end
    if isempty(y2_col)
        % y2_col = [255 204 204]./255;
        y2_col = [173, 239, 209]./255;
    end
    ax1_col   = [0 0 0];
    ax2_col   = [0 0 0];

    if isempty(fontsize)
        fontsize  = 12;
    end

    if isempty(linewidth)
        linewidth = 1.5;
    end
    
    h = {}; h_cdf = {};
    
    figure
    set(gcf,'color','w'); %/ set figure bg color to be white
    
    if draw_cdf || draw_pdf
        yyaxis left
    end
    if isempty(nbins)
        nbins = fix(length(y1)/5); %/ By default
    end

    %/ Convert into column vectors for consistency
    if size(y1, 1) == 1
        y1 = y1';  
    end
    if size(y2, 1) == 1
        y2 = y2';
    end

    if isempty(y1_binlimit) && isempty(y2_binlimit)
        y1_binlimit = [min([y1; y2], [], 'omitnan'), max([y1; y2], [], 'omitnan')];
        y2_binlimit = y1_binlimit;
    end

    if isempty(facealpha)
        facealpha = 0.6;
    end

    h{1} = histogram(y1, nbins, 'BinLimits', y1_binlimit, 'FaceAlpha', facealpha, 'Facecolor', y1_col, 'edgecolor', 'k', 'linewidth', linewidth*1.5);
    xlabel(var_label);  %/ e.g., K, m/s, W m^{-2}
    ylabel('Frequency'); 
    
    %/ whether to draw the 5th and 95th percentile lines
    if ~isempty(set_prcntl)
        if any(set_prcntl > 1)
            error('set_prcntl must be between 0 and 1!');
        end
        
        y1_lower_prcntl = quantile(y1, set_prcntl(1));
        y1_upper_prcntl = quantile(y1, set_prcntl(2));
        xline(y1_lower_prcntl, 'r--', 'linewidth', linewidth);    
        xline(y1_upper_prcntl, 'r--', 'linewidth', linewidth);    

        if ~isempty(y2)
            y2_lower_prcntl = quantile(y2, set_prcntl(1));
            y2_upper_prcntl = quantile(y2, set_prcntl(2));
            xline(y2_lower_prcntl, 'r--', 'linewidth', linewidth);    
            xline(y2_upper_prcntl, 'r--', 'linewidth', linewidth);   
        end
    end

    if ~isempty(y2)
        hold on;
        h{2} = histogram(y2, nbins, 'BinLimits', y2_binlimit, 'FaceAlpha', facealpha, 'Facecolor', y2_col, 'edgecolor', 'k', 'linewidth', linewidth*1.5);
        ylabel('Frequency');
        if ~isempty(y_lim)    
            ylim(y_lim);   
        end

        %/ Perform two-sample KS test
        if perform_KS_test
            % if isempty(asterisk)    
            %     asterisk    = {'(*)', '(**)', '(***)'};    
            % end
            % if isempty(alpha_level) 
            %     alpha_level = [0.05,    0.01,   0.001];    
            % end

            [~,pval] = kstest2(y1,y2);
            str_sig = sprintf(' (p=%.2f)', pval);
            
            % str_sig = {''};
            % for i = 1:length(alpha_level)
            %     ind = find(pval <= alpha_level(i));
            %     str_sig(ind) = repmat(asterisk(i), length(ind), 1);
            % end
            % 
            % if isempty(str_sig{:})
            %     str_sig = sprintf(' (p=%.2f)', pval);
            % else
            %     str_sig = sprintf(' %s', str_sig{:});  %/ Append a leading space
            % end
        else
            str_sig = '';
        end
    end
    hold on;
    if isempty(y_lim)   
        y_lim = [0, max([h{1}.BinCounts,h{2}.BinCounts])*1.5];
    end

    y_ticks = 1:1:20; %/ frequency, make it integer

    ylim(y_lim);
    yticks(y_ticks);
    yticklabels(y_ticks);


    if draw_pdf
        disp('draw also empirical pdf (by kernel distribution)');

        yyaxis right  %/ rmb this!
        x_intvl = 0.1;
        a1      = fitdist(y1,'Kernel','Width',4);
        y1_x    = y1_binlimit(1):x_intvl:y1_binlimit(end);
        y1_pdf  = pdf(a1,y1_x);
        plot(y1_x,y1_pdf,'Color', y1_col, 'LineStyle', '-', 'LineWidth', linewidth*2);
        hold on;
        ax = gca;
        ax.YAxis(1).Color = ax1_col;
        ax.YAxis(2).Color = ax2_col;

        if ~isempty(y2)
            yyaxis right
            y2_x   = y2_binlimit(1):x_intvl:y2_binlimit(end);
            a2     = fitdist(y2,'Kernel','Width',4);
            y2_pdf = pdf(a2,y2_x);
            plot(y2_x,y2_pdf,'Color', y2_col, 'LineStyle', '-', 'LineWidth', linewidth*2);
        end
        hold on;
        xlim([min([y1_binlimit,y2_binlimit]), max([y1_binlimit,y2_binlimit])])
        ylim([0, max([y1_pdf, y2_pdf])*1.5])
    end

    if draw_cdf
        disp('draw also cdf');
        yyaxis right %/ rmb this!
        [F, xi] = ksdensity(y1, 'Support', 'unbound','Function','cdf', 'NumPoints',length(y1));
        h_cdf{1} = plot(xi, F, 'LineWidth', linewidth*1.5, 'color', 'k', 'linestyle', '-');
        
        ax = gca;
        ax.YAxis(1).Color = ax1_col;
        ax.YAxis(2).Color = ax2_col;

        if ~isempty(y2)
            yyaxis right
            [F, xi2] = ksdensity(y2, 'Support', 'unbound','Function','cdf', 'NumPoints',length(y2));
            h_cdf{2} = plot(xi2, F, 'LineWidth',  linewidth*1.5, 'color', 'k', 'linestyle', '--');
        end
    
        ylabel('Cumulative densify function');
        xlim([min([y1_binlimit,y2_binlimit]), max([y1_binlimit,y2_binlimit])])
        ylim([0, 1])
    end
    box on;
    
    if ~isempty(y2)
        if draw_cdf
            str_ldg = {y1name, strcat(y2name, str_sig), sprintf('CDF - %s', y1name), sprintf('CDF - %s', y2name)};
        else
            str_ldg = {y1name, strcat(y2name, str_sig)};
        end
        legend(str_ldg, 'Location','northwest','edgecolor','none', 'fontsize', fontsize);
    end
    set(gca, 'fontsize', fontsize, 'linewidth', linewidth*2);
    

    %/ whether to show basic statistics (min, median, mean, max, sample size n)
    if basic_stat
        str_y1_min          = sprintf('y1: min    = %.1f', min(y1, [], 'omitnan'));
        str_y1_mean         = sprintf('y1: mean   = %.1f', mean(y1, 'omitnan'));
        str_y1_median       = sprintf('y1: median = %.1f', median(y1, 'omitnan'));
        str_y1_max          = sprintf('y1: max    = %.1f', max(y1, [], 'omitnan'));
        str_n               = sprintf('y1: n      = %d',   length(y1));
        
        if ~isempty(set_prcntl)
            str_y1_lower_prcntl  = sprintf('y1: %d-p   = %.1f', set_prcntl(1)*100, y1_lower_prcntl);
            str_y1_upper_prcntl  = sprintf('y1: %d-p   = %.1f', set_prcntl(2)*100, y1_upper_prcntl);
            str_y1_stat =  {str_y1_min, str_y1_lower_prcntl, str_y1_mean, str_y1_median, str_y1_upper_prcntl, str_y1_max, str_n};
        else
            str_y1_stat =  {str_y1_min, str_y1_mean, str_y1_median, str_y1_max, str_n};
        end

        %/ create a textbox [x_from_left, y_from_bottom, text_box_height, text_box_width]
        if ~isempty(y2)
            textbox_y1_pos = [0.50, 0.65, 0.25, 0.25];
        else
            textbox_y1_pos = [0.75, 0.65, 0.25, 0.25];
        end
        textbox_fontsize = fontsize*0.75;
        annotation('textbox', textbox_y1_pos, 'String', str_y1_stat, 'edgecolor', 'none', 'fontsize', textbox_fontsize)

        if ~isempty(y2)
            str_y2_min          = sprintf('y2: min    = %.1f', min(y2, [], 'omitnan'));
            str_y2_mean         = sprintf('y2: mean   = %.1f', mean(y2, 'omitnan'));
            str_y2_median       = sprintf('y2: median = %.1f', median(y2, 'omitnan'));
            str_y2_max          = sprintf('y2: max    = %.1f', max(y2, [], 'omitnan'));
            str_n               = sprintf('y2: n      = %d',   length(y2));
            
            if ~isempty(set_prcntl)
                str_y2_lower_prcntl  = sprintf('y2: %d-p   = %.1f', set_prcntl(1)*100, y2_lower_prcntl);
                str_y2_upper_prcntl  = sprintf('y2: %d-p   = %.1f', set_prcntl(2)*100, y2_upper_prcntl);
                str_y2_stat =  {str_y2_min, str_y2_lower_prcntl, str_y2_mean, str_y2_median, str_y2_upper_prcntl, str_y2_max, str_n};
            else
                str_y2_stat =  {str_y2_min, str_y2_mean, str_y2_median, str_y2_max, str_n};
            end

            textbox_y2_pos = [0.7, 0.65, 0.25, 0.25];
            annotation('textbox', textbox_y2_pos, 'String', str_y2_stat, 'edgecolor', 'none', 'fontsize', textbox_fontsize)
        end
    end
    
    titlename = strrep(titlename, '_', ' ');
    % title(titlename, 'fontsize', fontsize);

    if isempty(title_fontsize)
        title_fontsize = fontsize;
    end
    title_pos = [0.1, 1, 0.9, 0];
    annotation('textbox', title_pos, 'String', titlename, 'edgecolor', 'none', 'fontsize', title_fontsize)

    % hold on;
    % drawnow;
    pause(2);

    % set(gcf, 'visible', 'off')
    if savefig
        FigName_underscore = strrep(FigName_underscore, ' ', '_');
        if facealpha ~= 1
            png_dpi = 300;
            final_figname = char(fullfile(plotting_folder, strcat(FigName_underscore,'.png')));

            f = gcf;
            exportgraphics(f,final_figname,'Resolution',300, 'BackgroundColor','none','ContentType','vector');  %/ Use this; export_fig does not work--it outputs blank figures
            % saveas(gcf,final_figname)
            % export_fig(final_figname, '-r300','-png','-opengl', '-nocrop', '-transparent');
        else
            final_figname = char(fullfile(plotting_folder, strcat(FigName_underscore,'.pdf')));

            f = gcf;
            exportgraphics(f,final_figname,'Resolution',300, 'BackgroundColor','none','ContentType','vector');  %/ Use this; export_fig does not work--it outputs blank figures
            % export_fig(final_figname,'-pdf','-painters','-nocrop', '-transparent');
        end
        fprintf('!!! Figure saved into %s !!!\n', final_figname);
    end
    % set(gcf, 'visible', 'on')
    
    if zoomin_mode  %/ a zoom-in view of histogram 1
        figure
        set(gcf,'color','w');
        histogram(y1, fix(length(y1)/30), 'BinLimits',[min(y1), max(y1)], 'FaceAlpha', 1, 'Facecolor', ax1_col);
        xlabel(var_label);
        xlim([min(y1), max(y1)])
        ylabel('Frequency');
        ylim([0 10]); 
        
        if savefig
            FigName_underscore = strrep(FigName_underscore, ' ', '_');
            export_fig(char(fullfile(plotting_folder, strcat(FigName_underscore,'_zoomin.pdf'))),'-pdf','-painters','-nocrop', '-transparent');
        end
    end

end
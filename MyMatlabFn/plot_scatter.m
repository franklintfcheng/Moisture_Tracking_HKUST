%%
function plot_scatter(varargin)

    pnames = {'X', 'Y', 'X_label', 'Y_label', 'sc_colmap', 'titlename', 'figname', 'plotting_folder', 'savefig'};
          
    dflts  = cell(size(pnames));

    [          X,   Y,   X_label,   Y_label,   sc_colmap,  titlename,   figname,  plotting_folder,   savefig] ...
           = internal.stats.parseArgs(pnames, dflts, varargin{:});

%%
       
    figure
    set(gcf,'color','w'); 
    set(gcf,'position', [600, 50, 800, 650]) %/ x, y, width, height
%         axis equal
%         markersize = 60;
%         MarkerEdgeColor = 'k';
%         MarkerFaceColor = [51 153 204]./255;
%         scatter(X, Y, markersize, 'MarkerEdgeColor', MarkerEdgeColor, 'MarkerFaceColor', MarkerFaceColor);
    
    fontsize = 20;   
    linewidth = 2;
    X_range = max(X)-min(X);
    Y_range = max(Y)-min(Y);
    Xlim = [min(X)-X_range*0.1, max(X)+X_range*0.1];
    Ylim = [min(Y)-Y_range*0.1, max(Y)+Y_range*0.1];

    xline(0, 'k--', 'linewidth', linewidth*0.5);
    yline(0, 'k--', 'linewidth', linewidth*0.5);
    hold on;

    %/ Fit linear regression line
    mdl         = fitlm(X, Y, 'linear');
    intercept   = table2array(mdl.Coefficients(1,1));
    beta        = table2array(mdl.Coefficients(2,1));
    pval        = table2array(mdl.Coefficients(2,end));
    Rsquared    = mdl.Rsquared.Ordinary;
    Y_c         = X.*beta + intercept;
    X_regline   = [min(X), max(X)];
    Y_c_regline = X_regline.*beta + intercept;
    if pval < 0.05
        reg_linest = '-';  %/ sig.
    else
        reg_linest = '--'; %/ insig.
    end
    reg_label = sprintf('Slope=%.2G (p=%.1G), R-squared=%.2G', beta, pval, Rsquared);
    h = plot(X_regline, Y_c_regline, 'color', 'k', 'linewidth', linewidth*1.5);
    set(h, 'LineStyle', reg_linest);
    hold on;
    fprintf('*** Correlation betn. X and Y = %.2f ***\n', corr(X, Y, 'rows','complete')) %/ just for ref.
    %     plot(mdl, 'linewidth', linewidth)

    %/ Add scatters
    scatter(X, Y, [], sc_colmap, 'filled');
%     colorbar
% %         colormap hsv
%     sc_colmap = hsv(73);  %/ ensure there are 73 colors.
%     colormap(cmap);
%     caxis([0.5 73.5])  %/ ensure the ticks are always at the middle of the colors.
%     textfit(X, Y, num2cell(pentads), 'k' ,'VerticalAlignment','top','HorizontalAlignment','right',...
%             'fontsize', fontsize*0.7)
%     hold on;

    %/ Legend, labels, limits
    legend(h, reg_label, 'location', 'best');
    xlabel(strrep(X_label, '_', ' '));
    ylabel(strrep(Y_label, '_', ' '));
    xlim(Xlim);
    ylim(Ylim);
    set(gca, 'fontsize', fontsize)
    box on;
    grid off;
    hold on;
    title(titlename);
    
    if savefig
        export_fig(strcat(plotting_folder, figname),'-pdf','-painters', '-c[inf, inf, inf, inf]');
%             export_fig(strcat(plotting_folder, figname),'-png', '-r300', '-painters', '-c[inf, inf, inf, inf]');
    end

end
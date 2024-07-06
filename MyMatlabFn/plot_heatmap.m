%%
function plot_heatmap(varargin)

    pnames = {'X', 'Y', 'heatmap_data', 'texts', 'thres_white_color', 'contf_levels', 'colmap', 'X_ticks', 'X_ticklabels', 'Y_ticks', 'Y_ticklabels',... 
              'X_tickangle', 'Y_tickangle', 'xaxisLocation', 'add_diagline', 'scatter_X', 'scatter_Y', 'markersize', 'fontsize', 'linewidth',...
              'panel_wid', 'panel_hgt', 'titlename', 'FigName_underscore', 'plotting_folder', 'savefig', 'fig_mat'};
          
    dflts  = cell(size(pnames));

    [          X,   Y,   heatmap_data,   texts,   thres_white_color,   contf_levels,   colmap,   X_ticks,   X_ticklabels,   Y_ticks,   Y_ticklabels,...
               X_tickangle,   Y_tickangle,   xaxisLocation,    add_diagline,   scatter_X,   scatter_Y,   markersize,   fontsize,   linewidth,...
               panel_wid,   panel_hgt,   titlename,   FigName_underscore,  plotting_folder,   savefig,   fig_mat] ...
           = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    %/ NOTE:
    %/      This function uses surf() to mimic heatmap. -> greater flexibility.
    
    if isempty(fontsize)    fontsize = 12;    end
    if isempty(markersize)  markersize = 12;  end
    fprintf('Mean value of heatmap_data = %.2f \n', mean(heatmap_data, 'all'));
    
    figure
    set(gcf,'color','w'); 
    set(gcf,'position', [600, 50, panel_wid, panel_hgt]) %/ x, y, width, height
    axis equal
    
%     surf(X, Y, heatmap_data);  
%     view(0,90);  %/ view the XY-plane
%     shading faceted;                       %/ get black outlines on the cubes.
    
    %/ imagesc is better than surf when overlaying other plots.
    h = imagesc(X(1,:), Y(:,1), heatmap_data); 
    set(gca, 'YDir', 'normal');              %/ since imagesc has already used 'YDir' = 'reverse'!
    hold on;
    
    %/ Manually outline the grids (since shading faceted does not work for imagesc().
    nx = size(heatmap_data,2);
    ny = size(heatmap_data,1);
    edge_x = repmat((0:nx)+0.5, ny+1,1);
    edge_y = repmat((0:ny)+0.5, nx+1,1).';
    plot(edge_x ,  edge_y,   'k') % vertical lines
    plot(edge_x.', edge_y.', 'k') % horizontal lines

    xlim([min(X, [], 'all')-0.5, max(X, [], 'all')+0.5]);
    ylim([min(Y, [], 'all')-0.5, max(Y, [], 'all')+0.5]);
    
    if ~isempty(X_ticks)        xticks(X_ticks);            end
    if ~isempty(Y_ticks)        yticks(Y_ticks);            end
    if ~isempty(X_ticklabels)   xticklabels(X_ticklabels);  end
    if ~isempty(Y_ticklabels)   yticklabels(Y_ticklabels);  end
    if ~isempty(X_tickangle)    xtickangle(X_tickangle);    end
    if ~isempty(Y_tickangle)    ytickangle(Y_tickangle);    end
    
    colormap(colmap); 
    caxis([min(contf_levels), max(contf_levels)]);
    cb              = colorbar;
    cbar_YTick      = contf_levels(2:1:end-1);
    cbar_YTickLabel = round(cbar_YTick,2);
    set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', fontsize, 'Location', 'eastoutside')   
    set(get(cb,'Title'),'String','', 'Fontsize', fontsize)
    hold on;
    if add_diagline
        h_diag = refline(1,0);
        h_diag.Color = 'k';
        h_diag.LineWidth = linewidth;
        hold on;
    end

    %/ Write values as texts 
    [X_2D, Y_2D] = meshgrid(X,Y);
    X_2Dto1D = reshape(X_2D, [], 1);
    Y_2Dto1D = reshape(Y_2D, [], 1);
    
    heatmap_data_2Dto1D = reshape(heatmap_data, [], 1);
    if isempty(texts)
        texts_2Dto1D = string(heatmap_data_2Dto1D); %/ then just puts texts based on all non-nan values.
    else
    	texts_2Dto1D = reshape(texts, [], 1);
    end
    text_fontsize = fontsize*0.8;
%     a = 7*90/950;
%     text_fontsize = panel_hgt/ny;   %/ vary with the # of rows and the height of the panel.
    h_text = text(X_2Dto1D, Y_2Dto1D, texts_2Dto1D, 'HorizontalAlignment', 'center', 'fontsize', text_fontsize);
    

    %/ Draw texts in white for |corr| > 0.5 if the heatmap data is abt corr.
    if ~isempty(thres_white_color)
        ind = find(abs(heatmap_data_2Dto1D) > thres_white_color); 
        for i = 1:length(ind)
            h_text(ind(i)).Color = [1 1 1];
        end
    end
    
    %/ The easiest solution to hide nan cells.
    set(h, 'AlphaData', ~isnan(heatmap_data)) 
    
    if ~isempty(scatter_X) && ~isempty(scatter_Y)
        scatter(scatter_X, scatter_Y, markersize,'rs', 'linewidth', linewidth)
    end
%     if ~isempty(heatmap_marker)
%         scatter(X,Y, markersize,'ro', 'linewidth', linewidth)
% %         plot3(X,Y,heatmap_marker,'ro', 'markersize', markersize, 'linewidth', linewidth)
%     end
    set(gca, 'YDir','reverse')      %/ by default
    
    if isempty(xaxisLocation)  xaxisLocation = 'top';  end    
    set(gca,'xaxisLocation', xaxisLocation) 
    set(gca,'TickLength',[0 0])
    set(gca, 'FontSize', fontsize);
    title(titlename, 'FontSize', fontsize);

    if savefig
        if isequal(fig_mat, 'pdf')
            export_fig(char(strcat(plotting_folder, FigName_underscore, '.pdf')),'-pdf','-painters', '-c[inf, inf, inf, inf]'); %, '-nocrop', '-transparent');
        elseif isequal(fig_mat, 'png')
            export_fig(char(strcat(plotting_folder, FigName_underscore, '.png')),'-r300','-png','-opengl', '-c[inf, inf, inf, inf]');
        else
            error('wrong input of fig_mat!')
        end
        fprintf('!!! Plot is saved: %s !!! \n', string(FigName_underscore));
    end
end

% function quickplot_heatmap(varargin)
% 
%     pnames = {'heatmap_data', 'x_labels', 'y_labels', 'contf_levels', 'colmap',...
%               'fontsize', 'titlename', 'FigName_underscore', 'plotting_folder'};
%     dflts  = cell(size(pnames));
% 
%     [          heatmap_data, x_labels, y_labels, contf_levels, colmap,...
%                fontsize, titlename, FigName_underscore, plotting_folder] ...
%         = internal.stats.parseArgs(pnames, dflts, varargin{:});
% 
%     if isempty(fontsize)   fontsize = 12;    end
%     
%     figure
%     set(gcf,'color','w'); 
%     h = heatmap(x_labels, y_labels, heatmap_data);
%     colormap(colmap);
%     caxis([min(contf_levels), max(contf_levels)]);
%     %             set(h, 'linewidth', linewidth);           %/ change linewidth of heatmap is not possible..
%     
%     cb              = colorbar;
%     cbar_YTick      = contf_levels(2:1:end-1);
%     cbar_YTickLabel = round(cbar_YTick,2);
% 
%     set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', fontsize,...
%         'Location', 'eastoutside')   
%     set(get(cb,'Title'),'String','', 'Fontsize', fontsize)
% 
%     %/ A way to put x-label of heatmap on the top.
%     ax = gca;
%     axp = struct(ax);                           %/ you will get a warning
%     axp.Axes.XAxisLocation = 'top';
% 
% 
% %     s = struct(h);
% %     s.XAxis.TickLabelRotation = 90;  % vertical
%     axp.XAxis.TickLabelRotation = 20;  % angled
% 
%     set(h, 'FontSize', fontsize);
%     
%     title(titlename);
%     
%     if ~isempty(FigName_underscore)
%         export_fig(char(strcat(plotting_folder, FigName_underscore, '.pdf')),'-pdf','-painters', '-c[inf, inf, inf, inf]', '-nocrop', '-transparent');
%         fprintf('!!! Plot is saved: %s !!! \n', string(FigName_underscore));
%     end
% 
% end
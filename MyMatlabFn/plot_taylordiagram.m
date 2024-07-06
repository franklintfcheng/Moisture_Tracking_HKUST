function STATS = plot_taylordiagram(varargin)
        
    pnames = {'data', 'name_list', 'color_var', 'marker_var', 'lgd_location', 'show_text',...
              'fontsize', 'lgd_fontsize', 'linewidth', 'markersize', 'markerlinewi', 'show_legend', 'show_legend_only',...
              'figname', 'fig_fmt', 'png_dpi', 'savefig', 'plotting_folder'};

    dflts  = cell(1, length(pnames));

              [data,   name_list,   color_var,   marker_var,   lgd_location,   show_text,...
               fontsize,   lgd_fontsize,  linewidth, markersize,   markerlinewi,   show_legend,   show_legend_only,...
               figname,   fig_fmt,   png_dpi,   savefig,   plotting_folder]...
                           = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 5 Feb 2024
    %/
    %/ Description: This function is developed based on 'demo1_taylor' of the 
    %/              Matlab Exchange function by Zhaoxu Liu / slandarer.
    %/              See https://www.mathworks.com/matlabcentral/fileexchange/130889-taylor-diagram-class 
    %/=====================================================================

    % Calculate STD, RMSD and Correlation
    [~, ncolor] = size(data);  %/ row: samples, column: types
    STATS = zeros(4,ncolor);
    for i = 1:size(data,2)
        STATS(:,i) = SStats(data(:,1),data(:,i));
    end
    STATS(1,:) = [];
    
    % Set color and marker list
    if size(marker_var, 1)
        marker_var = marker_var'; %/ make it a column vector
    end

    ncolor  = size(color_var, 1);
    nmarker = length(marker_var);

    if ncolor*nmarker < ncolor
        error('The total # of combinations of color_var and marker_var is smaller than the data columns!');
    else
        color_list = nan(ncolor*nmarker, 3);
        for i = 1:size(color_var, 1)
            color_list((1:nmarker)+(i-1)*nmarker,:) = repmat(color_var(i,:), nmarker, 1);
        end
        facecolor_list = repmat({'none'}, ncolor*nmarker, 1);
        marker_list = repmat(marker_var, ncolor, 1);
    end
    color_list(1,:)     = [0 0 0];  %/ Use solid black to indicate the obs data
    facecolor_list{1,:} = [0 0 0];

    fontname = 'Times New Roman';
    name_list = strrep(name_list, '_', ' ');


    % Create taylor axes 
    % figure('Units','normalized','Position',[.2,.1,.52,.72]);
    % figure('Units','normalized','Position',[0, .1,.72,.72]);
    figure('Units','normalized','Position',[0, .1,1.0,.72]);
    set(gcf, 'Color', 'w');
    set(gca, 'Color', 'w');
    TD = STaylorDiag(STATS); %/ Open the code to modify the position of the labels
    
    % Scatter plots
    for i = 1:size(data,2)
        TD.SPlot(STATS(1,i),STATS(2,i),STATS(3,i),'Marker',marker_list{i},'MarkerSize', markersize,...
            'Linewidth', markerlinewi, 'Color',color_list(i,:),'MarkerFaceColor', facecolor_list{i,:});
    end
    
    % Show text by annotation
    if show_text
        TD.SText(STATS(1,1),STATS(2,1),STATS(3,1),{'reference';' '},'FontWeight','bold',...
            'FontSize',fontsize,'FontName',fontname,'Color',color_list(1,:),...
            'VerticalAlignment','bottom','HorizontalAlignment','center');
        
        for i = 1:size(data,2)
            TD.SText(STATS(1,i),STATS(2,i),STATS(3,i),"   "+string(name_list{i}),'FontWeight','bold',...
            'FontSize',fontsize,'FontName',fontname,...
            'VerticalAlignment','middle','HorizontalAlignment','left');
        end
    end
    
    if show_legend || show_legend_only
        if isempty(lgd_fontsize)
            lgd_fontsize = fontsize*0.8;
        end
        lgd_pos = [0.82    0.32    0.08    0.40];
        h = legend(name_list,'FontSize',lgd_fontsize,'FontName',fontname, 'Position', lgd_pos); %, 'location', lgd_location);
        % h.Position
        % h.Position(1) = h.Position(1) + 0.15;
        % h.Position(2) = h.Position(2) - 0.15;
        % h.Position
        % h = legend(name_list,'FontSize',lgd_fontsize, 'FontName',fontname);
        % h.Position(1) = h.Position(1) + 0.2;
        % h.Position(2) = h.Position(2) - 0.2;
        legend boxoff;
        if show_legend_only
            axis off %/ hide axis
            set(gca,'visible','off')
        end
    end
    % TD = [];
    % Create the plot
    % ax = axes(); 
    % hold on
    % h(1) = plot(linspace(1,5,25), rand(1,25), 'ro', 'DisplayName', 'foo');
    % h(2) = plot(1:5, rand(1,5), 'b-', 'DisplayName', 'bar');
    % % copy the objects
    % hCopy = copyobj(h, ax); 
    % % replace coordinates with NaN 
    % % Either all XData or all YData or both should be NaN.
    % set(hCopy(1),'XData', NaN', 'YData', NaN)
    % set(hCopy(2),'XData', NaN', 'YData', NaN)
    % % Note, these lines can be combined: set(hCopy,'XData', NaN', 'YData', NaN)
    % % To avoid "Data lengths must match" warning, assuming hCopy is a handle array, 
    % % use arrayfun(@(h)set(h,'XData',nan(size(h.XData))),hCopy)
    % % Alter the graphics properties
    % hCopy(1).MarkerSize = 15; 
    % hCopy(1).LineWidth = 2;
    % hCopy(2).LineWidth = 3; 
    % % Create legend using copied objects
    % legend(hCopy)

    %/ From demo2_taylor
    % Set other properties(设置其他属性)
    % TD.set('TickLength',[.015,.05])
    % TD.set('SLim',[0,300])
    % TD.set('RLim',[0,175])
    % TD.set('STickValues',0:50:300)
    % TD.set('SMinorTickValues',0:25:300)
    % TD.set('RTickValues',0:25:175)
    % TD.set('CTickValues',[.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.99])
    

        
    % Set Tick Label(修饰刻度标签)
    TD.set('STickLabelX','fontsize',fontsize)
    TD.set('STickLabelY','fontsize',fontsize)
    TD.set('RTickLabel', 'fontsize',fontsize)
    TD.set('CTickLabel', 'fontsize',fontsize)
        
    % Set Label(修饰标签)
    TD.set('SLabel','fontsize',fontsize*1.25)
    TD.set('CLabel','fontsize',fontsize*1.25)

    % Set Grid(修饰各个网格)
    if ~isempty(linewidth)
        grid_col = [.6 .6 .6];
        TD.set('SGrid','LineWidth',linewidth*0.5)
        TD.set('RGrid','LineWidth',linewidth*0.5)
        TD.set('CGrid','LineWidth',linewidth*0.5);
        
        % Set Axis(修饰各个轴)
        TD.set('SAxis', 'LineWidth',linewidth*1.5)
        TD.set('CAxis', 'LineWidth',linewidth*1.5)
        
        % Set Tick and MinorTick(修饰主次刻度)
        TD.set('STick', 'LineWidth',linewidth*1.5)
        TD.set('CTick', 'LineWidth',linewidth*1.5)
        TD.set('SMinorTick', 'LineWidth',linewidth)
        TD.set('CMinorTick', 'LineWidth',linewidth)
    end

    drawnow; pause(2);
    if savefig
        figname = strrep(figname, ' ', '_');
        final_FigName = char(fullfile(plotting_folder, strcat(figname, '.', fig_fmt)));
        % disp(final_FigName)
        
        if isequal(fig_fmt, 'png')
            if isempty(png_dpi)   
                png_dpi = 300;   
            end
            export_fig(final_FigName,['-r',num2str(png_dpi)],'-png','-opengl',  '-nocrop', '-transparent'); %, '-nocrop');
        else
            export_fig(final_FigName,'-pdf','-painters', '-nocrop', '-transparent'); %/ make background transparent (do not auto-crop; it may cut out part of the figure)
        end
        fprintf('!!! Plot is saved: %s !!! \n', string(final_FigName));
        fprintf('Time cost in saving figure: %f \n', toc);
    end
    set(gcf, 'Color', 'w');
    set(gca, 'Color', 'w');

end
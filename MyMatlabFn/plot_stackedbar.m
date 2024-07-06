function plot_stackedbar(varargin)

    pnames       = {'bar_mode', 'bar_data', 'bar_error_upper', 'bar_error_lower', 'bar_range', 'bar_rangeticks', 'sort_mode',...
                   'show_value', 'show_value_digit', 'text_on_stackedbar', 'text_next_to_stackedbar', 'group_labels', 'bar_labels_intvl', 'stack_labels',...
                    'bar_colmap', 'barplot_margin', 'edgecolor', 'y_labels', 'orientation', 'reverse_data_axis', 'box_on', 'panel_pos',...
                    'fontsize', 'show_lgd', 'lgd_fontsize', 'lgd_position', 'lgd_orient', 'text_fontsize', 'text_yshift', 'linewidth', 'markersize',...
                    'title_pos', 'titlename', 'title_fontsize', 'FigName_underscore', 'plotting_folder', 'savefig'};
    dflts        = cell(length(pnames), 1);
    
    [                bar_mode,   bar_data,   bar_error_upper,   bar_error_lower,   bar_range,   bar_rangeticks,   sort_mode,...
                     show_value,   show_value_digit,   text_on_stackedbar,   text_next_to_stackedbar,   group_labels,   bar_labels_intvl,   stack_labels,...
                     bar_colmap,   barplot_margin,   edgecolor,   y_labels,    orientation,  reverse_data_axis,   box_on,   panel_pos,...
                     fontsize,   show_lgd,  lgd_fontsize,  lgd_position,  lgd_orient, text_fontsize, text_yshift, linewidth, markersize,...
                     title_pos,   titlename,   title_fontsize,  FigName_underscore,   plotting_folder,   savefig] = ...
            internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

    %%
    %/===========================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Jun 19, 2024
    %/
    %/ bar_mode:
    %/     'grouped': display one group for each row of bar_data (group x stacks)
    %/     'stacked': display one bar   for each row of bar_data (bar   x stacks)
    %/===========================================================================

    %/ The dimension format of bar_data: (bars, stacks) => 2D 
    if isempty(bar_mode)    
        bar_mode = 'stacked';     
    elseif ~ismember(bar_mode, {'stacked', 'grouped'})
        error('''bar_mode'' can only take ''stacked'' or ''grouped''!');
    end

    if isempty(stack_labels)
        error('Input ''stack_labels''!');
    end
    
    %/ set bar_colmap
    if isempty(bar_colmap)
        N = length(stack_labels);
        if N > 12
            c1 = brewermap(12,   'paired');
            c2 = brewermap(N-12, 'Dark2');
            bar_colmap = [c1; c2];
        else
            bar_colmap = brewermap(N, 'paired');
        end
    end
    % if isempty(bar_labels_intvl)  bar_labels_intvl = 1;  end
    
    if length(size(bar_data)) > 2
        error('bar_data must be a vector or a 2D matrix only.')
    elseif length(size(bar_error_upper)) > 2
        error('bar_error_upper must be a vector or a 2D matrix only.')
    elseif length(size(bar_error_lower)) > 2
        error('bar_error_lower must be a vector or a 2D matrix only.')
    end
    
    %/ Check consistency between bar_data and bar_error_upper, bar_error_lower in dimensions
    if ~isempty(bar_error_upper) && ~isempty(bar_error_lower) && ...
            (~isequal(size(bar_data), size(bar_error_upper)) || ~isequal(size(bar_data), size(bar_error_lower)))
        error('Inconsistent dimensions between ''bar_data'' and ''bar_error_upper'' or ''bar_error_lower''!');
    end


    %/ Make sure that bar_data's dimension is in groups x stacks
    if ~isvector(bar_data)
        if size(bar_data,1) ~= length(group_labels)
            error('The 1st dim of ''bar_data'' must fit the length of ''group_labels''!');
        end
    elseif ~iscolumn(bar_data)
        bar_data        = bar_data'; 
        bar_error_upper = bar_error_upper';  %/ Trick: []' => [];
        bar_error_lower = bar_error_lower';
    end

    %/ broadcast vars
    bar_data_bc         = bar_data;
    bar_error_upper_bc  = bar_error_upper;
    bar_error_lower_bc  = bar_error_lower;
    group_labels_bc     = group_labels;
    stack_labels_bc     = stack_labels;
    bar_colmap_bc       = bar_colmap;

    %/ sort
    if ~isempty(sort_mode)
        if ~ismember(sort_mode, {'descend', 'ascend'})
            error('Wrong input of ''sort_mode''!');
        end
        [bar_data_bc, I] = sort(bar_data_bc, sort_mode);
        bar_error_upper_bc = bar_error_upper_bc(I);
        bar_error_lower_bc = bar_error_lower_bc(I);
        group_labels_bc    = group_labels_bc(I);
        stack_labels_bc  = stack_labels_bc(I);
        bar_colmap_bc    = bar_colmap_bc(I,:);
    end
    
    [nbars, nstacks] = size(bar_data_bc);
    ind_nve = bar_data_bc < 0;
    ind_pve = bar_data_bc >= 0;
    bar_data_pve = bar_data_bc;  bar_data_pve(ind_nve) = 0;
    bar_data_nve = bar_data_bc;  bar_data_nve(ind_pve) = 0;
    
    if ~isempty(bar_error_upper_bc) && ~isempty(bar_error_lower_bc)
        bar_data_pve_error_upper = bar_error_upper_bc; bar_data_pve_error_upper(ind_nve) = nan;
        bar_data_pve_error_lower = bar_error_lower_bc; bar_data_pve_error_lower(ind_nve) = nan;

        bar_data_nve_error_upper = bar_error_upper_bc; bar_data_nve_error_upper(ind_pve) = nan;
        bar_data_nve_error_lower = bar_error_lower_bc; bar_data_nve_error_lower(ind_pve) = nan;
    end
    
    figure
    if isempty(fontsize)     fontsize = 15;               end
    if isempty(lgd_fontsize) lgd_fontsize = fontsize*0.8; end
    if isempty(linewidth)    linewidth = 1;               end
    if isempty(markersize)   markersize = 20;             end
    if isempty(edgecolor)    edgecolor = 'none';          end
    if isempty(show_lgd)     show_lgd = 1;                end
    %/ barh()
    if isempty(bar_range)       
        bar_data_max = max(bar_data_bc, [], 'all');
        bar_data_min = min(bar_data_bc, [], 'all');
        gap          = (bar_data_max - bar_data_min)*0.1;
        bar_range    = [bar_data_min-gap, bar_data_max+gap];    
    end

    if isempty(bar_rangeticks)  bar_rangeticks = 'auto';        end
    if isempty(orientation)     orientation    = 'vertical';    end
    if isequal(orientation, 'horizontal')
        if isempty(panel_pos)  panel_pos = [400 100 650 1000];  end
        set(gcf, 'position',panel_pos) %/ x, y, w, h
        b_pve = barh(bar_data_pve, bar_mode, 'FaceColor','flat');
        hold on;
        b_nve = barh(bar_data_nve, bar_mode, 'FaceColor','flat');
        
        ytickangle(0)
        xlim(bar_range)
        xticks(bar_rangeticks)
        xlabel(y_labels)
        if reverse_data_axis
            %/ do nothing.
        else
            set(gca, 'YDir','reverse');
        end
    else
        if isempty(panel_pos)  panel_pos = [400 100 1200 650];   end
        set(gcf,'position',panel_pos)  %/ x, y, w, h
        
        if isequal(bar_mode, 'stacked')
            b_pve = bar(bar_data_pve, bar_mode, 'FaceColor','flat', 'edgecolor', edgecolor);
            hold on;
            b_nve = bar(bar_data_nve, bar_mode, 'FaceColor','flat', 'edgecolor', edgecolor);
            
            if ~isempty(bar_error_upper_bc)
                error('Error bar does not work with ''stacked'' bars!');
            end
            % xticks(1:bar_labels_intvl:length(group_labels_bc))
            % xticklabels(group_labels_bc(1:bar_labels_intvl:end))
            xticks(1:length(group_labels_bc))
            xticklabels(group_labels_bc)
            xtickangle(0)
            xlim([0, length(group_labels_bc)+1])  %/ add margins at two ends of the bars
            hold on;

        elseif isequal(bar_mode, 'grouped')
            %/ Loop to plot a bar to each value 
            %/ Why? Because only with this we can have the bar handle for a
            %/ correct legend!! Only for grouped bar, no need for stacked bar.

            % if length(stack_labels_bc) ~= numel(bar_data_pve)
            %     error('The number of stack_labels does not equal to number of values!');
            % end
            b_pve = cell(length(group_labels_bc),1); 
            b_nve = cell(length(group_labels_bc),1);
            for i = 1:length(group_labels_bc)
                b_pve{i} = bar(i, bar_data_pve(i,:), bar_mode, 'FaceColor','flat', 'edgecolor', edgecolor, 'BarWidth', .9);
                hold on;
                b_nve{i} = bar(i, bar_data_nve(i,:), bar_mode, 'FaceColor','flat', 'edgecolor', edgecolor, 'BarWidth', .9);
            end

            %/ Add error bar
            if ~isempty(bar_error_upper_bc)
                hold on;
                color_err      = [30 30 30]/255;
                linewidth_err  = linewidth*0.6;
                markersize_err = markersize;

                if ~isvector(bar_data_pve)  %/ If a 2D matrix, then there are nested cells in b_pve!
                    for i = 1:numel(b_pve)   
                        for j = 1:length(b_pve{i})
                            if ~isnan(bar_data_pve_error_upper(i,j))
                                % ctr_pve(i,:) = bsxfun(@plus, b_pve{i}(1,j).XData, [b_pve{i}(1,j).XOffset]');
                                ctr_pve(i,:) = b_pve{i}(1,j).XEndPoints;
                                ydt_pve(i,:) = b_pve{i}(1,j).YData;
        
                                % errorbar(x,y,yneg,ypos), where yneg and ypos define the lower and upper lengths of the vertical error bars
                                errorbar(ctr_pve(i,:), ydt_pve(i,:), bar_data_pve_error_lower(i,j), bar_data_pve_error_upper(i,j), 'color', color_err, 'marker', '.', 'MarkerSize', markersize_err, 'linewidth', linewidth_err, 'linestyle', 'none')
                            end
        
                            if ~isnan(bar_data_nve_error_upper(i,j))
                                % ctr_nve(i,:) = bsxfun(@plus, b_nve{i}(1,j).XData, [b_nve{i}(1,j).XOffset]');
                                ctr_nve(i,:) = b_nve{i}(1,j).XEndPoints;
                                ydt_nve(i,:) = b_nve{i}(1,j).YData;
                                errorbar(ctr_nve(i,:), ydt_nve(i,:), bar_data_nve_error_lower(i,j), bar_data_nve_error_upper(i,j), 'color', color_err, 'marker', '.', 'MarkerSize', markersize_err, 'linewidth', linewidth_err, 'linestyle', 'none')
                            end
                        end
                    end
                else
                    for i = 1:numel(b_pve)   
                        if ~isnan(bar_data_pve_error_upper(i))
                            % ctr_pve(i,:) = bsxfun(@plus, b_pve{i}.XData, [b_pve{i}.XOffset]');
                            ctr_pve(i,:) = b_pve{i}.XEndPoints;
                            ydt_pve(i,:) = b_pve{i}.YData;
    
                            % errorbar(x,y,yneg,ypos), where yneg and ypos define the lower and upper lengths of the vertical error bars
                            errorbar(ctr_pve(i,:), ydt_pve(i,:), bar_data_pve_error_lower(i), bar_data_pve_error_upper(i), 'color', color_err, 'marker', '.', 'MarkerSize', markersize_err, 'linewidth', linewidth_err, 'linestyle', 'none')
                        end
    
                        if ~isnan(bar_data_nve_error_upper(i))
                            % ctr_nve(i,:) = bsxfun(@plus, b_nve{i}.XData, [b_nve{i}.XOffset]');
                            ctr_nve(i,:) = b_nve{i}.XEndPoints;
                            ydt_nve(i,:) = b_nve{i}.YData;
                            errorbar(ctr_nve(i,:), ydt_nve(i,:), bar_data_nve_error_lower(i), bar_data_nve_error_upper(i), 'color', color_err, 'marker', '.', 'MarkerSize', markersize_err, 'linewidth', linewidth_err, 'linestyle', 'none')
                        end
                    end
                end
                hold off
            end
            % xticks(1:bar_labels_intvl:length(group_labels_bc))
            % xticklabels(group_labels_bc(1:bar_labels_intvl:end))
            xticks(1:length(group_labels_bc))
            xticklabels(group_labels_bc)
            xtickangle(0)
            if isempty(barplot_margin)
                barplot_margin = 1;
            end
            xlim([1-barplot_margin, length(group_labels_bc)+barplot_margin])  %/ add margins at two ends of the bars
            hold on;


            %/ ------- old code ---------- %/
            % b_pve = cell(length(stack_labels_bc),1); 
            % b_nve = cell(length(stack_labels_bc),1);
            % for i = 1:length(stack_labels_bc)
            %     b_pve{i} = bar(i, bar_data_pve(i), bar_mode, 'FaceColor','flat', 'edgecolor', edgecolor, 'BarWidth', .9);
            %     hold on;
            %     b_nve{i} = bar(i, bar_data_nve(i), bar_mode, 'FaceColor','flat', 'edgecolor', edgecolor, 'BarWidth', .9);
            % end
            % 
            % %/ Add error bar
            % if ~isempty(bar_error_upper_bc)
            %     hold on;
            %     color_err      = [30 30 30]/255;
            %     linewidth_err  = linewidth*0.6;
            %     markersize_err = markersize;
            %     for i = 1:numel(b_pve)   
            %         if ~isnan(bar_data_pve_error_upper(i))
            %             ctr_pve(i,:) = bsxfun(@plus, b_pve{i}.XData, [b_pve{i}.XOffset]');
            %             ydt_pve(i,:) = b_pve{i}.YData;
            % 
            %             % errorbar(x,y,yneg,ypos), where yneg and ypos define the lower and upper lengths of the vertical error bars
            %             errorbar(ctr_pve(i,:), ydt_pve(i,:), bar_data_pve_error_lower(i).', bar_data_pve_error_upper(i).', 'color', color_err, 'marker', '.', 'MarkerSize', markersize_err, 'linewidth', linewidth_err, 'linestyle', 'none')
            %         end
            % 
            %         if ~isnan(bar_data_nve_error_upper(i))
            %             ctr_nve(i,:) = bsxfun(@plus, b_nve{i}.XData, [b_nve{i}.XOffset]');
            %             ydt_nve(i,:) = b_nve{i}.YData;
            %             errorbar(ctr_nve(i,:), ydt_nve(i,:), bar_data_nve_error_lower(i).', bar_data_nve_error_upper(i).', 'color', color_err, 'marker', '.', 'MarkerSize', markersize_err, 'linewidth', linewidth_err, 'linestyle', 'none')
            %         end
            %     end
            %     hold off
            % end
            % xticks(1:bar_labels_intvl:length(group_labels_bc))
            % xticklabels(group_labels_bc(1:bar_labels_intvl:end))
            % xtickangle(0)
            % xlim([0, length(group_labels_bc)+1])  %/ add margins at two ends of the bars
            % hold on;
        end
        
        ylabel(y_labels)
        ylim(bar_range)
        yticks(bar_rangeticks)
        if reverse_data_axis
            set(gca, 'YDir','reverse'); 
        end
    end
    
    %/ Modify the bar colors 
    if isequal(bar_mode, 'stacked')
        for i = 1:numel(b_pve)
            b_pve(i).CData     = bar_colmap_bc(i,:);
            b_pve(i).FaceColor = bar_colmap_bc(i,:);
            b_nve(i).CData     = bar_colmap_bc(i,:);
            b_nve(i).FaceColor = bar_colmap_bc(i,:);
        end
        
    elseif isequal(bar_mode, 'grouped')
        if ~isvector(bar_data_pve)  %/ If a 2D matrix, then there are nested cells in b_pve!
            for i = 1:numel(b_pve)
                for j = 1:length(b_pve{i})
                    b_pve{i}(1,j).CData     = bar_colmap_bc(j,:);
                    b_pve{i}(1,j).FaceColor = bar_colmap_bc(j,:);
                    b_nve{i}(1,j).CData     = bar_colmap_bc(j,:);
                    b_nve{i}(1,j).FaceColor = bar_colmap_bc(j,:);
                end
            end
        else
            for i = 1:numel(b_pve)
                b_pve{i}.CData     = bar_colmap_bc(i,:);
                b_pve{i}.FaceColor = bar_colmap_bc(i,:);
                b_nve{i}.CData     = bar_colmap_bc(i,:);
                b_nve{i}.FaceColor = bar_colmap_bc(i,:);
            end
        end
    end

    %/===============================
    %/ Put text on the stackbars 
    if isempty(text_fontsize)
        text_fontsize = fontsize*0.9;
    end

    if isempty(text_yshift)
        text_yshift = 0;    %/ Shift the text in y-dir by text_yshift
    end

    if isequal(bar_mode, 'stacked')
        set(gca, 'FontSize', fontsize, 'linewidth', linewidth)

        ind = findismember_loop(stack_labels_bc, stack_labels); %/ keep the original order of stack_labels
        if show_lgd
            legend(b_pve(ind), strrep(stack_labels, '_', ' '), 'Location', lgd_position, 'Orientation', lgd_orient, 'edgecolor', 'none', 'color', 'none', 'AutoUpdate','off', 'fontsize', lgd_fontsize)
        end

        if ~isempty(text_on_stackedbar)
            yt = get(gca, 'YTick');
            barbase = cumsum([zeros(nbars, 1) bar_data_pve(:,1:end-1)],2);
            joblblpos = bar_data_pve/2 + barbase;
            for k = 1:size(bar_data_pve,1)
                text(joblblpos(k,:), yt(k)*ones(1,nstacks)+text_yshift, text_on_stackedbar(k,:), 'HorizontalAlignment','center', 'fontsize', text_fontsize)
            end
        end 

        %/ Put text next to stackbars
        if ~isempty(text_next_to_stackedbar)
            yt = get(gca, 'YTick');
            joblblpos = sum(bar_data_pve,2)*1.05;
            
            for k = 1:size(bar_data_pve,1)
                text(joblblpos(k,:), yt(k)+text_yshift, text_next_to_stackedbar(k), 'HorizontalAlignment','center', 'fontsize', text_fontsize)
            end
        end
        
    elseif isequal(bar_mode, 'grouped')
        %/ Put values on top / bottom of each bar (by default)
        if isempty(show_value_digit)
            show_value_digit = 1;
        end

        hb = {}; hT = []; 
        if ~isvector(bar_data_pve)  %/ If a 2D matrix, then there are nested cells in b_pve!
            for i=1:numel(b_pve)
                for j=1:length(b_pve{i})  
                    if ~isempty(bar_error_upper_bc)
                        if b_pve{i}(1,j).YData ~= 0
                            hT = [hT text(b_pve{i}(1,j).XEndPoints, b_pve{i}(1,j).YData+bar_data_pve_error_upper(i,j)+text_yshift, string(round(b_pve{i}(1,j).YData, show_value_digit).'), ...
                                                  'VerticalAlignment','bottom','horizontalalign','center', 'fontsize', text_fontsize)];
                            hb = [hb, b_pve];
                            
                        elseif b_nve{i}(1,j).YData ~= 0
                            hT = [hT text(b_nve{i}(1,j).XEndPoints, b_nve{i}(1,j).YData+bar_data_nve_error_lower(i,j)-text_yshift, string(round(b_nve{i}(1,j).YData, show_value_digit).'), ...
                                                  'VerticalAlignment','top','horizontalalign','center', 'fontsize', text_fontsize)];
                            hb = [hb, b_nve];
                        end
                    else
                        if b_pve{i}(1,j).YData ~= 0
                            %/ Checking which can indicate the x-pos of bars
                            % b_pve{i}(1,j).XOffset
                            % b_pve{i}(1,j).XEndPoints

                            hT = [hT text(b_pve{i}(1,j).XEndPoints, b_pve{i}(1,j).YData+text_yshift, string(round(b_pve{i}(1,j).YData, show_value_digit).'), ...
                                                  'VerticalAlignment','bottom','horizontalalign','center', 'fontsize', text_fontsize)];
                            hb = [hb, b_pve];
                        elseif b_nve{i}(1,j).YData ~= 0

                            hT = [hT text(b_nve{i}(1,j).XEndPoints, b_nve{i}(1,j).YData-text_yshift, string(round(b_nve{i}(1,j).YData, show_value_digit).'), ...
                                                  'VerticalAlignment','top','horizontalalign','center', 'fontsize', text_fontsize)];
                            hb = [hb, b_nve];
                        end
                    end
                end
            end
        else
            for i=1:numel(b_pve)  %/ iterate over number of bar objects
                if ~isempty(bar_error_upper_bc)
                    if b_pve{i}.YData ~= 0
                        hT = [hT text(b_pve{i}.XEndPoints, b_pve{i}.YData+bar_data_pve_error_upper(i)+text_yshift, string(round(b_pve{i}.YData, show_value_digit).'), ...
                                              'VerticalAlignment','bottom','horizontalalign','center', 'fontsize', text_fontsize)];
                        hb = [hb, b_pve];
                    elseif b_nve{i}.YData ~= 0
                        hT = [hT text(b_nve{i}.XEndPoints, b_nve{i}.YData+bar_data_nve_error_lower(i)-text_yshift, string(round(b_nve{i}.YData, show_value_digit).'), ...
                                              'VerticalAlignment','top','horizontalalign','center', 'fontsize', text_fontsize)];
                        hb = [hb, b_nve];
                    end
                else
                    if b_pve{i}.YData ~= 0
                        hT = [hT text(b_pve{i}.XEndPoints, b_pve{i}.YData+text_yshift, string(round(b_pve{i}.YData, show_value_digit).'), ...
                                              'VerticalAlignment','bottom','horizontalalign','center', 'fontsize', text_fontsize)];
                        hb = [hb, b_pve];
                    elseif b_nve{i}.YData ~= 0
                        hT = [hT text(b_nve{i}.XEndPoints, b_nve{i}.YData-text_yshift, string(round(b_nve{i}.YData, show_value_digit).'), ...
                                              'VerticalAlignment','top','horizontalalign','center', 'fontsize', text_fontsize)];
                        hb = [hb, b_nve];
                    end
                end
            end
            
        end
        set(gca, 'FontSize', fontsize, 'linewidth', linewidth)
        ind = findismember_loop(stack_labels_bc, stack_labels); %/ keep the original order of stack_labels
        
        if show_lgd
            if ~isvector(bar_data_pve)  %/ If a 2D matrix, then there are nested cells in b_pve!
                legend([b_pve{1}(ind)], strrep(stack_labels, '_', ' '), 'Location', lgd_position, 'Orientation', lgd_orient, 'edgecolor', 'none', 'color', 'none', 'AutoUpdate','off', 'fontsize', lgd_fontsize);
            else
                legend([b_pve{ind}], strrep(stack_labels, '_', ' '), 'Location', lgd_position, 'Orientation', lgd_orient, 'edgecolor', 'none', 'color', 'none', 'AutoUpdate','off', 'fontsize', lgd_fontsize);
            end
        end
    end
    if box_on
        box on  
    else
        box off  %/ by default
    end
    grid off  %/ will only use xline / yline. 
    
    if ~isempty(group_labels_bc)
        if ischar(group_labels_bc) group_labels_bc = {group_labels_bc}; end
        
        %/ Label the group only when there are > 1 groups
        if length(group_labels_bc) > 1
            if iscell(group_labels_bc)
                xticks(1:length(group_labels_bc))
                xticklabels(strrep(group_labels_bc, '_', ' '))
                xtickangle(0)
            elseif iscategorical(group_labels_bc)
                xticks(group_labels_bc)
                xtickangle(0)
                for i = 1:length(group_labels_bc)
                    xh = xline(group_labels_bc(i), 'color', 'k', 'linestyle', '-', 'linewidth', 1); %/ plot it before bar.

                    %/ A workaround to put xline at the back of the plot
                    %/ See https://www.mathworks.com/matlabcentral/answers/667763-constantline-always-on-the-top
                    warnState = warning('off','MATLAB:structOnObject');
                    cleanupObj = onCleanup(@()warning(warnState)); 
                    Sxh = struct(xh);        % Get undocumented properties (you'll get a warning)
                    clear('cleanupObj')      % Trigger warning reset
                    Sxh.Edge.Layer = 'back'; % Set ConstantLine uistack
                end
            end
        end
        hold on;
    end
    
    %/ Annotation is indep. of axes -> not distorting the aspect ratio.
    if isempty(title_fontsize)
        title_fontsize = fontsize*0.75;       
    end
    if isempty(title_pos)
        title_pos = [0.1, 1.02, 1, 0];     % x y w h,  set w = 0.8 to allow new lines
    end
    annotation( 'textbox',  'String', titlename, 'Color', 'k', ...
                'FontSize', title_fontsize, 'FitBoxToText','on', 'Units', 'normalized', 'EdgeColor', 'none',...
                'Position', title_pos)
    
%     th = title(strrep(titlename, '_', ' '), 'fontsize', fontsize*1.5);
%     titlePos = get(th, 'position');
%     titlePos(2) = titlePos(2);   %/ change y position of title
%     set(th, 'position', titlePos)
            
    %/ Remove panel bg color (Put these at the end, otherwise it may not work!)
    set(gcf, 'color', 'none');
    set(gca, 'color', 'none');     
    
    if savefig
        final_figname = char(strcat(plotting_folder, FigName_underscore,'.pdf'));
        export_fig(final_figname,'-pdf','-painters', '-nocrop', '-transparent'); % '-c[inf, inf, inf, inf]'); %, '-nocrop'); %'-transparent');
        fprintf('Figure saved to: %s\n', final_figname);
    end

    set(gcf, 'color', 'w');
    set(gca, 'color', 'w');  
end
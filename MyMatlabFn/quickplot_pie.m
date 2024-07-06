%%
function quickplot_pie(varargin)

pnames = {'pie_data', 'pie_label', 'pie_legend', 'pie_colmap', 'thres_for_noshow', 'sort_mode',...
          'label_pos', 'label_col', 'label_fontsize', 'fontsize', 'linewidth', 'pie_titlename', 'pie_figname', 'plotting_folder', 'show_lgd', 'plot_lgd_only',...
          'Orientation', 'NumColumns'};

dflts  = cell(size(pnames));

[          pie_data, pie_label, pie_legend, pie_colmap, thres_for_noshow, sort_mode,...
           label_pos, label_col, label_fontsize, fontsize, linewidth, pie_titlename, pie_figname, plotting_folder, show_lgd, plot_lgd_only,...
           Orientation, NumColumns] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});

%/ 'label_pos':  0 at center. 1 at the circumference. >1: outside the pie.

%%
if isempty(Orientation)   Orientation = 'vertical';   end
if isempty(NumColumns)    NumColumns  = 1;            end
if isempty(show_lgd)      show_lgd    = 1;            end

if plot_lgd_only
    figure
%     set(gcf,'Position', [0 300 750 750]);
    set(gcf,'Position', get(0, 'Screensize'));
    set(gcf,'color','w'); 

    h = pie(ones(size(pie_data)));
    colormap(gca,pie_colmap); %/ set to gca to allow different colormaps in different panels
    
    %/ move the entire pie chart upward. (more complex than I thought..)
    scale = 1;
    xpos = 0;
    ypos = 1;  
    for k = 1:numel(h)
        if strcmp(get(h(k),'Type'),'patch') % Patch graphics
            XData = get(h(k),'XData');  % Extract current data
            YData = get(h(k),'YData');
            set(h(k),'XData',XData*scale + xpos); % Insert modified data
            set(h(k),'YData',YData*scale + ypos);
            
        else    % Text labels
            Pos = get(h(k),'Position'); % Extract
            set(h(k),'Position',Pos*scale + [xpos ypos 0]); % Insert
            set(h(k),'FontSize',8);
        end
    end
    set(h, 'linewidth', linewidth);
    
    h_lgd = legend(pie_legend, 'Location', 'southoutside', 'Orientation', Orientation,...
                                'NumColumns', NumColumns, 'fontsize', fontsize, 'linewidth', linewidth*0.75);
    
    %/ NOTE: Do NOT write [h_lgd, h_lgd_Obj] = legend(...), for it will
    %create a new set of legend objects whose colors do not match the
    %colormap. Always begin from h_lgd.
    hL=findobj(h_lgd.Parent.CurrentAxes.Children,'type','patch');  % get the lines, not text
    set(hL,'linewidth', linewidth*0.75)                             % set their width property
    
    h_lgd.Position(1) = h_lgd.Position(1) + 0.15; %/ shift it to the right
    legend boxoff
%     set(gca, 'linewidth', linewidth*0.1);
%     hold on;
    
    if ~isempty(pie_figname)
        export_fig(char(strcat(plotting_folder, pie_figname, '_lgd.pdf')),'-pdf','-painters', '-c[inf, inf, inf, inf]', '-nocrop', '-transparent');
        fprintf('!!! Plot is saved: %s !!! \n', string(pie_figname));
    end

    %/ if set(hp, 'Visible', 'off'), the legend wont remain unchanged. FAILED.
    return
end


figure
set(gcf, 'Position', [0 300 750 500]);
% set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','w'); 

%/ [IMPORTANT] set non-positve values to a tiny dummy value
pie_data(pie_data <= 0) = 1e-10;    %/ set to a tiny dummy value to keep the items (otherwise the item will be skipped and outputs a wrongly-shifted legend!)


%/ Use black (white) label (when it is inside the pie) for light (dark) pie.
if isempty(label_col)
    label_col   = repmat([0 0 0], length(pie_data), 1);
    ind_usewhite = find(mean(pie_colmap,2)<0.35 & label_pos < 1);
    label_col(ind_usewhite,:) = repmat([1 1 1], length(ind_usewhite), 1);  
end

%/ Set fontsize of thelabels (e.g., the bigger the pie, the larger the label)
if isempty(label_fontsize)
    label_fontsize = fontsize*(1 + sqrt(pie_data)/3);
end

%/ Set the labels within large pies, otherwise outside the pie.
if isempty(label_pos)
    label_pos = repmat(1.1, length(pie_data), 1);    %/ outside (default)
    label_pos(pie_data >= 10 & pie_data < 20) = 0.6; %/ inside
    label_pos(pie_data >= 20) = 0.5;                 %/ more inside

elseif length(label_pos) == 1
    label_pos = repmat(label_pos, length(pie_data), 1);
end

%/ Sort the pies
if ~isempty(sort_mode)  %/ 'descend', 'ascend'   [] -> no sorting.
    [B, I] = sort(pie_data, sort_mode);
    pie_data        = B;
    pie_label       = pie_label(I);
    pie_colmap      = pie_colmap(I,:);
    label_col       = label_col(I,:);
    pie_legend      = pie_legend(I);
    label_fontsize  = label_fontsize(I);
    label_pos       = label_pos(I);
end

%/ do not show labels for values < thres_for_noshow.
if ~isempty(thres_for_noshow)
    ind = find(pie_data < thres_for_noshow); 
    if ~isempty(ind)
        pie_label(ind) = '';
    end
end

%==== plot pie chart for each hs ====%
h = pie(pie_data, pie_label);

for i = 1:numel(h)/2
%     disp(hp(i*2).Position)
    
    pos = h(i*2).Position;
    [theta, rho] = cart2pol(pos(1), pos(2));   %/ find the radius (rho) first
%     fprintf('rho = %.4f \n', rho);
    
    pos_new = [];
    [pos_new(1), pos_new(2)] = pol2cart(theta, label_pos(i)*rho); %/ change the rho, then convert back to x,y
    h(i*2).Position = [pos_new, 0]; %/ update the new position (x, y, z)
    
    h(i*2).HorizontalAlignment = 'center'; %/ important to get the label centered!
    h(i*2).VerticalAlignment   = 'middle'; %/ important to get the label centered!
    h(i*2).FontSize = label_fontsize(i);
    h(i*2).Color    = label_col(i,:);

%     hp(iHandle).BackgroundColor = [1 1 1 ];
end
% T = hp(strcmpi(get(hp,'Type'),'text'));
% P = cell2mat(get(T,'Position'));
% set(T,{'Position'},num2cell(P*0.25,2))

% delete(findobj(hp,'Type','text')) %option 2
hold on;
% set(findobj(hp,'type','text'),'fontsize',fontsize*4);
set(h, 'linewidth', linewidth);

title(pie_titlename, 'fontsize', fontsize*2);
colormap(gca,pie_colmap); %/ set to gca to allow different colormaps in different panels

if show_lgd
    h_lgd = legend(pie_legend, 'Location', 'eastoutside', 'fontsize', fontsize, 'linewidth', linewidth*0.75);
    legend boxoff
    hold on;
    %/ NOTE: Do NOT write [h_lgd, h_lgd_Obj] = legend(...), for it will
    %create a new set of legend objects whose colors do not match the
    %colormap. Always begin from h_lgd.
    hL=findobj(h_lgd.Parent.CurrentAxes.Children,'type','patch');  % get the lines, not text
    set(hL,'linewidth', linewidth*0.75)                             % set their width property
   
end

if ~isempty(pie_figname)
    %/ CAVEAT: auto-cropping will cut off the labels outside the pie!
    export_fig(char(strcat(plotting_folder, pie_figname, '.pdf')),'-pdf','-painters', '-transparent', '-nocrop'); % '-c[inf, inf, inf, inf]'
    fprintf('!!! Plot is saved: %s !!! \n', string(pie_figname));
end

%/ Save legend (by readjusting the linewidth - a workaround) - no need anymore
% set(hp, 'linewidth', 2);
% if ~isempty(pie_figname)
%     export_fig(char(strcat(plotting_folder, 'pie_ldg.pdf')),'-pdf','-painters', '-c[inf, inf, inf, inf]', '-nocrop', '-transparent');
%     fprintf('!!! Plot is saved: %s !!! \n', string(pie_figname));
% end

% close all
% figure
% x = [0:100]/100;
% y = exp(-5*x);
% plot(x, y)

% if i == 1
%     hp(i*2).Position = 0.05*hp(i*2).Position;
% elseif i == 3
%     hp(i*2).Position = 0.4*hp(i*2).Position;
%     hp(i*2).Position(1) = 0.6*hp(i*2).Position(1) + 0.1;  %/ nudge x pos
% else
%     hp(i*2).Position = 0.1*hp(i*2).Position;
% end

end
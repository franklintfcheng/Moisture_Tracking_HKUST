%%
function plot_pointmap(varargin)

% create a set of valid parameters and their default value
pnames = {'bndry_data', 'point_data', 'marker', 'markersize', 'markerfacecolor', 'linest', 'edgecolor', 'linewi',...
          'map_lon_lower', 'map_lon_upper', 'map_lat_lower', 'map_lat_upper', 'fontsize', 'titlename', 'savepath'};

dflts  = {[], [], 'none', 1, 'none', 'none', 'none', 1, -179, 179, -90, 90, 20, {}, {}};

%/ parse function arguments
[ bndry_data, point_data, marker, markersize, markerfacecolor, linest, edgecolor, linewi,...
  map_lon_lower, map_lon_upper, map_lat_lower, map_lat_upper, fontsize, titlename, savepath] ...
             = internal.stats.parseArgs(pnames, dflts, varargin{:});

%===============
figure
set(gcf,'Position', get(0, 'Screensize'));
set(gcf,'color','w');

% 'Miller Cylindrical'
m_proj('robin','longitudes',[map_lon_lower map_lon_upper], ...
               'latitudes', [map_lat_lower map_lat_upper]);

%/ Plot points
if ~isempty(point_data)
    m_line(point_data(:,1), point_data(:,2), 'marker', marker, 'markersize', markersize, 'markerfacecolor', markerfacecolor,...
          'linest', 'none', 'color', 'none', 'linewi', 0.001);
    hold on;
end

%/ Plot boundaries
if ~isempty(bndry_data)
    if iscell(bndry_data) %/ if point_data is a cell array, then loop over it to plot.
        for i = 1:length(bndry_data)
            m_line(bndry_data{i}(:,1),bndry_data{i}(:,2), 'marker', 'none', 'markersize', 0.001, 'markerfacecolor', 'none',...
                      'linest', linest, 'color', edgecolor, 'linewi', linewi);
            hold on;
        end
    else
        m_line(bndry_data(:,1),bndry_data(:,2), 'marker', 'none', 'markersize', 0.001, 'markerfacecolor', 'none',...
                  'linest', linest, 'color', edgecolor, 'linewi', linewi);
        hold on;
    end
end

m_coast('linewidth',1.5,'color',[.15 .15 .15]);
% m_gshhs_i('speckle','color',[0.2 0.2 0.2]);    % with speckle added
% m_gshhs('fr1','Color','b','LineWidth',1.2);
m_grid('xtick', -180:30:360,'ytick', map_lat_lower:20:map_lat_upper,...
       'linewi',3,'linest',':','gridcolor', [.2 .2 .2], 'tickdir','in', 'fontsize',fontsize-6);
        
hold on;
title(titlename, 'fontsize', fontsize+5);

if ~isempty(savepath)
    FigName_underscore = strrep(savepath, ' ', '_');
%     export_fig(char(strcat(savepath,FigName_underscore,'.pdf')),'-pdf','-painters',...
%                         '-transparent', '-nocrop');
    
    export_fig(char(strcat(FigName_underscore,'.png')),'-png','-r300');      
                   
end

end
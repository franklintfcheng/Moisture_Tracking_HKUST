%%
function quickplot_basin_topo(varargin)
    
    pnames = {'bndry_data', 'color', 'draw_ZZ_county', 'map_lon_lower', 'map_lon_upper', 'map_lat_lower', 'map_lat_upper',...
              'topo_intvl', 'topo_land_only', 'plateau_hgt', 'plateau_col', 'mask_mountains', 'pcolor_mode',  'shadedrelief_mode',...
              'map_proj', 'fontsize', 'linewi', 'coast_wi', 'coast_col', 'grid_mode',...
              'plotting_folder', 'FigName_underscore', 'savefig', 'fig_fmt'};
          
    dflts = cell(length(pnames), 1);
    
    [          bndry_data,  color, draw_ZZ_county, map_lon_lower, map_lon_upper, map_lat_lower, map_lat_upper,...
               topo_intvl,   topo_land_only,   plateau_hgt,   plateau_col,   mask_mountains, pcolor_mode,   shadedrelief_mode,...
               map_proj,   fontsize,   linewi,   coast_wi,   coast_col,   grid_mode,...
               plotting_folder,   FigName_underscore,   savefig,   fig_fmt] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    if isempty(shadedrelief_mode) shadedrelief_mode = 1;             end %/ turn it on by default
    if isempty(topo_intvl)        topo_intvl  = 250;                 end
    if isempty(plateau_col)       plateau_col = [255 51 204]./255;   end
    if isempty(linewi)            linewi    = 2.5;                   end     
    if isempty(map_proj)          map_proj  = 'Miller Cylindrical';  end
    if isempty(fontsize)          fontsize  = 20;                    end
    if isempty(coast_wi)          coast_wi  = 2.5;                   end
    if isempty(coast_col)         coast_col = [.4 .4 .4];            end
    if isempty(fig_fmt)           fig_fmt   = 'pdf';                 end
    if isempty(color)             color = 'k';                       end
    
    
    if savefig   %/ do not plot title when saving fig.
        titlename = [];
    else
        titlename = strrep(FigName_underscore, '_', ' ');
    end
    shift_title_y = 0.03;
    draw_bndry_onebyone = 0;
    
    %/ Parameters
    create_fig       = 1;
    glb_data_mode    = 0;
    if ~isempty(plateau_hgt)
        glb_plateau_mode = 1;
    else
        glb_plateau_mode = 0;
    end
    
    backcolor        = 'none';
    marker = 'none'; markersize = 0.0001; markerfacecolor = 'k';
    
    %----- contf -----%
    %/ Extract high-res topography
    m_proj(map_proj,'longitudes',[map_lon_lower map_lon_upper], 'latitudes', [map_lat_lower map_lat_upper]);
    tic;
    [contf_data, lon_2D_topo, lat_2D_topo]=m_tbase([map_lon_lower, map_lon_upper, map_lat_lower, map_lat_upper]); %REGION =[west east south north];
    toc;
    if topo_land_only
        contf_data(contf_data < 0) = -100;
        contf_levels = [-topo_intvl:topo_intvl:topo_intvl*26];
        colmap = my_colormap(length(contf_levels)-1,'topo_land_27lev');
        
        %/ set colors grey if below plateau_hgt
%         ind = find(contf_levels >= 0 & contf_levels < plateau_hgt);
%         colmap(ind,:) = repmat([.7 .7 .7], length(ind), 1);
        
    else
        contf_levels    = [-4000:topo_intvl:7000];
        colmap          = [m_colmap('blues',length(find(contf_levels < 0)));m_colmap('gland',length(find(contf_levels > 0)))];
    end
    
    contf_lon       = lon_2D_topo(1,:);
    contf_lat       = lat_2D_topo(:,1);
    contf_unit      = '';
    cbar_mode       = 1;
    cbar_location   = 'eastoutside';
    cbar_interval   = 2;
    cbar_YTick      = [-4000:1000:7000];
    cbar_YTickLabel = cbar_YTick;
    
    savepath = [];  
    plot_contfmap('contf_data', contf_data, 'contf_lon', contf_lon, 'contf_lat', contf_lat, 'contf_levels', contf_levels, 'shadedrelief_mode', shadedrelief_mode,...
                  'contf_unit', contf_unit, 'colmap', colmap, 'cbar_interval', cbar_interval, 'pcolor_mode', pcolor_mode,...
                  'bndry_data', bndry_data, 'draw_bndry_onebyone', draw_bndry_onebyone, 'color', color, 'marker', marker, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'linewi', linewi, 'color', color,...
                  'titlename', titlename, 'shift_title_y', shift_title_y, 'savepath', savepath, ...
                  'glb_data_mode', glb_data_mode, 'glb_plateau_mode', glb_plateau_mode, 'plateau_hgt', plateau_hgt, 'plateau_col', plateau_col, 'mask_mountains', mask_mountains,...
                  'map_proj', map_proj, 'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper, 'coast_col', coast_col, 'coast_wi', coast_wi, 'backcolor', backcolor,...
                  'fontsize', fontsize,  'create_fig', create_fig, 'grid_mode', grid_mode, 'cbar_mode', cbar_mode, 'cbar_location', cbar_location, 'cbar_YTick', cbar_YTick, 'cbar_YTickLabel', cbar_YTickLabel)
    fprintf('done \n')
    
    if draw_ZZ_county
        %/ load Chinese city borders
        %/ http://gaohr.win/site/blogs/2017/2017-04-18-GIS-basic-data-of-China.html
        filename =  strcat('/disk/r059/tfchengac/henan_floods_2021/prcssd_data_4plotting/Henan/Henan_county.shp');
        [S,A] = shaperead(filename,'UseGeoCoords',true);

        ind_Zhengzhou = find([A.CENTROID_Y] > 34.4 & [A.CENTROID_Y] < 34.9 & ...
                             [A.CENTROID_X] > 113 & [A.CENTROID_X] < 114.1);

        henan_county_bndry_data = [cat(2, S(ind_Zhengzhou).Lon)', cat(2, S(ind_Zhengzhou).Lat)'];

        m_line(henan_county_bndry_data(:,1), henan_county_bndry_data(:,2), 'marker', 'none', 'markersize', 0.001, 'markerfacecolor', 'none',...
                              'linest', '-', 'color', 'r', 'linewi', 2);
        drawnow; pause(0.05);
    end

    if savefig
        if isequal(fig_fmt, 'pdf')
            final_figname = char(strcat(plotting_folder, FigName_underscore,'.pdf'));
            export_fig(final_figname,'-pdf','-painters', '-c[inf, inf, inf, inf]'); %, '-nocrop'); %'-transparent');
            
        elseif isequal(fig_fmt, 'png')
            final_figname = char(strcat(plotting_folder, FigName_underscore,'.png'));
            export_fig(final_figname,'-r300','-png','-opengl'); %, '-nocrop');
%             export_fig(final_figname,'-png', '-r300','-painters', '-c[inf, inf, inf, inf]'); %, '-nocrop'); %'-transparent');
        else
            error('only ''pdf'' or ''png'' is available for input ''fig_fmt''!') 
        end
        fprintf('*** Fig is saved: %s ***\n', final_figname);
    end
end
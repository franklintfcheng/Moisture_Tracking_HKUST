%%
function plot_hgt_lonlat(varargin)

    pnames       = {'contf_data', 'contf_hori', 'contf_hgt', 'colmap', 'contf_levels', 'pcolor_mode', 'cbar_mode', 'cbar_location', 'cbar_interval', 'cbar_YTick', 'cbar_YTickLabel', 'draw_cbar_only',...
                    'cont_data',  'cont_hori',  'cont_hgt',  'cont_colmap', 'cont_levels', 'cont_linewi', 'skip_zero_cont',...
                    'vec_hori_data','vec_hgt_data', 'vec_hori', 'vec_hgt', 'vec_hori_factor', 'vec_hgt_factor', 'qscale', 'vector_color', 'vector_linewidth', 'vector_step_hori', 'vector_step_hgt', 'show_refvec_only', 'refvec_hori_mag', 'refvec_hgt_mag',...
                    'line_data',  'line_hori',  'add_topo_to_line_data', 'scatter_data', 'scatter_hgt', 'scatter_size', 'border_lines', 'border_color', 'border_linestyle', 'border_linewidth',...
                    'fontsize',  'linewidth', 'grid_mode', 'map_xlim', 'map_ylim', 'map_xticks', 'map_yticks', 'map_xticklabels', 'map_yticklabels', 'map_ylabel',...
                    'zm_or_mm', 'terrain_res', 'terrain_mode', 'terrain_bndry', 'topo_color',... 
                    'create_fig', 'gcf_position', 'plotting_folder', 'titlename', 'title_fontsize', 'title_pos', 'FigName_underscore', 'savefig'};
    dflts        = cell(length(pnames), 1);
    
    [                contf_data,   contf_hori,   contf_hgt,   colmap,   contf_levels,   pcolor_mode,   cbar_mode,   cbar_location,   cbar_interval,   cbar_YTick,   cbar_YTickLabel,   draw_cbar_only,...
                     cont_data,    cont_hori,    cont_hgt,    cont_colmap,   cont_levels,   cont_linewi,    skip_zero_cont,...
                     vec_hori_data, vec_hgt_data, vec_hori,     vec_hgt, vec_hori_factor, vec_hgt_factor,  qscale,   vector_color,   vector_linewidth,   vector_step_hori,   vector_step_hgt,   show_refvec_only,   refvec_hori_mag, refvec_hgt_mag,...
                     line_data,    line_hori,   add_topo_to_line_data,  scatter_data,   scatter_hgt,    scatter_size,  border_lines,   border_color,   border_linestyle,  border_linewidth,...
                     fontsize,    linewidth,   grid_mode,   map_xlim,    map_ylim,   map_xticks,  map_yticks,  map_xticklabels,   map_yticklabels, map_ylabel,...
                     zm_or_mm,   terrain_res,   terrain_mode,   terrain_bndry,   topo_color,...
                     create_fig, gcf_position, plotting_folder,   titlename,  title_fontsize, title_pos, FigName_underscore,   savefig] = ...
            internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

    %======================================================================
    %/ Author: Fandy Cheng
    %/ Date of creation: Jan 5, 2022
    %/
    %/ This function is designed for overlaying contourf, contour, vector and line plots
    %/ using MATLAB built-in functions.
    %/
    %/ E.g. Height vs lon/lat,  Homvoller plot, etc.
        
    %======================================================================
    
    %/ IMPORTANT: Assume contf_hgt, cont_hgt, etc.  must be in unit of km 
    %/          -> in order to have a good aspect ratio of grid and vector! 
    if max(contf_hgt) > 1000
       error('Check your contf_hgt! It must be in km!') 
    end
    %===============================================================
    
    if isempty(draw_cbar_only)       cbar_mode  = 1;                end
    if isempty(create_fig)           create_fig = 1;                end
    str_cbar    = []; str_qscale = [];  str_refvec = []; res_x = []; res_y = []; flag_skip_to_last = 0;
    
    %=============
    if create_fig
        figure
        set(gcf, 'color','w');
        if isempty(gcf_position)   
            gcf_position = [600 100 900 750];   % sets figure siz, landscape if always_time_as_y == 0 
        end  
        set(gcf,'position',gcf_position) % sets figure size  -> landscape
%         set(gcf, 'Position', get(0, 'Screensize'));
    end
    ax_ori = gca;

    %======== contf  ========%    
    if ~isempty(contf_data)
        contf_unit      = '';
        res_x = round(abs(diff(contf_hori(1:2))), 3);  %/ round up to 3 decimal points to avoid numerical errors
        res_y = round(abs(diff(contf_hgt(1:2))),  3);
            
        %/ Modify contf_data if contf_levels is uneven
        [contf_data, contf_levels, cbar_YTickLabel, cbar_YTick] = contfdata4uneven('contf_data', contf_data, 'contf_levels', contf_levels, 'cbar_YTickLabel', cbar_YTickLabel);
        
        if draw_cbar_only
            str_cbar    = '_cbar';
            set(ax_ori, 'visible', 'off');
            
            caxis([min(contf_levels) max(contf_levels)]);
            colormap(ax_ori, colmap)
            
            if isempty(cbar_YTick)       cbar_YTick      = contf_levels(2:cbar_interval:end-1);   end
            if isempty(cbar_YTickLabel)  cbar_YTickLabel = cbar_YTick;                            end
            
            if isnumeric(cbar_YTickLabel) %/ if is a numeric array, use sprintf to get the desired rounding.
                cbar_YTickLabel = arrayfun(@(x) sprintf('%.4g',x), cbar_YTickLabel, 'UniformOutput', false);
            end
            
            cb          = colorbar;
            cb_fontsize = fontsize*1.5;
            set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', cb_fontsize,...
                'Location', cbar_location)   
            set(get(cb,'Title'),'String',contf_unit, 'Fontsize', cb_fontsize)
            pos = get(cb, 'Position');
            
            if isequal(cbar_location, 'southoutside')
%                 pos(2) = pos(2)*5;
%                 set(cb, 'Position', pos);
                
            elseif isequal(cbar_location, 'eastoutside')
                pos(1) = pos(1)*0.8;
                set(cb, 'Position', pos);
            end
            flag_skip_to_last = 1;    %/ we then skip all the following lines until savefig.
        
        elseif show_refvec_only == 0         
            if show_refvec_only    %/ if show_refvec_only, we dont plot contourf.
                contf_data = zeros(length(contf_hori), length(contf_hgt));
            end
                
            if pcolor_mode
                h = pcolor(contf_hori-res_x/2, contf_hgt-res_y/2, contf_data');     %/ shift by half grid to correct the bug in pcolor.
                set(h, 'edgecolor','none');
                shading flat
            else
                warning('contourf may not show the discontinued grids.');
                contourf(contf_hori, contf_hgt, contf_data', [-99999999, contf_levels], 'edgecolor', 'none') %/ IMPORTANT: input contf_levels, otherwise the results are completely wrong!!
                shading interp
            end
            
            if cbar_mode
                if isempty(cbar_YTick)       cbar_YTick      = contf_levels(2:cbar_interval:end-1);   end
                if isempty(cbar_YTickLabel)  cbar_YTickLabel = cbar_YTick;                            end

                if isnumeric(cbar_YTickLabel) %/ if is a numeric array, use sprintf to get the desired rounding.
                    cbar_YTickLabel = arrayfun(@(x) sprintf('%.4g',x), cbar_YTickLabel, 'UniformOutput', false);
                end
                cb          = colorbar;
                cb_fontsize = fontsize*1.5;
                
                set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', cb_fontsize, 'Location', cbar_location)   
                set(get(cb,'Title'),'String',contf_unit, 'Fontsize', cb_fontsize)
                
                cb.Position = cb.Position + 1e-10;  %/ IMPORTANT: First stop auto resizing by setting cb pos. 
                
                %/ Fine-tune cb pos
                if isequal(cbar_location, 'eastoutside')
                    %/ Shift it eastward by a distance 3 times its width.
                    cb.Position(1) = cb.Position(1) + cb.Position(3)*3;  %/ x, y, width, height
                    
                elseif isequal(cbar_location, 'southoutside')
                    %/ Shift it southward by a distance 3 times its height.
                    cb.Position(2) = cb.Position(2) - cb.Position(4)*3;  %/ x, y, width, height
                end
            end
            
            caxis([min(contf_levels) max(contf_levels)]);
            colormap(ax_ori, colmap)
%             drawnow; pause(0.05);
            fprintf('*** Spatial avg of time-mean contf_data = %.2f ***\n', my_nanmean(contf_data, 'all'));
        end
    end
    hold on;
%     %/ Drop contf grids w/ sample sizes < x%. Only for Tb (ie., contf_dataname = '') or MCS_prcp!
%     if isempty(contf_dataname) || isequal(contf_dataname, 'MCS_prcp')
%         contf_data = drop_grids('field_data_3D', contf_data_3D, 'n_prct_to_drop', drop_MCS_grids_if_prct);
%         str_drop = sprintf('_drop%.1f%%', drop_MCS_grids_if_prct);
%     else
%         contf_data = my_nanmean(contf_data_3D, 3);  %/ 3D -> 2D
%         str_drop = '';
%     end
%     contf_data_3D = [];
%     hold on;
    
    if flag_skip_to_last == 0
        %======== vector ========%
        if ~isempty(vec_hori_data) && ~isempty(vec_hgt_data)
            res_x = round(abs(diff(vec_hori(1:2))), 3);  %/ round up to 3 decimal points to avoid numerical errors
            res_y = round(abs(diff(vec_hgt(1:2))),  3);
        
            [vec_hori_2D, vec_hgt_2D] = meshgrid(vec_hori, vec_hgt);  %/ lat x lon, as required by quiver().
            if size(vec_hori_data, 1) == length(vec_hori)    vec_hori_data = vec_hori_data';   end
            if size(vec_hgt_data, 1)  == length(vec_hori)    vec_hgt_data  = vec_hgt_data';    end

            %/ Multiply by given factors -> to output a good aspect ratio of vector.
            if isempty(vec_hori_factor)  vec_hori_factor = 1;  end
            if isempty(vec_hgt_factor)   vec_hgt_factor  = 1;  end
            
            if show_refvec_only
                str_refvec = '_refvec';
                
                %/ Horizontal reference vector
                X = vec_hori_2D;
                Y = vec_hgt_2D;
                %/ [IMPORTANT]: Make a U / V matrix with the reference vector and zeros elsewhere to get the correct quiver scale!
                U = zeros(size(vec_hori_data)); 
                V = zeros(size(vec_hgt_data)); 
                ind_hori_mid = round(length(vec_hori)/2);
                ind_hgt_mid  = round(length(vec_hgt)/2);
                
                U(ind_hgt_mid, ind_hori_mid) = refvec_hori_mag*vec_hori_factor;
                
                h_quiver = quiver(X, Y, U, V, 0, 'color', vector_color, 'linewidth', vector_linewidth);
                hU = get(h_quiver,'UData') ;
                hV = get(h_quiver,'VData') ;
                set(h_quiver,'UData',qscale*hU,'VData',qscale*hV) 
                hold on;
                text(X(ind_hgt_mid+1,ind_hori_mid+1), Y(ind_hgt_mid+1,ind_hori_mid+1),...
                         string(refvec_hori_mag), 'FontSize', fontsize,...
                         'verticalalignment','top',...
                         'horizontalalignment','left');
                hold on;
                
%               %/ Vertical reference vector
                X = vec_hori_2D;
                Y = vec_hgt_2D;
                %/ [IMPORTANT]: Make a U / V matrix with the reference vector and zeros elsewhere to get the correct quiver scale!
                U = zeros(size(vec_hori_data)); 
                V = zeros(size(vec_hgt_data)); 
                ind_hori_mid = round(length(vec_hori)/2);
                ind_hgt_mid  = round(length(vec_hgt)/2);
                
                V(ind_hgt_mid, ind_hori_mid) = refvec_hgt_mag*vec_hgt_factor;
                
                h_quiver = quiver(X, Y, U, V, 0, 'color', vector_color, 'linewidth', vector_linewidth);
                hU = get(h_quiver,'UData') ;
                hV = get(h_quiver,'VData') ;
                set(h_quiver,'UData',qscale*hU,'VData',qscale*hV) 
                hold on;
                text(X(ind_hgt_mid-3,ind_hori_mid-3), Y(ind_hgt_mid-3,ind_hori_mid-3),...
                         string(refvec_hgt_mag), 'FontSize', fontsize,...
                         'verticalalignment','bottom',...
                         'horizontalalignment','right');
                hold on;
                
                %/ Combined reference vector (Seems problematic.)
%                 X = vec_hori_2D;
%                 Y = vec_hgt_2D;
%                 %/ [IMPORTANT]: Make a U / V matrix with the reference vector and zeros elsewhere to get the correct quiver scale!
%                 U = zeros(size(vec_hori_data)); 
%                 V = zeros(size(vec_hgt_data)); 
%                 ind_hori_mid = round(length(vec_hori)/2);
%                 ind_hgt_mid  = round(length(vec_hgt)/2);
%                 
%                 U(ind_hgt_mid, ind_hori_mid) = refvec_hori_mag;
%                 V(ind_hgt_mid, ind_hori_mid) = refvec_hgt_mag;
%                 
%                 h_quiver = quiver(X, Y, U, V, 0, 'color', vector_color, 'linewidth', vector_linewidth);
%                 hU = get(h_quiver,'UData') ;
%                 hV = get(h_quiver,'VData') ;
%                 set(h_quiver,'UData',qscale*hU,'VData',qscale*hV) 
%                 hold on;

%                 vec_hori_data_ref = zeros(size(vec_hori_data));
%                 vec_hgt_data_ref  = zeros(size(vec_hgt_data));
%                 ind_hori_mid      = round(length(vec_hori)/2);
%                 ind_hgt_mid       = round(length(vec_hgt)/2);
%                 vec_hori_data_ref(ind_hgt_mid, ind_hori_mid) = refvec_hori_mag;
%                 %/ Draw an eastward vector at the center.
%                 h_quiver = quiver(vec_hori_2D, vec_hgt_2D, vec_hori_data_ref, vec_hgt_data_ref, 0, 'color', vector_color, 'linewidth', vector_linewidth); 
%             
            else
                h_quiver = quiver(vec_hori_2D  (2:vector_step_hgt:end-1, 2:vector_step_hori:end-1),...
                                  vec_hgt_2D   (2:vector_step_hgt:end-1, 2:vector_step_hori:end-1),...
                                  vec_hori_data(2:vector_step_hgt:end-1, 2:vector_step_hori:end-1)*vec_hori_factor,...
                                  vec_hgt_data (2:vector_step_hgt:end-1, 2:vector_step_hori:end-1)*vec_hgt_factor, 0, 'color', vector_color, 'linewidth', vector_linewidth);
            end

            hU = get(h_quiver,'UData') ;
            hV = get(h_quiver,'VData') ;
            set(h_quiver,'UData',qscale*hU,'VData',qscale*hV) 
%             drawnow; pause(0.05);
%             hold on;
            fprintf('*** Spatial avg of time-mean vec_hori_data = %.2f ***\n', my_nanmean(vec_hori_data, 'all'));
            fprintf('*** Spatial avg of time-mean vec_hgt_data = %.2f ***\n', my_nanmean(vec_hgt_data, 'all'));
        end
       
        if show_refvec_only == 0
            %======== cont ========% %/ avoid being overlapped by vectors.
            if ~isempty(cont_data)
                cont_labelsize = fontsize*0.9;
                for cc = 1:length(cont_levels)
                    if cont_levels(cc) == 0 && skip_zero_cont
                        continue; 
                    end

                    cont_color = cont_colmap(cc,:);
                    [CS_cont, h_cont, ] = contour(cont_hori, cont_hgt, cont_data', [cont_levels(cc), cont_levels(cc)],...
                                                  'LineWidth', cont_linewi, 'LineStyle','-', 'Color',cont_color);
                    if ~isempty(cont_labelsize)
                        clabel(CS_cont, h_cont, 'Color', cont_colmap(cc,:), 'FontSize',cont_labelsize,'FontWeight','bold', 'labelspacing',2000); %, 'BackgroundColor',[1 1 1]);
                        hold on; %/ if not hold on, then changes in clabel() will not be retained.
                    end
                end
    %             drawnow; pause(0.05);
                fprintf('*** Spatial avg of time-mean cont_data = %.2f ***\n', my_nanmean(cont_data, 'all'));
            end

            %======== Terrain / Mountain / Topography (mean / median) ========%
            if ~isempty(terrain_mode)
                if isempty(terrain_bndry)             error('terrain_bndry not provided!');                                                     end
                if zm_or_mm ~= 1 && zm_or_mm ~= 2     error('terrain_mode codes not set for the input zm_or_mm!');                              end
                if cbar_mode                          warning('Do not turn on cbar_mode! Since the cbar will shift the new axes used by topo!');  end
                if isempty(topo_color)                topo_color = [.2 .2 .2];   end
                if isempty(terrain_res)               terrain_res = 'low';      end
                %/ select the range to take zonal/meridonal average
                if zm_or_mm == 1
                    reg_topo_lon = [min(terrain_bndry(:,1)), max(terrain_bndry(:,1))];  %/ take the zonal mean of topo in the lon range of the basin.
                    reg_topo_lat = contf_hori;
                else
                    reg_topo_lon = contf_hori;
                    reg_topo_lat = [min(terrain_bndry(:,2)), max(terrain_bndry(:,2))];  %/ take the meridional mean of topo in the lat range of the basin.
                end

                map_proj = 'Miller Cylindrical'; 
                m_proj(map_proj);
                
                %/ Use high or low-res topo
                if isequal(terrain_res, 'high')
                    [topo, topo_lon_2D, topo_lat_2D] = m_tbase([min(reg_topo_lon), max(reg_topo_lon), min(reg_topo_lat), max(reg_topo_lat)]); %REGION =[west east south north];
                elseif isequal(terrain_res, 'low') 
                    [topo, topo_lon_2D, topo_lat_2D] = m_elev([min(reg_topo_lon), max(reg_topo_lon), min(reg_topo_lat), max(reg_topo_lat)]); 
                elseif isequal(terrain_res, 'follow_contf') %/ will follow the grid res of contf_hori.
                    topo_res = abs(diff(contf_hori(1:2)));
                    lon_grids = min(reg_topo_lon):topo_res:max(reg_topo_lon);
                    lat_grids = min(reg_topo_lat):topo_res:max(reg_topo_lat);
                    
                    topo = interp_topo('lon_grids', lon_grids, 'lat_grids', lat_grids)';
                    [topo_lon_2D, topo_lat_2D] = meshgrid(lon_grids, lat_grids);
                else
                    error('Invalid input of ''terrain_res''!')
                end
                
                topo        = flip(topo, 1);
                topo_lat_2D = flip(topo_lat_2D, 1);
                topo        = topo/1e3;            %/ m to km.
                topo        = topo';               %/ make sure in (lon,lat) dim
                topo_lon    = topo_lon_2D(1,:)';
                topo_lat    = topo_lat_2D(:,1);

                %/ First, set negative value (below sealevel) as nan (avoid taking mean with it).
                cond = (topo < 0);
                topo(cond) = nan;              
                %/ Second, take true mean terrain (oceans = 0)
                if terrain_mode == 1
                    topo = my_nanmean(topo, zm_or_mm); 

                %/ Or take mean terrain (oceans = nan)
                elseif terrain_mode == 2
                    topo = my_nanmean(topo, zm_or_mm); 
                else
                    error('Wrong input of terrain_mode!'); 
                end
                 %/ Third, restore NaN to zero in the zm / mm topo.
                topo(isnan(topo)) = 0; 

                if     zm_or_mm == 1     topo_hori = topo_lat;  
                elseif zm_or_mm == 2     topo_hori = topo_lon;  end

                ax_topo = axes;                                             %/ IMPORTANT: axes is to create a new axes!
                set(ax_topo, 'Visible', 'off')                              %/ IMPORTANT: hide the 2nd axes 
                hold on;
                bar(topo_hori', topo, 'Barwidth', 1, 'edgecolor', 'none', 'facecolor', topo_color, 'FaceAlpha', 1); %'facecolor', [.3 .3 .3], 'FaceAlpha', 0.9);

    %             xline(bndry_edge(1), 'r--', 'linewidth', linewidth*0.5);  %/ Optional
    %             xline(bndry_edge(2), 'r--', 'linewidth', linewidth*0.5);
                linkaxes([ax_ori ax_topo], 'xy');                                 %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
            end

            %======== border lines =======%
            %/ Do NOT put this after set(gca), or will cause annoying margins of contf!
            if ~isempty(border_lines) 
    %             bar(topo_hori, topo, 'Barwidth', 1, 'edgecolor', 'none', 'facecolor', [.6 .6 .6], 'FaceAlpha', 1); %'facecolor', [.3 .3 .3], 'FaceAlpha', 0.9);

                if isempty(border_color)       border_color = 'k';          end
                if isempty(border_linestyle)   border_linestyle = '--';     end
                if isempty(border_linewidth)
                    if isequal(border_linestyle, ':')
                        border_linewidth = linewidth*2;
                    else
                        border_linewidth = linewidth*1.5;
                    end
                end
                
                for i = 1:length(border_lines)
                    if size(border_color, 1)    ~= 1   border_color_bc     = border_color(i,:);    else  border_color_bc     = border_color;       end
                    if iscell(border_linestyle) && length(border_linestyle) ~= 1   border_linestyle_bc = border_linestyle{i};  else  border_linestyle_bc = border_linestyle;   end
                    if length(border_linewidth) ~= 1   border_linewidth_bc = border_linewidth(i);  else  border_linewidth_bc = border_linewidth;   end

                    xline(border_lines(i), 'color', border_color_bc, 'linestyle', border_linestyle_bc, 'linewidth', border_linewidth_bc);
                end
                hold on;

                ax_border = axes;                                             %/ IMPORTANT: axes is to create a new axes!
                set(ax_border, 'Visible', 'off')                              %/ IMPORTANT: hide the 2nd axes 
                hold on;
                linkaxes([ax_ori ax_border], 'xy');         %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
            end

            %======== line / scatter plot ========%
            if ~isempty(line_data) || ~isempty(scatter_data)
                if cbar_mode     warning('Since cbar_mode is on, always double check if the cbar will shift the new axes used by line_data.');  end

%                 ax_line = axes;                             %/ IMPORTANT: axes is to create a new axes!
%                 ax_line.Clipping = 'off';                   %/ IMPORTANT: turn clipping off. (allow to draw a point outside the axes!
%                 set(ax_line, 'Visible', 'off')              %/ IMPORTANT: hide the 2nd axes 
%                 hold on;

                if ~isempty(line_data)
                    if add_topo_to_line_data  %/ Assume all in km.
                        topo_res = abs(diff(line_hori(1:2)));
                        topo_lon = reg_topo_lon(1):topo_res:reg_topo_lon(end);
                        topo_lat = reg_topo_lat(1):topo_res:reg_topo_lat(end);
                        topo = interp_topo('lon_grids', topo_lon, 'lat_grids', topo_lat)';
                        topo        = topo/1e3;            %/ m to km.
                        topo        = topo';               %/ make sure in (lon,lat) dim
                        
                        %/ First, set negative value (below sealevel) as nan (avoid taking mean with it).
                        cond = (topo < 0);
                        topo(cond) = nan;              
                        %/ Second, take true mean terrain (oceans = 0)
                        if terrain_mode == 1
                            topo = my_nanmean(topo, zm_or_mm); 

                        %/ Or take mean terrain (oceans = nan)
                        elseif terrain_mode == 2
                            topo = my_nanmean(topo, zm_or_mm); 
                        else
                            error('Wrong input of terrain_mode!'); 
                        end
                         %/ Thrid, restore NaN to zero in the zm / mm topo.
                        topo(isnan(topo)) = 0; 
                        
                        
                        %/ Make sure the hori dim of line_data and topo
                        %/ matches together. Sometimes their orders may be the reverse.
                        if zm_or_mm == 1
                            ind = findismember_loop(topo_lat, line_hori);
                            line_data_add_topo = line_data+topo(ind)';
                        else
                            ind = findismember_loop(topo_lon, line_hori);
                            line_data_add_topo = line_data+topo(ind);
                        end
                        plot(line_hori, line_data_add_topo, 'color', [255 0 204]./255, 'linewidth', linewidth*2/3, 'linestyle', '-');
                    else
                        plot(line_hori, line_data, 'color', [255 0 204]./255, 'linewidth', linewidth*2/3, 'linestyle', '-');
                    end
                end

                if ~isempty(scatter_data)
                    scatter_facecolor  = [255 242  0]./255;
                    scatter_edgecolor  = [  0   0  0]./255;
                    color_alpha        = 1;

                    if ~isempty(scatter_size)
                        %/ Normalize the size
                        scatter_size = (scatter_size - min(scatter_size))/(max(scatter_size) - min(scatter_size));
                        scatter_size = 100 + round(600*scatter_size);
                        scatter_size(isnan(scatter_size)) = 1.e-15;
                        scatter_size(scatter_size == 0)   = 1.e-15;  %/ since size must not be zero or contain NaN.
                    else
                        scatter_size = 300;
                    end

                    if isequal(scatter_hgt, 'top')       scatter_hgt_pos = max(contf_hgt)*1.045;   end
                    if isequal(scatter_hgt, 'bottom')    scatter_hgt_pos = min(contf_hgt);   end

                    sc = scatter(scatter_data, scatter_hgt_pos, scatter_size, 'Marker', 'v',...
                                'MarkerFaceColor', scatter_facecolor, 'MarkerEdgeColor', scatter_edgecolor, 'linewidth', linewidth); 
                    sc.MarkerFaceAlpha = color_alpha;
                end

                ax_line = axes;                             %/ IMPORTANT: axes is to create a new axes!
                ax_line.Clipping = 'off';                   %/ IMPORTANT: turn clipping off. (allow to draw a point outside the axes!
                set(ax_line, 'Visible', 'off')              %/ IMPORTANT: hide the 2nd axes 
                hold on;
                linkaxes([ax_ori, ax_line], 'xy');                        %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
            end
        end
        %======== axis/title/figname ========%
        if isempty(map_xticks)          
            if zm_or_mm == 1
                map_xticks = [-90:5:90]; %/ For lat
            elseif zm_or_mm == 2
                map_xticks = [0:10:360]; %/ For lon
            end
        end
%         if isempty(map_xticks)          map_xticks = contf_hori(mod(contf_hori, 5) == 0); end
        if isempty(map_yticks)          map_yticks = contf_hgt(1:5:end);                  end 
        if isempty(map_xticklabels)     map_xticklabels = map_xticks;   end
        if isempty(map_yticklabels)     map_yticklabels = map_yticks;   end
            
        if grid_mode == 0               %/ default setting
            if isempty(map_xlim)   map_xlim        = [min(contf_hori) max(contf_hori)]; end
            if isempty(map_ylim)   map_ylim        = [min(contf_hgt)  max(contf_hgt)];  end  %/ IMPORTANT: Assume contf_hgt is in km. Manually cut it down to 8km.
            
            %/ flip it if not strictly increasing, otherwise will encounter a bug.
            if diff(map_xticks(1:2)) < 0    map_xticks = flip(map_xticks);   map_xticklabels = flip(map_xticklabels);  end
            if diff(map_yticks(1:2)) < 0    map_yticks = flip(map_yticks);   map_yticklabels = flip(map_yticklabels);  end
            
        elseif grid_mode == 1     %/ this is tailor-made for FLEXPART Lagrangian plot
            if isempty(map_xlim)   map_xlim        = [min(contf_hori)-res_x/2 max(contf_hori)-res_x/2];  end
            if isempty(map_ylim)   map_ylim        = [min(contf_hgt)-res_y/2  max(contf_hgt)-res_y/2];   end  %/ IMPORTANT: Assume contf_hgt is in km. Manually cut it down to 8km.
            
            map_xticks      = contf_hori(1:10:end);
            map_yticks      = contf_hgt(1:5:end)-res_y/2;
            map_xticklabels = map_xticks;   
            map_yticklabels = map_yticks;  
            
        elseif grid_mode == 2  %/ this is tailor-made for FLEXPART Lagrangian plot (if not plotting topo and BLH).
            res_y_label = diff(map_yticklabels(1:2));
            ind_8km = min(find(map_yticklabels >= 8));
            
            if isempty(map_xlim)   map_xlim        = [min(contf_hori)-res_x/2 max(contf_hori)-res_x/2];  end
            if isempty(map_ylim)   map_ylim        = [min(contf_hgt)-res_y/2                  ind_8km];  end %/ IMPORTANT: Assume contf_hgt is in km. Manually cut it down to 8km.
            
            map_xticks      = contf_hori(1:10:end);
            map_yticks      = contf_hgt(1:5:end)-res_y/2;
            map_xticklabels = map_xticklabels(1:10:end);  
            map_yticklabels = map_yticklabels(1:5:end)-res_y_label/2;      
        end  
        map_xlabel = ''; 
        map_ylabel = ''; 
        %/ Adding map_ylabel will will distort figure if line_data is not empty.
%         if isempty(line_data)
%             if isempty(map_ylabel)              map_ylabel = 'z (km)';          end
%         end

        %/ Adding map_xlabel will will distort the border lines (maybe not?).
%         if isempty(border_lines)
%             if     zm_or_mm == 1       map_xlabel = sprintf('Lat (%s)', char(176));         
%             elseif zm_or_mm == 2       map_xlabel = sprintf('Lon (%s)', char(176));  end  
%         end
        
        xlabel(ax_ori,  map_xlabel);
        ylabel(ax_ori, map_ylabel);     
        xlim(ax_ori,  map_xlim)
        ylim(ax_ori,  map_ylim)
        xticks(ax_ori,  map_xticks)     
        yticks(ax_ori,  map_yticks)
        xticklabels(ax_ori,  map_xticklabels)
        yticklabels(ax_ori,  map_yticklabels)

        set(ax_ori, 'linewidth', linewidth, 'fontsize', fontsize)
%         drawnow; pause(0.05);
        hold on;
        
        if isempty(title_fontsize) title_fontsize = fontsize*1.5;       end
        if isempty(title_pos)      title_pos      = [0.91,0.98,0.2,0];  end
        
        %/ Annotation is indep. of axes -> not distorting the aspect ratio.
        annotation( 'textbox', 'String', titlename, 'Color', 'k', ...
                    'FontSize', title_fontsize, 'Units', 'normalized', 'EdgeColor', 'none', ...
                    'Position', title_pos)

        
        %/ Set title after set gca!!
%         title(ax_ori,  titlename, 'Fontsize', fontsize*1.5);
    %     titlePos = get( th , 'position');
    % %     titlePos(2) = titlePos(2) + 0.2;
    %     set(th , 'Fontsize', fontsize - 10, 'position' , titlePos);
%         drawnow; pause(0.05);
        
        %/ grid off only works without set() function after it!
        box on
        grid off
    end
    
    FigName_underscore = strcat(FigName_underscore, str_qscale, str_cbar, str_refvec); %/ append vars to FigName_underscore.
    if savefig
        final_figname = char(strcat(plotting_folder, FigName_underscore,'.pdf'));
        export_fig(final_figname,'-pdf','-painters', '-c[inf, inf, inf, inf]', '-nocrop', '-transparent'); %, '-nocrop'); %
        fprintf('*** Figure saved into %s. ***\n', final_figname)
    end

end
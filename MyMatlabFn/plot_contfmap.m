function plot_contfmap(varargin)
    
    % create a set of valid parameters and their default value
    pnames = {'contf_data', 'contf_lon', 'contf_lat', 'contf_levels', 'contf_unit', 'colmap', 'cbar_interval', 'pcolor_mode', 'shadedrelief_mode',...
              'cont_data', 'cont_data_raw', 'cont_lon', 'cont_lat', 'cont_levels', 'cont_colmap', 'cont_linewi', 'cont_labelsize', 'cont_label_col', 'skip_zero_cont',...
              'point_data', 'marker', 'markersize', 'markerfacecolor', 'markeredgecolor', 'linewi', 'color', 'patch_color', 'bndry_data', 'bndry_patch_mode', 'patch_alpha', 'draw_bndry_onebyone', 'text_data', 'text_backgroundcolor', 'text_edgecolor', 'text_fontsize', 'titlename', 'title_pos', 'savepath', 'fig_fmt',...
              'point2_data', 'marker2', 'marker2size', 'marker2facecolor', 'marker2edgecolor',...
              'Udata',   'Vdata' , 'uv_lon', 'uv_lat',  'vec_step_lon',  'vec_step_lat',  'vector_levels',  'vector_color',  'vector_edgecolor',  'vecscale',  'vecscale2',  'shaftwidth',  'headlength',  'vec_lbs',  'vec_mag_ref',  'vec_ref_fontsize',  'vec_ref_lat_shift', 'curvature', 'show_refvec', 'draw_refvec_only',...
              'U2data', 'V2data', 'uv2_lon', 'uv2_lat', 'vec2_step_lon', 'vec2_step_lat', 'vector2_levels', 'vector2_color', 'vector2_edgecolor', 'vec2scale', 'vec2scale2', 'shaft2width', 'head2length', 'vec2_lbs', 'vec2_mag_ref', 'vec2_ref_fontsize', 'vec2_ref_lat_shift',...
              'glb_data_mode', 'glb_plateau_mode', 'plateau_hgt', 'plateau_col',...
              'gcf_position', 'map_proj', 'map_lon_lower', 'map_lon_upper', 'map_lat_lower', 'map_lat_upper', 'map_center', 'coast_col', 'coast_wi', 'coast_patch_col', 'fontsize', 'title_fontsize', 'create_fig', 'grid_mode', 'grid_linewi', 'bndry_linewi',...
              'cbar_mode', 'cbar_YTick', 'cbar_YTickLabel', 'xaxisloc', 'backcolor', 'cbar_location', 'draw_cbar_only', 'cbar_fontsize', 'cbar_position',...
              'traj_data', 'traj_time', 'traj_levels', 'traj_colmap', 'traj_linewi', 'traj_unit', 'hatch_data', 'hatch_lon', 'hatch_lat', 'hatch_thres_pve', 'hatch_thres_nve', 'color_hatch_pve', 'color_hatch_nve', 'hatch_mode', 'hatch_linewi', 'hatch_intvl',...
              'mask_mountains', 'draw_province', 'draw_river', 'draw_country', 'ax_panel', 'trans_bg', 'png_dpi'};
    
    dflts  = {[] [] [] [] [] [] [] 0  0 ...
              [] [] [] [] [] [] 1  1 'k' 1 ...
              [] 'o' 1 'r' [] 2 [] [] [] 0 0.4 [] [] 'none' 'none' 14 [] [] [] 'pdf' ...
              [] 'o' 1 'r' [] ...
              [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] 0 0 0 ...
              [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] ...
              0  [] [] [] ...
              [] [] [] [] [] [] [] [] 1.5 [] 14 [] 1 1 [] []...
              1  [] [] 'bottom' 'none' [] 0 [] []...
              [] [] [] [] [] [] [] [] [] [] [] [] [] 1  1.3  8 ...
              0  0  0  0  []  0 []};
          
    %/ parse function arguments
    [ contf_data, contf_lon, contf_lat, contf_levels, contf_unit, colmap, cbar_interval, pcolor_mode, shadedrelief_mode,...
      cont_data, cont_data_raw, cont_lon,  cont_lat,  cont_levels,  cont_colmap, cont_linewi, cont_labelsize, cont_label_col, skip_zero_cont,...
      point_data,  marker,  markersize,  markerfacecolor,  markeredgecolor, linewi, color, patch_color, bndry_data, bndry_patch_mode, patch_alpha, draw_bndry_onebyone, text_data, text_backgroundcolor, text_edgecolor, text_fontsize, titlename, title_pos, savepath, fig_fmt,...
      point2_data, marker2, marker2size, marker2facecolor, marker2edgecolor,...
      Udata,  Vdata,  uv_lon,  uv_lat,  vec_step_lon,  vec_step_lat,  vector_levels,  vector_color,  vector_edgecolor,  vecscale,  vecscale2,  shaftwidth,  headlength,  vec_lbs,  vec_mag_ref,  vec_ref_fontsize,  vec_ref_lat_shift, curvature, show_refvec, draw_refvec_only,...
      U2data,V2data, uv2_lon, uv2_lat, vec2_step_lon, vec2_step_lat, vector2_levels, vector2_color, vector2_edgecolor, vec2scale, vec2scale2, shaft2width, head2length, vec2_lbs, vec2_mag_ref, vec2_ref_fontsize, vec2_ref_lat_shift,...
      glb_data_mode, glb_plateau_mode, plateau_hgt, plateau_col,...
      gcf_position, map_proj, map_lon_lower, map_lon_upper, map_lat_lower, map_lat_upper, map_center, coast_col, coast_wi, coast_patch_col, fontsize, title_fontsize, create_fig, grid_mode, grid_linewi, bndry_linewi,...
      cbar_mode, cbar_YTick, cbar_YTickLabel, xaxisloc, backcolor, cbar_location, draw_cbar_only, cbar_fontsize, cbar_position,...
      traj_data, traj_time, traj_levels, traj_colmap, traj_linewi, traj_unit, hatch_data, hatch_lon, hatch_lat, hatch_thres_pve, hatch_thres_nve, color_hatch_pve, color_hatch_nve, hatch_mode, hatch_linewi, hatch_intvl,...
      mask_mountains, draw_province, draw_river, draw_country, ax_panel, trans_bg, png_dpi] ...
                 = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    if length(pnames) ~= length(dflts)
        error('Inconsistent length of pnames and dflts.');
    end

    %%
    %======================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Mar 22, 2024
    %======================================================================
    % create_fig = 1; gcf_position = []; shadedrelief_mode = [];
    % draw_province = []; bndry_patch_mode = []; markeredgecolor = [];
    % draw_river = []; draw_country = []; traj_data = []; point2_data = [];
    % text_data = []; trans_bg = []; savepath = [];

    %/ Make sure they are all in double
    contf_data    = double(contf_data);
    contf_lon     = double(contf_lon);
    contf_lat     = double(contf_lat);
    cont_data     = double(cont_data);
    cont_data_raw = double(cont_data_raw);
    cont_lon      = double(cont_lon);
    cont_lat      = double(cont_lat);
    Udata         = double(Udata);
    Vdata         = double(Vdata);
    uv_lon        = double(uv_lon);
    uv_lat        = double(uv_lat);
    U2data        = double(U2data);
    V2data        = double(V2data);
    uv2_lon       = double(uv2_lon);
    uv2_lat       = double(uv2_lat);

    %====== Pre-processing ======%
    if isempty(map_lon_lower) && isempty(map_lon_upper) && isempty(map_lat_lower) && isempty(map_lat_upper) 
        map_lon_lower = min(contf_lon);
        map_lon_upper = max(contf_lon);
        map_lat_lower = min(contf_lat);
        map_lat_upper = max(contf_lat);
    end
    if isempty(cbar_interval)   cbar_interval = 2;      end
    % if isempty(curvature)       curvature = 0;          end
    
    if glb_data_mode || (map_lon_lower < 0)
        flag_conv_to_circular = 1;   %/ Then shift the lon dim of data to [-179 180] 
    else
        flag_conv_to_circular = 0;
    end
    % flag_conv_to_circular
    if ~isempty(contf_data)     
        if isvector(contf_lon) && isvector(contf_lat)
            %/ transpose to column vector
            if size(contf_lon, 1) == 1                    contf_lon = contf_lon';   end
            if size(contf_lat, 1) == 1                    contf_lat = contf_lat';   end
            if size(contf_data, 1) == length(contf_lon)   contf_data = contf_data'; end
            if flag_conv_to_circular
                [contf_data, contf_lon_2D, contf_lat_2D] = conv_to_circular_data(contf_data, contf_lon, contf_lat, map_lon_upper); 
            else
                [contf_lon_2D, contf_lat_2D] = meshgrid(contf_lon, contf_lat);
            end
        else  %/ otherwise do nothing, since the input lon,lat are already matrices.
            contf_lon_2D = contf_lon;
            contf_lat_2D = contf_lat;
        end
    end
    if ~isempty(cont_data)     
        if isvector(cont_lon) && isvector(cont_lat)
            %/ transpose to column vector
            if size(cont_lon, 1) == 1                     cont_lon = cont_lon';     end
            if size(cont_lat, 1) == 1                     cont_lat = cont_lat';     end
            if size(cont_data, 1) == length(cont_lon)     cont_data = cont_data';   end
            if flag_conv_to_circular
                [cont_data, cont_lon_2D, cont_lat_2D] = conv_to_circular_data(cont_data, cont_lon, cont_lat, map_lon_upper); 
            else
                [cont_lon_2D, cont_lat_2D] = meshgrid(cont_lon, cont_lat);
            end
        else  %/ otherwise do nothing, since the input lon,lat are already matrices.
            cont_lon_2D = cont_lon;
            cont_lat_2D = cont_lat;
        end
    end
    if ~isempty(cont_data_raw)     
        if isvector(cont_lon) && isvector(cont_lat)
            %/ transpose to column vector
            if size(cont_lon, 1) == 1                     cont_lon = cont_lon';             end
            if size(cont_lat, 1) == 1                     cont_lat = cont_lat';             end
            if size(cont_data_raw, 1) == length(cont_lon) cont_data_raw = cont_data_raw';   end
            if flag_conv_to_circular
                [cont_data_raw, cont_lon_2D, cont_lat_2D] = conv_to_circular_data(cont_data_raw, cont_lon, cont_lat, map_lon_upper); 
            else
                [cont_lon_2D, cont_lat_2D] = meshgrid(cont_lon, cont_lat);
            end
        else  %/ otherwise do nothing, since the input lon,lat are already matrices.
            cont_lon_2D = cont_lon;
            cont_lat_2D = cont_lat;
        end
    end
    if ~isempty(Udata)     
        if isvector(uv_lon) && isvector(uv_lat)
            %/ transpose to column vector
            if size(uv_lon, 1) == 1               uv_lon = uv_lon'; end
            if size(uv_lat, 1) == 1               uv_lat = uv_lat'; end
            if size(Udata, 1) == length(uv_lon)   Udata  = Udata';  end
            if flag_conv_to_circular
                [Udata, x, y] = conv_to_circular_data(Udata, uv_lon, uv_lat, map_lon_upper); 
            else
                [x, y] = meshgrid(uv_lon, uv_lat);
            end
        else  %/ otherwise do nothing, since the input lon,lat are already matrices.
            x = uv_lon;
            y = uv_lat;
        end
    end
    if ~isempty(Vdata)     
        if isvector(uv_lon) && isvector(uv_lat)
            %/ transpose to column vector
            if size(uv_lon, 1) == 1               uv_lon = uv_lon'; end
            if size(uv_lat, 1) == 1               uv_lat = uv_lat'; end
            if size(Vdata, 1) == length(uv_lon)   Vdata  = Vdata';  end
            if flag_conv_to_circular
                [Vdata, x, y] = conv_to_circular_data(Vdata, uv_lon, uv_lat, map_lon_upper); 
            else
                [x, y] = meshgrid(uv_lon, uv_lat);
            end
        else  %/ otherwise do nothing, since the input lon,lat are already matrices.
            x = uv_lon;
            y = uv_lat;
        end
    end
    if ~isempty(U2data)     
        if isempty(uv2_lon)            uv2_lon = uv_lon;  end
        if isempty(uv2_lat)            uv2_lat = uv_lat;  end

        if isvector(uv2_lon) && isvector(uv2_lat)
            %/ transpose to column vector
            if size(uv2_lon, 1) == 1                uv2_lon = uv2_lon'; end
            if size(uv2_lat, 1) == 1                uv2_lat = uv2_lat'; end
            if size(U2data, 1) == length(uv2_lon)   U2data  = U2data';  end
            if flag_conv_to_circular
                [U2data, x2, y2] = conv_to_circular_data(U2data, uv2_lon, uv2_lat, map_lon_upper); 
            else
                [x2, y2] = meshgrid(uv2_lon, uv2_lat);
            end
        else  %/ otherwise do nothing, since the input lon,lat are already matrices.
            x2 = uv2_lon;
            y2 = uv2_lat;
        end
    end
    if ~isempty(V2data)    
        if isempty(uv2_lon)            uv2_lon = uv_lon;  end
        if isempty(uv2_lat)            uv2_lat = uv_lat;  end

        if isvector(uv2_lon) && isvector(uv2_lat)
            %/ transpose to column vector
            if size(uv2_lon, 1) == 1                uv2_lon = uv2_lon'; end
            if size(uv2_lat, 1) == 1                uv2_lat = uv2_lat'; end
            if size(V2data, 1) == length(uv2_lon)   V2data  = V2data';  end
            if flag_conv_to_circular
                [V2data, x2, y2] = conv_to_circular_data(V2data, uv2_lon, uv2_lat, map_lon_upper); 
            else
                [x2, y2] = meshgrid(uv2_lon, uv2_lat);
            end
        else  %/ otherwise do nothing, since the input lon,lat are already matrices.
            x2 = uv2_lon;
            y2 = uv2_lat;
        end
    end
    
    %/ If the contf_levels are uneven, replace values of contf_data with indices
    if ~isempty(contf_data)     
        if isempty(contf_levels)
            if length(unique(contf_data)) == 1   %/ If the matrix contains only one unique number (e.g., 0)
                contf_levels = [unique(contf_data)-1, 0, unique(contf_data)+1];
            else
                contf_levels = linspace(double(min(contf_data, [], 'all', 'omitnan')), double(max(contf_data, [], 'all', 'omitnan')), 13);
            end
        end
        if isempty(colmap)
            colmap = nclCM('BlueWhiteOrangeRed', length(contf_levels)-1);
        end
        if length(unique(string(diff(contf_levels)))) > 1
            flag_uneven    = 1;
            contf_data_ori = contf_data;  %/ Keep this original data for computing min, mean and max.
            [contf_data, contf_levels, cbar_YTickLabel, cbar_YTick] = contfdata4uneven('contf_data', contf_data, 'contf_levels', contf_levels, 'cbar_YTickLabel', cbar_YTickLabel);
        else
            flag_uneven    = 0;
            contf_data_ori = [];
        end
    end
    %===============
    if isempty(map_proj)
        if glb_data_mode
            map_proj = 'robin'; 
        else
            map_proj = 'Miller Cylindrical'; 
        end
    end
    str_relief = '';
    
    if draw_cbar_only
        figure
        set(gcf, 'Color', 'None');
        set(gca, 'Color', 'None');
    %     set(gcf, 'color','w');
    %     set(gcf, 'Position', get(0, 'Screensize'));
        
        if isempty(cbar_location)  cbar_location = 'eastoutside';    end
        if isequal(cbar_location, 'eastoutside')
            if isempty(cbar_position)  cbar_position = [100 100 900 200]; end  %/ x y w h
    %         if isempty(cbar_position)  cbar_position = [100 100 1100 650]; end
        elseif isequal(cbar_location, 'southoutside')
            if isempty(cbar_position)  cbar_position = [100 100 350 500];  end 
        else
            error('code not set for cbar_position == %s!', cbar_position);  
        end
    %     cbar_location
        set(gcf,'position',cbar_position) % sets figure size, x y w h 
        
        fprintf('*** Drawing cbar only... *** \n');
        m_proj('Miller Cylindrical','longitudes',[map_lon_lower map_lon_upper], 'latitudes', [map_lat_lower map_lat_upper]);
        hold on;
    
        m_contourf([1, 1; 1, 1], [1, 1; 1, 1], [nan,nan ; nan, nan], [-999999999 1],'edgecolor','none');  %/ just an empty contf_map to control the position of cbar.
    %     m_contourf([1, 1; 1, 1], [1, 1; 1, 1], [1, 1; 1, 1], [-999999999 1],'edgecolor','none');  %/ just an empty contf_map to control the position of cbar.
        shading interp
        
        %/ If contf_levels is empty but not traj_levels, plot the cbar for traj_levels
        if isempty(contf_levels) && ~isempty(traj_levels)
            contf_levels = traj_levels;
            colmap       = traj_colmap;
        end
        
        % min(contf_levels) 
        % max(contf_levels)
        clim([min(contf_levels) max(contf_levels)]);
        colormap(gca, colmap)
                            
        if isempty(cbar_YTick) && isempty(cbar_YTickLabel)
            cbar_YTick = contf_levels(2:cbar_interval:end-1);
            cbar_YTickLabel = cbar_YTick;
        end
        if isempty(cbar_fontsize)                   cbar_fontsize = fontsize*0.5;       end
    %     if isequal(cbar_location, 'southoutside')   cbar_fontsize = cbar_fontsize*0.5;  end
        cb = colorbar(cbar_location);
    %     cbar_YTick
    %     cbar_YTickLabel
        set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', cbar_fontsize) %/ cbar Ytick for diverging colormap
        set(get(cb,'Title'),'String',contf_unit, 'Fontsize', cbar_fontsize)
        
        %/ axis off does not work; hide axis by setting color to none in m_grid
        m_grid('xtick', [], 'ytick', [], 'linewi', 0.001, 'linest','none', 'gridcolor', 'none', 'color', 'none');  %<-- important to call this for a correct size of the plot! 
    
        if ~isempty(savepath)
            savepath = strcat(savepath, '_cbar');   %/ automatically add the suffix to distinguish the cbar-only figure
        end
    end
    
    %/ m_vec wind reference
    if draw_refvec_only
        if isempty(vec_mag_ref)         error('Empty vec_mag_ref! Set ''draw_refvec_only = 0'' if you don''t want to draw a reference vector!');   end
        if isempty(vec_ref_fontsize)    vec_ref_fontsize = fontsize*1.25;   end
        
        if create_fig
            close all
            figure
            set(gcf, 'color','w');
            set(gcf, 'Position', get(0, 'Screensize'));
        end
        
        fprintf('*** Drawing reference vector only... *** \n');
        if isequal(map_proj, 'ortho')
            if isempty(map_center)  map_center = [mean([map_lon_lower map_lon_upper]), mean([map_lat_lower map_lat_upper])]; end
    	    m_proj(map_proj,'long', map_center(1), 'lat', map_center(2));  %/ , 'rad', 10, 'rec', 'off' 
            vec_ref_lon = map_center(1);
            vec_ref_lat = map_center(2);
        else
            m_proj(map_proj,'longitudes',[map_lon_lower map_lon_upper], 'latitudes', [map_lat_lower map_lat_upper]);
            vec_ref_lon = mean([map_lon_upper, map_lon_lower]);
            vec_ref_lat = mean([map_lat_upper, map_lat_lower]);
        end
        hold on;
        
        %%%%%%% IMPORTANT: Fixing the fig frame for a correct drawing of vectors! (do NOT remove it!) %%%%%%%
        m_grid('xtick', [], 'ytick', [], 'linewi', 0.001, 'linest','none', 'gridcolor', 'none');
    %     m_coast('linewidth', coast_wi, 'color',[.4 .4 .4]); %/ [247 247 163]./255
    %     Udata_ref = zeros(size(Udata));
    %     Vdata_ref = zeros(size(Vdata));
    %     
    %     cen_ind = [floor(size(y,1)/2), floor(size(x,2)/2)]; 
    %     Udata_ref(cen_ind(1), cen_ind(2)) = vec_mag_ref;
    %     Vdata_ref(cen_ind(1), cen_ind(2)) = 0;
    %     
    %     m_vec(vecscale, x, y, Udata_ref*vecscale2, Vdata_ref*vecscale2,...
    %                    vector_color, 'shaftwidth',shaftwidth ,'headlength',headlength,...
    %                    'EdgeColor', vector_color, 'linewidth', 0.5);    %/ do not use edgecolor.
        
    
    
        [~, ht] = m_vec(vecscale, vec_ref_lon, vec_ref_lat, vec_mag_ref*vecscale2, 0,...
                        vector_color, 'shaftwidth',shaftwidth ,'headlength',headlength, 'key', vec_lbs);
    %                             'EdgeColor', vector_edgecolor, 'linewidth', 0.5);      %/ do not use edgecolor.
        set(ht,'FontSize',vec_ref_fontsize, 'color', vector_color);  %, 'HorizontalAlignment', 'left');
        hold on;
    %     set(hp,'EdgeColor', vector_edgecolor);                             %/ have to set the EdgeColor using set() if using 'key' in m_vec!
    
        %/ Add the reference vector for the 2nd vector data
        if ~isempty(U2data) && ~isempty(V2data)
            [~, ht2] = m_vec(vec2scale, vec_ref_lon, vec_ref_lat-15, vec2_mag_ref*vec2scale2, 0,...
                             vector2_color, 'shaftwidth',shaft2width ,'headlength',head2length, 'key', vec2_lbs);
            set(ht2,'FontSize',vec_ref_fontsize, 'color', vector2_color);  %, 'HorizontalAlignment', 'left');
        end
        hold on
        m_grid('xtick', [], 'ytick', [], 'xticklabel', [], 'yticklabel', [], 'linewi',2, 'linest','none', 'gridcolor', 'k', 'tickdir','in',...
                   'xaxisloc', xaxisloc, 'backcolor', 'none', 'fontsize',fontsize);  %/ call it before m_vec to be able to draw vectors!!
               
        %/ axis off does not work; hide axis by setting color to none in m_grid
    %     m_grid('xtick', [], 'ytick', [], 'linewi', 0.001, 'linest','none', 'gridcolor', 'none');  %<-- important to call this for a correct size of the plot! 
               
        if ~isempty(savepath)
            savepath = strcat(savepath, '_vecref');   %/ automatically add the suffix to distinguish the cbar-only figure
        end
        fig_fmt = 'pdf'; %/ set to pdf automatically
    end
        
    if draw_cbar_only == 0 && draw_refvec_only == 0
        if create_fig
            figure
            set(gcf, 'color','w');
            if ~isempty(gcf_position)
                set(gcf,'position',gcf_position) % sets figure size
            else
                set(gcf, 'Position', get(0, 'Screensize'));
            end
        end
        
        if shadedrelief_mode 
    %         warning('we have to use ''equidistant cylindrical'' to enable m_shadedrelief!');
    %         m_proj('equidistant cylindrical', 'longitudes',[map_lon_lower map_lon_upper], 'latitudes', [map_lat_lower map_lat_upper]);
            m_proj(map_proj,'longitudes',[map_lon_lower map_lon_upper], 'latitudes', [map_lat_lower map_lat_upper], 'rectbox','on');
            str_relief = '_relief';
            
        elseif isequal(map_proj, 'ortho')
            if isempty(map_center)  map_center = [mean([map_lon_lower map_lon_upper]), mean([map_lat_lower map_lat_upper])]; end
    %         map_center
            m_proj(map_proj,'long', map_center(1), 'lat', map_center(2));  %/ , 'rad', 10, 'rec', 'off' 
        else
            m_proj(map_proj,'longitudes',[map_lon_lower map_lon_upper], 'latitudes', [map_lat_lower map_lat_upper]);
        end
        hold on;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Contourf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(contf_data)
            colormap(gca, colmap)

            % contf_levels
            % [min(contf_levels) max(contf_levels)]
            clim([min(contf_levels) max(contf_levels)]); %/ must be put before m_shaderelief.
            
            %/ Check if lat and lon were created by [lat,lon] = meshgrid(lat_array,lon_array) 
            %/ or [lon,lat] = meshgrid(lon_array,lat_array)\
            [dlon1,dlon2] = gradient(contf_lon_2D); 
            [dlat1,dlat2] = gradient(contf_lat_2D); 
            flag_interp = 0;
            if isequal(dlat1,zeros(size(contf_lat_2D))) 
                res_lon = mode(dlon1, 'all');  %/ Find the most common resolution (do NOT take abs, otherwise pcolor may shift in the wrong direction!)
                res_lat = mode(dlat2, 'all');  %/ Find the most common resolution (do NOT take abs, otherwise pcolor may shift in the wrong direction!)
                if ~isequal(dlon2,zeros(size(contf_lon_2D)))
                    flag_interp = 1;
                    warning('Detected that contf_lon_2D, contf_lat_2D are non-monotonic, which may distort contourf! Hence, by default, interplating contf_data to %.2f%s x %.2f%s for best visualization...', res_lon, char(176), res_lat, char(176))
                end
            else
                res_lon = mode(dlon2, 'all');  %/ Find the most common resolution (do NOT take abs, otherwise pcolor may shift in the wrong direction!)
                res_lat = mode(dlat1, 'all');  %/ Find the most common resolution (do NOT take abs, otherwise pcolor may shift in the wrong direction!)
                if ~isequal(dlon2,zeros(size(contf_lon_2D)))
                    flag_interp = 1;
                    warning('Detected that contf_lon_2D, contf_lat_2D are non-monotonic, which may distort contourf! Hence, by default, interplating contf_data to %.2f%s x %.2f%s for best visualization...', res_lon, char(176), res_lat, char(176))
                end
            end
            % flag_interp
            if flag_interp
                lon_new = 0:res_lon:360-res_lon;
                lat_new = -90:res_lat:90;
                
                %/ Update
                contf_data = my_interp('lon_old', contf_lon_2D, 'lat_old', contf_lat_2D, 'data', contf_data, 'lon_new', lon_new, 'lat_new', lat_new, 'is_global', 0, 'lon_dim', 1);
                contf_lon = lon_new;
                contf_lat = lat_new;

                [contf_lon_2D, contf_lat_2D] = meshgrid(contf_lon, contf_lat);
                contf_lon_2D = contf_lon_2D'; contf_lat_2D = contf_lat_2D';
            end
            % res_lon = abs(unique(round(diff(contf_lon(1:2)), 3)));  %/ round up to 3 decimal points to avoid numerical errors
            % res_lat = unique(round(diff(contf_lat(1:2)), 3));  %/ Do NOT take absolute value, this should be correct. 
                
            if shadedrelief_mode %/ for drawing topography
                disp('m_shadedrelief is used');
                find(isnan(contf_data))
                m_shadedrelief(contf_lon, contf_lat, contf_data, 'lightangle',-45, 'gradient', 8, 'coord','geog');
                hold on;
                
            elseif pcolor_mode
                % res_lon
                % res_lat

                %/ https://www.mathworks.com/matlabcentral/fileexchange/50706-offsets-and-missing-data-via-pcolor-and-surf'
                warning('Since pcolor has offsets and edge problem, here we solve the offset problem by shifting by a half grid.\nFor more info, see https://www.mathworks.com/matlabcentral/fileexchange/50706-offsets-and-missing-data-via-pcolor-and-surf');
                m_pcolor(contf_lon_2D-res_lon/2, contf_lat_2D-res_lat/2, contf_data, 'edgecolor','none');
                % m_pcolor(contf_lon_2D, contf_lat_2D, contf_data, 'edgecolor','none');
                shading flat
                
                %/ There will be a half-grid blank near the lon boundaries after pcolor grid correction.
                %/ Zoom in the map by half grid in lon range to "remove" the blank (except for map_proj == 'ortho').
                if ~isequal(map_proj, 'ortho')
                    m_proj(map_proj,'longitudes',[map_lon_lower+res_lon/2 map_lon_upper-res_lon/2], 'latitudes', [map_lat_lower map_lat_upper]);
                end
                hold on;
            else
                % m_contourf(contf_lon_2D-res_lon/2, contf_lat_2D-res_lat/2, contf_data, [-999999999 contf_levels],'edgecolor', 'none');
                m_contourf(contf_lon_2D, contf_lat_2D, contf_data, [-999999999 contf_levels],'edgecolor','none');
            end
            
            if cbar_mode
                if isempty(cbar_location)      cbar_location = 'eastoutside';    end 
                if isempty(cbar_fontsize)      cbar_fontsize = fontsize;         end
                cb = colorbar(cbar_location);
                if isempty(cbar_YTick) && isempty(cbar_YTickLabel)
                    cbar_YTick = contf_levels(2:cbar_interval:end-1);
                    cbar_YTickLabel = cbar_YTick;
                end
                set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', cbar_fontsize) %/ cbar Ytick for diverging colormap
                set(get(cb,'Title'),'String',contf_unit, 'Fontsize', cbar_fontsize)
            end
            % set(gca, 'Clipping', 'off');
            % ax = gca;
            % ax.Clipping = 'off';
            % set(cb, 'YTick', contflevels,'Fontsize', 12) %/ cbar Ytick for diverging colormap
            % cbarrow;
            % drawnow; pause(0.05);

            %/ Show the min, mean, and max data based on the domain shown
            logi_lon_2D = (contf_lon_2D >= map_lon_lower) & (contf_lon_2D <= map_lon_upper);
            logi_lat_2D = (contf_lat_2D >= map_lat_lower) & (contf_lat_2D <= map_lat_upper);
            contf_data_bc = contf_data;
            contf_data_bc(~(logi_lon_2D & logi_lat_2D)) = nan;

            if flag_uneven
                fprintf('*** min(contf_data)  = %.2g ***\n',  min(contf_data_ori, [], 'all', 'omitnan'));
                fprintf('*** mean(contf_data) = %.2g ***\n',  mean(contf_data_ori,    'all', 'omitnan'));
                fprintf('*** max(contf_data)  = %.2g ***\n',  max(contf_data_ori, [], 'all', 'omitnan'));
            else
                fprintf('*** min(contf_data)  = %.2g ***\n',  min(contf_data_bc, [], 'all', 'omitnan'));
                fprintf('*** mean(contf_data) = %.2g ***\n',  mean(contf_data_bc,    'all', 'omitnan'));
                fprintf('*** max(contf_data)  = %.2g ***\n',  max(contf_data_bc, [], 'all', 'omitnan'));
            end
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hatch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(hatch_data)
            if isempty(color_hatch_pve)   color_hatch_pve = [0 204 102]./255;  end  %/ green  = +ve
            if isempty(color_hatch_nve)   color_hatch_nve = [153 0 255]./255;  end  %/ purple = -ve
            
            %/ CAVEAT: Using 'edgebased' will occassionally fill some small 'holes', which may be misleading.
            %/         However, not using 'edgebased' may miss out on some sponradic points that are significant.
            edgebased = 1; 
            draw_rings = 0;
            
            if hatch_mode == 1
                if isempty(hatch_thres_pve) && isempty(hatch_thres_nve)
                    %/ no distinguishment
                    hatch_data(~isnan(hatch_data)) = 1;   %/ first turn nonnan to 1
                    hatch_data(isnan(hatch_data)) = 0;    %/ then turn nan to 0  (do not change the order!)
                    hatch_bndry = get_bndry_from_logi('logical_2D', hatch_data, 'bndry_lon', hatch_lon, 'bndry_lat', hatch_lat, 'outputcell', 0,...
                                                      'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper,...
                                                      'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper,...
                                                      'glb_data_mode', glb_data_mode, 'CONN', 4, 'draw_rings', draw_rings, 'edgebased', edgebased);       
                    if ~isempty(hatch_bndry)
                        m_line(hatch_bndry(:,1),  hatch_bndry(:,2), 'color', color_hatch_pve,  'linest', 'none');  %/ remove the outline
                        m_hatch(hatch_bndry(:,1), hatch_bndry(:,2), 'single', 30, hatch_intvl, 'linest', '-', 'linewi', hatch_linewi, 'color', color_hatch_pve); % ...with hatching added
        %                 drawnow; pause(0.05);
                    end
                    
                else
                    %/ +ve hatching
                    hatch_pve = hatch_data;
                    hatch_pve(hatch_pve < hatch_thres_pve) = nan;
                    hatch_pve(~isnan(hatch_pve)) = 1;   %/ first turn nonnan to 1
                    hatch_pve(isnan(hatch_pve)) = 0;    %/ then turn nan to 0  (do not change the order!)
                    hatch_pve_bndry = get_bndry_from_logi('logical_2D', hatch_pve, 'bndry_lon', hatch_lon, 'bndry_lat', hatch_lat, 'outputcell', 0,...
                                                          'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper,...
                                                          'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper,...
                                                          'glb_data_mode', glb_data_mode, 'CONN', 4, 'draw_rings', draw_rings, 'edgebased', edgebased);       
                    if ~isempty(hatch_pve_bndry)
                        m_line(hatch_pve_bndry(:,1),  hatch_pve_bndry(:,2), 'color', color_hatch_pve,  'linest', 'none');  %/ remove the outline
                        m_hatch(hatch_pve_bndry(:,1), hatch_pve_bndry(:,2), 'single', 30, hatch_intvl, 'linest', '-', 'linewi', hatch_linewi, 'color', color_hatch_pve); % ...with hatching added
        %                 drawnow; pause(0.05);
                    end
                    
                    %/ -ve hatching
                    hatch_nve = hatch_data;
                    hatch_nve(hatch_nve > hatch_thres_nve) = nan;
        %             hatch_nve(hatch_nve > -1*hatch_thres) = nan;
                    hatch_nve(~isnan(hatch_nve)) = 1;
                    hatch_nve(isnan(hatch_nve)) = 0;
                    hatch_nve_bndry = get_bndry_from_logi('logical_2D', hatch_nve, 'bndry_lon', hatch_lon, 'bndry_lat', hatch_lat, 'outputcell', 0,...
                                                          'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper,...
                                                          'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper,...
                                                          'glb_data_mode', glb_data_mode, 'CONN', 4, 'draw_rings', draw_rings, 'edgebased', edgebased);
                    if ~isempty(hatch_nve_bndry)
                        m_line(hatch_nve_bndry(:,1),  hatch_nve_bndry(:,2), 'color', color_hatch_nve,  'linest', 'none');  %/ remove the outline
                        m_hatch(hatch_nve_bndry(:,1), hatch_nve_bndry(:,2), 'single', 30, hatch_intvl, 'linest', '-', 'linewi', hatch_linewi, 'color', color_hatch_nve); % ...with hatching added
        %                 drawnow; pause(0.05);
                    end
                end
            end
            
    %         if hatch_mode == 0
    %             %/ +ve hatch 
    %             [ind_lon_sig, ind_lat_sig] = find(~isnan(hatch_data) & hatch_data > hatch_thres);
    %             hatch_pve = nan(length(ind_lon_sig), 2);
    %             for k = 1:length(ind_lon_sig)
    %                 hatch_pve(k,1) = hatch_lon(ind_lon_sig(k));
    %                 hatch_pve(k,2) = hatch_lat(ind_lat_sig(k));
    %             end
    % 
    %             %/ -ve hatch 
    %             [ind_lon_sig, ind_lat_sig] = find(~isnan(hatch_data) & hatch_data < -1*hatch_thres);
    %             hatch_nve = nan(length(ind_lon_sig), 2);
    %             for k = 1:length(ind_lon_sig)
    %                 hatch_nve(k,1) = hatch_lon(ind_lon_sig(k));
    %                 hatch_nve(k,2) = hatch_lat(ind_lat_sig(k));
    %             end
    % 
    %             %/ if plotting in [-179, 180], we convert hatch data's lon to fit that range. (no need to use conv_to_circular_data beforehand)
    %             if map_lon_lower < 0        
    %                 a = hatch_pve(:,1);
    %                 a(a > 180) = a(a > 180) - 360;
    %                 hatch_pve(:,1) = a;
    % 
    %                 a = hatch_nve(:,1);
    %                 a(a > 180) = a(a > 180) - 360;
    %                 hatch_nve(:,1) = a;
    %             end
    %             
    %             %/ default setting 1
    %             hatch_marker = 'o'; hatch_markersize = 1.5; 
    % 
    %             m_line(hatch_pve(:,1), hatch_pve(:,2), 'marker', hatch_marker, 'markersize', hatch_markersize, 'markerfacecolor', color_hatch_pve,...
    %                   'linest', 'none', 'color', 'none', 'linewi', linewi, 'clip', 'point');
    %             drawnow; pause(0.05);
    % 
    %             m_line(hatch_nve(:,1), hatch_nve(:,2), 'marker', hatch_marker, 'markersize', hatch_markersize, 'markerfacecolor', color_hatch_nve,...
    %                   'linest', 'none', 'color', 'none', 'linewi', linewi, 'clip', 'point');
    %             drawnow; pause(0.05);
    %         end
            % fprintf('Time cost in drawing hatched pattern/points: %f \n', toc);
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Coast (patch; plot before Plateau) %%%%%
        if ~isempty(coast_patch_col) && ~isequal(coast_patch_col, 'none')
            m_coast('patch', coast_patch_col, 'edgecolor','none'); %/ [247 247 163]./255
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Global Plateau %%%%%%%%%%%%%%%%%%%%%%%%%
        if glb_plateau_mode
            if isempty(plateau_col)   plateau_col      = [255 51 204]./255;   end
            if isempty(plateau_hgt)   plateau_hgt      = 1500;                end
    
    %         m_coord('geographic');                                      %/ required before calling m_tbase.
    %         [topo, topo_lon_2D, topo_lat_2D] = m_tbase([0, 360, -90, 90]); %REGION =[west east south north];
    %         topo     = topo';               %/ make sure in (lon,lat) dim
    %         topo_lon = topo_lon_2D(1,:)';
    %         topo_lat = topo_lat_2D(:,1);
    % 
    %         topo = topo(1:end-1,2:end);
    %         topo_lon = topo_lon(1:end-1);
    %         topo_lat = topo_lat(2:end);
    % 
    %         [topo_lon_2D, topo_lat_2D] = meshgrid(topo_lon, topo_lat);
    %         [contf_lon_2D, contf_lat_2D] = meshgrid(contf_lon, contf_lat);
    % 
    %         %/ Interpolate high-res topo to the desired resolution (1 deg).
    %         topo_new = interp2(topo_lon_2D,topo_lat_2D,topo',contf_lon_2D,contf_lat_2D,'nearest');
    %         topo_new = topo_new';
    
            %/ output an slighlty coarser averaged topo (from 5-min high topo data)
            res = 0.25; 
            lon_grids = map_lon_lower:res:map_lon_upper;
            lat_grids = map_lat_lower:res:map_lat_upper;
            [lon_grids_2D, lat_grids_2D] = meshgrid(lon_grids,lat_grids);
    
            avg_topo_map = interp_topo('lon_grids', lon_grids, 'lat_grids', lat_grids);
            m_contour(lon_grids_2D, lat_grids_2D, avg_topo_map',...
                              [plateau_hgt, plateau_hgt],'LineWidth', coast_wi, 'Color',plateau_col);
    
    %         m_elev('contour', [plateau_hgt plateau_hgt], 'Color', plateau_col,'LineWidth', coast_wi);   %/ m_elev is more accurate!!! 
    %         drawnow; pause(0.05);
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Coast (outline) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(coast_col)    coast_col = [.4 .4 .4];      end
    
        if isempty(coast_patch_col) || isequal(coast_patch_col, 'none')   %/ Outline coasts if coast patch is empty
            m_coast('linewidth',coast_wi,'color', coast_col);
        end
        if draw_province   m_gshhs('fb2', 'linewidth',coast_wi, 'color', coast_col);    end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Coast (overlay the patch boundary (if any)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if any(bndry_patch_mode == 1)
            m_coast('linewidth',coast_wi,'color', coast_col);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Point (serve as stippling of contf) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(point_data)
            % tic
            if map_lon_lower < 0        %/ if plotting in [-179, 180], we convert point_data's lon to fit that range.
                a = point_data(:,1);
                a(a > 180) = a(a > 180) - 360;
                point_data(:,1) = a;
            end
            if isempty(markerfacecolor)  markerfacecolor = 'none';  end
            if isempty(markeredgecolor)  markeredgecolor = 'none';  end
            
            if size(markerfacecolor, 1) > 1 
                if size(markeredgecolor, 1) == 1
                    markeredgecolor = repmat(markeredgecolor, size(markerfacecolor,1), 1); %/ then replicate it
                end
                %/ Loop to draw points with different facecolors
                for k = 1:size(point_data, 1)
                    m_line(point_data(k,1), point_data(k,2), 'marker', marker, 'markersize', markersize, 'markerfacecolor', markerfacecolor(k,:),...
                            'linest', 'none', 'color', markeredgecolor(k,:), 'linewi', linewi, 'clip', 'point');
                    hold on;
                end
            else
                m_line(point_data(:,1), point_data(:,2), 'marker', marker, 'markersize', markersize, 'markerfacecolor', markerfacecolor,...
                        'linest', 'none', 'color', markeredgecolor, 'linewi', linewi, 'clip', 'point');
                hold on;
            end
    %         alpha(h, 0.5);
            % fprintf('Time cost in drawing points: %f \n', toc);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Contour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(cont_data)
            %/ First, check if lat and lon were created by [lat,lon] = meshgrid(lat_array,lon_array) 
            %/ or [lon,lat] = meshgrid(lon_array,lat_array)\
            [dlon1,dlon2] = gradient(cont_lon_2D); 
            [dlat1,dlat2] = gradient(cont_lat_2D); 
            flag_interp = 0;
            if isequal(dlat1,zeros(size(cont_lat_2D))) 
                res_lon = abs(mode(dlon1, 'all'));  %/ Find the most common resolution (in case of non-monotonic matrix, abs to handle decreasing or increasing lon)
                res_lat = abs(mode(dlat2, 'all'));  %/ Find the most common resolution (in case of non-monotonic matrix, abs to handle decreasing or increasing lat)
                if ~isequal(dlon2,zeros(size(cont_lon_2D)))
                    flag_interp = 1;
                    warning('Detected that cont_lon_2D, cont_lat_2D are non-monotonic, which may distort contourf! Hence, by default, interplating contf_data to %.2f%s x %.2f%s for best visualization...', res_lon, char(176), res_lat, char(176))
                end
            else
                res_lon = abs(mode(dlon2, 'all'));  %/ Find the most common resolution (in case of non-monotonic matrix, abs to handle decreasing or increasing lon)
                res_lat = abs(mode(dlat1, 'all'));  %/ Find the most common resolution (in case of non-monotonic matrix, abs to handle decreasing or increasing lat)
                if ~isequal(dlon2,zeros(size(cont_lon_2D)))
                    flag_interp = 1;
                    warning('Detected that cont_lon_2D, cont_lat_2D are non-monotonic, which may distort contourf! Hence, by default, interplating contf_data to %.2f%s x %.2f%s for best visualization...', res_lon, char(176), res_lat, char(176))
                end
            end

            if flag_interp
                lon_new = 0:res_lon:360-res_lon;
                lat_new = -90:res_lat:90;
                
                %/ Update
                cont_data     = my_interp('lon_old', cont_lon_2D, 'lat_old', cont_lat_2D, 'data', cont_data,     'lon_new', lon_new, 'lat_new', lat_new, 'is_global', 0, 'lon_dim', 1);
                if ~isempty(cont_data_raw)
                    cont_data_raw = my_interp('lon_old', cont_lon_2D, 'lat_old', cont_lat_2D, 'data', cont_data_raw, 'lon_new', lon_new, 'lat_new', lat_new, 'is_global', 0, 'lon_dim', 1);
                end
                cont_lon = lon_new;
                cont_lat = lat_new;

                [cont_lon_2D, cont_lat_2D] = meshgrid(cont_lon, cont_lat);
                cont_lon_2D = cont_lon_2D'; cont_lat_2D = cont_lat_2D';
            end

            cont_colmap_uni = unique(cont_colmap,'rows'); %/ Check how many contour colors are specified

            for cc = 1:length(cont_levels)
                cont_color_bc     = cont_colmap(cc,:);
                cont_linewi_bc    = cont_linewi;
                
                if size(cont_colmap_uni, 1) == 1
                    if cont_levels(cc) >= 0
                        cont_linestyle_bc = '-';
                    else
                        cont_linestyle_bc = '--';    %/ For one contour color, use dashed contour to indicate negative values
                    end
                else
                    cont_linestyle_bc = '-';         %/ by default
                end
                
                if cont_levels(cc) == 0
                    if skip_zero_cont
                        continue;
                    else
                        cont_color_bc     = 'k';
                        cont_linestyle_bc = '--';
                        cont_linewi_bc    = cont_linewi; 
                    end
                end

                if ~isempty(cont_data_raw)   %/ then we plot this raw cont data with dashed contour.
                    m_contour(cont_lon_2D, cont_lat_2D, cont_data_raw,...
                              [cont_levels(cc), cont_levels(cc)],'LineWidth', cont_linewi_bc*.5, 'LineStyle',':','Color',cont_color_bc);
                    hold on;
                end
                
                [CS_cont, h_cont] = m_contour(cont_lon_2D, cont_lat_2D, cont_data,...
                                              [cont_levels(cc), cont_levels(cc)], 'LineWidth',cont_linewi_bc, 'LineStyle',cont_linestyle_bc,'Color',cont_color_bc);
                hold on;
    
                if ~isempty(cont_labelsize) && mod(cc, 2) == 0
                    if isempty(cont_label_col)   cont_label_col = 'k';    end
                    
    %                 if ~isempty(cont_data_raw) %/ since raw cont is plot first, the label will be overlapped by sig. cont!
    %                     clabel(CS_cont_raw, h_cont_raw, 'Color', cont_label_col, 'FontSize',cont_labelsize,'FontWeight','bold', 'labelspacing',2000); %, 'BackgroundColor',[1 1 1]);
    %                 else
                    clabel(CS_cont, h_cont, 'Color', cont_label_col, 'FontSize',cont_labelsize,'FontWeight','bold', 'labelspacing',2000); %, 'BackgroundColor',[1 1 1]);
    %                 end
                    hold on; %/ if not hold on, then changes in clabel() will not be retained.
                end
            end

            %/ Show the min, mean, and max data based on the domain shown
            logi_lon_2D = (cont_lon_2D >= map_lon_lower) & (cont_lon_2D <= map_lon_upper);
            logi_lat_2D = (cont_lat_2D >= map_lat_lower) & (cont_lat_2D <= map_lat_upper);
            cont_data_bc = cont_data;
            cont_data_bc(~(logi_lon_2D & logi_lat_2D)) = nan;
            fprintf('*** min(cont_data)  = %.2g ***\n',  min(cont_data_bc, [], 'all', 'omitnan'));
            fprintf('*** mean(cont_data) = %.2g ***\n',  mean(cont_data_bc,    'all', 'omitnan'));
            fprintf('*** max(cont_data)  = %.2g ***\n',  max(cont_data_bc, [], 'all', 'omitnan'));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(bndry_data)
            if iscell(bndry_data) %/ if bndry_data is a cell array, then loop over it to plot.
                for i = 1:length(bndry_data)
    %                 ind = find(bndry_data{i}(:,1) > 180);
    %                 bndry_data{i}(ind,1) = bndry_data{i}(ind,1) - 360;  %/ convert to [-179 180] if ind is not empty.
                    bndry_data_each = bndry_data{i};
                    if isempty(bndry_data_each)
                        warning('empty bndry_data is detected!')
                        continue;
                    end
                    %/ Check if the lon range is global and starts from negative values (e.g. [-179 180])
    %                 if map_lon_lower < 0 && (map_lon_upper - map_lon_lower >= 359) 
                        %/ Then, use the function 'correct_bndry_at_dateline'
                        %/ to avoid "transverse lines" across the dateline
    %                     bndry_data_each = correct_bndry_at_dateline(bndry_data_each, map_lon_upper, map_lon_lower);
    %                 end
                    
                    %/ check if various colors have been assigned to each boundary
                    if size(color, 1)  ~= 1            bndry_color  = color(i,:);                    else  bndry_color  = color;                       end
                    if length(bndry_patch_mode) ~= 1   bndry_patch_mode_each = bndry_patch_mode(i);  else  bndry_patch_mode_each = bndry_patch_mode;   end
                    
                    if ~isempty(bndry_linewi)
                        if length(bndry_linewi) > 1
                            bndry_linewi_bc = bndry_linewi(i);
                        else
                            bndry_linewi_bc = bndry_linewi;
                        end
                    else
                        if size(linewi, 1) ~= 1            
                            bndry_linewi_bc = linewi(i,:);                   
                        else  
                            bndry_linewi_bc = linewi;                      
                        end
                    end
    
                    if map_lon_lower < 0 || map_lon_upper < 0
                        ind = find(bndry_data_each(:,1) > 180);
                        bndry_data_each(ind,1) = bndry_data_each(ind,1) - 360;  %/ convert to [-179 180] if ind is not empty.
                    elseif map_lon_lower >= 0 && map_lon_upper >= 0
                        ind = find(bndry_data_each(:,1) < 0);
                        bndry_data_each(ind,1) = bndry_data_each(ind,1) + 360;  %/ convert to [-179 180] if ind is not empty.
                    end

                    %/ Whether to fill color of the boundary (visualization)
                    if bndry_patch_mode_each
                        if isempty(patch_color)        
                            patch_color = color;
                            warning('''patch_color'' is missing. Set it to be the same as ''color'' by default.');             
                        end
                        
                        if min(bndry_data_each(:,1)) < map_lon_lower || max(bndry_data_each(:,1)) > map_lon_upper || ...
                           min(bndry_data_each(:,2)) < map_lat_lower || max(bndry_data_each(:,2)) > map_lat_upper
                            disp([min(bndry_data_each(:,1)), max(bndry_data_each(:,1))])
                            disp([min(bndry_data_each(:,2)), max(bndry_data_each(:,2))])
                            error('Enlarge your domain of plot, otherwise patch color won''t be drawn correctly!');
                        end
                        if size(patch_color, 1) ~= 1       
                            bndry_facecolor = patch_color(i,:);
                        else
                            bndry_facecolor = patch_color;  
                        end
                        
                        m_patch(bndry_data_each(:,1), bndry_data_each(:,2), bndry_facecolor,...
                                        'edgecolor', bndry_color, 'linewidth', bndry_linewi_bc*0.5, 'facealpha', patch_alpha);
    %                     hpatch
    %                     hpatch.FaceAlpha = patch_alpha;
    %                     alpha(patch_alpha) %/ this will affect all
                    else
                        m_line(bndry_data_each(:,1), bndry_data_each(:,2), 'marker', 'none', 'markersize', 0.001, 'markerfacecolor', 'none',...
                              'linest', '-', 'color', bndry_color, 'linewi', bndry_linewi_bc);
                    end
                    hold on;
                    if draw_bndry_onebyone %/ for identifying the basin name.
                        drawnow; pause(2);
                    end
                end
            else
                if map_lon_lower < 0 || map_lon_upper < 0
                    ind = find(bndry_data(:,1) > 180);
                    bndry_data(ind,1) = bndry_data(ind,1) - 360;  %/ convert to [-179 180] if ind is not empty.
                elseif map_lon_lower >= 0 && map_lon_upper >= 0
                    ind = find(bndry_data(:,1) < 0);
                    bndry_data(ind,1) = bndry_data(ind,1) + 360;  %/ convert to [-179 180] if ind is not empty.
                end

                if bndry_patch_mode
                    if isempty(patch_color) 
                        error('specify ''patch_color''!'); 
                    end
                    m_patch(bndry_data(:,1), bndry_data(:,2), patch_color);
                else
                    m_line(bndry_data(:,1),bndry_data(:,2), 'marker', 'none', 'markersize', 0.001, 'markerfacecolor', 'none',...
                              'linest', '-', 'color', color, 'linewi', linewi);
                end
                hold on;
            end
            % drawnow; pause(0.05); %/ don't put drawnow function in a for-loop, or it will get extremely slow!! 
            % if ~isempty(bndry_data)
            %     m_line(bndry_data(:,1), bndry_data(:,2), 'linewi',1.5,'color','k');
            % %     m_hatch(bndry_data(:,1), bndry_data(:,2), 'single',30, 5,'color', color);
            % %     drawnow; pause(0.05);
            % end
        end
    
        if draw_river      
    %         river_color = [42  109 181]./255;
            river_color = [30 160 244]./255;
            m_gshhs('fr2', 'linewidth', coast_wi, 'color', river_color, 'linewi', linewi*0.75); 
    %         m_gshhs('ir', 'linewidth', coast_wi, 'color', 'b');          
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(traj_data)
            %/ NOTE: Assume traj_data is in (noOftraj, trajtime, [x,y,z])
            %/       where variables is in (lon, lat, x)
            % tic
            [ntraj, ntrajtime, ~] = size(traj_data);
            fprintf('*** min(traj_data(:,:,3), [], ''all'') = %.2f ***\n', min(traj_data(:,:,3), [], 'all'));  
            fprintf('*** max(traj_data(:,:,3), [], ''all'') = %.2f ***\n', max(traj_data(:,:,3), [], 'all'));

            if isempty(traj_linewi)
                traj_linewi = linewi;
            end
            % NumWorkers = 40;
            % if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
            %     parpool('Threads', NumWorkers)  %/ much faster
            %     % parpool('Processes', NumWorkers_HPC)
            % end

            %/ Convert traj lon to [0, 360] for plotting trajectories across dateline 
            %/ (since FLEXPART by default writes the traj lon in [-179 180])
            if map_lon_upper > 180
                traj_lon = squeeze(traj_data(:,:,1));
                traj_lon(traj_lon < 0) = traj_lon(traj_lon < 0) + 360;
                traj_data(:,:,1) = traj_lon;  %/ Update
            end

            for n = 1:ntraj
                fprintf('*** Drawing trajs (%d/%d)... ***\n', n, ntraj);
                traj_data_each = squeeze(traj_data(n,:,:));

                %/ Broadcast vars
                traj_time_bc   = traj_time;
                traj_levels_bc = traj_levels;
                traj_colmap_bc = traj_colmap;

                %/ Avoid "transverse lines" across the dateline
                % traj_data_each = correct_bndry_at_dateline(traj_data_each, map_lon_upper, map_lon_lower);

                %/ Avoid "transverse lines" across the map domain by removing the overflown segment
                traj_data_each = rm_overflown_traj_seg('traj_data', traj_data_each, 'dim_lon', 1, 'dim_lat', 2, 'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper);

                nan_cnt_per_traj = 0;
                for tj = 1:ntrajtime-1
                    %/ First, mark the start point of the trajectory (in chronological order)
                    if (diff(traj_time_bc(1:2)) < 0 && tj == ntrajtime-1) ...  %/ if bwd traj is detected, mark the last  (i.e., t = -XX) point.
                        || (diff(traj_time_bc(1:2)) > 0 && tj == 1)                %/ if fwd traj is detected, mark the first (i.e., t = 0) point.
                        
                        m_line(traj_data_each(tj,1), traj_data_each(tj,2), 'marker','o', 'color', 'none', 'linewi',traj_linewi,...
                                 'linest','none','markersize',markersize*4, 'markerfacecolor', 'k');
                    end
                    hold on;
                    
                    %/ NOTE: Since we need two points to draw a line, here we
                    %/       average the intensity at two timesteps to define the color.
    
                    %/ Define the color based on the value of the 3rd variable
                    mean_value = mean(traj_data_each(tj:tj+1, 3));
                    
                    if isnan(mean_value)
                        nan_cnt_per_traj = nan_cnt_per_traj + 1;
    %                     warning('NaN detected in a traj (%d/%d). Skip drawing this traj segment!', nan_cnt_per_traj, ntrajtime-1);
                        continue;
                    else
                        ind = find(traj_levels_bc(1:end-1) < mean_value, 1, 'last');
                        if isempty(ind)   %/ If mean_value is lower than the lowest traj_levels_bc, ind will be empty. Then we replace it with 1 to use the minimum color.
                            ind = 1;  
                            % error('No index color for the value %f!', mean_value);   
                        end
                    end
    
                    m_line(traj_data_each(tj:tj+1,1), traj_data_each(tj:tj+1,2), 'marker','none', 'color', traj_colmap_bc(ind,:), 'linewi',traj_linewi,...
                                    'linest','-','markersize', markersize, 'markerfacecolor', 'none', 'clip', 'off');
                    hold on;
                end
    %             drawnow;
            end
            hold on;
            clim([min(traj_levels_bc) max(traj_levels_bc)]);
            colormap(gca, traj_colmap_bc)
    
            if cbar_mode
                cb = colorbar('eastoutside');
                cbar_YTick = traj_levels_bc(2:cbar_interval:end-1);
                cbar_YTickLabel = cbar_YTick;
    
                set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', fontsize) %/ cbar Ytick for diverging colormap
                set(get(cb,'Title'),'String', traj_unit, 'Fontsize', fontsize)
            end
            % set(cb, 'YTick', contflevels,'Fontsize', 12) %/ cbar Ytick for diverging colormap
            % cbarrow;
    %         drawnow; pause(0.05);
            % fprintf('Time cost in plotting trajs: %f \n', toc);
        end
        
        %%%%%%% IMPORTANT: Fixing the fig frame for a correct drawing of vectors! (do NOT remove it!) %%%%%%%
        m_grid('xtick', [], 'ytick', [], 'linewi', 0.001, 'linest','none', 'gridcolor', 'none');
    %     drawnow; pause(0.05)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isequal(map_proj, 'ortho')   clip = 'off';   else  clip = 'on';  end
        if ~isempty(Udata) && ~isempty(Vdata)
    
            %/ [IMPORTANT] To avoid cumulative vectors at the edges, 
            %/ we manually set vectors outside the plot to NaNs.
            cond_x = (x >= map_lon_lower & x <= map_lon_upper); 
            cond_y = (y >= map_lat_lower & y <= map_lat_upper); 
            cond = (cond_x & cond_y);
            Udata(~cond) = nan;
            Vdata(~cond) = nan;
            x(~cond) = nan;
            y(~cond) = nan;
    
            if ~isempty(vector_levels)  %/ then we draw different vector colors following vector_levels.   
                if length(vector_levels) ~= size(vector_color, 1)
                    error('The length of vector_levels has to be consistent with the rows of vector_color!');
                end
                
                UVmag = sqrt(Udata.^2 + Vdata.^2);
                
                for k = 1:length(vector_levels)
                    %/ a faster way (should be)
                    if k == length(vector_levels)
                        cond = (vector_levels(k) <= UVmag); %/ UVmag >= upper bound.
                    else
                        cond = (vector_levels(k) <= UVmag & UVmag < vector_levels(k+1)); %/ eg. 0 <= UVmag < 20.
                    end
                    x_bc = x;          x_bc(~cond) = nan;
                    y_bc = y;          y_bc(~cond) = nan;  
                    Udata_bc = Udata;  Udata_bc(~cond) = nan;  
                    Vdata_bc = Vdata;  Vdata_bc(~cond) = nan;  
                    
                    m_vec(vecscale, x_bc(2:vec_step_lat:end-1, 2:vec_step_lon:end-1),...
                                    y_bc(2:vec_step_lat:end-1, 2:vec_step_lon:end-1),...
                                    Udata_bc(2:vec_step_lat:end-1, 2:vec_step_lon:end-1)*vecscale2,...
                                    Vdata_bc(2:vec_step_lat:end-1, 2:vec_step_lon:end-1)*vecscale2,...
                                    vector_color(k,:), 'shaftwidth',shaftwidth ,'headlength',headlength,...
                                    'EdgeColor', vector_edgecolor, 'linewidth', 0.5,  'clip', clip);  %/ Do NOT use 'edgeclip' but 'clip'
    %                 hold on; drawnow;
                end
            else
                m_vec(vecscale, x(2:vec_step_lat:end-1, 2:vec_step_lon:end-1),...
                                y(2:vec_step_lat:end-1, 2:vec_step_lon:end-1),...
                                Udata(2:vec_step_lat:end-1, 2:vec_step_lon:end-1)*vecscale2,...
                                Vdata(2:vec_step_lat:end-1, 2:vec_step_lon:end-1)*vecscale2,...
                                vector_color, 'shaftwidth',shaftwidth ,'headlength',headlength,...
                                'EdgeColor', vector_edgecolor, 'linewidth', 0.5, 'clip', clip);    %/ Do NOT use 'edgeclip' but 'clip'
                            
    %             m_vec(vecscale, x(2:vec_step_lat:end-1, 2:vec_step_lon:end-1),...
    %                             y(2:vec_step_lat:end-1, 2:vec_step_lon:end-1),...
    %                             Udata(2:vec_step_lat:end-1, 2:vec_step_lon:end-1)*vecscale2,...
    %                             Vdata(2:vec_step_lat:end-1, 2:vec_step_lon:end-1)*vecscale2,...
    %                             vector_color, 'shaftwidth',shaftwidth ,'headlength',headlength,...
    %                             'edgeclip', 'on', 'curvature', curvature);    %/ Do NOT use 'edgeclip' but 'clip'
    %                         hpv6 = m_vec(100, -124, 53, 30, -40, 'r', 'curvature',30,'centered','no')
            end
            fprintf('*** Max abs(Udata) = %.2g, Max abs(Vdata) = %.2g ***\n', max(abs(Udata), [], 'all', 'omitnan'),...
                                                                              max(abs(Vdata), [], 'all', 'omitnan'));
            
            if show_refvec
                if isempty(vec_mag_ref)
                    error('''vec_mag_ref'' is undefined when show_refvec is on!');
                end
                lon_range = map_lon_upper - map_lon_lower;
                lat_range = map_lat_upper - map_lat_lower;
    
                if isempty(vec_ref_lat_shift)   vec_ref_lat_shift = 0.1;   end
    
                vec_ref_lon = map_lon_upper - lon_range*0.15;
                vec_ref_lat = map_lat_lower - lat_range*vec_ref_lat_shift;
    
                if vec_ref_lat < -90     warning('The reference vector is outside the latitude of -90, it thus may not be shown on the map!!');    end
    
                if isequal(map_proj, 'robin')
                    %/ Draw vectors outside the map
                    %/ IMPORTANT: if the reference vector is outside the latitude of -90, it may not be shown on the map!!
                    [hp, ht] = m_vec(vecscale, vec_ref_lon, vec_ref_lat, vec_mag_ref*vecscale2, 0,...
                                 vector_color, 'shaftwidth', shaftwidth, 'headlength', headlength, 'edgeclip', 'off', 'key', vec_lbs); %/ Do NOT turn 'edgeclip' on! It doesn't help, and may cause bug for map_prof = 'ortho'.
                elseif isequal(map_proj, 'Miller Cylindrical')
                    [hp, ht] = m_vec(vecscale, vec_ref_lon, vec_ref_lat, vec_mag_ref*vecscale2, 0,...
                                 vector_color, 'shaftwidth', shaftwidth, 'headlength', headlength, 'edgeclip', 'off', 'key', vec_lbs); %/ Do NOT turn 'edgeclip' on! It doesn't help, and may cause bug for map_prof = 'ortho'.
                end
    
                if isempty(vec_ref_fontsize)        vec_ref_fontsize = fontsize;     end
                set(ht,'FontSize',vec_ref_fontsize, 'color', vector_color, 'HorizontalAlignment', 'center');
                set(hp,'EdgeColor', vector_edgecolor);                             %/ have to set the EdgeColor using set() if using 'key' in m_vec!
            end
            % fprintf('Time cost in vector: %f \n', toc);
            
            %/ M_STREAMLINE (useless, many bugs and cannot use in all projection) 
    %         ind_lon = find(map_lon_lower+1 <= x(1,:) &  x(1,:) <= map_lon_upper-1);  size(ind_lon)
    %         ind_lat = find(map_lat_lower+1 <= y(:,1) &  y(:,1) <= map_lat_upper-1);  size(ind_lat)
    %          
    %         m_streamline(    x(ind_lat, ind_lon),     y(ind_lat, ind_lon),...
    %                      Udata(ind_lat, ind_lon), Vdata(ind_lat, ind_lon)); 
    
        end
        
        %%%%%%%%%%%%%%%%% Mask mountains (using 2nd colmap!) %%%%%%%%%%%%%%%%%%
        if mask_mountains
            %/ output an 1-deg averaged topo (from 5-min high topo data)
            lon_grids = map_lon_lower:1:map_lon_upper;
            lat_grids = map_lat_lower:1:map_lat_upper;
            [topo_lon_2D, topo_lat_2D] = meshgrid(lon_grids,lat_grids);
            topo = interp_topo('lon_grids', lon_grids, 'lat_grids', lat_grids);
                          
            %/ NOTE: m_elev outputs 1 deg topo; m_tbase outputs much higher res. 
            %/       Prefer the coarser one to do masking.
            
    %         [topo, topo_lon_2D, topo_lat_2D] = m_elev([map_lon_lower, map_lon_upper, map_lat_lower, map_lat_upper]); %REGION =[west east south north];
    %         topo     = topo';               %/ make sure in (lon,lat) dim
            topo_lon = topo_lon_2D(1,:)';
            topo_lat = topo_lat_2D(:,1);
            topo(topo <= plateau_hgt) = nan;
            topo(~isnan(topo))        = plateau_hgt;  % so now it's only a single-value field.
            
            contf_levels_topo = [0, plateau_hgt];
            
            colmap2 = [255 255  255;
                        40 14 14;]./255;  %/ masked by grey boxes
            hold on;
            
            if ~isempty(ax_panel)
                ax_topo = ax_panel;
            else
                ax_topo = axes;  %/ axes is to create a new axes!
            end
            set(ax_topo, 'Visible', 'off')                            %/ IMPORTANT: hide the 2nd axes 
            
            %/ NOTE: 'Since pcolor has offsets and edge problem, here we solve the offset problem by shifting by a half grid.\nFor more info, see https://www.mathworks.com/matlabcentral/fileexchange/50706-offsets-and-missing-data-via-pcolor-and-surf'
            res_lon = abs(unique(round(diff(topo_lon), 3)));  %/ round up to 3 decimal points to avoid numerical errors
            res_lat = unique(round(diff(topo_lat), 3));  %/ do not take absolute value, seems to be correct. Can check afterwards.
            
            m_pcolor(topo_lon_2D-res_lon/2, topo_lat_2D-res_lat/2, topo', 'edgecolor','none');
            shading flat
            
    %         m_contourf(topo_lon_2D-res_lon/2, topo_lat_2D-res_lat/2, topo', 'edgecolor','none');
    %         shading interp
            
            colormap(ax_topo, colmap2);                                  %/ IMPORTANT: set colormap for the 2nd axes
            
            if cbar_mode   %/ IMPORTANT: This is to keep cbar_mode the same for two axes. Otherwise the fig will shift.
                cb_topo = colorbar(ax_topo, 'Location', cbar_location);      
                cbar_YTick = contf_levels_topo;
                cbar_YTickLabel = cbar_YTick;
                set(cb_topo, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', fontsize) %/ cbar Ytick for diverging colormap
                set(get(cb_topo,'Title'),'String', '', 'Fontsize', fontsize)
            end
            clim(ax_topo,[min(contf_levels_topo) max(contf_levels_topo)]);
            
            if ~isempty(ax_panel) 
                linkaxes([ax_panel ax_topo]);                            %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
            else
                linkaxes([gca ax_topo]);                            %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vector2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(U2data) && ~isempty(V2data)
            
            %/ Copy the settings from Vector 1 if not given
            if isempty(vec2_step_lon)      vec2_step_lon = vec_step_lon;  end
            if isempty(vec2_step_lat)      vec2_step_lat = vec_step_lat;  end
            if isempty(vector2_levels)     vector2_levels = vector_levels;  end
            if isempty(vector2_color)      vector2_color = vector_color;  end
            if isempty(vector2_edgecolor)  vector2_edgecolor = vector_edgecolor; end
            if isempty(vec2scale)          vec2scale = vecscale; end
            if isempty(vec2scale2)         vec2scale2 = vecscale2; end
            if isempty(shaft2width)        shaft2width = shaftwidth; end
            if isempty(head2length)        head2length = headlength; end
            if isempty(vec2_lbs)           vec2_lbs = vec_lbs; end
            if isempty(vec2_mag_ref)       vec2_mag_ref = vec_mag_ref; end
            if isempty(vec2_ref_fontsize)  vec2_ref_fontsize = vec_ref_fontsize; end
            if isempty(vec2_ref_lat_shift) vec2_ref_lat_shift = 0.3;   end
    
            %/ [IMPORTANT] To avoid cumulative vectors at the edges, 
            %/ we manually set vectors outside the plot to NaNs.
            cond_x = (x2 >= map_lon_lower & x2 <= map_lon_upper); 
            cond_y = (y2 >= map_lat_lower & y2 <= map_lat_upper); 
            cond = (cond_x & cond_y);
            U2data(~cond) = nan;
            V2data(~cond) = nan;
            x2(~cond) = nan;
            y2(~cond) = nan;
    
            if ~isempty(vector2_levels)  %/ then we draw different vector colors following vector_levels.   
                if length(vector2_levels) ~= size(vector2_color, 1)
                    error('The length of vector2_levels has to be consistent with the rows of vector2_color!');
                end
                UV2mag = sqrt(U2data.^2 + V2data.^2);
                
                for k = 1:length(vector2_levels)
                    
                    %/ a faster way (should be)
                    if k == length(vector2_levels)
                        cond = (vector2_levels(k) <= UV2mag); %/ UVmag >= upper bound.
                    else
                        cond = (vector2_levels(k) <= UV2mag & UV2mag < vector2_levels(k+1)); %/ eg. 0 <= UVmag < 20.
                    end
                    x2_bc = x2;          x2_bc(~cond) = nan;
                    y2_bc = y2;          y2_bc(~cond) = nan;  
                    Udata2_bc = U2data;  Udata2_bc(~cond) = nan;  
                    Vdata2_bc = V2data;  Vdata2_bc(~cond) = nan;  
                    
                    m_vec(vec2scale, x2_bc(2:vec2_step_lat:end-1, 2:vec2_step_lon:end-1),...
                                    y2_bc(2:vec2_step_lat:end-1, 2:vec2_step_lon:end-1),...
                                    Udata2_bc(2:vec2_step_lat:end-1, 2:vec2_step_lon:end-1)*vec2scale2,...
                                    Vdata2_bc(2:vec2_step_lat:end-1, 2:vec2_step_lon:end-1)*vec2scale2,...
                                    vector2_color(k,:), 'shaftwidth',shaft2width ,'headlength',head2length,...
                                    'EdgeColor', vector2_edgecolor, 'linewidth', 0.5, 'clip', clip);  %/ Do NOT use 'edgeclip' but 'clip'
    %                 hold on; drawnow;
                end
            else
                m_vec(vec2scale, x2(2:vec2_step_lat:end-1, 2:vec2_step_lon:end-1),...
                                y2(2:vec2_step_lat:end-1, 2:vec2_step_lon:end-1),...
                                U2data(2:vec2_step_lat:end-1, 2:vec2_step_lon:end-1)*vec2scale2,...
                                V2data(2:vec2_step_lat:end-1, 2:vec2_step_lon:end-1)*vec2scale2,...
                                vector2_color, 'shaftwidth',shaft2width ,'headlength',head2length,...
                                'EdgeColor', vector2_edgecolor, 'linewidth', 0.5, 'clip', clip);    %/ Do NOT use 'edgeclip' but 'clip'
            end
            fprintf('*** Max abs(U2data) = %.2g, Max abs(V2data) = %.2g ***\n', max(abs(U2data), [], 'all', 'omitnan'),...
                                                                                max(abs(V2data), [], 'all', 'omitnan'));
            
            if show_refvec
                if isempty(vec2_mag_ref)
                    error('''vec2_mag_ref'' is undefined when show_refvec is on!');
                end
                lon_range = map_lon_upper - map_lon_lower;
                lat_range = map_lat_upper - map_lat_lower;
    
                
    
                vec2_ref_lon = map_lon_upper - lon_range*0.15;
                vec2_ref_lat = map_lat_lower - lat_range*vec2_ref_lat_shift;
    
                if vec2_ref_lat < -90     warning('The reference vector2 is outside the latitude of -90, it thus may not be shown on the map!!');    end
    
                if isequal(map_proj, 'robin')
                    %/ Draw vectors outside the map
                    %/ IMPORTANT: if the reference vector is outside the latitude of -90, it may not be shown on the map!!
                    [hp, ht] = m_vec(vec2scale, vec2_ref_lon, vec2_ref_lat, vec2_mag_ref*vec2scale2, 0,...
                                 vector2_color, 'shaftwidth', shaft2width, 'headlength', head2length, 'edgeclip', 'off', 'key', vec2_lbs);
                elseif isequal(map_proj, 'Miller Cylindrical')
                    [hp, ht] = m_vec(vec2scale, vec2_ref_lon, vec2_ref_lat, vec2_mag_ref*vec2scale2, 0,...
                                 vector2_color, 'shaftwidth', shaft2width, 'headlength', head2length, 'edgeclip', 'off', 'key', vec2_lbs);
                end
    
                set(ht,'FontSize',  vec2_ref_fontsize, 'color', vector2_color, 'HorizontalAlignment', 'center');
                set(hp,'EdgeColor', vector2_edgecolor);                             %/ have to set the EdgeColor using set() if using 'key' in m_vec!
            end
            % fprintf('Time cost in vector2: %f \n', toc);
        end
        
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Point2 (serve as marking of gauges/grids) %%%
        if ~isempty(point2_data)
            % tic
            if map_lon_lower < 0        %/ if plotting in [-179, 180], we convert point_data's lon to fit that range.
                a = point2_data(:,1);
                a(a > 180) = a(a > 180) - 360;
                point2_data(:,1) = a;
            end
            if isempty(marker2facecolor)  marker2facecolor = 'none';  end
            if isempty(marker2edgecolor)  marker2edgecolor = 'none';  end
    
            if size(marker2facecolor, 1) > 1 
                if size(marker2edgecolor, 1) == 1
                    marker2edgecolor = repmat(marker2edgecolor, size(marker2facecolor,1), 1); %/ then replicate it
                end
                %/ Loop to draw points with different facecolors
                for k = 1:size(point2_data, 1)
                    m_line(point2_data(k,1), point2_data(k,2), 'marker', marker2, 'markersize', marker2size, 'markerfacecolor', marker2facecolor(k,:),...
                            'linest', 'none', 'color', marker2edgecolor(k,:), 'linewi', linewi, 'clip', 'point');
                    hold on;
                end
            else
                m_line(point2_data(:,1), point2_data(:,2), 'marker', marker2, 'markersize', marker2size, 'markerfacecolor', marker2facecolor,...
                        'linest', 'none', 'color', marker2edgecolor, 'linewi', linewi, 'clip', 'point');
                
            end
    %         alpha(h, 0.5);
            % fprintf('Time cost in drawing points2: %f \n', toc);
        end
        %%%%%%%%%%%%%%%%% m_grid (ticks, box linewidth, etc.) %%%%%%%%%%%%%%%%%
        fprintf('grid_mode = %d\n', grid_mode)
        if grid_mode == -4               %/ global
            if isempty(grid_linewi)  grid_linewi = 2.5;   end
            m_grid('xtick', [], 'ytick', -90:30:90, 'LineWidth', grid_linewi, ...
                   'linest','none','gridcolor', [.1 .1 .1], 'xticklabels',[], 'tickdir', 'in', 'backcolor', backcolor, 'fontsize',fontsize);   
    
        elseif grid_mode == -3                %/ not drawing lon lat but EQ
            if isempty(grid_linewi)  grid_linewi = 2.5;   end
            m_grid('xtick', [], 'ytick', -90:30:90, 'LineWidth', grid_linewi, ...
                   'linest','-','gridcolor', [.1 .1 .1], 'xticklabels',[],'yticklabels',[], 'tickdir', 'in', 'backcolor', backcolor, 'fontsize',fontsize);   
    
        elseif grid_mode == -2             %/ not drawing lon lat and grids
            if isempty(grid_linewi)  grid_linewi = 1;   end
            m_grid('xtick', [], 'ytick', [], 'linewidth', grid_linewi, 'linest','none', 'color', 'none', 'box', 'off', 'gridcolor', 'none', 'tickdir','none',...
                   'xaxisloc', xaxisloc, 'backcolor', backcolor, 'fontsize',fontsize);
               
        elseif grid_mode == -1             %/ not drawing lon lat
            if isempty(grid_linewi)  grid_linewi = 2.5;   end
            m_grid('xtick', [], 'ytick', [],...
                   'linewidth',grid_linewi,'linest','none','gridcolor', 'w', 'tickdir','in',  'backcolor', backcolor, 'fontsize',fontsize);    
        
        elseif grid_mode == 0              %/ not drawing lon
            if isempty(grid_linewi)  grid_linewi = 2;   end
            m_grid('xtick', [], 'ytick', -60:30:60, 'linewidth', grid_linewi, 'linest','none', 'gridcolor', 'k', 'tickdir','in',...
                   'xaxisloc', xaxisloc, 'backcolor', backcolor, 'fontsize',fontsize);
        
        elseif grid_mode == 1         %/ normal
            xtick = -180:30:360;    ytick = -90:20:90;
            if isempty(grid_linewi)  grid_linewi = 2;   end
            m_grid('xtick', xtick, 'ytick', ytick, 'linewidth', grid_linewi, 'linest','none','gridcolor', 'k', 'tickdir','in',  'backcolor', backcolor, 'fontsize',fontsize);
    
        elseif grid_mode == 2         %/ draw grids on each 2 degree
            xtick = -180:2:360;    ytick = -90:1:90;
            if isempty(grid_linewi)  grid_linewi = 2.5;   end
            m_grid('xtick', xtick, 'ytick', ytick,... % 'yticklabels', xtick(x_st_ind:5:x_ed_ind), 'yticklabels',  ytick(y_st_ind:5:y_ed_ind),...
                   'linewidth', grid_linewi,'linest',':','gridcolor', 'w', 'tickdir','in', 'xaxisloc', xaxisloc,  'backcolor', backcolor, 'fontsize',fontsize);
    
        elseif grid_mode == 3         %/ for regional study
            xtick = -180:30:360;    ytick = -90:10:90;
            if isempty(grid_linewi)  grid_linewi = 2.5;   end
            m_grid('xtick', xtick, 'ytick', ytick,...
                   'linewidth',grid_linewi,'linest','none','gridcolor', 'w', 'tickdir','in',  'backcolor', backcolor, 'fontsize',fontsize);    
        
        elseif grid_mode == 4         %/ for local study
            xtick = -180:2:360;    ytick = -90:2:90;
            if isempty(grid_linewi)  grid_linewi = 2;   end
            m_grid('xtick', xtick, 'ytick', ytick,...
                   'linewidth', grid_linewi,'linest','none','gridcolor', 'w', 'tickdir','in',  'backcolor', backcolor, 'fontsize',fontsize);  
        
        elseif grid_mode == 5         %/ for regional study 2
            xtick = -180:20:360;    ytick = -85:15:90;
            if isempty(grid_linewi)  grid_linewi = 2.5;   end
            m_grid('xtick', xtick, 'ytick', ytick,...
                   'linewidth', grid_linewi, 'linest','none','gridcolor', 'w', 'tickdir','in',  'backcolor', backcolor, 'fontsize',fontsize);  
               
        elseif grid_mode == 6         %/ for regional study 3
            xtick = -180:10:360;    ytick = -90:10:90;
            if isempty(grid_linewi)  grid_linewi = 2.5;   end
            m_grid('xtick', xtick, 'ytick', ytick,...
                   'linewidth', grid_linewi, 'linest','none','gridcolor', 'w', 'tickdir','in',  'backcolor', backcolor, 'fontsize',fontsize);  

        elseif grid_mode == 7         %/ for tropical study 
            xtick = -180:30:360;    ytick = -80:20:80;
            if isempty(grid_linewi)  grid_linewi = 1.5;   end
            m_grid('xtick', xtick, 'ytick', ytick,...
                   'linewidth', grid_linewi, 'linest','none','gridcolor', 'w', 'tickdir','in',  'backcolor', backcolor, 'fontsize',fontsize);  
        
        elseif grid_mode == 8         %/ for tropical study 
            xtick = -180:60:360;    ytick = -80:20:80;
            if isempty(grid_linewi)  grid_linewi = 1.5;   end
            m_grid('xtick', xtick, 'ytick', ytick,...
                   'linewidth', grid_linewi, 'linest','none','gridcolor', 'w', 'tickdir','in',  'backcolor', backcolor, 'fontsize',fontsize);  
        
        end
        
        if draw_country     m_gshhs('fb1', 'linewidth', coast_wi*0.5, 'color', [.1 .1 .1]);     end
    
        %%%%%%%%%%%%%%%%% remove 2nd colorbar (if any) %%%%%%%%%%%%%%%%%%%%%%%%
        if mask_mountains && cbar_mode 
            set(cb_topo, 'Visible', 'off');     %/ IMPORTANT: remove it AFTER setting m_grid!
    %         drawnow; pause(1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Text on the map %%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(text_data)
            if size(text_data, 1) == 1          %/ e.g., a box at the bottom left
                verticalalignment = 'bottom';
                horizontalalignment = 'left';
            else
                verticalalignment = 'middle';
                horizontalalignment = 'center';
            end
    %         text_fontsize
            for i = 1:size(text_data, 1)
                m_text(text_data{i,1}, text_data{i,2}, text_data{i,3}, 'color', [text_data{i,4}], 'fontweight', 'normal',...
                        'verticalalignment', verticalalignment, 'horizontalalignment', horizontalalignment, 'fontsize', text_fontsize,...
                        'backgroundcolor', text_backgroundcolor, 'edgecolor', text_edgecolor);
            end
    %         drawnow; pause(0.05)  %/ don't put drawnow function in a for-loop, or it will get extremely slow!! 
        end
        
        if isempty(title_fontsize) title_fontsize = fontsize*0.8;  end
        titlename = strrep(titlename, '_', ' ');

        if isempty(title_pos)
            title_pos = [0.12,1.01,1.0,0];
        end
        
        %/ Use ylim to fine-tune the y-pos of the title
        % y_lim = ylim;
        % title_pos(2) = y_lim(2) + 0.1;
        % title_pos
        % x_lim = xlim;
        % title_pos(1) = mean(x_lim);

        %/ Annotation is indep. of axes -> not distorting the aspect ratio.
        annotation( 'textbox', 'String', titlename, 'Color', 'k', ...
                    'FontSize', title_fontsize, 'Units', 'normalized', 'EdgeColor', 'none', ...
                    'Position', title_pos)

        hold on;
        drawnow; 
    end

    if ~isempty(savepath)
        % tic
        FigName_underscore = strrep(savepath, ' ', '_');
        final_FigName = char(strcat(FigName_underscore, str_relief, '.', fig_fmt));
        disp(final_FigName)
        
        if isequal(fig_fmt, 'png')
            if isempty(png_dpi)   png_dpi = 300;   end
            if trans_bg
                export_fig(final_FigName,['-r',num2str(png_dpi)],'-png', '-opengl', '-nocrop', '-transparent'); %, '-nocrop');
            else
                export_fig(final_FigName,['-r',num2str(png_dpi)],'-png', '-opengl', '-nocrop');
            end
        elseif isequal(fig_fmt, 'svg')
            export_fig(final_FigName,'-svg','-painters', '-nocrop', '-transparent'); %/ This svg output can output transparent object with transparent background!
        else
            fig = get(0, 'CurrentFigure');
            hasTransparency = ~isempty(findall(fig,'-property','FaceAlpha','-and','-not','FaceAlpha',1)); %/ code script borrowed from 'export_fig'!
            if hasTransparency
    %             set(fig, 'InvertHardCopy', 'off'); %/ transparent background
                set(gcf, 'Color', 'None');
                set(gca, 'Color', 'None');
                print(gcf, '-dpdf',final_FigName, '-fillpage'); 
            else
                %/ Cropping often causes unnecessary inconsistency in fontsize. Disable it.
                if trans_bg
                    export_fig(final_FigName,'-pdf','-painters', '-nocrop', '-transparent'); %/ make background transparent (do not auto-crop; it may cut out part of the figure)
                else
                    export_fig(final_FigName,'-pdf','-painters', '-c[inf, inf, inf, inf]'); %, '-nocrop'); 
                end
            end
        end
        fprintf('!!! Plot is saved: %s !!! \n', string(final_FigName));
        % fprintf('Time cost in saving figure: %f \n', toc);
    end
    set(gcf, 'Color', 'w');
    set(gca, 'Color', 'w');

% 
% 
% %/    
%     S = [];
%     S.contf_data_events = contf_data_events; %/ for other use, e.g., ttest2_sig_fn
%     S.cont_data_events  = cont_data_events;  %/ for other use, e.g., ttest2_sig_fn
%     S.Udata_events      = Udata_events;      %/ for other use, e.g., ttest2_sig_fn
%     S.Vdata_events      = Vdata_events;      %/ for other use, e.g., ttest2_sig_fn
%     S.contf_data        = contf_data;
%     S.contf_lon         = contf_lon;
%     S.contf_lat         = contf_lat;
%     S.contf_unit        = contf_unit;
%     S.point_data        = point_data;
%     S.cont_data         = cont_data;
%     S.cont_lon          = cont_lon;
%     S.cont_lat          = cont_lat;
%     S.cont_unit         = cont_unit;
%     S.Udata             = Udata;
%     S.Vdata             = Vdata;
%     S.U2data            = U2data;
%     S.V2data            = V2data;
%     S.uv_lon            = uv_lon;
%     S.uv_lat            = uv_lat;
%     S.Rossby_westerly   = Rossby_westerly;
%     S.Kelvin_easterly   = Kelvin_easterly;
%     S.RK_ratio          = RK_ratio;
%     S.contf_AWM         = contf_AWM;
%     S.AWM_lon_extent = compute_AWM_lon_extent;
%     S.AWM_lat_extent = compute_AWM_lat_extent;
end
 
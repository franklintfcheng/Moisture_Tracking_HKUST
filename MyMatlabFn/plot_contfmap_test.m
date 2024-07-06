%%
function plot_contfmap_test(varargin)

% create a set of valid parameters and their default value
pnames = {'contf_data', 'contf_lon', 'contf_lat', 'contf_levels', 'contf_unit', 'colmap', 'cbar_interval', 'pcolor_mode',...
          'cont_data', 'cont_data_raw', 'cont_lon', 'cont_lat', 'cont_levels', 'cont_colmap', 'cont_linewi', 'cont_labelsize',...
          'point_data', 'marker', 'markersize', 'markerfacecolor', 'linewi', 'color', 'bndry_data', 'text_data', 'titlename', 'savepath', 'fig_fmt',...
          'Udata', 'Vdata', 'uv_lon', 'uv_lat', 'vec_step_lon', 'vec_step_lat', 'vector_color', 'vector_edgecolor', 'vecscale', 'vecscale2', 'shaftwidth', 'headlength', 'vec_lbs', 'vec_mag_ref',...
          'glb_data_mode', 'glb_plateau_mode', 'plateau_hgt', 'plateau_col',...
          'map_proj', 'map_lon_lower', 'map_lon_upper', 'map_lat_lower', 'map_lat_upper', 'coast_col', 'coast_wi', 'fontsize', 'create_fig', 'grid_mode',...
          'cbar_mode', 'cbar_YTick', 'cbar_YTickLabel', 'xaxisloc', 'backcolor', 'cbar_location', 'draw_cbar_only',...
          'traj_data', 'traj_levels', 'c_traj', 'traj_unit', 'hatch_data', 'hatch_lon', 'hatch_lat', 'hatch_thres', 'hatch_mode', 'hatch_linewi', 'hatch_density'};

dflts  = {[] [] [] [] [] [] [] 0 ...
          [] [] [] [] [] [] 1  1 ...
          [] 'o' 1 'r' 0.0001 [] [] [] [] [] [] ...
          [] [] [] [] [] [] [] 'w' [] [] [] [] [] [] ...
          [] [] [] [] ...
          [] [] [] [] [] [] 1.5 14 1 1 ...
          1 [] [] 'bottom' 'none' [] 0 ...
          [] [] [] [] [] [] [] [] 1 1.3 8};

%/ parse function arguments
[ contf_data, contf_lon, contf_lat, contf_levels, contf_unit, colmap, cbar_interval, pcolor_mode,...
  cont_data, cont_data_raw, cont_lon,  cont_lat,  cont_levels,  cont_colmap, cont_linewi, cont_labelsize,...
  point_data, marker, markersize, markerfacecolor, linewi, color, bndry_data, text_data, titlename, savepath, fig_fmt,...
  Udata, Vdata, uv_lon, uv_lat, vec_step_lon, vec_step_lat, vector_color, vector_edgecolor, vecscale, vecscale2, shaftwidth, headlength, vec_lbs, vec_mag_ref,...
  glb_data_mode, glb_plateau_mode, plateau_hgt, plateau_col,...
  map_proj, map_lon_lower, map_lon_upper, map_lat_lower, map_lat_upper, coast_col, coast_wi, fontsize, create_fig, grid_mode,...
  cbar_mode, cbar_YTick, cbar_YTickLabel, xaxisloc, backcolor, cbar_location, draw_cbar_only,...
  traj_data, traj_levels, c_traj, traj_unit, hatch_data, hatch_lon, hatch_lat, hatch_thres, hatch_mode, hatch_linewi, hatch_density] ...
             = internal.stats.parseArgs(pnames, dflts, varargin{:});

%===============

if draw_cbar_only
    fprintf('*** Drawing cbar only... *** \n');
    figure
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf, 'color','w');

    caxis([min(contf_levels) max(contf_levels)]);
    colormap(gca, colmap)

    if isempty(cbar_location)      cbar_location = 'eastoutside';    end 

    cb = colorbar(cbar_location);
    if isempty(cbar_YTick) && isempty(cbar_YTickLabel)
        cbar_YTick = contf_levels(2:cbar_interval:end-1);
        cbar_YTickLabel = cbar_YTick;
    end


    axis off
    if ismember(cbar_location, {'southoutside'})
        set(cb, 'position', [.1 .3 .4 .05]);    %/ [xposition yposition width height]
        set(cb, 'YAxisLocation','right')
        cb_fontsize = 20;
    elseif ismember(cbar_location, {'eastoutside'})
        set(cb, 'position', [.1 .3 .04 .4]);    %/ [xposition yposition width height]
        set(cb, 'YAxisLocation','right')
        cb_fontsize = 30;
    end
    set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', cb_fontsize) %/ cbar Ytick for diverging colormap
    set(get(cb,'Title'),'String',contf_unit, 'Fontsize', cb_fontsize)
    drawnow; pause(0.05);
    if ~isempty(savepath)
        savepath = strcat(savepath, '_cbar');   %/ automatically add the suffix to distinguish the cbar-only figure
    end
else
    if isempty(map_lon_lower) && isempty(map_lon_upper)  && isempty(map_lat_lower) && isempty(map_lat_upper) 
        map_lon_lower = min(contf_lon);
        map_lon_upper = max(contf_lon);
        map_lat_lower = min(contf_lat);
        map_lat_upper = max(contf_lat);
    end

    %/ transpose to column vector
    if size(contf_lon, 1) == 1                    contf_lon = contf_lon';   end
    if size(contf_lat, 1) == 1                    contf_lat = contf_lat';   end
    if size(contf_data, 1) == length(contf_lon)   contf_data = contf_data'; end

    if size(cont_lon, 1) == 1                     cont_lon = cont_lon';     end
    if size(cont_lat, 1) == 1                     cont_lat = cont_lat';     end
    if size(cont_data, 1) == length(cont_lon)     cont_data = cont_data';   end
    if size(cont_data_raw, 1) == length(cont_lon) cont_data_raw = cont_data_raw';   end

    if size(uv_lon, 1) == 1               uv_lon = uv_lon'; end
    if size(uv_lat, 1) == 1               uv_lat = uv_lat'; end
    if size(Udata, 1) == length(uv_lon)   Udata = Udata'; end
    if size(Vdata, 1) == length(uv_lon)   Vdata = Vdata'; end

    %/ check if it is indicated as global data
    if glb_data_mode
        if ~isempty(contf_data)     [contf_data, contf_lon_2D, contf_lat_2D]     = conv_to_circular_data(contf_data,    contf_lon,  contf_lat); end
        if ~isempty(cont_data)      [ cont_data,  cont_lon_2D,  cont_lat_2D]     = conv_to_circular_data(cont_data,     cont_lon,   cont_lat);  end
        if ~isempty(cont_data_raw)  [ cont_data_raw,  cont_lon_2D,  cont_lat_2D] = conv_to_circular_data(cont_data_raw, cont_lon,   cont_lat);  end
        if ~isempty(Udata)          [Udata, x, y] = conv_to_circular_data(Udata, uv_lon, uv_lat); end
        if ~isempty(Vdata)          [Vdata, x, y] = conv_to_circular_data(Vdata, uv_lon, uv_lat); end

    else
        if ~isempty(contf_data)                 [contf_lon_2D, contf_lat_2D] = meshgrid(contf_lon, contf_lat);  end
        if ~isempty(cont_data)                  [cont_lon_2D, cont_lat_2D]   = meshgrid(cont_lon, cont_lat);    end
        if ~isempty(Udata) || ~isempty(Vdata)   [x, y]                       = meshgrid(uv_lon, uv_lat);        end
    end

    if isempty(cbar_interval) cbar_interval = 2;  end

    if create_fig
        figure
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf, 'color','w');
    end

    if isempty(map_proj)
        if glb_data_mode
            map_proj = 'robin'; 
        else
            map_proj = 'Miller Cylindrical'; 
        end
    end

    m_proj(map_proj,'longitudes',[map_lon_lower map_lon_upper], 'latitudes', [map_lat_lower map_lat_upper]);
    hold on;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Contourf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(contf_data)
        tic
        if pcolor_mode
            %/ https://www.mathworks.com/matlabcentral/fileexchange/50706-offsets-and-missing-data-via-pcolor-and-surf'
            warning(sprintf('Since pcolor has offsets and edge problem, here we solve the offset problem by shifting by a half grid.\nFor more info, see https://www.mathworks.com/matlabcentral/fileexchange/50706-offsets-and-missing-data-via-pcolor-and-surf'));
            res_lon = abs(unique(round(diff(contf_lon), 3)));  %/ round up to 3 decimal points to avoid numerical errors
            res_lat = abs(unique(round(diff(contf_lat), 3)));

            m_pcolor(contf_lon_2D-res_lon/2, contf_lat_2D-res_lat/2, contf_data, 'edgecolor','none');
            shading flat
        else
            m_contourf(contf_lon_2D, contf_lat_2D, contf_data, [-999999999 contf_levels],'edgecolor','none');
        end
        caxis([min(contf_levels) max(contf_levels)]);
        colormap(gca, colmap)

        if cbar_mode
            if isempty(cbar_location)      cbar_location = 'eastoutside';    end 

            cb = colorbar(cbar_location);
            if isempty(cbar_YTick) && isempty(cbar_YTickLabel)
                cbar_YTick = contf_levels(2:cbar_interval:end-1);
                cbar_YTickLabel = cbar_YTick;
            end
            set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', fontsize) %/ cbar Ytick for diverging colormap
            set(get(cb,'Title'),'String',contf_unit, 'Fontsize', fontsize)
        end
        % set(cb, 'YTick', contflevels,'Fontsize', 12) %/ cbar Ytick for diverging colormap
        % cbarrow;
        drawnow; pause(0.05);

        fprintf('Time cost in m_contourf: %f \n', toc);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(point_data)
        tic
        if map_lon_lower < 0        %/ if plotting in [-179, 180], we convert point_data's lon to fit that range.
            a = point_data(:,1);
            a(a > 180) = a(a > 180) - 360;
            point_data(:,1) = a;
        end

        m_line(point_data(:,1), point_data(:,2), 'marker', marker, 'markersize', markersize, 'markerfacecolor', markerfacecolor,...
              'linest', 'none', 'color', 'none', 'linewi', linewi, 'clip', 'on');
        drawnow; pause(0.05);
        fprintf('Time cost in drawing points: %f \n', toc);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hatch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(hatch_data)
    
        %/ https://www.mathworks.com/matlabcentral/fileexchange/50706-offsets-and-missing-data-via-pcolor-and-surf'
        warning(sprintf('Since pcolor has offsets and edge problem, here we solve the offset problem by shifting by a half grid.\nFor more info, see https://www.mathworks.com/matlabcentral/fileexchange/50706-offsets-and-missing-data-via-pcolor-and-surf'));
        res_lon = abs(unique(round(diff(hatch_lon), 3)));  %/ round up to 3 decimal points to avoid numerical errors
        res_lat = abs(unique(round(diff(hatch_lat), 3)));

        hatch_lon = hatch_lon - res_lon/2;
        hatch_lat = hatch_lat - res_lat/2;
        
        %/ NOTE: For global hatch data, shifting a half grid will cause the first lon to be -ve.
        %/       Need to move the first -ve lon (e.g., -0.5) to the last to avoid bug
        if glb_data_mode && hatch_lon(1) < 0         
            hatch_lon = hatch_lon([2:end, 1]);  
            hatch_data = hatch_data([2:end, 1], :);
        end

        color_hatch_pve = [0 204 102]./255;  %/ green  = +ve
        color_hatch_nve = [153 0 255]./255;  %/ purple = -ve
        
        if hatch_mode == 1

            hatch_pve = hatch_data;
            hatch_pve(hatch_pve <= hatch_thres) = nan;
            hatch_pve(~isnan(hatch_pve)) = 1;
            hatch_pve(isnan(hatch_pve)) = 0;
            hatch_pve_bndry = get_bndry_from_logi('logical_2D', hatch_pve, 'bndry_lon', hatch_lon, 'bndry_lat', hatch_lat, 'outputcell', 0);

            hatch_nve = hatch_data;
            hatch_nve(hatch_nve >= -1*hatch_thres) = nan;
            hatch_nve(~isnan(hatch_nve)) = 1;
            hatch_nve(isnan(hatch_nve)) = 0;
            hatch_nve_bndry = get_bndry_from_logi('logical_2D', hatch_nve, 'bndry_lon', hatch_lon, 'bndry_lat', hatch_lat, 'outputcell', 0);
            
            
            m_line(hatch_pve_bndry(:,1), hatch_pve_bndry(:,2), 'color', color_hatch_pve, 'linest', 'none');  %/ remove the outline
            m_hatch(hatch_pve_bndry(:,1), hatch_pve_bndry(:,2), 'single', 30, hatch_density, 'linest', '-', 'linewi', hatch_linewi, 'color', color_hatch_pve); % ...with hatching added
            drawnow; pause(0.05);

            m_line(hatch_nve_bndry(:,1), hatch_nve_bndry(:,2), 'color', color_hatch_nve, 'linest', 'none');  %/ remove the outline
            m_hatch(hatch_nve_bndry(:,1), hatch_nve_bndry(:,2), 'single', 30, hatch_density, 'linest', '-', 'linewi', hatch_linewi, 'color', color_hatch_nve); % ...with hatching added
            drawnow; pause(0.05);
        end
        if hatch_mode == 0
            %/ +ve hatch 
            [ind_lon_sig, ind_lat_sig] = find(~isnan(hatch_data) & hatch_data > hatch_thres);
            hatch_pve = nan(length(ind_lon_sig), 2);
            for k = 1:length(ind_lon_sig)
                hatch_pve(k,1) = hatch_lon(ind_lon_sig(k));
                hatch_pve(k,2) = hatch_lat(ind_lat_sig(k));
            end

            %/ -ve hatch 
            [ind_lon_sig, ind_lat_sig] = find(~isnan(hatch_data) & hatch_data < -1*hatch_thres);
            hatch_nve = nan(length(ind_lon_sig), 2);
            for k = 1:length(ind_lon_sig)
                hatch_nve(k,1) = hatch_lon(ind_lon_sig(k));
                hatch_nve(k,2) = hatch_lat(ind_lat_sig(k));
            end

            %/ if plotting in [-179, 180], we convert hatch data's lon to fit that range. (no need to use conv_to_circular_data beforehand)
            if map_lon_lower < 0        
                a = hatch_pve(:,1);
                a(a > 180) = a(a > 180) - 360;
                hatch_pve(:,1) = a;

                a = hatch_nve(:,1);
                a(a > 180) = a(a > 180) - 360;
                hatch_nve(:,1) = a;
            end
            
            %/ default setting 1
            hatch_marker = 'o'; hatch_markersize = 1.5; 

            m_line(hatch_pve(:,1), hatch_pve(:,2), 'marker', hatch_marker, 'markersize', hatch_markersize, 'markerfacecolor', color_hatch_pve,...
                  'linest', 'none', 'color', 'none', 'linewi', linewi, 'clip', 'on');
            drawnow; pause(0.05);

            m_line(hatch_nve(:,1), hatch_nve(:,2), 'marker', hatch_marker, 'markersize', hatch_markersize, 'markerfacecolor', color_hatch_nve,...
                  'linest', 'none', 'color', 'none', 'linewi', linewi, 'clip', 'on');
            drawnow; pause(0.05);
        end
        fprintf('Time cost in drawing hatched pattern/points: %f \n', toc);
    end
    
    %/ testing whether the coast is correctly drawn...
    % hold on;
    % [~,LONGS, LATS] = m_lldist([81.8825 81.8825], [-89 89], 40);
    % m_line(LONGS, LATS, 'marker', 'none', 'markersize', 0.1, 'markerfacecolor', 'none',...
    %           'linest', '-', 'color', 'b', 'linewi', 2);
    % drawnow; pause(0.05);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Global Plateau %%%%%%%%%%%%%%%%%%%%%%%%%
    % tic
    if glb_plateau_mode
        m_elev('contour', [plateau_hgt plateau_hgt], 'Color', plateau_col,'LineWidth', coast_wi);   %/ m_elev is more accurate!!!
        drawnow; pause(0.05);
    end
    % fprintf('Time cost in topo: %f \n', toc);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Contour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(cont_data)

        for cc = 1:length(cont_levels)
            if cont_levels(cc) == 0 continue; end
            cont_color = cont_colmap(cc,:);

            if ~isempty(cont_data_raw)   %/ then we plot this raw cont data with dashed contour.
                m_contour(cont_lon_2D, cont_lat_2D, cont_data_raw,...
                          [cont_levels(cc), cont_levels(cc)],'LineWidth',0.5,...
                          'LineStyle','--','Color',cont_color);
                hold on;
            end

            [CS_cont, h_cont] = m_contour(cont_lon_2D, cont_lat_2D, cont_data,...
                                          [cont_levels(cc), cont_levels(cc)],'LineWidth',cont_linewi,...
                                          'LineStyle','-','Color',cont_color);
            if ~isempty(cont_labelsize)
    %             if mod(cc,2)== 0 % only label some contours
                clabel(CS_cont, h_cont, 'Color','k','FontSize',cont_labelsize,'FontWeight','bold', 'labelspacing',2000); %, 'BackgroundColor',[1 1 1]);
    %             end
                hold on; %/ if not hold on, then changes in clabel() will not be retained.
            end

        end
        drawnow; pause(0.05);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(traj_data)

        %/ NOTE: Assume traj_data is in (noOftraj, trajtime, variables)
        %/       where variables is in (lon, lat, x)
        [ntraj, ntrajtime, ~] = size(traj_data);

        ind = find(traj_data(:,:,1) > 180);
        if ~isempty(ind)            error('traj_data lon must in [-179 180]!');  end

        tic
        for n = 1:ntraj
            traj_data_each = squeeze(traj_data(n,:,:));

            nan_cnt_per_traj = 0;
            for tj = 1:ntrajtime-1
                %/ NOTE: Since we need two points to draw a line, here we
                %/       average the intensity at two timesteps to define the color.

                %/ Define the color based on the value of the 3rd variable
                mean_value = mean(traj_data_each(tj:tj+1, 3));

                if isnan(mean_value)
                    nan_cnt_per_traj = nan_cnt_per_traj + 1;
                    warning('NaN detected in a traj (%d/%d). Skip drawing this traj segment!', nan_cnt_per_traj, ntrajtime-1);

                else
                    [ind, ~] = max(find(traj_levels(1:end-1) < mean_value ));
                    if isempty(ind)    error('No index color for the value %f!', mean_value);   end
                end

                m_line(traj_data_each(tj:tj+1,1), traj_data_each(tj:tj+1,2), 'marker','none', 'color', c_traj(ind,:), 'linewi',linewi,...
                                'linest','-','markersize', markersize, 'markerfacecolor', 'none');
                hold on;

                %/ add a black point to the end of trajectory
                if tj == ntrajtime-1
                    m_line(traj_data_each(tj+1,1), traj_data_each(tj+1,2), 'marker','o', 'color', 'none', 'linewi',linewi,...
                             'linest','none','markersize',markersize+4, 'markerfacecolor', [0 0 0]);
                end
                hold on;
            end
            drawnow;
        end
        hold on;
        caxis([min(traj_levels) max(traj_levels)]);
        colormap(gca, c_traj)

        if cbar_mode
            cb = colorbar('eastoutside');
            cbar_YTick = traj_levels(2:cbar_interval:end-1);
            cbar_YTickLabel = cbar_YTick;

            set(cb, 'YTick', cbar_YTick, 'YTickLabel', cbar_YTickLabel, 'Fontsize', fontsize) %/ cbar Ytick for diverging colormap
            set(get(cb,'Title'),'String', traj_unit, 'Fontsize', fontsize)
        end
        % set(cb, 'YTick', contflevels,'Fontsize', 12) %/ cbar Ytick for diverging colormap
        % cbarrow;
        drawnow; pause(0.05);
        fprintf('Time cost in plotting trajs: %f \n', toc);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% coast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(coast_col)
        m_coast('linewidth',coast_wi,'color', [.4 .4 .4]);
    else
        m_coast('linewidth',coast_wi,'color', coast_col);
    end
    drawnow; pause(0.05);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(bndry_data)
        if iscell(bndry_data) %/ if point_data is a cell array, then loop over it to plot.
            for i = 1:length(bndry_data)

                ind = find(bndry_data{i}(:,1) > 180);
                bndry_data{i}(ind,1) = bndry_data{i}(ind,1) - 360;  %/ convert to [-179 180] if ind is not empty.

                m_line(bndry_data{i}(:,1),bndry_data{i}(:,2), 'marker', 'none', 'markersize', 0.001, 'markerfacecolor', 'none',...
                          'linest', '-', 'color', color, 'linewi', linewi);
            end
        else
            ind = find(bndry_data(:,1) > 180);
            bndry_data(ind,1) = bndry_data(ind,1) - 360;  %/ convert to [-179 180] if ind is not empty.

            m_line(bndry_data(:,1),bndry_data(:,2), 'marker', 'none', 'markersize', 0.001, 'markerfacecolor', 'none',...
                      'linest', '-', 'color', color, 'linewi', linewi);
        end
        drawnow; pause(0.05); %/ don't put drawnow function in a for-loop, or it will get extremely slow!! 
        % if ~isempty(bndry_data)
        %     m_line(bndry_data(:,1), bndry_data(:,2), 'linewi',1.5,'color','k');
        % %     m_hatch(bndry_data(:,1), bndry_data(:,2), 'single',30, 5,'color', color);
        % %     drawnow; pause(0.05);
        % end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(Udata) && ~isempty(Vdata)
        tic

    %     if corrUV                    %/ in this mode, Udata is actually Udata_sig.
    %         Udata(isnan(Udata)) = 0; %/ set NaN to zeros, since I need to compute wind magn.
    %         Vdata(isnan(Vdata)) = 0; 
    %         size(Udata)
    %         size(Vdata)
    %         vector_colormap = jet(9);
    %         vector_colormap(1,:) = [0.4 0.4 0.4]; %set r=0.1~0.2 be in grey
    %         for xx = 2:vec_step_lon:length(uv_lon)-1
    %             disp(strcat('drawing vectors', num2str(xx), {' of '}, num2str(length(uv_lon)-1)))
    %             
    %             for yy = 2:vec_step_lat:length(uv_lat)-1
    %                 UVCorr_Magn = max([abs(Udata(yy,xx)), abs(Vdata(yy,xx))]);
    %                 if UVCorr_Magn < 0.1
    %                     continue
    %                 else
    %                     ind = floor(UVCorr_Magn*10);
    %                     m_vec(vecscale, x(yy, xx), y(yy, xx),...
    %                           Udata(yy, xx)*vecscale2, Vdata(yy, xx)*vecscale2, vector_colormap(ind,:),...
    %                           'shaftwidth',shaftwidth ,'headlength',headlength, 'centered','yes', 'edgeclip','yes',...
    %                           'EdgeColor','k', 'linewidth', 0.5);
    %                     hold on;
    %                 end
    %             end
    %         end
    %         %/ corr wind ref
    %         i = 0;
    %         for vecmagn = 0.1:0.1:0.9
    %             i = i + 1;
    %             lbs = strcat(num2str(vecmagn));
    %             [~, htv5] = m_vec(vecscale, uv_lon(1) + (i-1)*17, -60,...
    %                             vecmagn*vecscale2, 0, vector_colormap(i,:),...
    %                             'shaftwidth',shaftwidth ,'headlength',headlength, 'centered','yes', 'edgeclip','no', 'key', lbs);
    %             hold on;
    %             [~, htv5] = m_vec(vecscale, uv_lon(1) + (i-1)*17, -60,...
    %                             vecmagn*vecscale2, 0, vector_colormap(i,:),...
    %                             'shaftwidth',shaftwidth ,'headlength',headlength, 'centered','yes', 'edgeclip','no', 'EdgeColor','k', 'linewidth', 0.5);
    %             set(htv5,'FontSize',fontsize, 'color', 'k', 'HorizontalAlignment', 'center');
    %             hold on
    %         end

    % else
        %/ Normal mode for vector drawing
        %/ NOTE: x, y, Udata and Vdata have the 1st dim = lat

        m_vec(vecscale, x(2:vec_step_lat:end-1, 2:vec_step_lon:end-1),...
                y(2:vec_step_lat:end-1, 2:vec_step_lon:end-1),...
                Udata(2:vec_step_lat:end-1, 2:vec_step_lon:end-1)*vecscale2,...
                Vdata(2:vec_step_lat:end-1, 2:vec_step_lon:end-1)*vecscale2,...
                vector_color, 'shaftwidth',shaftwidth ,'headlength',headlength, 'centered','yes', 'edgeclip','yes', 'EdgeColor', vector_edgecolor, 'linewidth', 0.5);
        hold on;

        %/ m_vec wind reference
        if ~isempty(vec_mag_ref)
            lon_range = map_lon_upper - map_lon_lower;
            lat_range = map_lat_upper - map_lat_lower;

            vec_ref_lon = map_lon_upper - lon_range*0.1;
            vec_ref_lat = map_lat_lower - lat_range*0.1;

            [hp, ht] = m_vec(vecscale, vec_ref_lon, vec_ref_lat, vec_mag_ref*vecscale2, 0,...
                            vector_color, 'shaftwidth',shaftwidth ,'headlength',headlength, 'centered','yes', 'edgeclip', 'off', 'key', vec_lbs);
            hold on;
            set(ht,'FontSize',fontsize, 'color', vector_color, 'HorizontalAlignment', 'center');
            set(hp,'EdgeColor', vector_edgecolor);                             %/ have to set the EdgeColor using set() if using 'key' in m_vec!
        end

        hold on;
        drawnow; pause(0.05);
        fprintf('Time cost in vectors: %f \n', toc);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% m_grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('grid_mode = %d\n', grid_mode)
    if grid_mode == 0              %/ not drawing girds
        m_grid('xtick', [], 'ytick', [-60:30:60], 'linewi',2, 'linest','none', 'gridcolor', 'k', 'tickdir','in',...
               'xaxisloc', xaxisloc, 'backcolor', backcolor, 'fontsize',fontsize);

    elseif grid_mode == 1         %/ normal
        xtick = -180:30:360;    ytick = -90:20:90;
        m_grid('xtick', xtick, 'ytick', ytick,...
               'linewi',2,'linest','none','gridcolor', 'k', 'tickdir','in',  'backcolor', backcolor, 'fontsize',fontsize);

    elseif grid_mode == 2         %/ draw grids on each 2 degree
        xtick = -180:1:360;    ytick = -90:1:90;
        m_grid('xtick', xtick, 'ytick', ytick,... % 'yticklabels', xtick(x_st_ind:5:x_ed_ind), 'yticklabels',  ytick(y_st_ind:5:y_ed_ind),...
               'linewi',2,'linest',':','gridcolor', 'k', 'tickdir','in', 'xaxisloc', xaxisloc,  'backcolor', backcolor, 'fontsize',fontsize-8);

    elseif grid_mode == 3         %/ for regional study
        xtick = -180:30:360;    ytick = -90:10:90;
        m_grid('xtick', xtick, 'ytick', ytick,...
               'linewi',2,'linest','none','gridcolor', 'k', 'tickdir','in',  'backcolor', backcolor, 'fontsize',fontsize);    

    end
    drawnow; pause(0.05) 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Text on the map %%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(text_data)
        for i = 1:length(text_data)
            h = m_text(text_data{i,1}, text_data{i,2}, text_data(i,3), 'color', [text_data{i,4}], 'fontweight', 'bold',...
                    'fontsize', fontsize-3, 'horizontalalignment', 'center');
        end
        drawnow; pause(0.05)  %/ don't put drawnow function in a for-loop, or it will get extremely slow!! 
    end

    tic
    th = title(titlename, 'fontsize', fontsize+5, 'fontweight', 'bold');
    % titlePos = get(th, 'position');
    % titlePos(2) = titlePos(2) + 0.06;   %/ change y position of title
    % set(th, 'position', titlePos)
    fprintf('Time cost in title: %f \n', toc);
    drawnow; pause(0.05)
end

if ~isempty(savepath)
    tic
    FigName_underscore = strrep(savepath, ' ', '_');

    if isempty(fig_fmt)
        fig_fmt = 'pdf';
    end

    if isequal(fig_fmt, 'png')
        export_fig(char(strcat(FigName_underscore, '.', fig_fmt)),'-r300','-png','-opengl', '-nocrop');
    else
        export_fig(char(strcat(FigName_underscore, '.', fig_fmt)),'-pdf','-painters','-nocrop');
    end

    fprintf('!!! Plot is saved: %s !!! \n', string(FigName_underscore));
    fprintf('Time cost in saving figure: %f \n', toc);
end


end
 
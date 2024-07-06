%%
function add_mask_mountains(west, east, south, north, plateau_hgt)
       
%     if mask_mountains
        m_grid('xtick', [], 'ytick', [], 'linewi', 0.0001, 'linest', 'none', 'gridcolor', 'none'); %/ do not delete this. otherwise figures will mismatch.
       
        %/ NOTE: m_elev outputs 1 deg topo; m_tbase outputs much higher res.
        %/       Prefer the coarser one to do masking.
       
        [topo, topo_lon_2D, topo_lat_2D] = m_elev([west, east, south, north]); %REGION =[west east south north];
%         [topo, topo_lon_2D, topo_lat_2D] = m_tbase([map_lon_lower, map_lon_upper, map_lat_lower, map_lat_upper]); %REGION =[west east south north];
        topo     = topo';               %/ make sure in (lon,lat) dim
        topo_lon = topo_lon_2D(1,:)';
        topo_lat = topo_lat_2D(:,1);
        topo(topo <= plateau_hgt) = nan;
        topo(~isnan(topo))        = plateau_hgt;  % so now it's only a single-value field.
       
        contf_levels_topo = [0, plateau_hgt];
       
        colmap2 = [  1   1   1;
                    .7 .7 .7;];  %/ masked by grey boxes
        hold on;
        %=============================
        ax_topo = axes;                %/ axes is to create a new axes!                  
        set(ax_topo, 'Visible', 'off') %/ IMPORTANT: hide the 2nd axes
        %=============================
       
        %/ NOTE: 'Since pcolor has offsets and edge problem, here we solve the offset problem by shifting by a half grid.\nFor more info, see https://www.mathworks.com/matlabcentral/fileexchange/50706-offsets-and-missing-data-via-pcolor-and-surf'
        res_lon = abs(unique(round(diff(topo_lon), 3)));  %/ round up to 3 decimal points to avoid numerical errors
        res_lat = unique(round(diff(topo_lat), 3));  %/ do not take absolute value, seems to be correct. Can check afterwards.
       
        m_pcolor(topo_lon_2D-res_lon/2, topo_lat_2D-res_lat/2, topo', 'edgecolor','none');
        shading flat
        colormap(ax_topo, colmap2);                                  %/ IMPORTANT: set colormap for the 2nd axes
        caxis(ax_topo,[min(contf_levels_topo) max(contf_levels_topo)]);
       
        %=============================
        set(ax_topo, 'Visible', 'off')  
        linkaxes([gca ax_topo]);                            %/ IMPORTANT: link the two overlaying axes so they match at all times to remain accurate
        hold on;
        %=============================
%     end
end
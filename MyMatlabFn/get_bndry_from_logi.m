%%
function bndry_data = get_bndry_from_logi(varargin)

    pnames = {'logical_2D', 'bndry_lon', 'bndry_lat', 'outputcell',...
              'map_lon_lower', 'map_lon_upper', 'map_lat_lower', 'map_lat_upper', 'glb_data_mode',...
              'CONN', 'hole_mode', 'draw_rings', 'edgebased'};
    dflts  = cell(1, length(pnames));
    [          logical_2D,   bndry_lon,   bndry_lat,   outputcell,...
               map_lon_lower, map_lon_upper, map_lat_lower, map_lat_upper, glb_data_mode,...
               CONN,   hole_mode,   draw_rings,   edgebased]...
            = internal.stats.parseArgs(pnames, dflts, varargin{:});

%%
    if     isempty(draw_rings)          draw_rings = 1;              end   %/ assume to draw rings (i.e. not filling holes)
    if     draw_rings                   hole_mode = 'holes';               %/ A must
    elseif isempty(hole_mode)           hole_mode = 'noholes';       end    
    if     isempty(CONN)                CONN      = 4;               end   %/ directions to connection (4 or 8)
    if     isempty(edgebased)           edgebased = 1;               end   %/ outline grids by their edges
    
    %/ Convert lon, lat into row vectors, Check the gridding
    bndry_lon_bc = bndry_lon;
    bndry_lat_bc = bndry_lat;
    flag_2Dlonlat = 0;
    if isvector(bndry_lon_bc) && isvector(bndry_lat_bc)
        if size(bndry_lon_bc, 2) == 1       bndry_lon_bc = bndry_lon_bc';    end 
        if size(bndry_lat_bc, 2) == 1       bndry_lat_bc = bndry_lat_bc';    end 
        if diff(bndry_lat_bc(1:2)) < 0      flag_lat_descending = 1;   else   flag_lat_descending = 0;  end
        res_lon = unique(abs(diff(bndry_lon_bc))/2);
        res_lat = unique(abs(diff(bndry_lat_bc))/2);
        
        if length(res_lon) ~= 1 || length(res_lat) ~= 1
            warning('Detect that the given lon, lat are not evenly distributed. Stop the transformation from centered-based to edge-based gridding.')
            edgebased = 0;   %/ overwrite
        end  
    else
        warning('Detect that the given lon, lat are 2D. Stop the transformation from centered-based to edge-based gridding.')
        flag_2Dlonlat = 1;
        edgebased = 0;
    end
    
    %/ [IMPORTANT] Re-gridding logical_2D from center-based to edge-based
    if edgebased
        logical_2D_edge = logical_2D;

        %/ Append one more row and column to fill the void after shifting
        bndry_lon_edge = [bndry_lon_bc + res_lon, bndry_lon_bc(end) + 3*res_lon];

        %/ Handle lon dim
        if bndry_lon_edge(end)-360 == bndry_lon_edge(1)  %/ then we need not to append extra grids (i.e., lon is circular)
            flag_circular_lon   = 1;
            bndry_lon_edge(end) = [];
        else
            flag_circular_lon = 0;
            logical_2D_edge   = [zeros(1, length(bndry_lat_bc)); logical_2D_edge];
        end

        %/ Handle lat dim
        if flag_lat_descending  %/ e.g. [90 89 .... -90]  -> [90.5 89.5 .... -90.5]
            bndry_lat_edge  = [bndry_lat_bc + res_lat,    bndry_lat_bc(end) - res_lat];        %/ it's fine if the edge exceeds 90N, and we append a more grid 
            logical_2D_edge(:,end+1) = zeros(size(logical_2D_edge, 1), 1);           %/ Append one more lat dim at the end

        else                    %/ e.g. [-90 -89 .... 90] -> [-90.5 -89.5 ....90.5]
            bndry_lat_edge  = [bndry_lat_bc(1) - res_lat, bndry_lat_bc + res_lat];
            logical_2D_edge = [zeros(size(logical_2D_edge, 1), 1), logical_2D_edge]; %/ Append one more lat dim at the beginning
        end

        %/ Add the three other vertex points for each target pixel 
        %/ Shift the gridding northwestward
        [ind_row, ind_col] = find(logical_2D_edge);
        for i = 1:length(ind_row)
            if ind_row(i) == 1 && flag_circular_lon   %/ then 0 means the last index of lon dim
                if flag_lat_descending
                    logical_2D_edge(ind_row(i),  ind_col(i)+1) = true;
                    logical_2D_edge(end,         ind_col(i))   = true;
                    logical_2D_edge(end,         ind_col(i)+1) = true;
                else
                    logical_2D_edge(ind_row(i),  ind_col(i)-1) = true;
                    logical_2D_edge(end,         ind_col(i))   = true;
                    logical_2D_edge(end,         ind_col(i)-1) = true;
                end
            else
                if flag_lat_descending
                    logical_2D_edge(ind_row(i),  ind_col(i)+1) = true;
                    logical_2D_edge(ind_row(i)-1,ind_col(i))   = true;
                    logical_2D_edge(ind_row(i)-1,ind_col(i)+1) = true;
                else
                    logical_2D_edge(ind_row(i),  ind_col(i)-1) = true;
                    logical_2D_edge(ind_row(i)-1,ind_col(i))   = true;
                    logical_2D_edge(ind_row(i)-1,ind_col(i)-1) = true;
                end
            end
        end
        logical_2D = logical_2D_edge;
        bndry_lon_bc  = bndry_lon_edge;
        bndry_lat_bc  = bndry_lat_edge;
    end
    %/=====================================================================
    if glb_data_mode %/ assume to plot in a lon range of [-179, 180) or [-360 + map_lon_upper, map_lon_upper)
        if isempty(map_lon_lower)        map_lon_lower = -179;   end 
        if isempty(map_lon_upper)        map_lon_upper = 180;    end 
        if isempty(map_lat_lower)        map_lat_lower = -89;    end 
        if isempty(map_lat_upper)        map_lat_upper = 89;     end 

        contains_nve = find(bndry_lon_bc < 0);   %/ check if contains -ve
        ind          = find(bndry_lon_bc > map_lon_upper);
        if ~isempty(ind) && ~isempty(contains_nve)
            error('The input bndry_lon_bc contains both -ve values and values > map_lon_upper (%.2f)! Suggest cleaning the -ve values first!', map_lon_upper);
        end

        %>>>>> The MOST IMPORTANT BLOCK TO PLOT HATCH IN LON RANGE [-179 180] or up to map_lon_upper <<<<%
        if ~isempty(ind)
            bndry_lon_bc(ind) = bndry_lon_bc(ind) - 360;                 
            ind = find(bndry_lon_bc < 0); 
            ind_trans = [ind, 1:ind(1)-1]';   
            logical_2D_trans = logical_2D(ind_trans, :);                           %/ put the left to the right: from [0,...,180,-179,...,-1] to [-179,...,-1,0,...,180]
            bndry_lon_bc = bndry_lon_bc(ind_trans);                                %/ also do so to lon vector.
        else
            %/ Then no need to transform lon and logical_2D.
            logical_2D_trans = logical_2D;
        end
    else
        if isempty(map_lon_lower)        map_lon_lower = 0;      end 
        if isempty(map_lon_upper)        map_lon_upper = 359;    end 
        if isempty(map_lat_lower)        map_lat_lower = -89;    end 
        if isempty(map_lat_upper)        map_lat_upper = 89;     end 
        
        logical_2D_trans = logical_2D;    %/ assume to plot in a lon range of [0, 360)
    end
    
    fprintf('*** mode: ''%s'' ***\n', hole_mode)
%     logical_2D_trans = bwareafilt(logical_2D_trans, 1, 'largest', CONN); % Take largest object only.
    [bndry_data, ~, ~, A] = bwboundaries(logical_2D_trans, CONN, hole_mode); %/ 8 directions to connect. clockwise.
                      %/ Do NOT output the label matrix since the logical
                      %/ matrix has been translated.
                      
    %/ [IMPORTANT] Remove unnecessary mid-points at straight lines
    for k = 1:length(bndry_data)
        bndry_data_nomidpoint = bndry_data{k};
        D = diff(bndry_data_nomidpoint, 1);
        ind_important_row = []; prev_row = [];
        for i = 1:size(D,1)
            if ~isequal(D(i,:), prev_row)
                ind_important_row = [ind_important_row, i];
            end
            prev_row = D(i,:);  %/ update
        end
        %/ Finally, append the first index at last to 'wrap up' the region [IMPORTANT]!
        ind_important_row = [ind_important_row, 1];
        bndry_data{k} = bndry_data_nomidpoint(ind_important_row,:);
    end
    
    if diff(bndry_lat_bc(1:2)) < 0
        %/  Since lat is in descending order, the boundary points from the bwboundaries function actually goes in counter-clockwise on the map! 
        %/  Now flip its orientation to clockwise on the map!
        bndry_data = cellfun(@(x) flip(x, 1), bndry_data, 'UniformOutput', false);
    end                 

    %/ Draw rings without filling holes for the use of m_hatch (hints from Rich Pawlowicz)
    %/ *** CAUTION: This function can NOT handle a child in a parent! ***
    if draw_rings
        remove_hole_list = []; remove_par_list = [];
        for par = 1:length(A)    
            bndry_par  = bndry_data{par};

            %/ handle parents (i.e., outer boundaries)
            if ~isempty(map_lon_lower) && ~isempty(map_lon_upper) && ~isempty(map_lat_lower) && ~isempty(map_lat_upper)
                par_xx = bndry_lon_bc(bndry_par(:,1));
                par_yy = bndry_lat_bc(bndry_par(:,2));
                ind = find(par_xx > map_lon_lower & par_xx < map_lon_upper & ...
                           par_yy > map_lat_lower & par_yy < map_lat_upper);
                if isempty(ind)
        %             fprintf('par = %d, skipped \n', par)
                    remove_par_list = [remove_par_list; par];
                    continue;                                                      %/ skip this parent object as it is out of the map!
                else
                    %/ IMPORTANT: if ANY bndry nodes found within the domain, 
                    %/            MAKE an in-map point be the start point 
                    %/            Why? Cos if using the original start point that is out of the map, 
                    %/            then it will mess up the hatched pattern!

%                     bndry_par(end,:) = [];                                %/ remove the original end point (necessary?)
                    bndry_par = bndry_par([ind(1):end, 1:ind(1)-1], :);     %/ change the start point
                    bndry_par(end+1,:) = bndry_par(1,:);                    %/ add the end point by copying the new start point
                end
            end

            %/ handle childs (i.e. holes)
            hole_list = find(A(:, par));                                           %/ return the label of holes of a parent
            for i = 1:length(hole_list)                                            %/ if hole_list is empty, for loop return nothing.
                hole = hole_list(i);
                bndry_hole = bndry_data{hole};
                bndry_hole = flip(bndry_hole, 1);                                  %/ reverse the hole's boundary to make it counter-clockwise on the map

                %/ append hole's boundary after the parent's, then append the starting pt of the parent to complete a 'ring'!
                bndry_par = [bndry_par; bndry_hole; bndry_par(1,:)]; 
            end
            bndry_data{par} = bndry_par;                                           %/ update the parent's boundary
            remove_hole_list = [remove_hole_list; hole_list];                      %/ record labels of all holes for deletion later
        end

        remove = unique([remove_par_list; remove_hole_list]);
        bndry_data(remove) = [];                                                   %/ remove the child objects from the boundary data!
    end
    
    if flag_2Dlonlat   %/ then assume their dims to be in (lon,lat)
        if outputcell
            for k = 1:length(bndry_data)
                ind_lon = bndry_data{k}(:,1);
                ind_lat = bndry_data{k}(:,2);
                a = nan(length(ind_lon), 2);
                for i = 1:length(ind_lon)
                    a(i,:) = [bndry_lon_bc(ind_lon(i), ind_lat(i))', bndry_lat_bc(ind_lon(i),ind_lat(i))']; %/ convert indices to lon,lat
                end
                bndry_data{k} = a;
            end
        else  %/ then numerical array using [nan, nan] to separate boundaries.
            for k = 1:length(bndry_data)
                ind_lon = bndry_data{k}(:,1);
                ind_lat = bndry_data{k}(:,2);
                a = nan(length(ind_lon), 2);
                for i = 1:length(ind_lon)
                    a(i,:) = [bndry_lon_bc(ind_lon(i), ind_lat(i))', bndry_lat_bc(ind_lon(i),ind_lat(i))']; %/ convert indices to lon,lat
                end
                bndry_data{k} = a;
                bndry_data{k}(end+1,:) = [nan, nan];    %/ append nan rows (as a delimiter) to each cell.
            end
            bndry_data = cat(1, bndry_data{:});  %/ concatenate all bndry cells.   
        end
    else
        if outputcell
            for k = 1:length(bndry_data)
                bndry_data{k} = [bndry_lon_bc(bndry_data{k}(:,1))', bndry_lat_bc(bndry_data{k}(:,2))']; %/ convert indices to lon,lat
            end

        else  %/ then numerical array using [nan, nan] to separate boundaries.
            for k = 1:length(bndry_data)
                bndry_data{k} = [bndry_lon_bc(bndry_data{k}(:,1))', bndry_lat_bc(bndry_data{k}(:,2))']; %/ convert indices to lon,lat
                bndry_data{k}(end+1,:) = [nan, nan];    %/ append nan rows (as a delimiter) to each cell.
            end
            bndry_data = cat(1, bndry_data{:});  %/ concatenate all bndry cells.   
        end
    end
end
            
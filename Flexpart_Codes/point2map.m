function [ind_lon, ind_lat] = point2map(varargin)

    pnames = {'pt_x', 'pt_y', 'lon', 'lat', 'NumWorkers'};  
    dflts  = cell(1, length(pnames));
    [pt_x, pt_y, lon, lat, NumWorkers] = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
%%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 28 Jun 2024
    %/
    %/ Description: To output indices of lon, lat for the given point
    %/              (pt_x, pt_y)
    %/ 
    %/              pt_x, pt_y can be scalar or vectors, but not matrices.
    %/=====================================================================
    
    %/ Make sure they are all of the same double type
    pt_x = double(pt_x);
    pt_y = double(pt_y);
    lon = double(lon);
    lat = double(lat);

    res_lon = unique(diff(lon));
    res_lat = unique(diff(lat));
    
    if length(res_lon) ~= 1 || length(res_lat) ~= 1
        error('The function cannot handle uneven lon/lat!!')
    end

    % if ~isempty(find(lon < 0, 1))       
    %     error('lon must be [0, 359]!!');         
    % end

    pt_x(pt_x < 0) = pt_x(pt_x < 0) + 360;   

%     a_lon = 1/res_lon;
%     a_lat = 1/res_lat;

    % pt_x
    % res_lon
    x_nearest = round(pt_x/res_lon, 0)*res_lon; %/ Feel free to validate this trick!
    y_nearest = round(pt_y/res_lat, 0)*res_lat;
    x_nearest(x_nearest == 360) = 0;  %/ Restore 360E (if any) to 0E after round-up. 

    ind_lon = nan(length(x_nearest),1);
    ind_lat = nan(length(y_nearest),1); 
    
    if ~isempty(NumWorkers) 
        if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
            parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
        end

        parfor i = 1:length(x_nearest)   %/ much faster if the input pt_x pt_y arrays are long
            if isempty(find(lon - x_nearest(i) == 0, 1)) || isempty(find(lat - y_nearest(i) == 0, 1))
                %/ No grid cell is found for the given point. set ind NaN;
                ind_lon(i) = nan;
                ind_lat(i) = nan;
            else
            
                ind_lon(i) = find(lon - x_nearest(i) == 0);    %/ Then find the nearest 1 deg grid.
                ind_lat(i) = find(lat - y_nearest(i) == 0);
            end
        end

    else
        for i = 1:length(x_nearest)  
            if isempty(find(lon - x_nearest(i) == 0, 1)) || isempty(find(lat - y_nearest(i) == 0, 1))
                %/ No grid cell is found for the given point. set ind NaN;
                ind_lon(i) = nan;
                ind_lat(i) = nan;
            else
                ind_lon(i) = find(lon - x_nearest(i) == 0);    %/ Then find the nearest 1 deg grid.
                ind_lat(i) = find(lat - y_nearest(i) == 0);
            end
        end
    end
%     if find(isnan(ind_lon))   error('ind_lon contains nan!! Check the code!!'); end
%     if find(isnan(ind_lat))   error('ind_lat contains nan!! Check the code!!'); end
        
end
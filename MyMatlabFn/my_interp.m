function [data_new, lon_2D_new, lat_2D_new] = my_interp(varargin)

    pnames = {'lon_old', 'lat_old', 'data', 'lon_new', 'lat_new', 'is_global', 'lon_dim', 'interp_method', 'NumWorkers'};
    dflts  = cell(1, length(pnames));
              [lon_old,   lat_old,   data,   lon_new,   lat_new,   is_global,   lon_dim,   interp_method,   NumWorkers]...
                           = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: Mar 22, 2024
    %/
    %/ Description: This function was designed to interpolate data of
    %/              3D or higher dimension using simple parfor/for-loops.           
    %/              
    %/              Why? Because interp2 request 3D meshgrid, which can 
    %/              oftentimes blow up the memory.
    %/=====================================================================
    
    %/ Check the input data
    if isempty(data)                    error('data is empty! Please check.');      end
    if length(size(data)) > 4           error('This function currently only handles data of dimensions up to 4D!');    end
    if isempty(interp_method)           interp_method = 'linear';  end
    if isempty(lon_dim)    
        lon_dim = 1;
        % warning('lon_dim is missing. Assuming lon_dim == 1.');  
    end
    time_dim = length(size(data));
    [~, ~, n3, n4] = size(data);
    
    %/ Broadcast vars
    lon_old_bc = lon_old;
    lat_old_bc = lat_old;
    lon_old_res = abs(diff(lon_old_bc(1:2)));

    %/ Check the input lon lat
    if ~isvector(lon_old_bc) && ~isvector(lat_old_bc)   
        fprintf('*** Detected lon_old and lat_old are matrices. Transpose them if they are in (lon, lat) -> be consistent with the dim ordering of meshgrid().*** \n');
        lon_2D_old = lon_old_bc; 
        lat_2D_old = lat_old_bc;
        if lon_dim == 1   
            lon_2D_old = lon_2D_old';
            lat_2D_old = lat_2D_old';
        end
    else
        [lon_2D_old, lat_2D_old] = meshgrid(lon_old_bc, lat_old_bc); %/ NOTE: the meshgrid is in (nlat, nlon)
    end

    if is_global    %/ For global data, append longitudes at BOTH ends (for avoiding NaNs after interpolation)
        if lon_dim == 1
            lon_2D_old = [lon_2D_old(:,1)-lon_old_res, lon_2D_old, lon_2D_old(:,end)+lon_old_res]; %/ e.g, -0.625(APPENDED), 0, ..., 359.375, 360.625(APPENDED)
            lat_2D_old = [lat_2D_old(:,end), lat_2D_old, lat_2D_old(:,1)];
        elseif lon_dim == 2
            lon_2D_old = [lon_2D_old(1,:)-lon_old_res, lon_2D_old, lon_2D_old(end,:)+lon_old_res];
            lat_2D_old = [lat_2D_old(end,:), lat_2D_old, lat_2D_old(1,:)];
        else
            error('lon_dim can only be 1 or 2!')
        end

        % size(data)
        if time_dim == 3
            if lon_dim == 1
                data = cat(1, data(end,:,:), data, data(1,:,:));
            elseif lon_dim == 2
                data = cat(2, data(:,end,:), data, data(:,1,:));
            else
                error('lon_dim can only be 1 or 2!')
            end
        elseif time_dim == 4
            if lon_dim == 1
                data = cat(1, data(end,:,:,:), data, data(1,:,:,:));
            elseif lon_dim == 2
                data = cat(2, data(:,end,:,:), data, data(:,1,:,:));
            else
                error('lon_dim can only be 1 or 2!')
            end
        else
            error('code not set!')
        end
        % size(data)
    end


    % if is_global    %/ For global data, append longitudes at BOTH ends (for avoiding NaNs after interpolation)
    %     if lon_dim == 1
    %         lon_old_bc = [lon_old_bc(1,:)-lon_old_res; lon_old_bc; lon_old_bc(end,:)+lon_old_res]; %/ e.g, -0.625(APPENDED), 0, ..., 359.375, 360.625(APPENDED)
    %     elseif lon_dim == 2
    %         lon_old_bc = [lon_old_bc(1,:)-lon_old_res, lon_old_bc, lon_old_bc(end,:)+lon_old_res];
    %     else
    %         error('lon_dim can only be 1 or 2!')
    %     end
    % 
    %     % size(data)
    %     if lon_dim == 1
    %         data = cat(1, data(end,:,:), data, data(1,:,:));
    %     elseif lon_dim == 2
    %         data = cat(2, data(:,end,:), data, data(:,1,:));
    %     else
    %         error('lon_dim can only be 1 or 2!')
    %     end
    %     % size(data)
    % end
      

    if ~isvector(lon_new) && ~isvector(lat_new)
        lon_2D_new = lon_new;
        lat_2D_new = lat_new;
        fprintf('*** Detected lon_new and lat_new are matrices. ***\n'); 
    else
        [lon_2D_new, lat_2D_new] = meshgrid(lon_new, lat_new); %/ NOTE: the meshgrid is in (nlat, nlon)
    end
    
    if isequal(lon_2D_old, lon_2D_new) && isequal(lat_2D_old, lat_2D_new)
        fprintf('!!! The requested new gridding is the same as the old gridding! No interpolation is performed. !!!\n')
        data_new = data;
        return
    end

    nlon_new = size(lon_2D_new', 1);
    nlat_new = size(lon_2D_new', 2);
    
    ndim = length(size(data));
    data_new = nan(nlon_new, nlat_new, n3, n4);
    
    if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
        parpool('Threads', NumWorkers); %/ use process-based parpool 
        % parpool('Threads');
    end
 
    fprintf('*** Start interpolating (%s) ... ***\n', interp_method); tic;
    if isequal(interp_method, 'grid2grid') %/ my code block to "interpolate" data that is non-monotonically increasing 
                                          %/ But alternatively we can directly use griddata() instead!
        error('This ''grid2grid'' method is not very recommended. Are you sure to use it?');
        
        %/ Check if the lon, lat dim order is flipped; Transpose them if so. 
        if isequal(size(lon_2D_old), [size(data,2), size(data,1)])
            lon_2D_old = lon_2D_old';
            lat_2D_old = lat_2D_old';
        end

        d1 = abs(diff(lon_new(1:2)));
        d2 = abs(diff(lat_new(1:2)));

        lon_2Dto1D_old = reshape(lon_2D_old, [], 1);
        lat_2Dto1D_old = reshape(lat_2D_old, [], 1);
        data           = reshape(data, n1*n2, n3, n4);  %/ Trick
        
        for i = 1:nlon_new
            for j = 1:nlat_new
                %/ Find all the grid cells that fall within the new grid box 
                ind = find(lon_2Dto1D_old >= lon_new(i)-d1/2 & lon_2Dto1D_old < lon_new(i)+d1/2 &...
                           lat_2Dto1D_old >= lat_new(j)-d2/2 & lat_2Dto1D_old < lat_new(j)+d2/2);

                %/ Take the mean (by default)
                if ndim == 2
                    data_new(i,j) = mean(data(ind), 1, 'omitnan');
                elseif ndim == 3
                    data_new(i,j,:) = mean(data(ind,:), 1, 'omitnanm');
                elseif ndim == 4
                    data_new(i,j,:,:) = mean(data(ind,:,:), 1, 'omitnan');
                end
            end
        end
    else
        %================================= Tips ==================================%
        %         % Fast to create interpolant F and evaluate multiple times
        %         F = scatteredInterpolant(X,Y,V)
        %         v1 = F(Xq1,Yq1)
        %         v2 = F(Xq2,Yq2)
        % 
        %         % Slower to compute interpolations separately using griddata
        %         v1 = griddata(X,Y,V,Xq1,Yq1)
        %         v2 = griddata(X,Y,V,Xq2,Yq2)
        %=========================================================================%
        % size(lon_2D_old)
        % size(lat_2D_old)
        % size(data')
        % size(lon_2D_new)
        % size(lat_2D_new)
        if ndim == 2
            data_new = griddata(lon_2D_old,  lat_2D_old, data', lon_2D_new, lat_2D_new, interp_method)'; %/ Note the transpose to "correct" the seq of dims.
            
        elseif ndim == 3
            parfor t = 1:n3   %/ parfor loop is faster, for-loop may not initial multi-tread computation
            % for t = 1:n3   %/ parfor loop is faster, for-loop may not initial multi-tread computation
                data_new(:,:,t) = griddata(lon_2D_old,  lat_2D_old, data(:,:,t)', lon_2D_new, lat_2D_new, interp_method)'; %/ Note the transpose to "correct" the seq of dims.
            end
        
        elseif ndim == 4
            parfor t = 1:n3  %/ parfor loop is faster, for-loop may not initial multi-tread computation
            % for t = 1:n3  %/ parfor loop is faster, for-loop may not initial multi-tread computation
                for k = 1:n4
                    data_new(:,:,t,k) = griddata(lon_2D_old,  lat_2D_old, data(:,:,t,k)', lon_2D_new, lat_2D_new, interp_method)'; %/ Note the transpose to "correct" the seq of dims.
                end
            end
        else
            error('The input data is neither 2D, 3D nor 4D. Check the data dimension!')
        end
        
%         if flag_flip
%             data_new   = flip(data_new, 2); %/ restore the original lat ordering
%             lat_2D_new = flip(lat_2D_new', 2);
%             lon_2D_new = lon_2D_new';
%         end
    end
    % fprintf('*** Interpolation Done. Time cost: %.2f sec ***\n', toc);
    
    %/ Close parpool when not using
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj);
    end

end
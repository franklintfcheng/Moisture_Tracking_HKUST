function s = vertinte(varargin)

    pnames = {'data', 'level', 'level_unit', 'lnP_coor', 'sp_mode', 'sp', 'topo_lon', 'topo_lat',  'interp_mode', 'take_vertsum', 'take_vertmean', 'zm_or_mm', 'NumWorkers'}; 
    dflts  = {    [],      [],        'hPa',          0,        [],   [],         [],         [],             [],              0,               0,         [],           30};
    [           data,   level,   level_unit,   lnP_coor,   sp_mode,   sp,   topo_lon,   topo_lat,    interp_mode,   take_vertsum,   take_vertmean,   zm_or_mm,   NumWorkers] = ...
        internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/=====================================================================
    %/      Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: 25 Apr 2024
    %/
    %/ Description: The function ouputs a mass-weighted vertical integral
    %/              or a layer mean if 'take_vertmean' is toggled on.
    %/
    %/  'data':                     It must be 3D (lon, lat, level) or
    %/                              4D (lon, lat, level, time). 
    %/                              The 3rd dim is always the level dim.
    %/
    %/  'sp':                       Surface pressure (useful to determine
    %/                              data below terrain).
    %/
    %/  'lnP_coor':                 Trenberth & Guillemot (1995) said using dlnP in the
    %/                              vertical integration can reach a slightly 
    %/                              higher accuracy. See if needed.
    %/
    %/  'sp_mode':            'sp': Remove pressure-level data based on
    %                               surface pressure.
    %/                      'topo': Remove pressure-level data 
    %/                              below terrain prior to the
    %/                              calculation, using z2p_hydrostatic().
    %/                          []: Do nothing.
    %/
    %/  'interp_mode':    'linear': Perform linear interpolation to the
    %/                              level of sp.
    %/=====================================================================

    % fprintf('*** Running vertinte (vertical integration) with take_vertmean = %d ... ***\n', take_vertmean)
    %/ Broadcast vars
    level_bc = double(level); %/ Always ensure it is double (not integer), otherwise interp1 would report a bug
    
    %/ Check how many levels are there.
    if length(level_bc) == 1   
        warning('Since there is only one level, no vertical integral is performed. Returning the data...')  
        s = data;
        return
    end
    
    if isequal(level_unit, 'hPa')     %/ by default
        CONV_p = 100;               %/ convert from hPa to Pa
    elseif isequal(level_unit, 'Pa')
        CONV_p = 1;
    else
        error('P-levels are only in hPa or Pa. Modify the function code if needed.')
    end

    level_bc = level_bc*CONV_p;
    if size(level_bc,2) == 1   %/ switch to a row array
        level_bc = level_bc';    
    end   

    [nlon, nlat, nlev, ntime] = size(data);      %/ if 3D, then ntime = 1
    if nlev ~= length(level_bc)  
        error('nlev dim inconsistent!'); 
    end

    if ~isempty(sp_mode) && ~isequal(sp_mode, 'sp') && ~isequal(sp_mode, 'topo')
        error('The only valid input of ''sp_mode''     is [], ''sp'' or ''topo''!'); 
    end
    
    if ~isempty(interp_mode) && ~isequal(interp_mode, 'linear')
        error('The only valid input of ''interp_mode'' is   []   or ''linear''!');
    end
    %======================================================================
    % NOTE: The logic of 'sp' and 'topo' mode is different. 
    %         the former replaces the lowest level with sp (above terrain), 
    %         the latter simply removes levels below terrain
    %======================================================================
    
    %/ [IMPORTANT] Remove values below terrain before integration
    if isempty(sp_mode)
        data_bc_new  = data;  %/ (nlon, nlat, nlev, ntime) or [nlon, nlat, nlev)
        level_bc_new = level_bc; %/ (   1, nlev)

    else
        if isequal(sp_mode, 'sp')     %/ Using surface pressure
            if isempty(sp)  
                error('Input ''sp'' data when setting ''sp_mode'' == ''sp''!'); 
            end
            sp_bc = sp;
            sp_bc = reshape(sp_bc, nlon, nlat, 1, ntime);  %/ lon lat time -> lon lat level(sfc) time  

        elseif isequal(sp_mode, 'topo')    %/ Using topography and hydrostatic balance to estimate surface pressure (sp)
            if isempty(topo_lon) || isempty(topo_lat)  
                error('topo_lon and topo_lat are not set for sp_mode == ''%s''!', sp_mode); 
            end
            [lon_2D_new, lat_2D_new] = meshgrid(topo_lon, topo_lat); 
            
            %/ high-res topo (see if necessary) - a more accurate result.
            m_proj('Miller Cylindrical'); %/ NOTE: put m_proj before m_tbase() to initialize coordinate system, otherwise will encounter a bug.
            [ELEV_highres, LONG, LAT] = m_tbase([min(topo_lon), max(topo_lon), min(topo_lat), max(topo_lat)]); 
            ELEV_highres_interp = interp2(LONG, LAT, ELEV_highres, lon_2D_new, lat_2D_new, 'nearest');  
            ELEV_highres_interp = ELEV_highres_interp';  %/ convert to lon x lat.
            ELEV_highres_interp(ELEV_highres_interp < 0) = 0;
            
            %/ Take zonal/merid mean topo when requested
            if ~isempty(zm_or_mm)
                ELEV_highres_interp = mean(ELEV_highres_interp, zm_or_mm);
            end

            sp_bc = z2p_hydrostatic('z', ELEV_highres_interp, 'unit', 'm'); %/ It convert altitude into pressure (hPa)
            if isequal(level_unit, 'Pa')
                sp_bc = sp_bc*100;  %/ hPa -> Pa;
            end
        
            %/ Since topo is unchanged, for coding convenience, we use repmat
            %/ to enlarge its dim to 4D
            sp_bc = repmat(sp_bc, 1, 1, 1, ntime);
        end
                
        [nlon_sp, nlat_sp, ~] = size(sp_bc);
        if ~isequal(nlon, nlon_sp) || ~isequal(nlat, nlat_sp) 
            error('Inconsistent lon/lat dimension between ''data'' and ''sp''!')
        end
        % if max(sp_bc, [], 'all') < 1e4
        %     error('sp appears to be in hPa, make sure it is in Pa!')
        % end
        
        level_bc = reshape(level_bc, 1, 1, nlev);
        level_bc = repmat(level_bc, nlon, nlat, 1, ntime);
    
        %/ Set all pressure levels greater than sp to be nan
        cond = (level_bc > sp_bc);
        data_bc_nan  = data;  data_bc_nan(cond)  = nan;      
        level_bc_nan = level_bc; level_bc_nan(cond) = nan;
        
        %/ Replace the lowest p-level with sp 
        if isempty(gcp('nocreate')) && ~isempty(NumWorkers)
            parpool('Threads', NumWorkers) 
        end
        
        level_bc_reshape     = reshape(level_bc,     nlon*nlat, nlev, ntime);
        level_bc_nan_reshape = reshape(level_bc_nan, nlon*nlat, nlev, ntime);
        sp_reshape           = reshape(sp_bc,        nlon*nlat,    1, ntime);
        data_bc_reshape      = reshape(data,      nlon*nlat, nlev, ntime);
        data_bc_nan_reshape  = reshape(data_bc_nan,  nlon*nlat, nlev, ntime);
        level_bc_new         = cell(nlon*nlat, 1);
        data_bc_new          = cell(nlon*nlat, 1);
        [~, k_lowest_lv]        = max(level);      %/ Get the index of the lowest p-level; For ECMWF data, p-level is always from small (TOA) to large (sfc).
        % ntime
        % size(level_bc)
        % size(level_bc_old_reshape)
        % size(data)
        % size(data_bc_reshape)
                                    
        if isequal(interp_mode, 'linear') && ~isempty(find(isnan(data), 1))
            warning('The original data contains NaNs! Likely that the data have considered topo. Linear extrapolation may yield NaN at some grids.')
        end

        % for ij = 1:nlon*nlat  %/ For debugging
        parfor ij = 1:nlon*nlat
            level_bc_reshape_ij     = squeeze(level_bc_reshape(ij,:,:));           %/ (nlev, ntime)
            level_bc_nan_reshape_ij = squeeze(level_bc_nan_reshape(ij,:,:));       %/ (nlev, ntime)
            sp_reshape_ij           = squeeze(sp_reshape(ij,:,:));                 %/ (   1, ntime)
            data_bc_reshape_ij      = squeeze(data_bc_reshape(ij,:,:));            %/ (nlev, ntime)
            data_bc_nan_reshape_ij  = squeeze(data_bc_nan_reshape(ij,:,:));        %/ (nlev, ntime)
            %/ If ntime = 1, transpose it to make sure level is always in
            %/ the 1st dim of 'level_bc_reshape_ij' and 'sp_reshape_ij'
            if ntime == 1
                level_bc_reshape_ij     = level_bc_reshape_ij';
                level_bc_nan_reshape_ij = level_bc_nan_reshape_ij';
                data_bc_reshape_ij      = data_bc_reshape_ij';
                data_bc_nan_reshape_ij  = data_bc_nan_reshape_ij';
            end
            sp_reshape_ij = sp_reshape_ij';
            
            % size(level_bc_reshape_ij)
            for t = 1:ntime
                k = find(~isnan(squeeze(level_bc_nan_reshape_ij(:,t))), 1, 'first');  %/ Because the pressure levels greater than sp have been set nan
                if isempty(k)                                                       %/ Meaning that sp must be > 1e5 Pa.
                    k = k_lowest_lv;  %/ Replace the lowest p-level with sp
                end
                level_bc_nan_reshape_ij(k,t) = sp_reshape_ij(t);
                
                %/ IF perform interpolation of data to the level of sp
                if isequal(interp_mode, 'linear')

                    ind_k = k:k+1;                                      %/ For ECMWF data, p-level is always from small (TOA) to large (sfc)
                    x     = level_bc_reshape_ij(ind_k,t);               %/ The original level
                    xq    = level_bc_nan_reshape_ij(ind_k,t);           %/ The quired level
                    v     = data_bc_reshape_ij(ind_k,t);

                    % x
                    % xq
                    % v
                    data_bc_nan_reshape_ij(ind_k,t) = interp1(x,v,xq,'linear','extrap');  %/ Using 1D interp (with extrapolation!) and update 
                    
                    % %/ Double-check
                    % if ij == 1 && t == 1
                    %     x
                    %     xq
                    %     v
                    %     data_bc_nan_reshape_ij(ind_k,t)
                    % end
                end
            end
            level_bc_new{ij} = level_bc_nan_reshape_ij;
            if isequal(interp_mode, 'linear')
                data_bc_new{ij} = data_bc_nan_reshape_ij;
            end
        end
        
        %/ Update the data after sp_mode
        level_bc_new = cat(3, level_bc_new{:});                           %/ Concatenate the grid dimension at the 3rd dim
        level_bc_new = permute(level_bc_new, [3 1 2]);                    %/ Permute the grid dimension from the 3rd to the 1st dim
        level_bc_new = reshape(level_bc_new, nlon, nlat, nlev, ntime);    %/ Reshape to restore (lon,lat,lev,time)
        if ~isempty(interp_mode)
            data_bc_new = cat(3, data_bc_new{:});                         %/ Concatenate the grid dimension at the 3rd dim
            data_bc_new = permute(data_bc_new, [3 1 2]);                  %/ Permute the grid dimension from the 3rd to the 1st dim
            data_bc_new = reshape(data_bc_new, nlon, nlat, nlev, ntime);  %/ Reshape to restore (lon,lat,lev,time)
        else
            data_bc_new = data;  
        end
        clear data; 
    end

    %/ Make sure level_bc_new_2D is a double (for Matlab does not allow
    %/ operation between double and integer numbers)
    level_bc_new = double(level_bc_new); 
    
    %/ Vertical integration using the Trapesoidal Rule (it must be in
    %/ forward scheme in order to sum the area under the curve.
    if length(size(data_bc_new)) == 3                       %/ [lon lat level]
        data_2D         = reshape(data_bc_new,  [], nlev);  %/ Reshape data to 2D
        level_bc_new_2D = reshape(level_bc_new, [], nlev);  %/ Reshape level to 2D
        clear data_bc_new;
        
    elseif length(size(data_bc_new)) == 4                               %/ [lon lat level time]
        data_bc_new    = permute(data_bc_new,  [1,2,4,3]);             %/ Permute data to [lon lat time level]
        data_2D        = reshape(data_bc_new,  [], nlev);              %/ Reshape data to 2D

        clear data_bc_new;
        if ~isempty(sp_mode)
            level_bc_new_permute = permute(level_bc_new, [1,2,4,3]);         %/ Permute level to [lon lat time level]
            level_bc_new_2D      = reshape(level_bc_new_permute, [], nlev);  %/ Reshape level to 2D
        else
            level_bc_new_2D = level_bc_new;  %/ Otherwise, it's already a 2D row vector
        end
    end
    s    = zeros(nlon*nlat*ntime, 1);
    dP   = abs(diff(level_bc_new_2D,1,2));                                  %/ abs() to ensure the integral is from small to large
    % dP
    % size(data_2D)
    lnP  = log(level_bc_new_2D);  
    dlnP = abs(diff(lnP,1,2));                                              %/ For lnP_coor == 1
    
    %/ Perform Trapesoidal Rule
    for p = 1:nlev-1
        if lnP_coor             %/ in lnP coordinate
            layer_mean  = mean(data_2D(:,p:p+1).*level_bc_new_2D(:,p:p+1), 2, 'omitnan');  %/ The lower-level may contain nan due to levels below topo in some weather conditions
            s_EachLayer = layer_mean .* dlnP(:,p);
            s_EachLayer(isnan(s_EachLayer)) = 0;                            %/ In this case, we just omit the operation in the vertical integration. (not to change it unless you have strong reasons!)
            s = s + s_EachLayer;
        else                    %/ in P coordinate
            layer_mean  = mean(data_2D(:,p:p+1), 2, 'omitnan');             %/ The lower-level may contain nan due to levels below topo in some weather conditions
            s_EachLayer = layer_mean .* dP(:,p);
            s_EachLayer(isnan(s_EachLayer)) = 0;                            %/ In this case, we just omit the operation in the vertical integration. (not to change it unless you have strong reasons!)
            s = s + s_EachLayer;                                            %/ Aka. trapesoidal integration 
        end
    end
    clear data_2D;

    %/ Final Output
    if take_vertsum
        %/ do nothing
        fprintf('Output the ordinary vertical integral (without dividing it by g)...\n')

    elseif take_vertmean
        fprintf('Output the vertical mean...\n')
        s = s./sum(dP, 2, 'omitnan'); %/ Take a layer mean
        
    else
        fprintf('Output the mass weighted vertical integral...\n')
        g = 9.81;
        s = s./g;   %/ Mass-weighted vertical integral (See Eq. 2.5 in Neelin & Held 1987) 
    end

    %/ Reshape from 2D back to the original dimensions
    s = reshape(s, nlon, nlat, ntime); 

end


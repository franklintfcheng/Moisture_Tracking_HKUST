function [reg_2D, reg_bndry_list, reg_name_list, reg_id_list, reg_centr_list] = reg_extractor(varargin)

    pnames = {'lon', 'lat', 'slct_reg', 'mth', 'output_HR_ocean', 'outline_by_grid', 'draw_rings', 'regime_ver', 'output_logical', 'data_folder', 'saveload_cond_landocean', 'savemat',  'savenc', 'nc_remark', 'recompute', 'verbose'};
    dflts  = {   [],    [],         [],    [],                 0,                 0,            0,        'n/a',                0,            [],                         1,         0,         0,          '',           0,         0};
    [           lon,   lat,   slct_reg,   mth,   output_HR_ocean,   outline_by_grid,   draw_rings,   regime_ver,   output_logical,   data_folder,   saveload_cond_landocean,   savemat,    savenc,   nc_remark,   recompute,   verbose] = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %======================================================================
    %/      Author: Franklin CHENG (franklin.cheng@ust.hk)
    %/ Last update: 24 Apr 2024
    %/
    %/ Description: This function is to output all my previously prescribed
    %/              box regions or the global river basins from HydroSHEDS.
    %/
    %=====================================================================
    if isempty(data_folder)
        error('Specify ''data_folder''!');   
    end
    
    %/ Always UPDATE 'str_mth'!!
    str_mth = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec',....
               'MAM', 'JJA', 'SON', 'DJF', 'AMJJAS', 'ONDJFM', 'JFD', 'nonJJA', 'MJJASO', 'NDJFMA'};

    %/ If input lon lat are both a vector
    if isvector(lon) && isvector(lat)
        flag_2D_lonlat  = 0;
        lon_bc          = lon;
        lat_bc          = lat;
        if size(lon_bc, 1) == 1    lon_bc = lon_bc';   end  %/ convert it into a column 
        if size(lat_bc, 1) == 1    lat_bc = lat_bc';   end  %/ convert it into a column 
        nlon = length(lon_bc);
        nlat = length(lat_bc);
        lon_bc(lon_bc > 180) = lon_bc(lon_bc > 180) - 360;     %/ for inpolygon later.
        [lon_bc_2D, lat_bc_2D] = meshgrid(lon_bc, lat_bc);
        lon_bc_2D     = lon_bc_2D';  
        lat_bc_2D     = lat_bc_2D';
        lon_bc_2Dto1D = reshape(lon_bc_2D, [], 1);
        lat_bc_2Dto1D = reshape(lat_bc_2D, [], 1);
    else
        %/ If lon lat are matrices (e.g., uneven grid cells)
        flag_2D_lonlat          = 1;
        lon_bc                  = lon;  
        lat_bc                  = lat;
        [nlon, nlat]            = size(lon_bc);
        lon_bc(lon_bc > 180)    = lon_bc(lon_bc > 180) - 360;   %/ for inpolygon later.
        lon_bc_2D               = lon_bc;
        lat_bc_2D               = lat_bc;
        lon_bc_2Dto1D           = reshape(lon_bc, [], 1);
        lat_bc_2Dto1D           = reshape(lat_bc, [], 1);
    end
    
    %/ [IMPORTANT] Pre-process cond_land & cond_ocean [Save time and Map_Toolbox license!]
    lon_grids = lon_bc;
    lat_grids = lat_bc;
    res = abs(diff(lon_grids(1:2)));
    filename_cond_landocean = strcat(data_folder, sprintf('cond_landocean_%dx%d_lon%.2f_%.2f_from%.2fto%.2f_lat%.2f_%.2f_from%.2fto%.2f_res%.2fdeg',...
                                     length(lon_grids), length(lat_grids), min(lon_grids), max(lon_grids), lon_grids(1), lon_grids(end), min(lat_grids), max(lat_grids), lat_grids(1), lat_grids(end), res), '.mat');
    
    if saveload_cond_landocean == 0
        [~, cond_land, cond_ocean] = show_land_or_ocean_hydrosheds('lon_grids', lon_grids, 'lat_grids', lat_grids);
    else
        if ~isfile(filename_cond_landocean)
            [~, cond_land, cond_ocean] = show_land_or_ocean_hydrosheds('lon_grids', lon_grids, 'lat_grids', lat_grids);

            if savemat 
                if verbose
                    fprintf('\n*** To avoid using Map_Toolbox, saving cond_land & cond_ocean into %s ...***\n', filename_cond_landocean)
                end
                save(filename_cond_landocean, 'cond_land', 'cond_ocean', '-v7.3');
            end
        else
            if verbose
                fprintf('\n*** To avoid using Map_Toolbox, loading cond_land & cond_ocean from %s ...***\n', filename_cond_landocean)
            end
            load(filename_cond_landocean, 'cond_land', 'cond_ocean');
        end
    end

    %/ Set a dummy value -9999 for oceans, +9999 for land, 10000 for TP.
    value_land         = 9999;
    value_ocean        = -9999;
    reg_2D             = nan(nlon,nlat);
    reg_2D(cond_land)  = value_land;
    reg_2D(cond_ocean) = value_ocean;
    if ~isempty(find(isnan(reg_2D), 1))   
        error('Check the completeness of cond_land and cond_ocean!');   
    end
    
    %/ Load the TP boundaries
    TP_bndry_data = box_region('TP'); %/ Using {} to explicitly make it as a cell.
    try
        [cond_TP, ~] = inpoly2([lon_bc_2Dto1D, lat_bc_2Dto1D], TP_bndry_data);                          %/ inpoly2 is 600xx faster than inpolygon!! 
        cond_TP = reshape(cond_TP, nlon, nlat);
    catch 
%         if verbose warning(E.identifier, 'Error caught!: %s. Proceed with inpolygon().', E.message); end
        [cond_TP, ~] = inpolygon(lon_bc_2D, lat_bc_2D, TP_bndry_data(:,1), TP_bndry_data(:,2));    %/ sometimes inpoly2 could run into a bug, then we use inpolygon.
    end

    labels_customized_land  = {'Pearl', 'NewGuinea', 'NESA', 'WNMexico', 'CMexico', 'SMexico', 'Antarctica'};
    labels_hydrosheds_conti = {'Africa', 'Europe', 'Siberia', 'Asia', 'Australia', 'South America', 'North America', 'Arctic', 'Greenland'};
    basinnamelist           = hydrosheds_name('level', 3);
    labels_hydrosheds       = basinnamelist(:,1);
    where_are_hydrosheds    = basinnamelist(:,2);
    
    %/ TP basins 
    labels_TP_basins        = {'TP-Indus',   'TP-Ganges', 'TP-Irrawaddy', 'TP-Mekong',...
                               'TP-Yangtze', 'TP-Yellow', 'TP-Gobi1', 'TP-Tarim', 'TP-Syr Darya', 'TP-Amu Darya',...
                               'TP-CapsianSea East Coast', 'TP-Helmand',...
                               'S. Inner TP', 'N. Inner TP', 'Qaidam'};
    
    %/ Mid-lower basins in connection w/ the TP basins (e.g., TP-XXX)
    labels_ML_basins        = {'ML-Indus',   'ML-Ganges', 'ML-Irrawaddy', 'ML-Mekong',...
                               'ML-Yangtze', 'ML-Yellow', 'ML-Gobi1', 'ML-Tarim', 'ML-Syr Darya', 'ML-Amu Darya',...
                               'ML-CapsianSea East Coast', 'ML-Helmand'};
                           
    %/ TP basins and their Mid-lower reaches
    labels_TP_and_ML_basins = [labels_TP_basins, labels_ML_basins];
        
    %/ Original basins with headwater from TP
    labels_TP_basins_original = {'Indus', 'Ganges', 'Irrawaddy', 'Mekong', 'Yangtze',...
                                 'Yellow', 'Gobi1', 'Tarim', 'Syr Darya', 'Amu Darya', 'CapsianSea East Coast', 'Helmand',...
                                 'S. Inner TP', 'N. Inner TP', 'Qaidam'};
    

    %/ [IMPORTANT]: Ocean-only basins
    labels_ocean = {'GoG', 'NTAO', 'STAO', 'WSAO', 'ESAO', 'SO', 'GoM', 'CaribeanSea', 'WNAO', 'ENAO', 'HudsonBay', 'AO',...
                    'MeditSea', 'BlackSea', 'CaspianSea', 'AS', 'BoB', 'RedSea', 'PersianOmanGulf', 'WIO', 'EIO', 'SCS', 'GoT',...
                    'WTPO', 'WSPO', 'YSECS', 'SoJ', 'WNPO', 'ENPO', 'ETPO', 'ESPO'};    

    %/ [IMPORTANT]: Group all the river basins that must take land grids only
    %/              -> This will avoid mistakenly treat new box region as land only
    % size(labels_hydrosheds')
    % size(labels_TP_and_ML_basins)
    labels_land = [labels_hydrosheds', labels_TP_and_ML_basins, labels_TP_basins_original];


    %/ Random domains that contain both land and oceans (for TP study)
    labels_Pan_domain = {'Pan_IM',                  'Pan_Westerly',           'Pan_EAM',...
                         'Pan_IM_noTP',             'Pan_Westerly_noTP',      'Pan_EAM_noTP'};
                     
    labels_regime = {'NH_MidLat_Westerly',      'NH_Polar_Westerly',      'IM',           'EAM',      'Easterly',      'Overturning_elsewhere', ...
                     'NH_MidLat_Westerly_noTP', 'NH_Polar_Westerly_noTP', 'IM_noTP',      'EAM_noTP', 'Easterly_noTP', 'Overturning_elsewhere_noTP', 'SH_Westerly', 'Elsewhere_noTP'};


    labels_IPCC_PAK = {'IPCC-PAK-NEAF', 'IPCC-PAK-SEAF', 'IPCC-PAK-WCA', 'IPCC-PAK-TIB',...
                       'IPCC-PAK-ARP', 'IPCC-PAK-SAS', 'IPCC-PAK-ARS', 'IPCC-PAK-BOB',...
                       'IPCC-PAK-EIO', 'IPCC-PAK-SIO', 'IPCC-PAK-rest'}; 
    if isempty(slct_reg)
        warning('[reg_extractor]: slct_reg is empty. Return nothing!')
        return;
    elseif ischar(slct_reg)
        slct_reg = {slct_reg};
    end
    
    %/ Check whether the input 'slct_reg' refers to a collection of basins.
    if isequal(slct_reg, {'global'})        %/ Simply output a matrix with all ones
        if isvector(lon) && isvector(lat)
            reg_2D         = ones(length(lon), length(lat));
        else
            reg_2D         = ones(size(lon));  %/ when the input lon lat are 2D matrices
        end
        reg_bndry_list = [];
        reg_name_list  = slct_reg;   
        reg_id_list    = [];
        reg_centr_list = [];
        return;
        
    elseif any(ismember(slct_reg, {'oceans'}))    
        slct_reg_bc = labels_ocean;   
        remaining_reg = setdiff(slct_reg,{'oceans'});
        slct_reg_bc   = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)
        
    elseif any(ismember(slct_reg, {'hydrosheds'}))      
        slct_reg_bc = [labels_hydrosheds_conti, {'Antarctica'}];     
        remaining_reg = setdiff(slct_reg,{'hydrosheds'});
        slct_reg_bc   = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)
        
    elseif any(ismember(slct_reg, {'hydrosheds_oceans'}))      
        slct_reg_bc = [labels_hydrosheds_conti, {'Antarctica'}, labels_ocean];    
        remaining_reg = setdiff(slct_reg,{'hydrosheds_oceans'});
        slct_reg_bc   = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)
        
    elseif any(ismember(slct_reg, {'hydrosheds_TPbasins_oceans', 'hydrosheds_TPbasins'}))  %/ Consider TP basins and their Mid-lower reaches.
        labels_hydrosheds_bc = labels_hydrosheds;       %/ Make a new bc var, do NOT modify labels_hydrosheds itself!
        
        %/ First, remove the original TP basins to avoid redundancy
        ind = findismember_loop(labels_hydrosheds_bc, labels_TP_basins_original);
        if ~isequal(labels_hydrosheds(ind)', labels_TP_basins_original)
            error('Inconsistency in basin name! Check ''labels_hydrosheds'' or ''labels_TP_basins_ori''!');
        end
        labels_hydrosheds_bc(ind) = [];
        
        if isequal(slct_reg, {'hydrosheds_TPbasins_oceans'})
            slct_reg_bc = [labels_hydrosheds_bc', {'Antarctica'}, labels_TP_and_ML_basins, labels_ocean];  
        elseif isequal(slct_reg, {'hydrosheds_TPbasins'})
            slct_reg_bc = [labels_hydrosheds_bc', {'Antarctica'}, labels_TP_and_ML_basins]; 
        end
        remaining_reg = setdiff(slct_reg,{'hydrosheds_TPbasins_oceans', 'hydrosheds_TPbasins'});
        slct_reg_bc   = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)
        
    elseif any(ismember(slct_reg, {'hydrosheds_TP_oceans'}))  %/ Consider TP basins and their Mid-lower reaches.
        labels_hydrosheds_bc = labels_hydrosheds;             %/ Make a new bc var, do NOT modify labels_hydrosheds itself!
        
        %/ First, remove the original TP basins to avoid redundancy
        ind = findismember_loop(labels_hydrosheds_bc, labels_TP_basins_original);
        if ~isequal(labels_hydrosheds(ind)', labels_TP_basins_original)
            error('Inconsistency in basin name! Check ''labels_hydrosheds'' or ''labels_TP_basins_ori''!');
        end
        labels_hydrosheds_bc(ind) = [];
        slct_reg_bc = [labels_hydrosheds_bc', {'Antarctica', 'TP'}, labels_ML_basins, labels_ocean];  
        remaining_reg = setdiff(slct_reg,{'hydrosheds_TP_oceans'});
        slct_reg_bc   = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)
        
    elseif any(ismember(slct_reg, {'customized_land'}))  
        slct_reg_bc = labels_customized_land;     
        remaining_reg = setdiff(slct_reg,{'customized_land'});
        slct_reg_bc   = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)

    elseif any(ismember(slct_reg, {'TP_basins'}))   
        slct_reg_bc = labels_TP_basins;   
        remaining_reg = setdiff(slct_reg,{'TP_basins'});
        slct_reg_bc   = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)
        
    elseif any(ismember(slct_reg, {'ML_basins'}))  
        slct_reg_bc = labels_ML_basins;  
        remaining_reg = setdiff(slct_reg,{'ML_basins'});
        slct_reg_bc   = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)
        
    elseif any(ismember(slct_reg, {'TP_grids_1x1', 'TP_grids_0.5x0.5'}))
        res_lon = abs(unique(diff(lon)));
        res_lat = abs(unique(diff(lat)));
        if length(res_lon) ~= 1 || length(res_lat) ~= 1  
            error('Check if lon, lat are uniformly distributed!'); 
        end
        queried_res = split(slct_reg, 'x');
        queried_res = split(queried_res{2}, '_');
        queried_res = str2num(queried_res{1});
        
        if ~isequal(res_lon, queried_res) || ~isequal(res_lat, queried_res)
            disp(['res_lon = ', res_lon]);
            disp(['res_lat = ', res_lat]);
            disp(['queried_res = ', queried_res]);
            error('Check if lon, lat resolution is consistent with %s!', slct_reg{:});
        end
        cond_TP_2Dto1D = reshape(cond_TP, [], 1);
        TP_grids_lon   = lon_bc_2Dto1D(cond_TP_2Dto1D);
        TP_grids_lat   = lat_bc_2Dto1D(cond_TP_2Dto1D);
        slct_reg_bc    = cellstr(strcat('TP_grids_', string(res_lon), 'x', string(res_lat), '_', string(TP_grids_lon), 'E', string(TP_grids_lat), 'N')); %/ cellstr -> convert string to cell
        
        remaining_reg  = setdiff(slct_reg,{'TP_grids_1x1', 'TP_grids_0.5x0.5'});
        slct_reg_bc    = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)
        
    elseif any(ismember(slct_reg, {'TP_and_ML_basins'}))        
        slct_reg_bc   = labels_TP_and_ML_basins;  
        remaining_reg = setdiff(slct_reg,{'TP_and_ML_basins'});
        slct_reg_bc   = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)
        
    elseif any(ismember(slct_reg, {'TP_basins_ori'}))    
        slct_reg_bc   = labels_TP_basins_original;  
        remaining_reg = setdiff(slct_reg,{'TP_basins_ori'});
        slct_reg_bc   = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)

    elseif any(ismember(slct_reg, {'IPCC-PAK'}))    
        slct_reg_bc   = labels_IPCC_PAK;   
        remaining_reg = setdiff(slct_reg,{'IPCC-PAK'});
        slct_reg_bc   = [slct_reg_bc, remaining_reg];   %/ append the remaining regions (if any)
    else
        slct_reg_bc = slct_reg;
    end
    
    %/ Load IPCC-PAK regions if quried
    if any(ismember(slct_reg_bc, labels_IPCC_PAK))
        mask_filename = '/disk/r059/tfchengac/FLEXPART/domfill_EA_0.5deg_MPI_202207-08_PakCase/Masks/IPCCregions_Pakistancase.nc';
        % ncdisp(mask_filename)
        mask                = ncread(mask_filename, 'mask');  %/ lon x lat x regions
        mask_lon            = ncread(mask_filename, 'lon');
        mask_lat            = ncread(mask_filename, 'lat');
        mask_region_abbrevs = ncread(mask_filename, 'abbrevs');
        % mask_region_id      = ncread(mask_filename, 'region');
        % mask_region_fullnames = ncread(mask_filename, 'names');

        if size(mask_lon, 2) == 1   
            mask_lon = mask_lon';   %/ To row vector
        end
        if size(mask_lat, 2) == 1   
            mask_lat = mask_lat';   %/ To row vector
        end

        if sign(diff(mask_lat(1:2))) ~= sign(diff(lat(1:2)))  %/ Make sure the ordering of mask_lat and lat matches to each other!
            mask_lat = flip(mask_lat, 2);
            mask = flip(mask, 2);
        end

        %/ Remove the point at two poles. 1440 721 11 -> 1440 719 11
        % mask_lat = mask_lat(2:end-1);
        % mask     = mask(:,2:end-1,:);  
        ind_poles = (mask_lat == 90 | mask_lat == -90);
        mask_lat(ind_poles) = [];
        mask(:,ind_poles,:) = [];

        ind_poles = (lat == 90 | lat == -90);
        lat(ind_poles) = [];
        
        size(mask_lat)
        size(mask)
        size(lat)
        
        d_lon = mask_lon - lon;
        d_lat = mask_lat - lat;

        eps = 1e-8;
        if max(abs(d_lon),[],'all') > eps
            flag_bndry2reg2D = 1;
            warning('The input lon does not match the lon of IPCC-PAK mask! Using bndry_data to output reg_2D instead.');
        elseif max(abs(d_lat),[],'all') > eps
            flag_bndry2reg2D = 1;
            warning('The input lat does not match the lat of IPCC-PAK mask! Using bndry_data to output reg_2D instead.');
        else
            flag_bndry2reg2D = 0;
        end
    end

    %/ Get a list of basin boundaries by contenation
    
    reg_bndry_list    = cell(length(slct_reg_bc), 1); %/ pre-allocation
    reg_name_list     = cell(length(slct_reg_bc), 1); %/ pre-allocation
    reg_id_list       = cell(length(slct_reg_bc), 1); %/ pre-allocation
    reg_centr_list    = cell(length(slct_reg_bc), 1); %/ pre-allocation
    for k = 1:length(slct_reg_bc)
        if ismember(slct_reg_bc(k), labels_ocean) %/ if for oceans
            if output_HR_ocean
                str_HR      = '_HR';
            else
                str_HR      = '_LR';
            end
        else
            str_HR = '';
        end
        bndry_filename = strcat(data_folder, slct_reg_bc{k}, '_bndry_data', str_HR, '.mat');
        
        if isfile(bndry_filename) && recompute == 0
%             disp('yes')
            if verbose  
                fprintf('*** Loading bndry data from %s ...***\n', bndry_filename);   
            end
            load(bndry_filename, 'bndry_data', 'reg_name', 'reg_id', 'reg_centr');
        else
            %/ If input is 'Asia', then output all its river basins.
            if any(ismember(slct_reg_bc{k}, labels_hydrosheds_conti))
                %/ Load Level 3 river basin from hydrosheds
                [bndry_data, id_list, basin_name_list, ~] = retrieve_hydrosheds('slct_conti', slct_reg_bc{k}, 'simplify_bndry', 0);
                reg_name = basin_name_list;
                reg_id   = id_list;

            %/ If the name matches that of a hydrobasin (assumed to be level 3)
            elseif any(ismember(slct_reg_bc{k}, labels_hydrosheds))
                %/ Find the corresponding continent.
                I = find(ismember(labels_hydrosheds, slct_reg_bc{k}));
                if isempty(I)           error_status = 'missing';
                elseif length(I) > 1    error_status = 'redundant';
                else                    error_status = [];              end
                if ~isempty(error_status)
                    error('You are inquiring a %s name ''%s'' in ''hydroshed_name()''! Specify it uniquely in the function!', error_status, slct_reg_bc{k});
                end

                [bndry_data, id_list, basin_name_list, ~] = retrieve_hydrosheds('slct_conti', where_are_hydrosheds{I}, 'simplify_bndry', 0);

                %/ Extract the one that is enquried
                ind = find(ismember(basin_name_list, slct_reg_bc{k}));

                bndry_data = bndry_data(ind);
                reg_name   = basin_name_list(ind);
                reg_id     = id_list(ind);

            elseif contains(slct_reg_bc{k}, {'TP_grids_'})  %/ e.g., 'TP_grids_1x1_27N99E', 'TP_grids_1x1_27N100E', etc.
                %/ Outline the grid cell
                point_lon  = TP_grids_lon(k);
                point_lat  = TP_grids_lat(k);
                bndry_data = point2bndry('point_lon',point_lon,'point_lat',point_lat,'res_lon',res_lon,'res_lat',res_lat);
%                 bndry_data = {bndry_data};  %/ wrap the set of boundry cells into a big cell first (will be expanded later)
                reg_name   = slct_reg_bc{k};
                reg_id     = [];
                
            elseif ismember(slct_reg_bc{k}, {'Inner_TP', 'TP-Indus_TP-Ganges'})
                if isequal(slct_reg_bc{k}, 'Inner_TP')
                    slct_reg_comb = {'S. Inner TP', 'N. Inner TP'};
                elseif isequal(slct_reg_bc{k}, 'TP-Indus_TP-Ganges')
                    slct_reg_comb = {'TP-Indus', 'TP-Ganges'};
                end
                
                bndry_data_comb = cell(length(slct_reg_comb), 1);
                for ii = 1:length(slct_reg_comb)
                    %/ Trick: recursion (using the function from within the function)
                    [~, bndry_data_bc, ~, ~, ~] = reg_extractor('lon', lon, 'lat', lat, 'slct_reg', slct_reg_comb{ii}, 'mth', mth,...
                                                          'outline_by_grid', 0, 'draw_rings', 0, 'data_folder', data_folder,...
                                                          'savemat', savemat, 'recompute', 0); %<- do not recompute
                    bndry_data_comb{ii} = bndry_data_bc{:}; %/ avoid nested cells
                end

                if any(ismember(slct_reg_bc{k}, {'Inner_TP', 'TP-Indus_TP-Ganges'}))
                    poly1 = polyshape(bndry_data_comb{1});
                    poly2 = polyshape(bndry_data_comb{2});
                    polyout = union(poly1,poly2);
                    vertices_new = close_vertices(polyout.Vertices); %/ my function to close the vertices -> no "leaking"
                    bndry_data = {vertices_new};
                end
                % bndry_data  = {cat(1, bndry_data_comb{:})}; %/ concatenate, then make it a cell for coding convenience
                reg_name    = slct_reg_bc(k);
                reg_id      = nan;   
                
            elseif contains(slct_reg_bc{k}, 'TP-') || contains(slct_reg_bc{k}, 'ML-')  %/ TP_basins that intersect with TP boundaries.
%                 disp('yes')
                reg_name_in_part = strsplit(slct_reg_bc{k}, '-');

                %/ Load the river basin boundaries (if 'TP-Yangtze', then we load
                %/  the boundary of 'Yangtze').
                I = find(ismember(labels_hydrosheds, reg_name_in_part(2)));
                if isempty(I)           error_status = 'missing';
                elseif length(I) > 1    error_status = 'redundant';
                else                    error_status = [];              end
                if ~isempty(error_status)
                    error('You are inquiring a %s name ''%s'' in ''hydroshed_name()''! Specify it uniquely in the function!', error_status, slct_reg_bc{k});
                end

                [basin_bndry_data, ~, basin_name_list, ~] = retrieve_hydrosheds('slct_conti', where_are_hydrosheds{I}, 'simplify_bndry', 0);

                %/ Extract the one that is enquried
                ind = find(ismember(basin_name_list, reg_name_in_part(2)));

                %/ *Intersection* of two boundaries
                %/ e.g., 'TP' intersects 'Yangtze' = 'TP-Yangtze'
                poly1   = polyshape(TP_bndry_data);
                poly2   = polyshape(basin_bndry_data{ind});
                polyout = intersect(poly1, poly2);
                vertices_new = close_vertices(polyout.Vertices); %/ my function to close the vertices -> no "leaking"
                
                if contains(slct_reg_bc{k}, 'ML-')
                    %/ We then subtract 'TP-Yangtze' from 'Yangtze'
                    %/ e.g., 'Yangtze' - 'TP-Yangtze' = 'ML-Yangtze'
                    
                    %/ Do NOT append the start point in this case!!
                    TP_basin_bndry_bc  = {polyout.Vertices};  
                    poly3   = polyshape(TP_basin_bndry_bc{:});
                    polyout = subtract(poly2, poly3);
                    vertices_new = close_vertices(polyout.Vertices);  %/ my function to close the vertices -> no "leaking"
                end

                %/ Do NOT append the start point in this case!!
                bndry_data  = {vertices_new};
                reg_name    = slct_reg_bc(k);
                if contains(slct_reg_bc{k}, 'ML-')
                    ind = findismember_loop(labels_ML_basins, slct_reg_bc(k));
                    reg_id      = 20000+ind;
                else
                    ind = findismember_loop(labels_TP_basins, slct_reg_bc(k));
                    reg_id      = 10000+ind;
                end
                 
            elseif ismember(slct_reg_bc(k), labels_Pan_domain) 
 
                %/ See if we need to subtract TP from the domain
                if contains(slct_reg_bc(k), '_noTP')
                    Pan_domain = slct_reg_bc{k};
                    Pan_domain(end-4:end) = '';  %/ remove char '_noTP' 
                    Pan_domain_bndry_data = box_region(Pan_domain);

                    poly1   = polyshape(Pan_domain_bndry_data);
                    poly2   = polyshape(TP_bndry_data);
                    polyout = subtract(poly1, poly2);
                    vertices_new = close_vertices(polyout.Vertices); %/ my function to close the vertices -> no "leaking"
                    
                    %/ Supplementing points in any long straight line (to allow for curvature when map_proj = 'ortho'
                    vertices_new = add_pts_to_bndry('bndry_data', vertices_new, 'suppl_Dx', 2, 'suppl_Dy', 2);    
                    bndry_data   = {vertices_new};  
                else
                    bndry_data = {box_region(slct_reg_bc(k))};
                end
                reg_name    = slct_reg_bc(k);
                ind         = findismember_loop(labels_Pan_domain, slct_reg_bc(k));
                reg_id      = 30000+ind;            
                
            elseif ismember(slct_reg_bc(k), labels_regime) %/ NOTE: for robustness, reg_2D is readily modified without generating bndry_data.
                %/ First, load the regime matrix and check dimension
                regime_filename = strcat(data_folder, sprintf('uv300x850_regime_1971-2010_%s', str_mth{mth}), '.mat');
                load(regime_filename, 'regime', 'regime_name', 'regime_lon', 'regime_lat');
                if verbose
                    fprintf('!!! Regime matrix loaded from %s. !!!\n', regime_filename);
                end
                ind_lon     = findismember_loop(regime_lon, lon); %/ do not use lon_bc, since it is from -179 to 180
                ind_lat     = findismember_loop(regime_lat, lat);
                regime      = regime(ind_lon,ind_lat);
                regime_lon  = regime_lon(ind_lon)';
                regime_lat  = regime_lat(ind_lat)';
                if ~isequal(regime_lon, lon) || ~isequal(regime_lat, lat)  error('Dimension inconsistent!');  end
                
                if contains(slct_reg_bc(k), 'NH_MidLat_Westerly')
                    cond_1 = (lat_bc_2D > 15 & lat_bc_2D <= 60);
                    cond_2 = (regime == findismember_loop(regime_name, {'Westerly'}));
                    
                elseif contains(slct_reg_bc(k), 'NH_Polar_Westerly')
                    cond_1 = (lat_bc_2D > 60 & lat_bc_2D <= 90);
                    cond_2 = (regime == findismember_loop(regime_name, {'Westerly'}));
                    
                elseif contains(slct_reg_bc(k), 'SH_Westerly')
                    cond_1 = (lat_bc_2D >= -90 & lat_bc_2D <= -15);
                    cond_2 = (regime == findismember_loop(regime_name, {'Westerly'}));
                    
                elseif any(contains(slct_reg_bc(k), {'IM', 'EAM', 'Overturning_elsewhere'}))
                    if contains(slct_reg_bc(k), 'IM')
                        Monsoon_bndry_data = box_region('Pan_IM');
                        
                    elseif contains(slct_reg_bc(k), 'EAM')
                        Monsoon_bndry_data = box_region('Pan_EAM');
                        
                    elseif contains(slct_reg_bc(k), 'Overturning_elsewhere')
                        Monsoon_bndry_data = box_region({'Pan_IM', 'Pan_EAM'}); %/ Later will get overturning regime outside the IM and EAM domains 
                    end
                    try
                        [cond_1, ~] = inpoly2([lon_bc_2Dto1D, lat_bc_2Dto1D], Monsoon_bndry_data);                          %/ inpoly2 is 600xx faster than inpolygon!! 
                        cond_1 = reshape(cond_1, nlon, nlat);
                    catch 
%                         if verbose warning(E.identifier, 'Error caught!: %s. Proceed with inpolygon().', E.message); end
                        [cond_1, ~] = inpolygon(lon_bc_2D, lat_bc_2D, Monsoon_bndry_data(:,1), Monsoon_bndry_data(:,2));    %/ sometimes inpoly2 could run into a bug, then we use inpolygon.
                    end
                    
                    %/ Get overturning regime outside the IM and EAM domains
                    if contains(slct_reg_bc(k), 'Overturning_elsewhere')
                        cond_1 = ~cond_1;  
                    end
                    cond_2 = (regime == findismember_loop(regime_name, {'Overturning'}));
                    
                elseif contains(slct_reg_bc(k), 'Easterly')
                    cond_1 = (lat_bc_2D > -30 & lat_bc_2D < 30);
                    cond_2 = (regime == findismember_loop(regime_name, {'Easterly'}));
                    
                elseif contains(slct_reg_bc(k), 'Elsewhere')
                    %/ Trick: recursion (using the function from within the function)
%                     a = {'NH_MidLat_Westerly',      'NH_Polar_Westerly',      'IM',           'EAM',      'Easterly',      'Overturning_elsewhere',      'SH_Westerly'};
                    a = {'TP', 'NH_MidLat_Westerly_noTP', 'NH_Polar_Westerly_noTP', 'IM_noTP', 'EAM_noTP', 'Easterly_noTP', 'Overturning_elsewhere_noTP', 'SH_Westerly'};
                    if verbose fprintf('*** [reg_extractor]: Perform recursion... ***\n');  end
                    [cond_1, ~, ~, ~] = reg_extractor('lon', lon, 'lat', lat, 'slct_reg', a, 'mth', mth,...
                                                      'outline_by_grid', 0, 'draw_rings', 0, 'data_folder', data_folder,...
                                                      'savemat', savemat, 'recompute', 0); %<- do not recompute
                    cond_1(~isnan(cond_1)) = 1;
                    cond_1(isnan(cond_1))  = 0;
                    cond_1 = ~cond_1;            %/ the inverse will be 'Elsewhere'
                    cond_2 = ones(size(cond_1)); %/ dummy
                else
                    error('code not set!')
                end
                
                %/ Based on the conditions, do classification
                if contains(slct_reg_bc(k), '_noTP')
                    reg_2D(cond_1 & cond_2 & ~cond_TP) = k;
                else
                    reg_2D(cond_1 & cond_2) = k;
                end
                
                %/ Subtle modification (based on 'regime_ver')
                if ~isequal(regime_ver, 'n/a') 
                    if isempty(regime_ver)
                        %/ No modification is needed.
                    
                    elseif isequal(regime_ver, 'V2')  
                        if mth == 14 && contains(slct_reg_bc(k), 'Elsewhere')  %/ Only modify JJA classification. Do it after finalizing 'Elsewhere'
                            k_Easterly_noTP = find(ismember(slct_reg_bc, {'Easterly_noTP'}));
                            cond_3 = (reg_2D == k_Easterly_noTP);

                            North_India_bndry_data = box_region('Easterly_North_India');
                            try
                                [cond_4, ~] = inpoly2([lon_bc_2Dto1D, lat_bc_2Dto1D], North_India_bndry_data);                          %/ inpoly2 is 600xx faster than inpolygon!! 
                                cond_4 = reshape(cond_4, nlon, nlat);
                            catch 
        %                         if verbose warning(E.identifier, 'Error caught!: %s. Proceed with inpolygon().', E.message); end
                                [cond_4, ~] = inpolygon(lon_bc_2D, lat_bc_2D, North_India_bndry_data(:,1), North_India_bndry_data(:,2));    %/ sometimes inpoly2 could run into a bug, then we use inpolygon.
                            end
                            k_IM_noTP = find(ismember(slct_reg_bc, {'IM_noTP'}));
                            reg_2D(cond_3 & cond_4) = k_IM_noTP;
                        end
                    else
                        error('code not set!');
                    end
                end
                
                bndry_data  = [NaN NaN; NaN NaN;]; %/ make dummy bndry_data 
%                 bndry_data  = [];
                reg_name    = slct_reg_bc(k);
                ind         = findismember_loop(labels_regime, slct_reg_bc(k));
                reg_id      = 40000+ind;
                
            elseif ismember(slct_reg_bc(k), {'Antarctica'}) %/ Use 'get_bndry_from_logi' of the land-masked logical_2D to outline Antarctica
                if flag_2D_lonlat
                    error('''outline_by_grid'' currently does not work for uneven lon lat!');
                end
                
                %/ Get a rough outline of Antarctica (i.e., South Pole)
                bndry_data = box_region(slct_reg_bc(k));
                
                try
                    [mask_2D, ~] = inpoly2([lon_bc_2Dto1D, lat_bc_2Dto1D], bndry_data);                          %/ inpoly2 is 600xx faster than inpolygon!! 
                    mask_2D = reshape(mask_2D, nlon, nlat);
                catch 
                    % if verbose warning(E.identifier, 'Error caught!: %s. Proceed with inpolygon().', E.message); end
                    [mask_2D, ~] = inpolygon(lon_bc_2D, lat_bc_2D, bndry_data(:,1), bndry_data(:,2));    %/ sometimes inpoly2 could run into a bug, then we use inpolygon.
                end
                
                %/ Find the land grids that fall in the South Pole
                logical_2D = (cond_land & mask_2D);  
                
                %/ Create the final logical matrix for outlining boudaries
                hole_mode  = 'holes';
                bndry_data = get_bndry_from_logi('logical_2D', logical_2D, 'bndry_lon', lon_bc, 'bndry_lat', lat_bc,...
                                                 'glb_data_mode', 1, 'outputcell', 0, 'hole_mode', hole_mode, 'draw_rings', draw_rings);
                bndry_data = {bndry_data};
                reg_name   = slct_reg_bc(k);
                reg_id     = nan;
                
            elseif ismember(slct_reg_bc(k), labels_ocean) %/ if for oceans
                % disp('yes')
                %/ Load my prescribed ploygon regions
                LR_bndry_data = {box_region(slct_reg_bc(k))}; %/ Using {} to explicitly make it as a cell.
                
                if output_HR_ocean
                    %/ Trick: recursion (using the function from within the function)
                    [~, hydrosheds_bndry_list, ~, ~] = reg_extractor('lon', lon, 'lat', lat, 'slct_reg', 'hydrosheds',...
                                                                     'outline_by_grid', 0, 'draw_rings', 0, 'data_folder', data_folder,...
                                                                     'savemat', savemat, 'recompute', 0); %<- do not recompute

                    %/ First, find the nearby hydroshed basins that intersect
                    %/  with the coarse ocean basin. (slow)
                    poly1 = polyshape(LR_bndry_data{:});
                    ind_nearby_basin = [];
                    for ii = 1:length(hydrosheds_bndry_list)
%                         disp(ii)
                        poly2   = polyshape(hydrosheds_bndry_list{ii});
                        polyout = intersect(poly1, poly2);
                        
                        if ~isempty(polyout.Vertices)
                            ind_nearby_basin = [ind_nearby_basin; ii]; 
                        end
                    end
%                     disp(hydrosheds_name_list(ind_nearby_basin))

                    %/ Concatenate bndry data of the basins nearby the ocean basin
                    nearby_basin_bndry = cat(1, hydrosheds_bndry_list{ind_nearby_basin});

                    %/ Do a subtraction
                    poly3   = polyshape([nearby_basin_bndry; [nan nan]]);
                    polyout = subtract(poly1, poly3);
                    vertices_new = close_vertices(polyout.Vertices); %/ my function to close the vertices -> no "leaking"
                    
                    bndry_data  = {vertices_new};
                else
                    bndry_data = LR_bndry_data; 
                end
                reg_name = slct_reg_bc(k);
                ind      = findismember_loop(labels_ocean, slct_reg_bc(k));
                reg_id   = -1*ind;                         %/ Use negative id to denote ocean basins
                
            elseif ismember(slct_reg_bc(k), {'land', 'ocean'})
                bndry_data = {[nan nan; nan nan]};
                reg_name   = slct_reg_bc(k);
                reg_id     = nan;
                
            elseif ismember(slct_reg_bc(k), {'Pakistan_box', 'Australia_box', 'Scotland_box'})

                if ismember(slct_reg_bc(k), {'Pakistan_box'})
                    mask_filename = '/disk/r059/tfchengac/FLEXPART/domfill_EA_0.5deg_MPI_202207-08_PakCase/Masks/Mask_PakistanFlood_box_lon-180to180.nc';
                elseif ismember(slct_reg_bc(k), {'Australia_box'})
                    mask_filename = '/disk/r059/tfchengac/FLEXPART/domfill_EA_0.5deg_MPI_202202_AusCase/Masks/mask_australia_box_lon-180to180.nc';
                elseif ismember(slct_reg_bc(k), {'Scotland_box'}) 
                    mask_filename = '/disk/r059/tfchengac/FLEXPART/domfill_EA_0.5deg_MPI_202309-10_ScotCase/Masks/mask_scotland_box_lon-180to180.nc';
                end
                    
                % ncdisp(mask_filename)
                mask       = ncread(mask_filename, 'mask');
                mask_lon   = ncread(mask_filename, 'lon');
                mask_lat   = ncread(mask_filename, 'lat');
                
                if ismember(slct_reg_bc(k), {'Pakistan_box'})
                    mask(988,246) = 0;  %/ remove the straight pixel
                end
                logical_2D = logical(mask);
                
                %/ Outlining boudaries from logical matrix
                hole_mode  = 'holes';
                bndry_data = get_bndry_from_logi('logical_2D', logical_2D, 'bndry_lon', mask_lon, 'bndry_lat', mask_lat,...
                                                 'glb_data_mode', 1, 'outputcell', 0, 'hole_mode', hole_mode, 'draw_rings', 0, 'edgebased', 1);
                bndry_data = {bndry_data};
                reg_name   = slct_reg_bc(k);
                reg_id     = nan;
                
            elseif ismember(slct_reg_bc(k), labels_IPCC_PAK) %/ load IPCC AR6 regions (designed for Pakistan case)

                str = slct_reg_bc{k};
                queried_region = strrep(str, 'IPCC-PAK-', ''); %/ remove substring
                ind = findismember_loop(mask_region_abbrevs, queried_region);
                if isempty(ind)
                    error('The quried region %s is not found from the mask data!', queried_region);
                end
                logical_2D = logical(mask(:,:,ind));

                %/ Outlining boudaries from logical matrix
                hole_mode  = 'holes';
                bndry_data = get_bndry_from_logi('logical_2D', logical_2D, 'bndry_lon', mask_lon, 'bndry_lat', mask_lat,...
                                                 'glb_data_mode', 1, 'outputcell', 0, 'hole_mode', hole_mode, 'draw_rings', 0, 'edgebased', 1);

                bndry_data = {bndry_data};
                reg_name   = slct_reg_bc(k);
                reg_id     = nan;
                
            else
                %/ Load my prescribed ploygon regions
%                 disp('yes')
                bndry_data = {box_region(slct_reg_bc(k))}; %/ Using {} to explicitly make it as a cell.
                reg_name   = slct_reg_bc(k);
                reg_id     = k;
            end
            
            %/ Get the centroid of the boundary and save it if not 'land' or 'ocean'.
            reg_centr = [];
            if ~ismember(slct_reg_bc(k), [{'land', 'ocean'}, labels_regime])
                if any(ismember(slct_reg_bc{k}, labels_hydrosheds_conti)) %/ We have to do a looping over all basins in that continent
                    reg_centr = nan(length(bndry_data), 2);
                    for i = 1:length(bndry_data)
                        polyin = polyshape(bndry_data{i}(:,1), bndry_data{i}(:,2));
                        [x,y] = polycenter(polyin);  %/ this is designed by Chad, much better than centr().
                        reg_centr(i,:) = [x,y];
                    end
                else
                    polyin = polyshape(bndry_data{:}(:,1), bndry_data{:}(:,2));
                    [x,y] = polycenter(polyin);  %/ this is designed by Chad, much better than centr().
                    reg_centr = [x,y];
                end
                
                if savemat
                    if verbose
                        if isempty(bndry_data) 
                            warning('Detected that bndry data is empty. Will not save it. ***\n')
                        else
                            fprintf('Saving bndry data into %s ... ***\n', bndry_filename)
                        end
                    end
                    save(bndry_filename, 'bndry_data', 'reg_name', 'reg_id', 'reg_centr', '-v7.3');
                end

                if savenc
                    ncfilename = strrep(bndry_filename, '.mat', '.nc');
                    fprintf('*** Saving bndary vertices (lon,lat) into NetCDF file: %s *** \n', ncfilename);
                    
                    data              = bndry_data{:};  %/ Transpose it since the original DRM_src_region is in (lat lon)!
                    data_shortname    = 'vertices';
                    data_standardname = 'vertices';
                    data_units        = 'degree';
                    data_longname     = 'A list of vertices in (lon,lat) outlining the boundary of the region';

                    othervars = [];
                    othervars.reg_name = reg_name;

                    write_nc('ncfilename', ncfilename, 'data', data, 'data_shortname', data_shortname,...
                             'data_standardname', data_standardname, 'data_longname', data_longname, 'data_units', data_units,...
                             'lon', [], 'lat', [], 'plev', [], 'time', [], 'time_unit', [], ...
                             'date', [], 'date_format', [], 'remark', nc_remark, 'othervars', othervars);
                    ncdisp(ncfilename);
                    % ncread(ncfilename, 'vertices') 
                    % ncread(ncfilename, 'reg_name')
                end
            end
        end
        reg_bndry_list{k}    = bndry_data;
        reg_name_list{k}     = reg_name;
        reg_id_list{k}       = reg_id;
        reg_centr_list{k}    = reg_centr;
    end
    
    reg_bndry_list    = cat(1, reg_bndry_list{:}); %/ concatenate the nested cells -> cell or numeric
    if iscell(reg_name)
        reg_name_list = cat(1, reg_name_list{:});  %/ concatenate the nested cells (if any) 
    end
    reg_id_list       = cat(1, reg_id_list{:});    %/ concatenate the nested cells -> numeric
    reg_centr_list    = cat(1, reg_centr_list{:}); %/ concatenate the nested cells -> numeric
    
    %/ If only one region, then the boundary list will not be in cell
    %/ format. Correct this case for coding consistency.
    if ~iscell(reg_bndry_list)    reg_bndry_list = {reg_bndry_list};     end
    if ~iscell(reg_name_list)     reg_name_list  = {reg_name_list};      end
    
    %/ Check if there is repetitive basins
    if length(unique(reg_id_list)) ~= length(reg_id_list)
%         [B, I] = unique(reg_id_list)
        error('There exists repetitive basin ids! Check your code!');
    end
    
    %/ Loop over reg_bndry_list to produce 2D matrix with indices ('reg_2D')
    for i = 1:length(reg_name_list)
        if ismember(reg_name_list{i}, labels_regime)   %/ As we have already modified reg_2D for labels_regime
            continue;

        elseif ismember(reg_name_list{i}, labels_IPCC_PAK) && flag_bndry2reg2D == 0 %/ If flag_bndry2reg2D ==0, we can use the mask for IPCC-PAK regions to output reg_2D right away.
            str = reg_name_list{i};
            queried_region = strrep(str, 'IPCC-PAK-', ''); %/ remove substring
            ind = findismember_loop(mask_region_abbrevs, queried_region);
            if isempty(ind)
                error('The quried region %s is not found from the mask data!', queried_region);
            end
            logical_2D = logical(mask(:,:,ind));
            reg_2D(logical_2D) = i;

        else
            bndry_data = reg_bndry_list{i};
            if ismember(reg_name_list(i), {'land', 'ocean'})
                if ismember(reg_name_list(i), {'land'})
                    logical_2D = cond_land;
                elseif ismember(reg_name_list(i), {'ocean'})
                    logical_2D = cond_ocean; 
                end
            else
                cond_unassigned_land  = (reg_2D == value_land);
                cond_unassigned_ocean = (reg_2D == value_ocean);
                try
                    [mask_2D, ~] = inpoly2([lon_bc_2Dto1D, lat_bc_2Dto1D], bndry_data);                  %/ inpoly2 is 600xx faster than inpolygon!! 
                    mask_2D = reshape(mask_2D, nlon, nlat);
                catch
                    % if verbose warning(E.identifier, 'Error caught!: %s. Proceed with inpolygon().', E.message);  end
                    [mask_2D, ~] = inpolygon(lon_bc_2D, lat_bc_2D, bndry_data(:,1), bndry_data(:,2));    %/ sometimes inpoly2 could run into a bug, then we use inpolygon.
                end
                
                %/ [IMPORTANT]: Distinguish those take land/ocean only from those take both grids.
                if ismember(reg_name_list(i), labels_ocean)
                    logical_2D = (cond_unassigned_ocean & mask_2D);
                elseif ismember(reg_name_list(i),  labels_land)
                    logical_2D = (cond_unassigned_land & mask_2D);
                else
                    logical_2D = ((cond_unassigned_ocean | cond_unassigned_land) & mask_2D);
                end
            end
    
            %/ Outputs
            reg_2D(logical_2D) = i;
    
            %/ Outline boundaries by grids (based on lon_bc lat_bc) if not hydrosheds. 
            if outline_by_grid && ismember(reg_name_list(i), labels_ocean)
                if flag_2D_lonlat
                    error('''outline_by_grid'' currently does not work for uneven lon lat!');
                end
                %/ create the final logical matrix for outlining boudaries
                hole_mode = 'holes';
                reg_bndry_list{i} = get_bndry_from_logi('logical_2D', logical_2D, 'bndry_lon', lon_bc, 'bndry_lat', lat_bc,...
                                                     'glb_data_mode', 1, 'outputcell', 0, 'hole_mode', hole_mode, 'draw_rings', draw_rings);
            else
                reg_bndry_list{i} = bndry_data;
            end
        end
    end
    
%     %/ Double check if all boudnaries are well-defined [No use, polyshap() will modify the complex boundaries automatically to be well-defined.]
%     for i = 1:length(reg_name_list)
%         fprintf('Checking if boudnaries are well-defined: %s...\n', reg_name_list{i});
%         bndry_data = reg_bndry_list{i};
%         a = polyshape(bndry_data);    
%         if issimplified(a) == 0
%             error('The bndry data for region ''%s'' seems not well defined! Check again please.', reg_name_list{i});
%         end
%     end
    
    %/ Double check if all ocean/land/TP grids are assigned
    if contains(slct_reg, 'oceans') 
        if ~isempty(find(reg_2D == value_ocean, 1)) %/ check if all oceanic grid cells are assigned.
            [ind_lon, ind_lat] = find(reg_2D == value_ocean);
            if flag_2D_lonlat
                unassigned_lon_lat_grids = [lon_bc(ind_lon,ind_lat), lat_bc(ind_lon,ind_lat)];
            else
                unassigned_lon_lat_grids = [lon_bc(ind_lon), lat_bc(ind_lat)];
            end
    %         disp(unassigned_lon_lat_grids);
            if verbose 
                warning('There are %d oceanic grid cells left unassigned, given res == %.2f! Check your code!', length(ind_lon), res); 
            end
        end
    end
    
    if contains(slct_reg, 'hydrosheds') 
        if ~isempty(find(reg_2D == value_land, 1)) %/ check if all oceanic grid cells are assigned.
            [ind_lon, ind_lat] = find(reg_2D == value_land);
            if flag_2D_lonlat
                unassigned_lon_lat_grids = [lon_bc(ind_lon,ind_lat), lat_bc(ind_lon,ind_lat)];
            else
                unassigned_lon_lat_grids = [lon_bc(ind_lon), lat_bc(ind_lat)];
            end
    %         disp(unassigned_lon_lat_grids);
            if verbose 
                warning('There are %d land grid cells left unassigned, given res == %.2f! Check your code!', length(ind_lon), res); 
            end
        end
    end
    
    if any(contains(slct_reg, 'TP_basins'))
        cond_unassigned_TP = (cond_TP & (reg_2D == value_land));
        if any(reshape(cond_unassigned_TP, [], 1)) %/ check if all oceanic grid cells are assigned.
            [ind_lon, ind_lat] = find(cond_unassigned_TP);
            if flag_2D_lonlat
                unassigned_lon_lat_grids = [lon_bc(ind_lon,ind_lat), lat_bc(ind_lon,ind_lat)];
            else
                unassigned_lon_lat_grids = [lon_bc(ind_lon), lat_bc(ind_lat)];
            end
%             disp(unassigned_lon_lat_grids);
            if verbose 
                warning('There are %d TP grid cells left unassigned. See if you need to check your code!', length(ind_lon)); 
            end
        end
    end
    
    %/ Restore unassigned land / ocean grids to NaN
    reg_2D(reg_2D == value_land)  = NaN;
    reg_2D(reg_2D == value_ocean) = NaN;
    
    if verbose 
        [ind_lon, ~] = find(isnan(reg_2D));
        if ~isempty(ind_lon)
            warning('[reg_extractor] Detected NaN in reg_2D! Check if this is expected!')
            % disp([ind_lon(1:3), ind_lat(1:3)])
        end
    end

    %/ If outputing reg_2D in logical
    if output_logical
        reg_2D(~isnan(reg_2D)) = 1;
        reg_2D(isnan(reg_2D))  = 0;
        reg_2D = logical(reg_2D);
    end
end

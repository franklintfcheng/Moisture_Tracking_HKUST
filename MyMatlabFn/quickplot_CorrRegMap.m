%%
function [contf_data, hatch_data, cont_data, titlename] = quickplot_CorrRegMap(varargin)

    pnames = {'year', 'lag_year', 'slct_month', 'corr_or_reg', 'onept_or_local', 'CorrType', 'alpha',...
              'predictor', 'predictor_name', 'predictor_dates', 'save_CorrRegMap', 'recompute', 'prcssd_data_folder',...
              'bndry_data', 'contf_dataname', 'hatch_dataname', 'cont_dataname', 'vector_dataname',...
              'pcolor_mode', 'draw_cbar_only', 'cbar_mode', 'cbar_location', 'dataset', 'select_field',...
              'glb_data_mode', 'mask_mountains', 'map_proj', 'map_lon_lower', 'map_lon_upper', 'map_lat_lower', 'map_lat_upper', 'grid_mode',...
              'all_in_one_format', 'fontsize', 'markersize', 'linewi',  'plot_or_not', 'ax_panel', 'plotting_folder', 'shift_title_y', 'savefig'};
    
    dflts  = cell(length(pnames), 1);
    
    [          year,   lag_year,   slct_month,   corr_or_reg,   onept_or_local,   CorrType,   alpha,...
               predictor,   predictor_name,    predictor_dates,   save_CorrRegMap,  recompute,   prcssd_data_folder,...
               bndry_data,   contf_dataname,   hatch_dataname,  cont_dataname,   vector_dataname,...
               pcolor_mode,   draw_cbar_only,    cbar_mode,  cbar_location,   dataset,   select_field,...
               glb_data_mode, mask_mountains, map_proj,  map_lon_lower,   map_lon_upper,   map_lat_lower,   map_lat_upper, grid_mode,...
               all_in_one_format, fontsize,   markersize,   linewi,   plot_or_not, ax_panel,      plotting_folder,  shift_title_y, savefig] = ...
                                            internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    %/ initialization
    contf_data = []; contf_lon = []; contf_lat = []; contf_levels = []; contf_unit = []; colmap = []; cbar_interval = []; cbar_YTick = []; cbar_YTickLabel = [];
    hatch_data = []; hatch_lon = []; hatch_lat = []; hatch_thres_pve  = []; hatch_thres_nve  = []; color_hatch_pve = []; color_hatch_nve = [];   hatch_mode = [];    hatch_linewi  = []; hatch_intvl= [];
    cont_data = []; cont_data_raw = []; cont_lon = []; cont_lat =[]; cont_levels=[]; cont_colmap=[]; cont_linewi=[]; cont_labelsize=[]; cont_label_col=[]; skip_zero_cont = 0;
    Udata=[]; Vdata=[]; uv_lon=[]; uv_lat=[]; vec_step_lon=[]; vec_step_lat=[]; vector_levels = []; vector_color=[]; vector_edgecolor=[]; vecscale=[]; vecscale2=[]; shaftwidth=[];
    headlength=[]; vec_lbs=[]; vec_mag_ref=[]; vec_ref_fontsize=[]; vec_ref_lat_shift=[];
    contf_unitconv = 1;  cont_unitconv = 1;
    
    str_alpha = sprintf('a%.2f', alpha);   
    if corr_or_reg == 1 && onept_or_local == 1      str_plottype = 'OnePtCorrMap'; 
    elseif corr_or_reg == 2 && onept_or_local == 1  str_plottype = 'OnePtRegMap';       end
    
    %/ Compute lag-correlation / lag-regression.
    predictor_dates_dt = int2datetime(predictor_dates, 'yyyyMM');
    ind                = find(predictor_dates_dt.Month == slct_month);
    slct_lag_months    = datetime2int(predictor_dates_dt(ind) + calyears(lag_year), 'yyyyMM'); %/ datetime format is convienient.
    ind                = find(floor(slct_lag_months/1e2) >= year(1) & floor(slct_lag_months/1e2) <= year(end));
    slct_lag_months    = slct_lag_months(ind);
    
    if ismember(select_field, {'monthly_3MA'})
        str_mth_bc = {'DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND','NDJ'};
    elseif ismember(select_field, {'monthly'})
        str_mth_bc = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
    else
        error('You need to adapt the function for the given ''selectfield''!');
    end
    str_lagtime = sprintf('%s(%d)', str_mth_bc{slct_month}, lag_year); %/ Lag-Time arrays
    
    %/ Create var_list to loop from.
    if ~isempty(vector_dataname) vector_dataname_uv = strcat({'u', 'v'}, vector_dataname);
    else                         vector_dataname_uv = {[], []};            end
    var_list = [contf_dataname, cont_dataname, vector_dataname_uv];
    var_list = var_list(~cellfun('isempty',var_list)); %/ remove empty cells
    
    for k = 1:length(var_list)
        
        %/ to save and load the CorrRegMap data -> save time and cpu.
        CorrRegMap_filename     = strcat(prcssd_data_folder, str_plottype, '_',    str_alpha, '_', predictor_name, '_', var_list{k}, '_', str_lagtime, '.mat');
        CorrRegMap_sig_filename = strcat(prcssd_data_folder, str_plottype, '_sig_', str_alpha, '_', predictor_name, '_', var_list{k}, '_', str_lagtime, '.mat');
        
        if recompute == 0 && isfile(CorrRegMap_filename) && isfile(CorrRegMap_sig_filename)
            fprintf('*** CorrRegMap / CorrRegMap_sig data has been computed! Now loading from it...\n');
            load(CorrRegMap_filename,     'CorrRegMap');
            load(CorrRegMap_sig_filename, 'CorrRegMap_sig');
        else
            predictand     = dataset.(var_list{k}).(select_field);
            ind_lag_months = findismember_loop(dataset.(var_list{k}).date_yyyymm_AllYr, slct_lag_months);
            CorrdataX      = predictand(:,:,ind_lag_months);
            
            %/ double-check
            a = dataset.(var_list{k}).date_yyyymm_AllYr(ind_lag_months);
            disp(a)
            disp(length(a));
            
            if lag_year <= 0          %/ if lag_year = -1, predictor will be from 1980, predictand from 1979.
                CorrdataY = predictor(1-lag_year:end);     %/ CorrdataY -> predictor in the regression. If lag_year = -1, n = length(year)-1, and so on.
            else                      %/ if lag_year = 1, predictor will be from 1979, predictand from 1980.
                CorrdataY = predictor(1:end-lag_year);
            end
            
%             if length(ind_lag_months) == length(slct_lag_months) %/ that means predictand has data for all dates even at the time-lag.
%                 disp('yes')
%                 CorrdataX = CorrdataX(:,:,1:end+lag_year);       %/ remove the last one to keep the time length consistent.
%             end
            
            NumWorkers = 20;
    %         disp(size(CorrdataY))
    %         disp(size(CorrdataX))
            if corr_or_reg == 1 && onept_or_local == 1
                [CorrRegMap, CorrRegMap_sig, ~, ~, ~, ~, ~] = OnePtCorrMap('CorrdataX', CorrdataX, 'CorrdataY', CorrdataY,...
                                                                     'corr_or_reg', corr_or_reg, 'onept_or_local',onept_or_local,...
                                                                     'CorrType', CorrType, 'alpha', alpha, 'NumWorkers', NumWorkers);
                
            elseif corr_or_reg == 2 && onept_or_local == 1
                [~, ~, CorrRegMap, CorrRegMap_sig, ~, ~, ~] = OnePtCorrMap('CorrdataX', CorrdataX, 'CorrdataY', CorrdataY,...
                                                                     'corr_or_reg', corr_or_reg, 'onept_or_local',onept_or_local,...
                                                                     'CorrType', CorrType, 'alpha', alpha, 'NumWorkers', NumWorkers);
               
            end
                    
            if save_CorrRegMap
                fprintf('*** Saving into %s... ***\n', CorrRegMap_filename);
                save(CorrRegMap_filename, 'CorrRegMap');
                
                fprintf('*** Saving into %s... ***\n', CorrRegMap_sig_filename);
                save(CorrRegMap_sig_filename, 'CorrRegMap_sig');
            end
        end
        
        %/ assign to suitable variables for plotting.
        if ismember(var_list{k}, {contf_dataname})
            contf_data = CorrRegMap;
            hatch_data = CorrRegMap_sig;
        end

        if ismember(var_list{k}, {cont_dataname})
            cont_data = CorrRegMap_sig;
        end

        if ~isempty(vector_dataname)
            if ismember(var_list{k}, vector_dataname_uv(1))
                Udata = CorrRegMap_sig;
            end
            if ismember(var_list{k}, vector_dataname_uv(2))
                Vdata = CorrRegMap_sig;
            end
        end
    end
    
    %/ NOTE: if multiple casedates, then we take the corresp. time mean of data         
    if isempty(glb_data_mode)  glb_data_mode = 0;     end
    if isempty(fontsize)       fontsize    = 20;      end
    if isempty(markersize)     markersize  = 3;       end
    if isempty(linewi)         linewi      = 2.5;     end
    if isempty(pcolor_mode)    pcolor_mode = 0;       end
    if isempty(grid_mode)      grid_mode   = 2;       end
    if isempty(plot_or_not)    plot_or_not = 1;       end
    if mask_mountains          glb_plateau_mode = 0;  else  glb_plateau_mode = 1;  end  %/ if masking mountains, then no need to draw contour.
    if isempty(draw_cbar_only) draw_cbar_only = 0;               end
    if isempty(cbar_mode)      cbar_mode     = 1;                end 
    if isempty(cbar_location)  cbar_location = 'eastoutside';    end
        
    %/ parameters
    plateau_hgt      = 1500;
    plateau_col      = [255 51 204]./255;
    coast_wi         = 2.5;
    coast_col        = [.4 .4 .4];
    backcolor        = 'none';      %/ [.7 .7 .7]
    marker = 'none'; markerfacecolor = 'k'; color = 'k'; 
    
    %====== contf ======%
    if ~isempty(contf_dataname)
        contf_unit      = '';
        contf_lon       = dataset.(contf_dataname).lon;
        contf_lat       = dataset.(contf_dataname).lat;
        cbar_interval   = 2;
        cbar_YTick      = [];
        cbar_YTickLabel = [];
        
        size(contf_data)
        size(contf_lon)
        size(contf_lat)
            
        if corr_or_reg == 1 
            contf_levels   = -0.8:0.1:0.8;
            NoOfColor      = length(contf_levels)-1;
            colmap         = brewermap(NoOfColor, '*RdYlBu');
            colmap(NoOfColor/2:NoOfColor/2+1, :)  = [1 1 1; 1 1 1];
            
        elseif corr_or_reg == 2
            if ismember(contf_dataname, {'prcp'})
                contf_levels  = [-9:1:9];               %/ in mm/day
                NoOfColor     = length(contf_levels)-1;
                colmap          = brewermap(NoOfColor, '*BrBG');
                colmap(NoOfColor/2:NoOfColor/2+1,:) = [1 1 1; 1 1 1];

            elseif ismember(contf_dataname, {'pw'})     %/ pw == tcwv
                contf_levels  = [-21:3:21];               %/ in mm (instantaneous)
                NoOfColor     = length(contf_levels)-1;
                colmap          = brewermap(NoOfColor, '*BrBG');
                colmap(NoOfColor/2:NoOfColor/2+1,:) = [1 1 1; 1 1 1];

            elseif ismember(contf_dataname, {'omega500'})    
                contf_levels  = [-110:10:110];               %/ in hPa/day 
                colmap          = brewermap(length(contf_levels)-1, 'RdBu'); %/ looks better
                colmap((length(contf_levels)-1)/2:(length(contf_levels)-1)/2+1,:) = [1 1 1; 1 1 1;];

            elseif ismember(contf_dataname, {'IVT'})    
                contf_levels = [0:50:500];               %/ in kg/m/s 
                colmap       = jet(length(contf_levels)-1);

            elseif isequal(contf_dataname, 'Z850')
                contf_levels  = [-3:.3:3];
                colmap        = brewermap(length(contf_levels)-1, 'RdBu'); %/ looks better
                
            elseif ismember(contf_dataname, {'OLR'})    
                contf_levels  = [-50:5:50];               %/ in W m^-2 
                NoOfColor     = length(contf_levels)-1;
                colmap          = brewermap(NoOfColor, '*BrBG');
                colmap(NoOfColor/2:NoOfColor/2+1,:) = [1 1 1; 1 1 1];

            elseif ismember({contf_dataname}, {'Q1_vi', 'Q2_vi'})
                contf_levels  = [-440:40:440];              %/ in W m^-2 
                NoOfColors    = length(contf_levels)-1;
                colmap = brewermap(NoOfColors, '*RdYlBu');
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [ 1 1 1; 1 1 1];  

            elseif ismember(contf_dataname, {'SSHF', 'SLHF'})    
                contf_levels  = [-220:20:220];               %/ in W m^-2 
                colmap        = brewermap(length(contf_levels)-1, '*PiYG');
    %             colmap        = cmocean('balance', length(contf_levels)-1);
                NoOfColors    = length(contf_levels)-1;
                colmap(NoOfColors/2:NoOfColors/2+1,:) = [1 1 1; 1 1 1];

            elseif ismember(contf_dataname, {'T2m'})     %/ no need to subtract 273 from it when showing slope coef.
                contf_levels  = [-9:1:9];
                NoOfColor      = length(contf_levels)-1;
                colmap         = brewermap(NoOfColor, '*RdYlBu');
                colmap(NoOfColor/2:NoOfColor/2+1, :)  = [1 1 1; 1 1 1];

            elseif ismember(contf_dataname, {'SST'})     %/ no need to subtract 273 from it when showing slope coef.
                contf_levels  = [-4.5:0.5:4.5]/100;
                NoOfColor      = length(contf_levels)-1;
                colmap         = brewermap(NoOfColor, '*RdYlBu');
                colmap(NoOfColor/2:NoOfColor/2+1, :)  = [1 1 1; 1 1 1];
                
            elseif ismember(contf_dataname, {'psi200'})
                %/ NOTE:
                %/ In geostrophic balance (See Holton or notes in Notion):
                %/   High (low) stream function -> high (low) geopotential height in N.H.
                %/                              -> low (high) geopotential height in S.H.
                %/ https://www.cpc.ncep.noaa.gov/products/CDB/Tropics/figt22.shtml

                contf_unitconv = 10^-6;
                contf_unit     = '10^{6} m^{2} s^{-1}';
                contf_levels   = [-20:2:20]; %/ This makes sure level interval to be const. 
                NoOfColor      = length(contf_levels)-1; 
                colmap         = brewermap(NoOfColor,'*RdBu');
                colmap(NoOfColor/2:NoOfColor/2+1,:) = [1 1 1; 1 1 1;];

            elseif ismember(contf_dataname, {'S','S1','S2'})
                contf_levels   = [-50:5:50];  %/ +ve: source, -ve: sink.
                contf_unitconv = 10^11; %/ scaling
                contf_unit     = '10^{-11} s^{-2}';
                NoOfColor      = length(contf_levels)-1; 
                colmap         = brewermap(NoOfColor, '*RdGy');
                colmap(NoOfColor/2:NoOfColor/2+1,:) = [1 1 1; 1 1 1;];

            else
                error('The levels of the request data has not been defined!')
            end
        end
    end
    
    %====== hatching =====%
    if ~isempty(hatch_dataname) || ~isempty(hatch_data)
        hatch_mode = 1;
        
        if ~isempty(hatch_data)   
            hatch_lon       = dataset.(contf_dataname).lon;
            hatch_lat       = dataset.(contf_dataname).lat;
            hatch_thres_pve = 0;
            hatch_thres_nve = 0;
            color_hatch_pve = [0 0 0]; 
            color_hatch_nve = [0 0 0]; 
            hatch_linewi    = 3;
            hatch_intvl     = 10;
            
        else
            if isequal(casedate_fmt, 'yyyymmddHHMM')        
                X_dates     = datetime2int('X_datetime', dataset.(hatch_dataname).date_UTC_AllYr, 'format', 'yyyymmddHHMM');
                ind_date    = findismember_loop(X_dates, casedates);

            elseif isequal(casedate_fmt, 'yyyymmdd')        
                X_dates     = dataset.(hatch_dataname).date_yyyymmdd_AllYr; 
                ind_date    = findismember_loop(X_dates, casedates);

            elseif isequal(casedate_fmt, 'index')          
                ind_date    = casedates;  %/ assume casedates store pentad.
            end
            
            hatch_data      = my_nanmean(dataset.(hatch_dataname).(select_field)(:,:,ind_date), 3);
            hatch_lon       = dataset.(hatch_dataname).lon;
            hatch_lat       = dataset.(hatch_dataname).lat;
            hatch_thres_pve = 0;
            hatch_thres_nve = 0;
            color_hatch_pve = [0 0 0]; 
            color_hatch_nve = [0 0 0]; 
            hatch_linewi    = 1.5;
            hatch_intvl     = 14;
        end
    end
    
    %====== cont ======%
    if ~isempty(cont_dataname)
        fprintf('*** Mean cont data = %.2f ***\n', nanmean(cont_data, 'all'));
        cont_lon       = dataset.(cont_dataname).lon;
        cont_lat       = dataset.(cont_dataname).lat;
        cont_linewi    = 3;
        cont_labelsize = fontsize*0.8;
        cont_label_col = 'k';
        skip_zero_cont = 1;   %/ whether to skip the zero-contour.
        
        if corr_or_reg == 1 
            cont_levels = [-0.8:0.05:0.8];        
            cont_colmap = create_div_colmap([50 102 255]./255, [255 51 102]./255, cont_levels);
%             cont_colmap = create_div_colmap([103 34 160]./255, [255 102 51]./255, cont_levels);
            
        elseif corr_or_reg == 2
            if isequal(cont_dataname, 'Z850')
                cont_levels = [-3:0.3:3];
                cont_colmap = create_div_colmap([50 102 255]./255, [255 51 102]./255, cont_levels);
                
            elseif ismember({cont_dataname}, {'dMSE_dt_vi', 'MSE_Res'})
                cont_levels      = [-55:5:55];             %/ in W m^-2 
                skip_zero_cont   = 0;
                
            elseif ismember(cont_dataname, {'psi200'})
                %/ NOTE:
                %/ In geostrophic balance (See Holton or My Notion):
                %/   High (low) stream function -> high (low) geopotential height in N.H.
                %/                              -> low (high) geopotential height in S.H.
                %/ https://www.cpc.ncep.noaa.gov/products/CDB/Tropics/figt22.shtml
                cont_levels   = [-100:10:100];   %/ x 10^6
                cont_unitconv = 10^-6;
                
            elseif isequal(cont_dataname, 'Z200')
                cont_levels   = [-200:20:200];   
                cont_colmap = create_div_colmap([50 102 255]./255, [255 51 102]./255, cont_levels);
                
            elseif isequal(cont_dataname, 'slp')       
                cont_levels    = [1000:2:1020];       %/ in hPa
                cont_colmap    = jet(length(cont_levels));
                
            elseif isequal(cont_dataname, {'T2m', 'SST'})  %/ no need to substract 273 from it for slope coef.
                cont_levels    = [-10:2:30];
                cont_colmap = create_div_colmap([50 102 255]./255, [255 51 102]./255, cont_levels);
                
            elseif any(ismember(cont_dataname, {'uWS_LL','uWS_UL','vWS_LL','vWS_UL'}))
                cont_levels    = [-50:5:50];   %/ in m/s
                cont_colmap = create_div_colmap([50 102 255]./255, [255 51 102]./255, cont_levels);
                
            elseif isequal(cont_dataname, 'omega500') 
                cont_levels    = [flip([40:40:280])*-1, [40:40:280]];   %/ in hPa/day
                cont_colmap = create_div_colmap([50 102 255]./255, [255 51 102]./255, cont_levels);
                
            elseif ismember(cont_dataname, {'OLR'})
                cont_levels  = [-100:20:100];               %/ in W m^-2 
                cont_colmap = create_div_colmap([50 102 255]./255, [255 51 102]./255, cont_levels);
                
            elseif ismember({cont_dataname}, {'MSE_vi'})
                cont_data = cont_data./1e9;  
                cont_levels = linspace(2.925, 3.2, 12);             %/ in 10^9 J m^-2 
            else
                error('The levels of the request data has not been defined!')
            end
        end
    end
    
    %====== vector ======%
    if ~isempty(vector_dataname)
        uv_lon = dataset.(vector_dataname_uv{1}).lon;
        uv_lat = dataset.(vector_dataname_uv{1}).lat;
        
        if corr_or_reg == 1 
            vector_levels     = [0:0.1:1];   
            vector_color      = jet(length(vector_levels));
            vector_edgecolor  = 'k';
            vec_step_lon      = 5;
            vec_step_lat      = 3;
            vecscale          = 4;                       % the smaller value the bigger vector. for winds
            vecscale2         = 2;                       % control shaft length.
            shaftwidth        =  2.2;                    % control width of the shaft, the larger value the thicker
            headlength        =  shaftwidth*5;           % control length of the head, the larger value the longer
            vec_mag_ref       = 0.5;
            vec_lbs           = strcat(num2str(vec_mag_ref)); %/ vector reference
            
        elseif corr_or_reg == 2
            if ismember(vector_dataname, {'IVT'})
                vector_color      = [.2 .2 .2];
                vector_edgecolor  = 'w';
                vec_step_lon      = 5;
                vec_step_lat      = 3;
                vecscale          = 10;    % the smaller value the bigger vector. for winds
                vecscale2         = 2;           % control shaft length.
                shaftwidth        = 3.5;    % control width of the shaft, the larger value the thicker
                headlength        = shaftwidth*4; % control length of the head, the larger value the longer
                vec_ref_lat_shift = 0.08;   %/ the smaller, the less the ref vector shift to the south
                vec_mag_ref       = 500;
                vec_lbs           = strcat(num2str(vec_mag_ref), {' kg/m/s'});  %/ vector reference
                
            elseif ismember(vector_dataname, {'850'})
                vector_color      = [102 255 102]./255; %/ dark green
                vector_edgecolor  = 'k';    
                vec_step_lon      = 5;              %/ important difference may only be shown with high density vectors.
                vec_step_lat      = 3;   
                vecscale          = 0.5;            % the smaller value the bigger vector. for winds
                vecscale2         = 3;    
                shaftwidth        = 2;              % control width of the shaft, the larger value the thicker
                headlength        = shaftwidth*4;   % control length of the head, the larger value the longer  
                vec_ref_lat_shift = 0.08;           % the smaller, the less the ref vector shift to the south
                vec_mag_ref       = 10;
                vec_lbs           = strcat(num2str(vec_mag_ref), {' m/s'});  %/ vector reference
                
            elseif ismember(vector_dataname, {'200', '200chi'})
                vector_color      = [.2 .2 .2];
                vector_edgecolor  = 'w';    
                vec_step_lon      = 4;           %/ important difference may only be shown with high density vectors.
                vec_step_lat      = 2; 
                vecscale          = 1;             % the smaller value the bigger vector. for winds
                vecscale2         = 2;      
                shaftwidth        = 3.5;            % control width of the shaft, the larger value the thicker
                headlength        = shaftwidth*4;   % control length of the head, the larger value the longer  
                vec_ref_lat_shift = 0.08;           % the smaller, the less the ref vector shift to the south
                vec_mag_ref       = 10;
                vec_lbs           = strcat(num2str(vec_mag_ref), {' m/s'});  %/ vector reference
                
            elseif ismember(vector_dataname, {'WAF200'})  %/ order of magnitude 10^2 m^2 s^-1
                vector_color      = [147 255 50]./255; 
                vector_edgecolor  = 'w';    
                vec_step_lon      = 4;           %/ important difference may only be shown with high density vectors.
                vec_step_lat      = 2;  
                vecscale          = 10;             % the smaller value the bigger vector. for winds
                vecscale2         = 2;  
                shaftwidth        = 3.5;            % control width of the shaft, the larger value the thicker
                headlength        = shaftwidth*4;   % control length of the head, the larger value the longer  
                vec_ref_lat_shift = 0.08;           % the smaller, the less the ref vector shift to the south
                vec_mag_ref       = 10;
                vec_lbs           = strcat(num2str(vec_mag_ref), {' m**2 s**-1'});  %/ vector reference
            end
        end
    end
    
    %====== to plot or not to plot ======%
    if plot_or_not
        %/ axis/title/figname
        
        titlename = strcat(str_plottype, {' '}, str_alpha, {' '}, predictor_name, {' '},...
                            strjoin(var_list, ' '), {' '}, str_lagtime);
        titlename = strrep(titlename, '_', ' ');
        
%         titlename = strcat(titlename, str_var); %/ append vars to titlename.
    %     titlename = strcat(titlename, str_var, str_cbar, str_vecref); %/ append vars to titlename.
        
    %     th = title(titlename);
    %     titlePos = get( th , 'position');
    %     titlePos(2) = titlePos(2) + 0.2;
    %     set(th , 'Fontsize', fontsize - 10, 'position' , titlePos);
    %     drawnow; pause(0.05);
        
        if all_in_one_format == 1  %/ tentative formatting for a plot with 4x2 panels
            fontsize       = fontsize*1.5;
            markersize     = markersize*1.5;
            linewi         = linewi*1.5;
            coast_wi       = coast_wi*1.5;
            cont_linewi    = cont_linewi*1.5;
            cont_labelsize = cont_labelsize*1.5;
            
%             vec_step_lon   = 1;
%             vec_step_lat   = 1;
%             vecscale       = vecscale /2;         % the smaller value the bigger vector. for winds
%             vecscale2      = vecscale2*0.5;       % control shaft length.
            shaftwidth   = shaftwidth*1.5;    % control width of the shaft, the larger value the thicker
            headlength   = headlength*1.5; % control length of the head, the larger value the longer  
        end
        
        plot_contfmap('contf_data', contf_data*contf_unitconv, 'contf_lon', contf_lon, 'contf_lat', contf_lat, 'contf_levels', contf_levels,...
                      'contf_unit', contf_unit, 'colmap', colmap, 'cbar_interval', cbar_interval, 'pcolor_mode', pcolor_mode,...
                      'cont_data',  cont_data*cont_unitconv,  'cont_data_raw', cont_data_raw*cont_unitconv, 'cont_lon', cont_lon, 'cont_lat', cont_lat, 'cont_levels', cont_levels, 'cont_colmap', cont_colmap, 'cont_linewi', cont_linewi, 'cont_labelsize', cont_labelsize, 'cont_label_col', cont_label_col, 'skip_zero_cont', skip_zero_cont,...
                      'Udata', Udata, 'Vdata', Vdata, 'uv_lon', uv_lon, 'uv_lat', uv_lat, 'vec_step_lon', vec_step_lon, 'vec_step_lat', vec_step_lat,...
                      'vector_levels', vector_levels, 'vector_color', vector_color, 'vector_edgecolor', vector_edgecolor, 'vecscale', vecscale, 'vecscale2', vecscale2, 'shaftwidth', shaftwidth, 'headlength', headlength, 'vec_lbs', vec_lbs, 'vec_mag_ref', vec_mag_ref, 'vec_ref_fontsize', vec_ref_fontsize, 'vec_ref_lat_shift', vec_ref_lat_shift,...
                      'hatch_data', hatch_data, 'hatch_lon', hatch_lon, 'hatch_lat', hatch_lat, 'hatch_thres_pve', hatch_thres_pve, 'hatch_thres_nve', hatch_thres_nve,...
                      'color_hatch_pve', color_hatch_pve, 'color_hatch_nve', color_hatch_nve, 'hatch_mode', hatch_mode, 'hatch_linewi', hatch_linewi, 'hatch_intvl', hatch_intvl,...
                      'bndry_data', bndry_data, 'marker', marker, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'linewi', linewi, 'color', color,...
                      'titlename', titlename, 'shift_title_y', shift_title_y, 'savepath', [], ...
                      'glb_data_mode', glb_data_mode, 'mask_mountains', mask_mountains, 'glb_plateau_mode', glb_plateau_mode, 'plateau_hgt', plateau_hgt, 'plateau_col', plateau_col,...
                      'map_proj', map_proj, 'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper, 'coast_col', coast_col, 'coast_wi', coast_wi, 'backcolor', backcolor,...
                      'fontsize', fontsize,  'create_fig', 1, 'grid_mode', grid_mode, 'draw_cbar_only', draw_cbar_only, 'cbar_mode', cbar_mode, 'cbar_location', cbar_location, 'cbar_YTick', cbar_YTick, 'cbar_YTickLabel', cbar_YTickLabel,...
                      'ax_panel', ax_panel)
        fprintf('done \n')
        
        if savefig
    %         export_fig(char(strcat(plotting_folder, titlename,'.pdf')),'-pdf','-painters', '-nocrop'); %'-transparent');
            export_fig(char(strcat(plotting_folder, titlename,'.pdf')),'-pdf','-painters', '-c[inf, inf, inf, inf]'); %, '-nocrop'); %'-transparent');
        end
    end
end
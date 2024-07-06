function [dataset] = MSE_Q1_Q2_wrapper(varargin)
    
    % create a set of valid parameters and their default value
    pnames = {'dataset', 'select_field', 'stlevel',  'edlevel', 'year', 'prcssd_data_folder', 'recompute', 'save_4D_data', 'save_mat_or_nc'};
    dflts  = cell(length(pnames), 1);
    
    %/ parse function arguments
    [         dataset,    select_field,    stlevel,    edlevel,   year,    prcssd_data_folder, recompute,   save_4D_data, save_mat_or_nc] ...
           = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/==========================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 26 Apr 2024
    %/==========================================

    fprintf('*** Running MSE_Q1_Q2_wrapper... ***\n')
    st_mth = 1; ed_mth = 12;  %/ Extract data for all dates (by default)
    str_years = sprintf('%d-%d', year(1), year(end));
    
    if isequal(stlevel, 100) && isequal(edlevel, 1000)
        str_level = '';   
    else
        str_level = sprintf('%dto%d', stlevel, edlevel); 
    end
    
    date_append = []; 
    cnt = 0;
    for t = 1:length(year)
        cnt    = cnt + 1;
        stDate = datetime(year(t), st_mth, 1,                       'Format', 'yyyy-MM-dd');
        edDate = datetime(year(t), ed_mth, eomday(year(t), ed_mth), 'Format', 'yyyy-MM-dd');
        
        casedates_annual = stDate:edDate;
        casedates_annual = (casedates_annual.Year*1e4 + casedates_annual.Month*100 + casedates_annual.Day)';
        casedates        = casedates_annual;
%         casedates = cat(1, casedates, casedates_annual);
        
        %/ by default, dataname of 'q' means a q contains 100-1000 hPa.
        T_dataname = 'T';  %/ ERA5: T (K)
        q_dataname = 'q';  %/ ERA5: q (kg/kg)
        Z_dataname = 'Z';  %/ ERA5: Z (m)       Should've been divided by 9.81 in data pre-processing.
        U_dataname = 'U';  %/ ERA5: U (m/s)
        V_dataname = 'V';  %/ ERA5: V (m/s)
        W_dataname = 'W';  %/ ERA5: Omega (Pa/s)
        
        %/ Load data
        str_years_bc  = sprintf('%d', year(t));   %/ to correctly read the pre-processed annual mat data
        always_reload = 1;                        %/ since we loop over years, need to always reload.
        
        dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', T_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
        dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', q_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
        dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', Z_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
        dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', U_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
        dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', V_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
        dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', W_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
        
        T = dataset.(T_dataname).(select_field);
        q = dataset.(q_dataname).(select_field);
        Z = dataset.(Z_dataname).(select_field);
        U = dataset.(U_dataname).(select_field);
        V = dataset.(V_dataname).(select_field);
        W = dataset.(W_dataname).(select_field);
        
        %===== Compute MSE (m) (requisite: T, Z, q) ======%
        r   = 6371e3;
        lon = dataset.(q_dataname).lon;
        lat = dataset.(q_dataname).lat;
        hx  = diff(lon(1:2))*pi/180;   %/ rmb to convert degree into radian!!!!
        hy  = diff(lat(1:2))*pi/180;   %/ rmb to convert degree into radian!!!! NOTE: hy can be -ve. <- correct.
        
        %/ Get pressure data (hPa)
        P          = dataset.(q_dataname).level;
        level_unit = 'hPa';
        
        %/ Compute the 4D MSE
        MSE = MoistStaticEnergy('T', T, 'Z', Z, 'q', q);
%         size(MSE)
        
        %/ Compute Q1 Q2 [J kg-1 s-1]
        [Q1, Q2] = Yanai_Q1Q2('T', T, 'q', q, 'U', U, 'V', V, 'W', W, 'lon', lon, 'lat', lat, 'P', P);
        
        %  Notes about gradient()
        %       it computes central          difference for interoir data points.
        %       it computes forward/backward difference for data points at edges
        %
        %       [FX,FY,FZ,...,FN] = gradient(F) returns the N components of the numerical gradient of F, where F is an array with N dimensions.
        %       The first output FX is always the gradient along the *2nd* dimension of F, going across columns.
        %       The second output FY is always the gradient along the *1st* dimension of F, going across rows.
        %       For the third output FZ and the outputs that follow, the Nth output is the gradient along the Nth dimension of F.
        
        dt = 24*3600;                              %/ unit time step (e.g., one day = 24*3600 s)
        [~, ~, ~, dMSE_dt] = gradient(MSE, dt);    %/ J kg-1 s-1

        [~, dMSE_dx, ~, ~] = gradient(MSE, hx);    
        unit_conv          = (1./(r*cosd(lat)))';       %/ 1 x lat
        dMSE_dx            = unit_conv.*dMSE_dx;        %/ J kg-1 m-1.   div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))
        
        [dMSE_dy, ~, ~, ~] = gradient(MSE, hy);         
        unit_conv          = 1./r;                     
        dMSE_dy            = unit_conv.*dMSE_dy;        %/ J kg-1 m-1.   div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))
        
        [~, ~, dMSE_dA, ~] = gradient(MSE);                          %/ Chain rule to do derivate with *non-uniform* spacings.
        dP_dA              = gradient(P * 100);                      %/ get non-uniform gradient of pressure level first.
        dP_dA              = reshape(dP_dA, [1 1 length(dP_dA)]);    %/ 1 x 1 x plev
        dMSE_dP            = dMSE_dA./dP_dA;                         %/ J kg-1 Pa-1   elementwise division.  
        
        dataset.('MSE').(select_field)       = MSE;
        dataset.('dMSE_dt').(select_field)   = dMSE_dt;
        dataset.('U_dMSE_dx').(select_field) = U.*dMSE_dx;
        dataset.('V_dMSE_dy').(select_field) = V.*dMSE_dy;
        dataset.('W_dMSE_dP').(select_field) = W.*dMSE_dP;
        dataset.('Q1').(select_field)        = Q1;
        dataset.('Q2').(select_field)        = Q2;
        
%         a = squeeze(W(90, 30, :, 150))
%         b = squeeze(dMSE_dP(90, 30, :, 150))
%         c = squeeze(MSE(90, 30, :, 150))
%         d = squeeze(dataset.('W_dMSE_dP').(select_field)(90, 30, :, 150))
%         isequal( a.*b , d)
%         all = [plev', a, b, c, d]
        
        %/ Mass-weighted vertical integral (between stlevel and edlevel)
        ind_P    = find(stlevel <= P & P <= edlevel);  %/ This is correct, otherwise we'll miss the 650-700hPa layer
        P_subset = P(ind_P);
        
        MSE_4D_namelist = {'MSE', 'dMSE_dt', 'U_dMSE_dx', 'V_dMSE_dy', 'W_dMSE_dP', 'Q1', 'Q2'};
        MSE_vi_namelist = strcat(MSE_4D_namelist, '_vi', str_level);
        for k = 1:length(MSE_vi_namelist)
            %/ Subset the correct plev dimension from the 4D data.
            data_subset = dataset.(MSE_4D_namelist{k}).(select_field)(:,:,ind_P,:);
            
            %/ NOTE: topo_mode is on.
            X_vi = vertinte('data', data_subset, 'level', P_subset, 'level_unit', level_unit, 'take_vertmean', 0,...
                            'topo_mode', 1, 'topo_lon', lon, 'topo_lat', lat);   
            fprintf('*** mean %s = %.2f ***\n', MSE_vi_namelist{k}, mean(X_vi, 'all'));
            
            if cnt == 1 
                dataset.(MSE_vi_namelist{k}).(select_field) = [];
            end
            dataset.(MSE_vi_namelist{k}).(select_field) = cat(3, dataset.(MSE_vi_namelist{k}).(select_field), X_vi);
        end
        date_append = [date_append; dataset.(q_dataname).date_yyyymmdd_AllYr];
        
        %/ Save data as mat or nc
        if save_mat_or_nc == 1
            if save_4D_data
                %/ Save 4D data year by year
                for i = 1:length(MSE_4D_namelist)
                    filename = strcat(prcssd_data_folder, MSE_4D_namelist{i}, '_', select_field, '_', num2str(year(t)), '.mat');
                    if isfile(filename) && recompute == 0 
                        disp('Data has been processed. Skip writing.');
                    else
                        dataset.(MSE_4D_namelist{i}).lon                  = dataset.(q_dataname).lon;
                        dataset.(MSE_4D_namelist{i}).lat                  = dataset.(q_dataname).lat;
                        dataset.(MSE_4D_namelist{i}).level                = dataset.(q_dataname).level;
                        dataset.(MSE_4D_namelist{i}).date_yyyymmdd_AllYr  = dataset.(q_dataname).date_yyyymmdd_AllYr;

                        if isequal(MSE_4D_namelist{i}, 'MSE')
                            dataset.(MSE_4D_namelist{i}).unit = 'J kg**-1';                         % <-- mind the unit.
                        else
                            dataset.(MSE_4D_namelist{i}).unit = 'W kg**-1';                         % <-- mind the unit.
                        end

                        S = dataset.(MSE_4D_namelist{i});
                        save(filename, 'S', '-v7.3');
                        disp(['*** Written into ', filename, '... ***'])
                    end
                end
            end
            
            %/ Will always save column-integrated MSE (all years)
            if cnt == length(year)
                for i = 1:length(MSE_vi_namelist)
                    filename = strcat(prcssd_data_folder, MSE_vi_namelist{i}, '_', select_field, '_', str_years, '.mat');
                    
                    if isfile(filename) && recompute == 0 
                        disp('Data has been processed. Skip writing.');
                    else
                        dataset.(MSE_vi_namelist{i}).lon                  = dataset.(q_dataname).lon;
                        dataset.(MSE_vi_namelist{i}).lat                  = dataset.(q_dataname).lat;
                        dataset.(MSE_vi_namelist{i}).level                = P_subset;     %/ mind this!
                        dataset.(MSE_vi_namelist{i}).date_yyyymmdd_AllYr  = date_append;
                        
                        if isequal(MSE_vi_namelist{i}, 'MSE_vi')
                            dataset.(MSE_vi_namelist{i}).unit = 'J m**-2';                         % <-- mind the unit.
                        else
                            dataset.(MSE_vi_namelist{i}).unit = 'W m**-2';                         % <-- mind the unit.
                        end
                        
                        S = dataset.(MSE_vi_namelist{i});
                        save(filename, 'S', '-v7.3');
                        disp(['*** Written into ', filename, '... ***'])
                    end
                end
            end
        elseif save_mat_or_nc == 2
            error('code not ready for save_mat_or_nc == 2');
        end
    end
end

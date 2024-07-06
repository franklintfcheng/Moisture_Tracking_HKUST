%%
function [dataset] = wind_div_wrapper(varargin)
    
    % create a set of valid parameters and their default value
    pnames = {'dataset', 'select_field', 'stlevel',  'edlevel', 'year', 'smooth_method', 'smooth_pt',...
              'clear_UV', 'prcssd_data_folder', 'recompute', 'save_mat_or_nc'};
    dflts  = cell(length(pnames), 1);
    
    %/ parse function arguments
    [         dataset,    select_field,    stlevel,    edlevel,   year,  smooth_method,   smooth_pt,...
              clear_UV,   prcssd_data_folder, recompute,  save_mat_or_nc] ...
           = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    %/ NOTE:
    %/      This function is copied from MSE_Q1_Q2_wrapper(), exclusively
    %  designed for computing layer-mean wind divergence.
    %/      For more complicated derivatives, see MSE_Q1_Q2_wrapper().
    
    fprintf('*** Running wind_div_wrapper... ***\n')
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
        casedates_annual = [casedates_annual.Year*1e4 + casedates_annual.Month*100 + casedates_annual.Day]';
        casedates        = casedates_annual;
%         casedates = cat(1, casedates, casedates_annual);
        
        %/ by default, dataname of 'q' means a q contains 100-1000 hPa.
%         T_dataname = 'T';  %/ ERA5: T (K)
%         q_dataname = 'q';  %/ ERA5: q (kg/kg)
%         Z_dataname = 'Z';  %/ ERA5: Z (m)       Should've been divided by 9.81 in data pre-processing.
        U_dataname = 'U';  %/ ERA5: U (m/s)
        V_dataname = 'V';  %/ ERA5: V (m/s)
%         W_dataname = 'W';  %/ ERA5: Omega (Pa/s)
        
        %/ Load data
        str_years_bc  = sprintf('%d', year(t));   %/ to correctly read the pre-processed annual mat data
        always_reload = 1;                        %/ since we loop over years, need to always reload.
        
%         dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', T_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
%         dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', q_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
%         dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', Z_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
        dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', U_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
        dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', V_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
%         dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', W_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
        
%         T = dataset.(T_dataname).(select_field);
%         q = dataset.(q_dataname).(select_field);
%         Z = dataset.(Z_dataname).(select_field);
        U = dataset.(U_dataname).(select_field);
        V = dataset.(V_dataname).(select_field);
%         W = dataset.(W_dataname).(select_field);
        
        
        %===== Apply a 5-point smoothing on wind fields ====%. (Berry & Reeder 2014; Wodzicki & Rapp 2016) 
        str_sm = '';
        if isequal(smooth_method, 'movmean')
            %/ Use 'convn' to do a 5*5 moving mean smoothing on N-dim data.
            %https://www.mathworks.com/matlabcentral/answers/358411-i-need-a-code-that-produce-a-moving-average-matrix-with-a-5-5-window
            smooth_window = ones(smooth_pt,smooth_pt)/(smooth_pt^2);    
            U = convn(U, smooth_window, 'same'); %/ 'same' -> output the same size of matrix with zero-padded edges.
            V = convn(V, smooth_window, 'same'); %/ 'same' -> output the same size of matrix with zero-padded edges.
            str_sm = '_sm';
            fprintf('*** UV wind fields have been smoothed using a %d-pt %s. ***\n', smooth_pt, smooth_method);
            
        elseif ~isempty(smooth_method)
            error('The input smooth_method is not valid.');
        end

        %===== Compute MSE (m) (requisite: T, Z, q) ======%
        r   = 6371e3;
        lon = dataset.(U_dataname).lon;
        lat = dataset.(U_dataname).lat;
        hx  = diff(lon(1:2))*pi/180;      %/ rmb to change degree to radian!!!!
        hy  = diff(lat(1:2))*pi/180;      %/ rmb to change degree to radian!!!! NOTE: hy can be -ve. <- correct.
        P   = dataset.(U_dataname).level; %/ get the pressure data

        %  Notes about gradient()
        %       it computes central          difference for interoir data points.
        %       it computes forward/backward difference for data points at edges
        %
        %       [FX,FY,FZ,...,FN] = gradient(F) returns the N components of the numerical gradient of F, where F is an array with N dimensions.
        %       The first output FX is always the gradient along the *2nd* dimension of F, going across columns.
        %       The second output FY is always the gradient along the *1st* dimension of F, going across rows.
        %       For the third output FZ and the outputs that follow, the Nth output is the gradient along the Nth dimension of F.
        
        [~, dU_dx, ~, ~] = gradient(U, hx);    
        unit_conv        = (1./(r*cosd(lat)))';     %/ 1 x lat
        dU_dx            = unit_conv.*dU_dx;        %/ s-1.   div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))
        
        [dV_dy, ~, ~, ~] = gradient(V, hy);         
        unit_conv        = 1./r;                     
        dV_dy            = unit_conv.*dV_dy;        %/ s-1.   div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))

        D = dU_dx + dV_dy;   %/ wind divergence (D), 4D
        
        %/ Layer-mean D (between stlevel and edlevel)
        ind_P    = find(stlevel <= P & P <= edlevel);  %/ This is correct, otherwise we'll miss the 650-700hPa layer
        P_subset = P(ind_P);
        
        data_3D_vm_list = strcat({'D_vm'}, str_level, str_sm);
        for k = 1:length(data_3D_vm_list)
            %/ Subset the correct plev dimension from the 4D data.
            data_subset = D(:,:,ind_P,:);

            %/ Use vertinte to perform vertmean in topo_mode.
            % D_vm = vertmean('data', data_subset, 'level', P_subset);
            D_vm = vertinte('data', data_subset, 'level', P_subset, 'topo_mode', 1,...
                            'topo_lon', lon, 'topo_lat', lat, 'take_vertmean', 1);
                        
            fprintf('*** mean abs. %s = %.2g ***\n', data_3D_vm_list{k}, mean(abs(D_vm), 'all'));
            
            if cnt == 1 
                dataset.(data_3D_vm_list{k}).(select_field) = [];
            end
            dataset.(data_3D_vm_list{k}).(select_field) = cat(3, dataset.(data_3D_vm_list{k}).(select_field), D_vm);
        end
        date_append = [date_append; dataset.(U_dataname).date_yyyymmdd_AllYr];
        
        %/ Save as mat or nc
        if save_mat_or_nc == 1         
            %/ 3D data for ALL years (after fully appended).
            if cnt == length(year)
                for i = 1:length(data_3D_vm_list)
                    filename = strcat(prcssd_data_folder, data_3D_vm_list{i}, '_', select_field, '_', str_years, '.mat');
                    
                    if isfile(filename) && recompute == 0 
                        disp('Data has been processed. Skip writing.');
                    else
                        dataset.(data_3D_vm_list{i}).lon                  = dataset.(U_dataname).lon;
                        dataset.(data_3D_vm_list{i}).lat                  = dataset.(U_dataname).lat;
                        dataset.(data_3D_vm_list{i}).level                = P_subset;     %/ mind this!
                        dataset.(data_3D_vm_list{i}).date_yyyymmdd_AllYr  = date_append;
                        dataset.(data_3D_vm_list{i}).unit                 = 's**-1';      % <-- mind the unit.
                        
                        S = dataset.(data_3D_vm_list{i});
                        save(filename, 'S', '-v7.3');
                        disp(['*** Written into ', filename, '... ***'])
                    end
                end
            end
            %/ Save 4D data year by year
%             for i = 1:length(data_4D_namelist)
%                 filename = strcat(prcssd_data_folder, data_4D_namelist{i}, '_', select_field, '_', num2str(year(t)), '.mat');
%                 if isfile(filename) && recompute == 0 
%                     disp('Data has been processed. Skip writing.');
%                 else
%                     dataset.(data_4D_namelist{i}).lon                  = dataset.(U_dataname).lon;
%                     dataset.(data_4D_namelist{i}).lat                  = dataset.(U_dataname).lat;
%                     dataset.(data_4D_namelist{i}).level                = dataset.(U_dataname).level;
%                     dataset.(data_4D_namelist{i}).date_yyyymmdd_AllYr  = dataset.(U_dataname).date_yyyymmdd_AllYr;
%                     
%                     if isequal(data_4D_namelist{i}, 'MSE')
%                         dataset.(data_4D_namelist{i}).unit = 'J kg**-1';                         % <-- mind the unit.
%                     else
%                         dataset.(data_4D_namelist{i}).unit = 'W kg**-1';                         % <-- mind the unit.
%                     end
%                     
%                     S = dataset.(data_4D_namelist{i});
%                     save(filename, 'S', '-v7.3');
%                     disp(['*** Written into ', filename, '... ***'])
%                 end
%             end
        elseif save_mat_or_nc == 2
            error('code not ready for save_mat_or_nc == 2');
        end
    end
    
    if clear_UV
        dataset.(U_dataname) = [];
        dataset.(V_dataname) = [];
    end
    
end

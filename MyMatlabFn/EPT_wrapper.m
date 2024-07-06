%%
function [dataset] = EPT_wrapper(varargin)
    
    % create a set of valid parameters and their default value
    pnames = {'dataset', 'select_field', 'stlevel',  'edlevel', 'year', 'prcssd_data_folder', 'save_mat_and_nc'};
    dflts  = cell(length(pnames), 1);

    %/ parse function arguments
    [         dataset,    select_field,    stlevel,    edlevel,   year,    prcssd_data_folder,  save_mat_and_nc] ...
           = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    fprintf('*** Running EPT_wrapper... ***\n')
    
    if stlevel > edlevel 
        error('Input stlevel should be smaller than edlevel.');
    end
    
    str_years = sprintf('%d-%d', year(1), year(end));
    str_level = sprintf('%d_%d', stlevel, edlevel);  
    
%     T_dataname = strcat('T', str_level);
%     q_dataname = strcat('q', str_level);
    
    T_dataname = 'T';
    q_dataname = 'q';
%     U_dataname = 'U';
%     V_dataname = 'V';
    
    %/ Extract data for all dates
    st_mth = 1; ed_mth = 12;
    
%     casedates = [];
%     for t = 1:length(year)
%         stDate = datetime(year(t), st_mth, 1,                       'Format', 'yyyy-MM-dd');
%         edDate = datetime(year(t), ed_mth, eomday(year(t), ed_mth), 'Format', 'yyyy-MM-dd');
% 
%         casedates_annual = stDate:edDate;
%         casedates_annual = [casedates_annual.Year*1e4 + casedates_annual.Month*100 + casedates_annual.Day]';
%         casedates = cat(1, casedates, casedates_annual);
%     end
    date_append                       = [];
    dataset.('EPT').(select_field)    = [];
    dataset.('dEPTdy').(select_field) = [];
    cnt = 0;
    for t = 1:length(year)
        cnt    = cnt + 1;
        stDate = datetime(year(t), st_mth, 1,                       'Format', 'yyyy-MM-dd');
        edDate = datetime(year(t), ed_mth, eomday(year(t), ed_mth), 'Format', 'yyyy-MM-dd');
        
        casedates_annual = stDate:edDate;
        casedates_annual = [casedates_annual.Year*1e4 + casedates_annual.Month*100 + casedates_annual.Day]';
        casedates        = casedates_annual;
        
        %/ Load 4D data year by year
        str_years_bc = sprintf('%d', year(t));  %/ to correctly read the pre-processed annual mat data
        always_reload = 1;    
        dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', T_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
        dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', q_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
%         dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', U_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
%         dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', V_dataname, 'select_field', select_field, 'str_years', str_years_bc, 'data_folder', prcssd_data_folder, 'always_reload',  always_reload);
        
        %===== Compute EPT (only need q and T) ======%
        P          = dataset.(q_dataname).level;
        ind        = find(P >= stlevel & P <= edlevel);
        P_subset   = P(ind);
        q_subset   = dataset.(q_dataname).(select_field)(:,:,ind,:);
        T_subset   = dataset.(T_dataname).(select_field)(:,:,ind,:);
        
        level_unit = 'hPa';
        EPT        = EquiPoteTemp('P', P_subset, 'q', q_subset, 'T', T_subset);
        size(EPT)
        
        %/ vertical mean of EPT (<EPT>) in K
        EPT_vm = vertinte('data',  EPT, 'level', P_subset, 'level_unit', level_unit, 'take_vertmean', 1); %/ vertical mean EPT (700-925 hPa)      
        
        dataset.('EPT').lon                  = dataset.(q_dataname).lon;
        dataset.('EPT').lat                  = dataset.(q_dataname).lat;
        dataset.('EPT').level                = P_subset;
        dataset.('EPT').unit                 = 'K';                         % <-- mind the unit.
        
        %/ meridional gradient of <EPT> in K/km
        r   = 6371e3;
        lat = dataset.(q_dataname).lat;
        hy  = diff(lat(1:2))*pi/180;                      %/ IMPORTANT: rmb to change degree to radian!!!!
        [dEPTdy, ~, ~] = gradient(EPT_vm, hy);                 
        dEPTdy                                     = 1./r*dEPTdy*1000*100;        %/ K/m --> K/100km.  div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))
        dataset.('dEPTdy').lon                     = dataset.(q_dataname).lon;
        dataset.('dEPTdy').lat                     = dataset.(q_dataname).lat;
        dataset.('dEPTdy').level                   = P_subset;
        dataset.('dEPTdy').unit                    = 'K/km';                % <-- mind the unit.
        
        
        fprintf('mean EPT_vm = %.2f\n', mean(EPT_vm, 'all'));
        fprintf('mean dEPTdy = %.2f\n', mean(dEPTdy, 'all'));
    
        dataset.('EPT').(select_field)    = cat(3, dataset.('EPT').(select_field),    EPT_vm);
        dataset.('dEPTdy').(select_field) = cat(3, dataset.('dEPTdy').(select_field), dEPTdy);
        date_append                       = [date_append; dataset.(q_dataname).date_yyyymmdd_AllYr];
    
        %/ Save EPT and dEPTdy as mat or nc (after fully appended)
        if cnt == length(year) && save_mat_and_nc 
            EPT_shortname = {'EPT', 'dEPTdy'};
            
            EPT_longname  = [{sprintf('Vertically integrated (%s hPa) equivalent potential temperature', strrep(str_level, '_', '-'))},...
                             {sprintf('Meridional gradient of vertically integrated (%s hPa) equivalent potential temperature', strrep(str_level, '_', '-'))}];
            
            for i = 1:length(EPT_shortname)
                
                dataset.(EPT_shortname{i}).date_yyyymmdd_AllYr = date_append;
                
                S        = dataset.(EPT_shortname{i});
                filename = strcat(prcssd_data_folder, EPT_shortname{i},'_',select_field, '_', str_years);
            
                %/ save as mat
                filename_bc = strcat(filename, '.mat');
                save(filename_bc, 'S','-v7.3');
                disp(['*** Written into ', filename_bc, '... ***'])

                %/ save as nc
                filename_bc = strcat(filename, '.nc');
                write_nc('ncfilename',  filename_bc, 'data', dataset.(EPT_shortname{i}).(select_field),...
                         'data_shortname', EPT_shortname{i}, 'data_standardname', EPT_shortname{i}, 'data_longname', EPT_longname{i},...
                         'data_unit', dataset.(EPT_shortname{i}).unit, 'lon', dataset.(EPT_shortname{i}).lon, 'lat', dataset.(EPT_shortname{i}).lat,...
                         'time', [], 'date', dataset.(EPT_shortname{i}).date_yyyymmdd_AllYr, 'remark', []);
                ncdisp(filename_bc)
                disp(['*** Written into ', filename_bc, '... ***'])
            end
        end
    end
end

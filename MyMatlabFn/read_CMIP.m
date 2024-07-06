function CMIPdata = read_CMIP(varargin)

    pnames = {'project_name',      'CMIP_folder',      'model_list',   'var_list',      'exp_list',       'slct_year_list',...
              'noleap',            'slct_level_range', 'units_level',  'ori_field',     'output_field',...
              'compute_MME',       'MME_name',         'MME_lon_res',  'MME_lat_res',   'interp_method',  'NumWorkers',...
              'anom_hist_period',  'output_folder',    'savemat',      'recompute',     'recompute_MME',  'load_or_not',...
              'delete_prob_files'};

    dflts =  cell(1, length(pnames));
    [          project_name,        CMIP_folder,        model_list,     var_list,       exp_list,          slct_year_list,...
               noleap,              slct_level_range,   units_level,    ori_field,      output_field,...
               compute_MME,         MME_name,           MME_lon_res,    MME_lat_res,    interp_method,     NumWorkers,...
               ~,                   output_folder,      savemat,        recompute,      recompute_MME,     load_or_not,...
               delete_prob_files] ...
                            = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/==========================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 19 Mar 2024
    %/==========================================

    %/ A list of accumulated variables (Not using it, but keep it for now)
    % acc_var_list = {'pr', 'evspsbl'};

    %/ A list of ocean variables (Not using it, but keep it for now)
    %/ See https://pcmdi.llnl.gov/mips/cmip3/variableList.html#Table_O1a
    % ocean_list = {'hfogo', 'stfmmc', 'zos', 'zostoga', 'zosga', 'tos',...
    %               'sic', 'sit', 'usi', 'vsi', 'wfo', 'stfbarot', 'hfcorr',...
    %               'wfcorr', 'tauucorr', 'tauvcorr', 'zobt', 'qflux',...
    %               'so', 'thetao', 'rhopoto', 'uo', 'vo', 'wo',...
    %               'zmlo', 'htovdiff', 'htovgyre', 'htovovrt', 'sltovdiff',...
    %               'sltovgyre', 'sltovovrt', 'sbl', 'hfsib', 'sltfsib'};

    CMIPdata = [];
    if ischar(CMIP_folder)
        CMIP_folder = {CMIP_folder};
    elseif ~iscell(CMIP_folder)
        error('Please input ''CMIP_folder'' (as a cell/string array)!');
    end
    
    if compute_MME
        load_or_not = 1; %/ automatically turn load_or_not on.
    end

    if ischar(model_list) model_list = {model_list};  end
    if ischar(var_list)   var_list = {var_list};      end
    if ischar(exp_list)   exp_list = {exp_list};      end

    %/ Convert it into string from cell, otherwise isequal would not work.
    if iscell(output_field)
        output_field = output_field{:};
    end

    if isequal(units_level, 'hPa')
        unitconv_level = 100;   
    elseif isequal(units_level, 'Pa')
        unitconv_level = 1;
    else
        error('Set ''units_level'' to be ''hPa'' or ''Pa''!');
    end

    %/ Defaults
    dt_slct_mo = []; dt_slct_hr = []; dt_slct_min = [];  %/ time interval in mo/hr/min
    
    %/ Set dates (if slct_year is not empty)
    if isempty(slct_year_list)
        slct_year_list = 1850:2014;  
    end

    if noleap 
        if isequal(output_field, 'daily')
            str_noleap = '_noleap'; 
        elseif isequal(output_field, 'monthly')
            str_noleap = '';
            noleap = 0; 
            warning('For original field of monthly, resume ''noleap'' to 0.');
        else
            error('code not ready!');
        end
    else
        str_noleap = ''; 
    end

    problem_files = {}; slct_level = [];
    for k = 1:length(var_list)
        var = var_list{k};
        for n = 1:length(exp_list)
            if iscell(slct_year_list)
                slct_year = slct_year_list{n};
            else
                slct_year = slct_year_list;
            end
            str_slct_year = sprintf('_%d-%d', slct_year(1), slct_year(end));

            if isequal(output_field, 'monthly')
                date_fld           = 'date_yyyymm';
                output_date_format = 'yyyyMM';
                slct_output_dates  = date_array_gen('year_list', slct_year, 'dt_slct_mo', 1, 'output_date_format', output_date_format, 'noleap', noleap);   

            elseif isequal(output_field, 'daily')
                date_fld           = 'date_yyyymmdd_AllYr';
                output_date_format = 'yyyyMMdd';
                slct_output_dates  = date_array_gen('year_list', slct_year, 'dt_slct_hr', 24, 'output_date_format', output_date_format, 'noleap', noleap);
        
            elseif isequal(output_field, 'subdaily')
                date_fld           = 'date_yyyymmddHHMM';
                output_date_format = 'yyyyMMddHHmm';
                slct_output_dates  = date_array_gen('year_list', slct_year, 'dt_slct_mo', dt_slct_mo, 'dt_slct_hr', dt_slct_hr, 'dt_slct_min', dt_slct_min, 'output_date_format', output_date_format, 'noleap', noleap);          
            else
                error('code not set!'); 
            end
            if isempty(slct_output_dates)
                error('slct_output_dates is empty! Check your function!');
            end

            flag_to_compute_MME = 0; 
            for m = 1:length(model_list)  %/ compute MME (not ready)
                model        = model_list{m};
                model_strrep = strrep(model, '-', '_');  %/ Sadly, '-' is not allowed in the field name. Replacing it with '_'.
                exp          = exp_list{n};
                exp_strrep   = strrep(exp,   '-', '_');

                if isempty(ori_field)     
                    ori_field = read_CMIP_ori_field('model', model, 'var', var, 'output_field', output_field);
                end 
                slct_ens     = read_CMIP_ens('project_name',  project_name, 'model_list',  model_list,  'var', var, 'exp', exp, 'ori_field', ori_field);

                if isequal(ori_field, 'monthly')
                    LEN = 6;
                    ori_date_format = 'yyyyMM';
                    slct_ori_dates  = date_array_gen('year_list', slct_year, 'dt_slct_mo', 1, 'output_date_format', ori_date_format, 'noleap', noleap);   
        
                elseif isequal(ori_field, 'daily')
                    LEN = 8;
                    ori_date_format = 'yyyyMMdd';
                    slct_ori_dates  = date_array_gen('year_list', slct_year, 'dt_slct_hr', 24, 'output_date_format', ori_date_format, 'noleap', noleap);
        
                elseif contains(ori_field, 'hr')
                    LEN = 12;
                    ori_date_format = 'yyyyMMddHHmm';
                    slct_ori_dates  = date_array_gen('year_list', slct_year, 'dt_slct_mo', dt_slct_mo, 'dt_slct_hr', dt_slct_hr, 'dt_slct_min', dt_slct_min, 'output_date_format', ori_date_format, 'noleap', noleap);          
                else
                    error('unknown date format for ori_field = %s', ori_field)
                end
                if isempty(slct_ori_dates)
                    error('slct_ori_dates is empty! Check your function!');
                end

                flag_skip_the_model = 0;
                if ~isempty(slct_level_range)  
                    % var
                    % slct_level_range
                    slct_level = int64(slct_level_range{k});  
                end
                
                %/ Time scale/Frequency (day/Oday/Amon/Omon/etc.)
                ts = read_CMIP_ts('model', model, 'var', var, 'ori_field', ori_field);

                %/ Grid format (gn/gr/gr1/etc.)
                grid = read_CMIP_grid('model', model, 'exp', exp, 'var', var, 'ori_field', ori_field);
                
                if ismember(var, {'ua', 'va', 'wap', 'ta', 'hus', 'div'}) 
                    if isempty(slct_level)     
                        error('slct_level is empty! Check ''slct_level_range''!');  
                    end
                    if numel(slct_level) == 1
                        var_strrep = strcat(var, num2str(slct_level*unitconv_level/100));
                    else
                        var_strrep = strcat(var, num2str(min(slct_level)*unitconv_level/100),'to',num2str(max(slct_level)*unitconv_level/100));
                    end
                else
                    var_strrep = var;
                end
                
                %/ Check if MME data has been post-processed
                domain_str = '_global';
                if compute_MME 
                    [~,d]        = cellfun(@size,slct_ens);
                    max_ens_size = max(d);
                    MME_ens      = sprintf('%dEM', max_ens_size);
                    
                    MME_output_filename = fullfile(output_folder, strcat(MME_name,'_',exp,'_',var_strrep,'_',MME_ens,'_',output_field,domain_str,str_slct_year,str_noleap,'.mat'));
                    if isfile(MME_output_filename) && recompute_MME == 0
                        disp(MME_output_filename)
                        fprintf('*** Postprocessed MME file exists. Double check if the ensemble list is correct... ***\n');
                        load(MME_output_filename, 'S');
                        
                        if isequal(S.ens_list, slct_ens)
                            fprintf('*** Postprocessed MME has correct ens_list. No recomputation is required. ***\n');
                            fld = fieldnames(S);
                            for f = 1:length(fld)
                                CMIPdata.(MME_name).(fld{f}) = S.(fld{f});
                            end
                            CMIPdata.(MME_name).(strcat(var,'_lon')) = S.lon;  %/ Since different var may have different griddings, even for the same model and exp.
                            CMIPdata.(MME_name).(strcat(var,'_lat')) = S.lat;  %/ Since different var may have different griddings, even for the same model and exp.
                            
                            clear S;
                            data = CMIPdata.(MME_name).(strcat(var_strrep,'_',exp_strrep));
                            data_units = CMIPdata.(MME_name).(strcat(var_strrep,'_units'));
                            fprintf('Mean of data: %.3g %s\n', mean(data,'all','omitnan'),data_units);
                            fprintf('Min  of data: %.3g %s\n', min(data,[],'all','omitnan'),data_units);
                            fprintf('Max  of data: %.3g %s\n', max(data,[],'all','omitnan'),data_units);
                            break;    %/ Break the model loop
                        else
                            flag_to_compute_MME = 1;
                            warning('*** Postprocessed MME has inconsistent ens_list. Recomputation is required. ***\n');
                        end
                    else
                        flag_to_compute_MME = 1;
                    end
                else
                    flag_to_compute_MME = 0;
                end

                %/ Check if the individual model data has been post-processed
                slct_ens_indiv = slct_ens{m,:};
                slct_ens_indiv(cellfun(@isempty,slct_ens_indiv)) = [];  %/ remove empty cell from the query 
                str_ens_list    = strcat('_',strjoin(slct_ens_indiv,'_'));  %/ set dummy ens string for post-processing data
                output_filename = fullfile(output_folder, strcat(model,'_',exp,'_',var_strrep, str_ens_list,'_',output_field,domain_str,str_slct_year,str_noleap,'.mat'));
                
                if isfile(output_filename) && recompute == 0
                    flag_to_read = 0;
                    if load_or_not
                        fprintf('*** Postprocessed file %s exists. Loading from it... ***\n', output_filename);
                        load(output_filename, 'S');
                        fld = fieldnames(S);
                        for f = 1:length(fld)
                            CMIPdata.(model_strrep).(fld{f}) = S.(fld{f});
                        end
                        CMIPdata.(model_strrep).(strcat(var,'_lon')) = S.lon;  %/ Since different var may have different griddings, even for the same model and exp.
                        CMIPdata.(model_strrep).(strcat(var,'_lat')) = S.lat;  %/ Since different var may have different griddings, even for the same model and exp.
                        clear S;
                    else
                        fprintf('*** Postprocessed file %s exists. No loading as required. ***\n', output_filename);
                    end
                else
                    flag_to_read = 1;
                    fprintf('*** Postprocessed file %s not exists (OR recompute == 1). Loading the raw files... ***\n', output_filename);
                end

                %/ Read data from the raw files / processed data
                if flag_to_read
                    p_Wang1988 = [300, 500, 700];                      %/ The two-layer model of Wang 1988 
                    p2         = p_Wang1988(2)*100;                    %/ hPa to Pa
                    pb = 1000;
                    pt = 100;
                    p_unit = 'hPa';
                    if isequal(p_unit, 'hPa')
                        CONV_p = 100;
                    else
                        CONV_p = 1;
                    end
                    if ismember(var, {'S500', 'S2_Wang1988'}) %/ Static Stability at 500 hPa
                        var_list_nest = {'ta', 'ta', 'ta'}; 
                        slct_level_range_nest = num2cell(p_Wang1988);
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);

                        %/ Compute potential temperature 
                        T500  = temp_CMIPdata.(model_strrep).(strcat('ta500_',exp_strrep));
                        PT = [];
                        for i = 1:length(slct_level_range_nest)
                            p = slct_level_range_nest{i};
                            PT.(strcat('p',num2str(p))) = compute_PT('T', temp_CMIPdata.(model_strrep).(strcat('ta', num2str(p), '_', exp_strrep)),...
                                                                     'p', p, 'p_unit', p_unit);

                            %/ Spare memory
                            temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('ta', num2str(p), '_', exp_strrep));
                        end

                        dp = (p_Wang1988(3) - p_Wang1988(1))*CONV_p;  %/ hPa to Pa
                        
                        if ismember(var, {'S500'})
                            %/ Compute dry static stability at 500 hPa using conventional equation
                            S500 = (T500./PT.p500).*(PT.p300 - PT.p700)/dp;
                            if any(~isreal(S500))
                                error('C0 contains complex number!')
                            end

                            CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = S500;
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'K Pa^{-1}';  

                        elseif ismember(var, {'S2_Wang1988'})
                            %/ Compute dry static stability at 500 hPa using Eq. (3.2) in Wang (1988), 
                            %/ useful when computing the nondimensional I.

                            %/ Assume for dry air, the air density can be
                            %/ obtained from ideal gas law (rho = P/RT)
                            R      = 287;                %/ J kg-1 K-1
                            rho500 = 500*100./(R.*T500); %/ kg m^-3
                            S2 = (1./PT.p500./rho500).*(PT.p300 - PT.p700)/dp;

                            if any(~isreal(S2))
                                error('S2 contains complex number!')
                            end

                            CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = S2;
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'm^{2} s^{2} Pa^{-2}';  
                            clear rho500;
                        end

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;
                        clear PT;
                        clear T500;
                        % fprintf('Mean of %s: %.3g %s\n', var, mean(CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)),'all','omitnan'),  CMIPdata.(model_strrep).(strcat(var_strrep,'_units')));
                        % fprintf('Min  of %s: %.3g %s\n', var, min(CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)),[],'all','omitnan'),CMIPdata.(model_strrep).(strcat(var_strrep,'_units')));
                        % fprintf('Max  of %s: %.3g %s\n', var, max(CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)),[],'all','omitnan'),CMIPdata.(model_strrep).(strcat(var_strrep,'_units')));

                    elseif ismember(var, {'S5'}) %/ Static Stability at 500 hPa (using 400 and 600 hPa)
                        var_list_nest = {'ta', 'ta', 'ta'}; 
                        slct_level_range_nest = num2cell([400, 500, 600]);
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);

                        %/ Compute potential temperature 
                        T500  = temp_CMIPdata.(model_strrep).(strcat('ta500_',exp_strrep));
                        
                        
                        PT = [];
                        for i = 1:length(slct_level_range_nest)
                            p = slct_level_range_nest{i};
                            p_unit = 'hPa'; %<- mind this!!!
                            PT.(strcat('p',num2str(p))) = compute_PT('T', temp_CMIPdata.(model_strrep).(strcat('ta', num2str(p), '_', exp_strrep)),...
                                                                     'p', p, 'p_unit', p_unit);

                            %/ Spare memory
                            temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('ta', num2str(p), '_', exp_strrep));
                        end
  
                        %/ Compute dry static stability at 500 hPa using conventional equation
                        if isequal(p_unit, 'hPa')
                            CONV_p = 100;
                        elseif isequal(p_unit, 'Pa')
                            CONV_p = 1;
                        end

                        S5 = (T500./PT.p500).*(PT.p400 - PT.p600)/((600-400)*CONV_p);
                        if any(~isreal(S5))
                            error('S5 contains complex number!')
                        end

                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = S5;
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'K Pa^{-1}';  

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;
                        clear PT;
                        clear T500;
                        
                    elseif ismember(var, {'S400'}) %/ Static Stability at 400 hPa (using 300 and 500 hPa)
                        var_list_nest = {'ta', 'ta', 'ta'}; 
                        slct_level_range_nest = num2cell([300, 400, 500]);
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);

                        %/ Compute potential temperature 
                        T400  = temp_CMIPdata.(model_strrep).(strcat('ta400_',exp_strrep));
                        
                        PT = [];
                        for i = 1:length(slct_level_range_nest)
                            p = slct_level_range_nest{i};
                            p_unit = 'hPa'; %<- mind this!!!
                            PT.(strcat('p',num2str(p))) = compute_PT('T', temp_CMIPdata.(model_strrep).(strcat('ta', num2str(p), '_', exp_strrep)),...
                                                                     'p', p, 'p_unit', p_unit);

                            %/ Spare memory
                            temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('ta', num2str(p), '_', exp_strrep));
                        end
  
                        %/ Compute dry static stability at 500 hPa using conventional equation
                        if isequal(p_unit, 'hPa')
                            CONV_p = 100;
                        elseif isequal(p_unit, 'Pa')
                            CONV_p = 1;
                        end

                        S400 = (T400./PT.p400).*(PT.p300 - PT.p500)/((500-300)*CONV_p);
                        if any(~isreal(S400))
                            error('S400 contains complex number!')
                        end

                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = S400;
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'K Pa^{-1}';  

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;
                        clear PT T400 S400;

                    elseif ismember(var, {'S'}) %/ Static Stability 
                        var_list_nest = {'ta'}; 
                        slct_level_range_nest = {[pt, pb]};
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);

                        p = double(temp_CMIPdata.(model_strrep).level); %/ Make sure it's double, integer will cause bug in gradient()
                        disp(p)
                        p_unit = 'Pa'; CONV_p = 1;
                        p_dim  = 3;
                        T = temp_CMIPdata.(model_strrep).(strcat('ta', num2str(pt),'to',num2str(pb), '_', exp_strrep));
                        PT = compute_PT('T', T, 'p_dim', p_dim,  'p', p, 'p_unit', p_unit);
                        
                        %/ Static Stability
                        if p_dim == 3
                            [~, ~, dPT_dA, ~] = gradient(PT);                      %/ Chain rule to do derivate with *non-uniform* spacings.
                        else
                            error('code not set!');
                        end
                        dp_dA     = gradient(p * CONV_p);                        %/ get non-uniform gradient of pressure level first.
                        dp_dA     = reshape(dp_dA, [ones(1, p_dim-1), length(dp_dA)]);      %/ 1 x 1 x plev
                        dPT_dP    = dPT_dA./dp_dA;                         %/ J kg-1 Pa-1   elementwise division.  
                        SS        = -T./PT.*dPT_dP; 
                        
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = SS;
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'K Pa^{-1}';  

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end

                        %/ Spare memory
                        clear temp_CMIPdata;
                        clear T PT dPT_dP SS;

                    elseif ismember(var, {'uIVT', 'vIVT'}) 
                        var_list_nest         = {'ua',    'va',    'hus',  'ps'}; 
                        slct_level_range_nest = {[pt, pb],[pt, pb],[pt, pb],[]};
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);
                        sp      = temp_CMIPdata.(model_strrep).(strcat('ps_', exp_strrep));
                        sp_lon  = double(temp_CMIPdata.(model_strrep).(strcat('ps_','lon')));
                        sp_lat  = double(temp_CMIPdata.(model_strrep).(strcat('ps_','lat')));

                        hus     = temp_CMIPdata.(model_strrep).(strcat('hus', num2str(pt),'to',num2str(pb), '_', exp_strrep));
                        hus_lon = double(temp_CMIPdata.(model_strrep).(strcat('hus_','lon')));
                        hus_lat = double(temp_CMIPdata.(model_strrep).(strcat('hus_','lat')));

                        ua      = temp_CMIPdata.(model_strrep).(strcat('ua',  num2str(pt),'to',num2str(pb), '_', exp_strrep));
                        va      = temp_CMIPdata.(model_strrep).(strcat('va',  num2str(pt),'to',num2str(pb), '_', exp_strrep));
                        uv_lon  = double(temp_CMIPdata.(model_strrep).(strcat('ua_','lon'))); %/ Dimensions of ua should be equal to that of va.
                        uv_lat  = double(temp_CMIPdata.(model_strrep).(strcat('ua_','lat')));
                        
                        p_unit = 'Pa'; 
                        p = double(temp_CMIPdata.(model_strrep).level); %/ Make sure it's double, integer will cause bug in gradient()
                        disp(p)

                        if ~isequal(sp_lon,uv_lon) || ~isequal(sp_lat,uv_lat) 
                            warning('Detected that for some reason %s output ''sp'' has a different gridding (%d x %d) with ''ua'' (%d x %d). Interpolating ''sp'' to that of ''ua''...',...
                                     model, length(sp_lon), length(sp_lat), length(uv_lon), length(uv_lat));
                            NumWorkers_interp = 40;
                            sp = my_interp('lon_old', sp_lon, 'lat_old', sp_lat, 'data', sp, 'lon_new', uv_lon, 'lat_new', uv_lat,...
                                            'is_global', 1, 'lon_dim', 1, 'interp_method', interp_method, 'NumWorkers', NumWorkers_interp);
                            
                        end
                        if ~isequal(hus_lon,uv_lon) || ~isequal(hus_lat,uv_lat) 
                            warning('Detected that for some reason %s output ''hus'' has a different gridding (%d x %d) with ''ua'' (%d x %d). Interpolating ''hus'' to that of ''ua''...',...
                                     model, length(hus_lon), length(hus_lat), length(uv_lon), length(uv_lat));
                            NumWorkers_interp = 40;
                            hus = my_interp('lon_old', hus_lon, 'lat_old', hus_lat, 'data', hus, 'lon_new', uv_lon, 'lat_new', uv_lat,...
                                            'is_global', 1, 'lon_dim', 1, 'interp_method', interp_method, 'NumWorkers', NumWorkers_interp);
                            
                        end
                        temp_CMIPdata.(model_strrep).lon = uv_lon; %/ Update
                        temp_CMIPdata.(model_strrep).lat = uv_lat; %/ Update

                        % [uIVT, vIVT] = calc_IVT('U', ua, 'V', va, 'q', hus, 'level', p, 'level_unit', p_unit,...
                        %                         'topo_lon', [], 'topo_lat', [], 'sp_mode', [], 'interp_mode', []);

                        sp_mode = 'sp'; interp_mode = 'linear'; lnP_coor = 1;
                        % sp_mode = 'sp'; interp_mode = 'linear'; lnP_coor = 0;
                        [uIVT, vIVT] = calc_IVT('U', ua, 'V', va, 'q', hus, 'level', p, 'level_unit', p_unit,...
                                                'topo_lon', [], 'topo_lat', [], 'lnP_coor', lnP_coor, 'sp_mode', sp_mode, 'sp', sp, 'interp_mode', interp_mode);

                        if isequal(var, 'uIVT')
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = uIVT;
                        elseif isequal(var, 'vIVT')
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = vIVT;
                        else
                            error('check your var input!');
                        end
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'kg/m/s';  

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end

                        %/ Spare memory
                        clear temp_CMIPdata; 

                    elseif ismember(var, {'MC925'}) %/ MC925 == 925-hPa moisture convergence
                        var_list_nest         = {'ua', 'va', 'hus'}; 
                        p_BL                  = 925;
                        slct_level_range_nest = {p_BL, p_BL, p_BL};
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);
                        
                        hus     = temp_CMIPdata.(model_strrep).(strcat('hus', num2str(p_BL), '_', exp_strrep));
                        hus_lon = double(temp_CMIPdata.(model_strrep).(strcat('hus_','lon')));
                        hus_lat = double(temp_CMIPdata.(model_strrep).(strcat('hus_','lat')));

                        ua      = temp_CMIPdata.(model_strrep).(strcat('ua',  num2str(p_BL), '_', exp_strrep));
                        va      = temp_CMIPdata.(model_strrep).(strcat('va',  num2str(p_BL), '_', exp_strrep));
                        uv_lon  = double(temp_CMIPdata.(model_strrep).(strcat('ua_','lon'))); %/ Dimensions of ua should be equal to that of va.
                        uv_lat  = double(temp_CMIPdata.(model_strrep).(strcat('ua_','lat')));
                        
                        % p_unit = 'Pa'; 
                        % p = double(temp_CMIPdata.(model_strrep).level); %/ Make sure it's double, integer will cause bug in gradient()
                        % disp(p)
                        if ~isequal(hus_lon,uv_lon) || ~isequal(hus_lat,uv_lat) 
                            warning('Detected that for some reason %s output ''hus'' has a different gridding (%d x %d) with ''ua'' (%d x %d). Interpolating ''hus'' to that of ''ua''...',...
                                     model, length(hus_lon), length(hus_lat), length(uv_lon), length(uv_lat));
                            NumWorkers_interp = 40;
                            hus = my_interp('lon_old', hus_lon, 'lat_old', hus_lat, 'data', hus, 'lon_new', uv_lon, 'lat_new', uv_lat,...
                                            'is_global', 1, 'lon_dim', 1, 'interp_method', interp_method, 'NumWorkers', NumWorkers_interp);
                        end
                        temp_CMIPdata.(model_strrep).lon = uv_lon; %/ Update
                        temp_CMIPdata.(model_strrep).lat = uv_lat; %/ Update

                        hx = diff(uv_lon(1:2))*pi/180;   %/ rmb to change degree to radian!!!!
                        [~, dAdx, ~] = gradient(ua.*hus, hx);
                        clear ua;
                        
                        hy = diff(uv_lat(1:2))*pi/180;   %/ rmb to change degree to radian!!!!
                        [dAdy, ~, ~] = gradient(va.*hus, hy);
                        clear va hus;
                        
                        [~, lat_2D] = meshgrid(uv_lon, uv_lat);
                        lat_2D = lat_2D';
                        
                        r = 6371e3;   %/ in m
                        conv = -1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D));  %/ Convergence == -(du/dx + dv/dy) in cartesian coor.
                        
                        cond = isinf(conv);
                        conv(cond) = nan;  %/ set inf to nan, as it occurs at lat = 90 or -90 when multiplying with 1./(r*cosd(lat_2D))!
                        
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = conv;
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 's^{-1}';  

                        if isequal(var, 'uIVT')
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = uIVT;
                        elseif isequal(var, 'vIVT')
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = vIVT;
                        else
                            error('check your var input!');
                        end
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'kg/m/s';  

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end

                        %/ Spare memory
                        clear temp_CMIPdata; 

                    elseif ismember(var, {'MC1000to850'}) %/ 1000-850  == 925-hPa moisture convergence
                        pt_BL                 = 850;
                        pb_BL                 = 1000;
                        var_list_nest         = {'ua', 'va', 'hus'}; 
                        slct_level_range_nest = {[pt_BL, pb_BL],[pt_BL, pb_BL],[pt_BL, pb_BL]};
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);
                        
                        disp(temp_CMIPdata.(model_strrep)) %/ checking the fields
                        str_lv  = sprintf('%dto%d', min(slct_level_range_nest{1}), max(slct_level_range_nest{1}));
                        hus     = temp_CMIPdata.(model_strrep).(strcat('hus', str_lv, '_', exp_strrep));
                        hus_lon = double(temp_CMIPdata.(model_strrep).(strcat('hus', '_lon')));
                        hus_lat = double(temp_CMIPdata.(model_strrep).(strcat('hus', '_lat')));
                        temp_CMIPdata.(model_strrep).(strcat('hus', str_lv, '_', exp_strrep)) = [];  %/ spare memory;

                        ua      = temp_CMIPdata.(model_strrep).(strcat('ua',  str_lv, '_', exp_strrep));
                        temp_CMIPdata.(model_strrep).(strcat('ua', str_lv, '_', exp_strrep)) = [];  %/ spare memory;
                        va      = temp_CMIPdata.(model_strrep).(strcat('va',  str_lv, '_', exp_strrep)); 
                        temp_CMIPdata.(model_strrep).(strcat('va', str_lv, '_', exp_strrep)) = [];  %/ spare memory;
                        uv_lon  = double(temp_CMIPdata.(model_strrep).(strcat('ua', '_lon'))); %/ Dimensions of ua should be equal to that of va.
                        uv_lat  = double(temp_CMIPdata.(model_strrep).(strcat('ua', '_lat')));
                        

                        if ~isequal(hus_lon,uv_lon) || ~isequal(hus_lat,uv_lat) 
                            warning('Detected that for some reason %s output ''hus'' has a different gridding (%d x %d) with ''ua'' (%d x %d). Interpolating ''hus'' to that of ''ua''...',...
                                     model, length(hus_lon), length(hus_lat), length(uv_lon), length(uv_lat));
                            NumWorkers_interp = 40;
                            hus = my_interp('lon_old', hus_lon, 'lat_old', hus_lat, 'data', hus, 'lon_new', uv_lon, 'lat_new', uv_lat,...
                                            'is_global', 1, 'lon_dim', 1, 'interp_method', interp_method, 'NumWorkers', NumWorkers_interp);
                        end
                        temp_CMIPdata.(model_strrep).lon = uv_lon; %/ Update
                        temp_CMIPdata.(model_strrep).lat = uv_lat; %/ Update

                        p_unit = 'Pa'; 
                        p      = double(temp_CMIPdata.(model_strrep).level); %/ Make sure it's double, integer will cause bug in gradient()
                        disp(p)

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;  %/ spare memory

                        fprintf('Computing dAdx...\n')
                        hx = diff(uv_lon(1:2))*pi/180;   %/ rmb to change degree to radian!!!!
                        [~, dAdx, ~] = gradient(ua.*hus, hx);
                        clear ua;
                        
                        fprintf('Computing dAdy...\n')
                        hy = diff(uv_lat(1:2))*pi/180;   %/ rmb to change degree to radian!!!!
                        [dAdy, ~, ~] = gradient(va.*hus, hy);
                        clear va hus;
                        
                        [~, lat_2D] = meshgrid(uv_lon, uv_lat);
                        lat_2D = lat_2D';
                        
                        fprintf('Computing MC...\n')
                        r  = 6371e3;   %/ in m
                        MC = -1 * (1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D)));  %/ Convergence == -(du/dx + dv/dy) in cartesian coor.
                        clear dAdx dAdy;

                        cond     = isinf(MC);
                        MC(cond) = nan;  %/ set inf to nan, as it occurs at lat = 90 or -90 when multiplying with 1./(r*cosd(lat_2D))!
                        
                        %/ Vertical integartion of moisture convergence (1000-850 hPa)
                        fprintf('Performing vertinte...\n')
                        wd     = 997;  %/ density of water 
                        sp     = [];   sp_mode = []; interp_mode = []; lnP_coor = 0;
                        MC     = 1000/wd*vertinte('data', MC, 'level', p, 'level_unit', p_unit, 'lnP_coor', lnP_coor, 'sp_mode', sp_mode, 'sp', sp, 'interp_mode', interp_mode);
                        CONV   = 24*3600;

                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = MC.*CONV;  %/ convert mm s-1 -> mm day-1
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'mm day^{-1}';   %/ convert mm s-1 -> mm day-1 
                        
                        %/ Spare memory
                        clear MC; 

                    elseif ismember(var, {'C0_Wang1988', 'alpha_Wang1988'}) 
                        %/ The long gravity wave speed of the gravest baroclinic mode 
                        %/ based on Eq. (3.9) in Wang (1988)

                        var_list_nest = {'S2_Wang1988'}; 
                        slct_level_range_nest = [];
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);
                        
                        S2    = temp_CMIPdata.(model_strrep).(strcat('S2_Wang1988_',exp_strrep));
                        C0    = sqrt(1/2*S2*dp^2);             %/ m s-1
                        
                        if any(~isreal(C0))
                            error('C0 contains complex number!')
                        end
                        

                        Cp    = 1004;                          %/ J kg-1 K-1
                        R     = 287;                           %/ J kg-1 K-1
                        Lc    = 2.26e6;                        %/ J kg-1
                        b     = 0.9;                           %/ non-dimensional
                        alpha = (2*Cp*p2*C0.^2)/(R*b*Lc*dp);   %/ non-dimensional
                        
                        if ismember(var, {'C0_Wang1988'})
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = C0;
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'm s^{-1}';  

                        elseif ismember(var, {'alpha_Wang1988'})
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = alpha;
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = '';  
                        end

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;

                    elseif ismember(var, {'q3_bar_Wang1988', 'q1_bar_Wang1988'}) 
                        var_list_nest = {'hus'}; 
                        if isequal(var, 'q3_bar_Wang1988')
                            slct_level_range_nest = {[500, 925]}; %/ Vertical mean of q3 (500-925), for no 900 hPa in CMIP6 output!
                        
                        elseif isequal(var, 'q1_bar_Wang1988')
                            slct_level_range_nest = {[100, 500]}; %/ Vertical mean of q1 (100-500)
                        end
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);

                        str_lv = sprintf('%dto%d', min(slct_level_range_nest{1}), max(slct_level_range_nest{1}));
                        q = temp_CMIPdata.(model_strrep).(strcat('hus', str_lv, '_', exp_strrep));
                        level = temp_CMIPdata.(model_strrep).level;
                        disp(level)

                        take_vertmean = 1;  %/ Perform vertical mean
                        level_unit = 'Pa';  %/ The default unit of level in CMIP6 data 
                        
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = vertinte('data', q, 'level', level, 'level_unit', level_unit, 'take_vertmean', take_vertmean);
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'kg/kg';  

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;

                    elseif ismember(var, {'dq_bar_Wang1988'}) 
                        var_list_nest = {'q3_bar_Wang1988', 'q1_bar_Wang1988'}; 
                        slct_level_range_nest = [];
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);
                        
                        q3_bar = temp_CMIPdata.(model_strrep).(strcat('q3_bar_Wang1988', '_', exp_strrep));
                        temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('q3_bar_Wang1988', '_', exp_strrep)); %/ clear out 

                        q1_bar = temp_CMIPdata.(model_strrep).(strcat('q1_bar_Wang1988', '_', exp_strrep));
                        temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('q1_bar_Wang1988', '_', exp_strrep)); %/ clear out 
 
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = (q3_bar - q1_bar);
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'kg/kg';  

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;

                    elseif ismember(var, {'I_Wang1988'}) 
                        var_list_nest = {'q3_bar_Wang1988', 'q1_bar_Wang1988', 'alpha_Wang1988'}; 
                        slct_level_range_nest = [];
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);
                        
                        q3_bar = temp_CMIPdata.(model_strrep).(strcat('q3_bar_Wang1988', '_', exp_strrep));
                        temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('q3_bar_Wang1988', '_', exp_strrep)); %/ clear out 

                        q1_bar = temp_CMIPdata.(model_strrep).(strcat('q1_bar_Wang1988', '_', exp_strrep));
                        temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('q1_bar_Wang1988', '_', exp_strrep)); %/ clear out 
                        
                        alpha = temp_CMIPdata.(model_strrep).(strcat('alpha_Wang1988', '_', exp_strrep));
                        temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('alpha_Wang1988', '_', exp_strrep)); %/ clear out 

                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = (q3_bar - q1_bar)./alpha;
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = '';  

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;

                    elseif ismember(var, {'C1_Wang1988'}) 
                        var_list_nest = {'C0_Wang1988', 'I_Wang1988'}; 
                        slct_level_range_nest = [];
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);
                        
                        C0 = temp_CMIPdata.(model_strrep).(strcat('C0_Wang1988', '_', exp_strrep));
                        temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('C0_Wang1988', '_', exp_strrep)); %/ clear out 

                        I = temp_CMIPdata.(model_strrep).(strcat('I_Wang1988', '_', exp_strrep));
                        temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('I_Wang1988', '_', exp_strrep)); %/ clear out 
 
                        C1 = C0.*real(sqrt((1-I)));

                        if any(~isreal(C1))
                            error('C1 contains complex number!')
                        end

                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = C1; %/ Take the real part of sqrt((1-I))
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'm/s';  

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;

                    elseif ismember(var, {'Lcpr'}) 
                        var_list_nest = {'pr'}; 
                        slct_level_range_nest = [];
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);
                        
                        Lc = 2.26e6;   %/ J kg-1
                        if isequal(output_field, 'monthly')
                            CONV          = 1/30/24/3600;
                        elseif isequal(output_field, 'daily')
                            CONV          = 1/24/3600;
                        else
                            error('code not set!'); 
                        end
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = Lc * CONV * temp_CMIPdata.(model_strrep).(strcat('pr', '_', exp_strrep));
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'W m^{-2}';  

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;

                    % elseif ismember(var, {'h'})  %/ h == MSE
                    %     var_list_nest = {'hus', 'ta', 'zg'}; 
                    %     slct_level_range_nest = repmat({[100, 1000]}, 1, length(var_list_nest));
                    % 
                    %     temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                    %                          'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                    %                          'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                    %                          'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                    %                          'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);
                    % 
                    %     q3_bar = temp_CMIPdata.(model_strrep).(strcat('hus', '_', exp_strrep));
                    %     temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('q3_bar_Wang1988', '_', exp_strrep)); %/ clear out 
                    % 
                    %     q1_bar = temp_CMIPdata.(model_strrep).(strcat('tas', '_', exp_strrep));
                    %     temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('q1_bar_Wang1988', '_', exp_strrep)); %/ clear out 
                    % 
                    %     alpha = temp_CMIPdata.(model_strrep).(strcat('zg', '_', exp_strrep));
                    %     temp_CMIPdata.(model_strrep) = rmfield(temp_CMIPdata.(model_strrep), strcat('alpha_Wang1988', '_', exp_strrep)); %/ clear out 
                    % 
                    %     CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = (q3_bar - q1_bar)./alpha;
                    %     CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = '';  
                    % 
                    %     %/ Copy the rest of the fields
                    %     other_flds = {'lon', 'lat', 'ens', flds_date, 'institution', 'institution_id'};
                    %     for f = 1:length(other_flds)
                    %         CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                    %     end
                    %     clear temp_CMIPdata;

                    elseif isequal(var, 'SS') %/ Static Stability (ta500 - ta850)
                        var_list_nest        = {'ta', 'ta'}; 
                        slct_level_range_nest = num2cell([850, 500]);
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);

                        %/ ta850 - ta500
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = temp_CMIPdata.(model_strrep).(strcat('ta500_',exp_strrep)) - temp_CMIPdata.(model_strrep).(strcat('ta850_',exp_strrep));
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = temp_CMIPdata.(model_strrep).('ta850_units');  
                        
                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;

                    elseif isequal(var, 'SS2') %/ Static Stability (ta300 - ta500)
                        var_list_nest        = {'ta', 'ta'}; 
                        slct_level_range_nest = num2cell([500, 300]);
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);

                        %/ ta500 - ta300
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = temp_CMIPdata.(model_strrep).(strcat('ta300_',exp_strrep)) - temp_CMIPdata.(model_strrep).(strcat('ta500_',exp_strrep));
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = temp_CMIPdata.(model_strrep).('ta500_units');  
                        
                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;

                    elseif isequal(var, 'hus500wap300') %/ Diagnosed precipitation (1/g * q500 * Omega300)
                        var_list_nest        = {'hus', 'wap'}; 
                        slct_level_range_nest = num2cell([500, 300]);
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);

                        %/ 1/g * q500 * Omega300 [mm/s]
                        unit_conv = 24*3600;   %/ mm/s -> mm/day
                        g         = 9.81;
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = unit_conv .* 1/g .* temp_CMIPdata.(model_strrep).(strcat('hus500_',exp_strrep)) .* temp_CMIPdata.(model_strrep).(strcat('wap300_',exp_strrep));
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 'mm day^{-1}';
                        
                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;

                    elseif isequal(var, 'div') %/ wind divergence
                        if isempty(slct_level_range) 
                            error('Specify ''slct_level_range'' for %s!', var);
                        end
                        var_list_nest = {'ua', 'va'}; slct_level_range_nest = [slct_level_range, slct_level_range];
                        temp_CMIPdata = read_CMIP('project_name', project_name, 'CMIP_folder',  CMIP_folder, 'model_list',  model, 'var_list', var_list_nest,...
                                             'exp_list', exp, 'slct_year', slct_year, 'noleap', noleap, 'slct_level_range', slct_level_range_nest, 'units_level', units_level,...
                                             'ori_field', ori_field, 'output_field', output_field, 'output_folder', output_folder,...
                                             'compute_MME', 0, 'MME_name', [], 'MME_lon_res', MME_lon_res, 'MME_lat_res', MME_lat_res, 'interp_method', interp_method, 'NumWorkers', NumWorkers,...
                                             'savemat', 1, 'recompute', 0, 'recompute_MME', 0, 'load_or_not', 1);

                        %/ Divergence (du/dx + dv/dy)
                        lon = temp_CMIPdata.(model_strrep).lon;
                        lat = temp_CMIPdata.(model_strrep).lat;
                        
                        u = temp_CMIPdata.(model_strrep).(strcat('ua',num2str(slct_level_range{:}),'_',exp_strrep));
                        v = temp_CMIPdata.(model_strrep).(strcat('va',num2str(slct_level_range{:}),'_',exp_strrep));

                        hx = diff(lon(1:2))*pi/180;   %/ rmb to change degree to radian!!!!
                        [~, dAdx, ~] = gradient(u, hx);
                        clear u;

                        hy = diff(lat(1:2))*pi/180;   %/ rmb to change degree to radian!!!!
                        [dAdy, ~, ~] = gradient(v, hy);
                        clear v;

                        [~, lat_2D] = meshgrid(lon, lat);
                        lat_2D = lat_2D';
        
                        r = 6371e3;   %/ in m
                        div = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D));  %/ == du/dx + dv/dy in cartesian coor.

                        cond = isinf(div);
                        div(cond) = nan;  %/ set inf to nan, as it occurs at lat = 90 or -90 when multiplying with 1./(r*cosd(lat_2D))!

                        CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = div;
                        CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = 's^{-1}';  

                        %/ Copy the rest of the fields
                        other_flds = {'lon', 'lat', 'ens', date_fld, 'institution', 'institution_id'};
                        for f = 1:length(other_flds)
                            CMIPdata.(model_strrep).(other_flds{f}) = temp_CMIPdata.(model_strrep).(other_flds{f});
                        end
                        clear temp_CMIPdata;
                        clear u v dAdx dAdy div;

                    else
                        %/ Find the available ensembles in the path
                        raw_filename_ens_cell = cell(length(CMIP_folder), 1);
                        for ii = 1:length(CMIP_folder)
                            raw_file = dir(fullfile(CMIP_folder{ii}, strcat(var,'_',ts,'_',model,'_',exp,'_*_',grid,'_*')));
                            raw_filename_parts        = cellfun(@(x) strsplit(x, '_'), {raw_file.name}, 'UniformOutput', false);
                            raw_filename_ens_cell{ii} = cellfun(@(x) x{end-2}, raw_filename_parts, 'UniformOutput', false)'; %/ Transpose to column vector before storing into the cell array
                        end
                        raw_filename_ens = cat(1, raw_filename_ens_cell{:});
                        raw_filename_ens_unique = unique(raw_filename_ens);

                        %/ If slct_ens is not given, we determine it based on data availability on the server
                        if isempty(slct_ens)
                            ens_list = raw_filename_ens_unique; 
                        else
                            missing_ens = slct_ens_indiv;
                            missing_ens(ismember(slct_ens_indiv, raw_filename_ens_unique)) = [];
                            
                            if ~isempty(missing_ens) 
                                missing_ens_msg = sprintf('No ensemble %s for %s %s %s %s at %s on disk!', strjoin(missing_ens,' '), model, var, ts, exp, grid);
                                warning(missing_ens_msg);
                                problem_files = [problem_files; missing_ens_msg];
                                continue;
                            end
                            ens_list = raw_filename_ens_unique(ismember(raw_filename_ens_unique, slct_ens_indiv)); %/ select the desire ensemble(s) for the model (why? since some models output data with r1i1p1f2)
                        end
                        data_EM = []; 
                        CMIPdata.(model_strrep) = [];  %/ clear up the field 
                        flag_skip_the_model = 0;

                        %/ Loop over all available ensembles on disk
                        for e = 1:length(ens_list) %/ compute ME
                            if flag_skip_the_model   
                                break;  
                            end
                            ens = ens_list{e};
                            
                            %/ Get the paths to the target files
                            name_only_cell = cell(length(CMIP_folder), 1);
                            filename_list_cell = cell(length(CMIP_folder), 1);
                            for ii = 1:length(CMIP_folder)
                                listing       = dir(CMIP_folder{ii});
                                filename_list = {listing.name};
                                prefix        = strcat(var,'_',ts,'_',model,'_',exp,'_',ens,'_',grid,'_');
                                fprintf('*** Searching files in %s ... ***\n', fullfile(CMIP_folder{ii},prefix))
                                logical_1D    = contains(filename_list, prefix);
                                name_only_cell{ii} = filename_list(logical_1D)';
                                filename_list_cell{ii} = fullfile(CMIP_folder{ii}, filename_list(logical_1D))'; %/ transpose to column vector before storing into the cell array
                            end
                            name_only     = cat(1, name_only_cell{:});
                            filename_list = cat(1, filename_list_cell{:});
                            [~,ia,~]      = unique(name_only);    %/ [IMPORTANT] Check if there are duplicated files (it happens when both folders store the same data)
                            filename_list = filename_list(ia);        %/ [IMPORTANT] Then remove the duplicated ones

                            % A = {'aaa', 'bbb', 'ccc','aaa', 'bbb', 'ddd'};
                            % [C,ia,ic] = unique(A)
                            % A(ia)    % == C
                            % C(ic)    % == A

                            if isempty(filename_list)
                                error('No such files in the directory %s with a prefix ''%s''', strjoin(CMIP_folder, ' or '), prefix);  
                            end
                            
                            %/ File loop
                            data = []; dates = [];
                            for f = 1:length(filename_list)
                                filename = filename_list{f};
                                fprintf('*** Loading %s... ***\n', filename);
                                
                                %/ Some model files are overlapping in dates, skip them.
                                if isequal(ori_field, 'monthly')
                                    if  isequal(model, 'FGOALS-g3') && contains(filename, '201001-201412')   || ...
                                        isequal(model, 'GISS-E2-1-G') && contains(filename, '194101-199012') || ...
                                        isequal(model, 'NorESM2-LM') && contains(filename, '204001-205012')
                                        continue;  
                                    end
                                end
                                
                                %/ Get the basic info (NOTE: the reference period can be different for different files (e.g., BCC-CSM2-MR)
                                flag_serious_error = 0;
                                % try
                                    %/ First, convert time into dates
                                    raw_time             = double(ncread(filename,'time'));
                                    raw_time_units       = ncreadatt(filename,'time','units');
                                    raw_time_units_parts = strsplit(raw_time_units,' '); %/ split the string by ' '
                                    if contains(raw_time_units_parts{1}, 'day')
                                        time_units_conv = 1;
                                    elseif contains(raw_time_units_parts{1}, 'hour')
                                        time_units_conv = 1/24;
                                    elseif contains(raw_time_units_parts{1}, 'minute')
                                        time_units_conv = 1/24/60;
                                    elseif contains(raw_time_units_parts{1}, 'second')
                                        time_units_conv = 1/24/60/60;
                                    else
                                        error('code not set for %s!', raw_time_units_parts);     
                                    end
                                    raw_ref_date   = datetime(raw_time_units_parts{3},'InputFormat','yyyy-MM-dd');

                                    filename_parts = strsplit(filename, '_');
                                    filename_date  = strsplit(erase(filename_parts{end},'.nc'),'-');  %/ Since CMIP models use a 365-day calender, we first obtain the dates range: e.g., 245101-250012
                                    
                                    if numel(filename_date{1}) ~= LEN || numel(filename_date{2}) ~= LEN
                                        warning('filename contains dates that is not in %s. Likely this is NOT a RAW file. Skipping it.', ori_date_format);
                                        continue;
                                    end
                                    
                                    t1 = datetime(filename_date{1}, 'InputFormat',ori_date_format);
                                    t2 = datetime(filename_date{2}, 'InputFormat',ori_date_format);
                                    if isequal(ori_field, 'monthly')
                                        t0 = raw_ref_date;  %/ Should start from the reference dates (t0)
                                        tv = (t0:calmonths(1):t2)';
                                    elseif isequal(ori_field, 'daily')
                                        t0 = raw_ref_date;  %/ Should start from the reference dates (t0)
                                        tv = (t0:caldays(1):t2)';
                                    elseif isequal(ori_field, '3hr')
                                        HH = t1.Hour;
                                        mm = t1.Minute;
                                        t0 = raw_ref_date + hours(HH) + minutes(mm);  %/ For subdaily data, it may not start from 00:00; we have to check the time when t1 begins
                                        tv = (t0:hours(3):t2)';
                                    else
                                        error('code not set for %s!', ori_field);     
                                    end

                                    model_calendar = ncreadatt(filename,'time','calendar');
                                    if ismember(model_calendar, {'365_day', 'noleap'})
                                        flag_365day = 1;
                                    elseif ismember(model_calendar, {'proleptic_gregorian', 'gregorian', 'standard', 'julian'})
                                        flag_365day = 0;
                                    else
                                        flag_serious_error = 1; %/ Make it a serious error to stop the program
                                        error('Code not ready for the calendar ''%s''', model_calendar);
                                    end

                                    ind_extra_leap = [];
                                    if flag_365day
                                        if isequal(ori_field, 'monthly')
                                            leapday = double((leapyear(tv.Year) & tv.Month == 2));  %/ A column vector of 0/1 to shift the raw_time (on a 365-day calender) element-wise
                                            
                                        elseif contains(ori_field, {'daily', 'hr'})
                                            leapday = double((tv.Month == 2 & tv.Day == 29));       %/ A column vector of 0/1 to shift the raw_time (on a 365-day calender) element-wise
                                        end
                                        %/ [IMPORTANT] accumulate the leapday since the reference dates
                                        %/  so that it can be correctly shifted for all years in the period
                                        leapday_cumsum = cumsum(leapday);
        
                                        %/ Finally, retrieve the target period (t1 to t2)
                                        ind_t1 = find(tv == t1);
                                        leapday_cumsum = leapday_cumsum(ind_t1:end);

                                        if length(raw_time) == length(leapday_cumsum)
                                            raw_date_dt = int2datetime(str2num(datestr(raw_ref_date + (raw_time + leapday_cumsum)*time_units_conv,'yyyymmdd')), 'yyyyMMdd');
                                        else
                                            warning('For some reasons, length(raw_time) ~= length(leapday_cumsum). Now generating raw_dates based on the filename instead... ');
                                            if isequal(ori_field, 'monthly')
                                                raw_date_dt = (t1:calmonths(1):t2)';
                                                
                                            elseif isequal(ori_field, 'daily')
                                                raw_date_dt = (t1:caldays(1):t2)';
                                                if flag_365day
                                                    raw_date_dt(raw_date_dt.Month == 2 & raw_date_dt.Day == 29) = [];
                                                end
                                                
                                            elseif isequal(ori_field, '3hr')
                                                raw_date_dt = (t1:hours(3):t2)';
                                                if flag_365day
                                                    raw_date_dt(raw_date_dt.Month == 2 & raw_date_dt.Day == 29) = [];
                                                end
                                            end
                                        end
                                    elseif ismember(model_calendar, {'julian'}) 
                                        %/=========================================================================================
                                        %/ GOAL: To convert julian into standard calendar by moving the extra leap days from data.
                                        %/=========================================================================================
                                        years = (t1.Year:t2.Year)';
                                        julian_dates   = date_array_gen('year_list', years, 'calendar', model_calendar); %/ date in int
                                        raw_date_dt    = date_array_gen('year_list', years, 'calendar', []);             %/ date in int
                                        ind_extra_leap = find(~ismember(julian_dates, raw_date_dt));  %/ length(julian_dates) must >= length(standard_dates)
                                        raw_date_dt    = int2datetime(raw_date_dt, 'yyyyMMdd');       %/ Convert date (int) into datetime
                                    else
                                        leapday_cumsum = 0; % no change to raw_time
                                        raw_date_dt = int2datetime(str2num(datestr(raw_ref_date + (raw_time + leapday_cumsum)*time_units_conv,'yyyymmdd')), 'yyyyMMdd');
                                    end
                                    raw_dates = datetime2int(raw_date_dt, ori_date_format);
                                    disp(raw_dates([1:3, end-2:end]));
                                    
                                    %/ [IMPORTANT] Skip if the file if no data for the requested dates -> save RAM!
                                    if ~any(ismember(raw_dates, slct_ori_dates))
                                        fprintf('*** File (%d-%d) skipped due to no data for the requested dates (%d-%d). ***\n', raw_dates(1),raw_dates(end),slct_ori_dates(1),slct_ori_dates(end))
                                        continue;   
                                    end

                                    %/ Get the field names in NetCDF data
                                    info = ncinfo(filename);
                                    field_names = {info.Variables.Name}; 

                                    %/ [IMPORTANT] Loop over to check which name of lon lat the file is used.
                                    lon_names = {'lon', 'longitude', 'nav_lon'};
                                    lat_names = {'lat', 'latitude',  'nav_lat'};
                                    lon = []; 
                                    for i = 1:length(lon_names)
                                        if ismember(lon_names{i}, field_names)
                                            lon = double(ncread(filename, lon_names{i}));
                                            break;
                                        end
                                    end
                                    lat = [];
                                    for i = 1:length(lat_names)
                                        if ismember(lat_names{i}, field_names)
                                            lat = double(ncread(filename, lat_names{i}));
                                            break;
                                        end
                                    end
 
                                    if isempty(lon) || isempty(lat)
                                        ncdisp(filename);
                                        error('Missing lon, lat in %s! Check their exact names!', filename)
                                    end

                                    %/ Use ncinfo to auto check the dimension of the data (ncread is slow!)
                                    finfo   = ncinfo(filename);
                                    ind_var = findismember_loop({finfo.Variables.Name}, var);
                                    ndim    = length(finfo.Variables(ind_var).Size);
                                    if ndim == 4
                                        level = int64(ncread(filename, 'plev'));  %/ int64; in Pa, *in descending order*
                                        if numel(slct_level) == 1
                                            st_level    = findismember_loop(level, slct_level*unitconv_level);
                                            count_level = 1;
                                        else
                                            %/ First, double-check if the level range is valid
                                            if isempty(find(level == min(slct_level)*unitconv_level, 1))
                                                disp(level)
                                                error('No p-level corresponds to the queried min level of %d %s! Revise ''slct_level_range'' or check your data!', min(slct_level), units_level);
                                            end
                                            if isempty(find(level == max(slct_level)*unitconv_level, 1))
                                                disp(level)
                                                error('No p-level corresponds to the queried max level of %d %s! Revise ''slct_level_range'' or check your data!', max(slct_level), units_level);
                                            end
                                            ind         = find(level >= min(slct_level)*unitconv_level & level <= max(slct_level)*unitconv_level);
                                            st_level    = ind(1);
                                            count_level = length(ind);
                                        end
                                        if isempty(st_level)
                                            disp(level)
                                            error('No %.2f in the p-level data!', max(slct_level)*unitconv_level);
                                        end
                                        level    = int64(ncread(filename, 'plev', st_level, count_level)); %/ Updated; int64; in Pa; 
                                        raw_data = double(ncread(filename, var, [1, 1, st_level, 1], [Inf, Inf, count_level, Inf])); %/ Subset the selected level
                                        raw_data(:,:,:,ind_extra_leap) = [];   %/ if model_calendar == 'julian'
                                    else
                                        raw_data  = double(ncread(filename, var));
                                        raw_data(:,:,ind_extra_leap) = [];     %/ if model_calendar == 'julian'
                                    end 
                                    % size(raw_data)
                                    
                                    %/ Replace missing/unreasonably large values with NaN (for easy coding)
                                    missing_val = ncreadatt(filename,var,'missing_value');
                                    if ~isempty(find(raw_data == missing_val, 1))
                                        warning('Detected %d missing values. Replacing them with NaNs...', length(find(raw_data == missing_val)))
                                        raw_data(raw_data == missing_val) = nan;
                                    end

                                    if ~isempty(find(abs(raw_data) > 1e10, 1))
                                        warning('Detected %d unreasonably large values (abs > 1e10)! Likely they are also missing values. Replacing them with NaNs...', length(find(abs(raw_data) > 1e10)));
                                        raw_data(abs(raw_data) > 1e10) = nan;
                                    end

                                    %/ Double check if there is still any unreasonably large value
                                    a = max(abs(raw_data), [], 'all', 'omitnan');
                                    if a > 1e10   
                                        flag_serious_error = 1;
                                        error('*** DOUBLE-CHECK: max(abs(raw_data), [], ''all'', ''omitnan'') = %.3g ***\n', a);
                                    else
                                        fprintf('*** DOUBLE-CHECK: max(abs(raw_data), [], ''all'', ''omitnan'') = %.3g ***\n', a);
                                    end
                                    
                                    %/ Convert subdaily to daily (b4 setting -ve pr and evspsbl to zeros for consistency)
                                    if isequal(ori_field, '3hr')
                                        %/ NOTE: CMIP6 output pr (or the like) in unit of flux (kg m-2 s-1), 
                                        %/ So they are instantaneous rate -> ins
                                        ins_or_acc = 'ins';
                                        
                                        %/ raw_dates should be in yyyyMMddHHmm
                                        raw_data  = subdaily2daily('subdaily_data', raw_data, 'dates', raw_dates, 'ins_or_acc', ins_or_acc);
                                        raw_dates = unique(raw_dates./1e4);  %/ yyyyMMddHHmm -> yyyyMMdd, don't forget to take the unique values!
                                    end

                                    %/ Convert daily to monthly (b4 setting -ve pr and evspsbl to zeros for consistency)
                                    if isequal(ori_field, 'daily') && isequal(output_field, 'monthly')
                                        %/ NOTE: CMIP6 output pr (or the like) in unit of flux (kg m-2 s-1), 
                                        %/ So they are instantaneous rate -> ins
                                        ins_or_acc = 'ins';
                                        
                                        %/ raw_dates should be in yyyyMM
                                        [raw_data, raw_dates] = daily2monthly('daily_data', raw_data, 'dates', raw_dates, 'ins_or_acc', ins_or_acc, 'noleap', flag_365day);
                                        % raw_dates = unique(raw_dates./1e2);  %/ yyyyMMdd -> yyyyMM, don't forget to take the unique values!
                                    end
                                    
                                    %/ Set -ve pr and evspsbl to zeros
                                    if ismember(var, {'pr', 'evspsbl'})
                                        if ~isempty(find(raw_data < 0, 1)) 
                                            warning('Variable %s contains %d negative values (mean: %.2g)! Correcting them to zeros...',...
                                                var, numel(find(raw_data < 0)), mean(raw_data(raw_data < 0)));
                                            raw_data(raw_data < 0) = 0;  %/ Correction
                                        end
                                    end

                                    %/ Show the # of NaNs
                                    if ~isempty(find(isnan(raw_data), 1)) 
                                        % raw_data(:,:,1,10)
                                        warning('Variable %s contains %d NaN values!', var, numel(find(isnan(raw_data))))
                                    end
                                    
                                    data = cat(ndim,  data, raw_data);  %/ concatenate along the last dimension (assume to be the time dimension)
                                    dates = cat(1,  dates, raw_dates);

                                % catch ME
                                %    if flag_serious_error
                                %         error('Error from reading %s: %s', model, ME.message);
                                %    else
                                %         warning('Error from reading %s:', model);
                                %         warning('%s; skip reading this model.', filename, ME.message);
                                %    end
                                %    flag_skip_the_model = 1;
                                %    problem_files = [problem_files; filename];
                                %    break;
                                % end
                            end
                            
                            %/ Subset dates (IMPORTANT: it avoids any incomplete data in the last year)
                            if ~isempty(slct_year) 
                                if isempty(dates)
                                    flag_skip_the_model = 1;
                                    
                                    missing_ens_msg = sprintf('Missing years of %d-%d in ensemble %s for %s %s %s %s at %s on disk!', slct_year(1), slct_year(end), ens, model, var, ts, exp, grid);
                                    warning(missing_ens_msg);
                                    problem_files = [problem_files; missing_ens_msg];
                                    % problem_files = [problem_files; filename];
                                    break;
                                end

                                if isequal(output_field, 'daily')
                                    if flag_365day == 1 && noleap == 0
                                        %/ Remove leap days from slct_dates (for double-checking)
                                        slct_output_dates(mod(slct_output_dates,1e4) == 229) = [];  %/ for noleap calendar, remove Feb 29 from slct_dates 
                                        
                                    elseif flag_365day == 0 && noleap == 1 
                                        %/ Remove leap days from the raw data and dates
                                        ind_leaps = find(mod(dates,1e4) == 229);
                                        dates(ind_leaps) = [];
                                        if ndim == 3
                                            data(:,:,ind_leaps) = [];
                                        elseif ndim == 4
                                            data(:,:,:,ind_leaps) = [];
                                        else
                                            error('code not ready!')
                                        end
                                    end
                                end
                                ind_dates = findismember_loop(dates, slct_output_dates);

                                %/ Check if all the queried dates are available (for original daily data)
                                if numel(ind_dates) > numel(slct_output_dates)
                                    %/ find(diff(floor(dates/100)) == 0)
                                    error('Seems like there are overlapping dates among the raw files! Please check!');
                                    
                                elseif numel(ind_dates) < numel(slct_output_dates)
                                    warning('Seems like not all the queried dates are available! Please check! Skipping this model...');
                                    flag_skip_the_model = 1;  %/ break the loop and proceed to the next model
                                    problem_files = [problem_files; filename];
                                    break;
                                end

                                %/ Update the dates
                                dates = dates(ind_dates);  

                                %/ Update the data
                                if ndim == 3
                                    data = data(:,:,ind_dates);     
                                elseif ndim == 4
                                    data = data(:,:,:,ind_dates);  
                                else
                                    error('code not set!')
                                end
                            end
    
                            %/ Unit Conversion
                            unit_addition = 0;
                            unit_multiply = 1;
                            raw_data_units = ncreadatt(filename, var, 'units'); 
                            if ismember(var, {'pr', 'evspsbl'}) && isequal(raw_data_units, 'kg m-2 s-1')
                                if isequal(output_field, 'monthly') 
                                    eom           = eomday(2022,1:12)';      %/ Consider the diffrence in the days of month (NOTE: CMIP6 does not simulate the leap day)                     
                                    unit_multiply = repmat(3600*24.*eom, numel(dates)/12, 1);   %/ Replicate the 12-element array for all years.
                                    unit_multiply = reshape(unit_multiply, [ones(1,ndim-1), length(unit_multiply)]);  %/ Reshape for element-wise multiplication
                                    data_units    = 'mm month^{-1}';
                                elseif isequal(output_field, 'daily') 
                                    unit_multiply = 3600*24;       %/ multiply for a day 
                                    data_units    = 'mm day^{-1}';
                                else
                                    error('code not set!'); 
                                end
                            elseif ismember(var, {'tos'}) && isequal(raw_data_units, 'degC')
                                unit_addition = 273.15;  %/ Convert degC to Kelvin
                                data_units = 'K';
                            else
                                data_units = raw_data_units;
                            end
                            data = data.*unit_multiply + unit_addition;
                            fprintf('Mean of data: %.3g %s\n', mean(data,'all','omitnan'),data_units);
                            fprintf('Min  of data: %.3g %s\n', min(data,[],'all','omitnan'),data_units);
                            fprintf('Max  of data: %.3g %s\n', max(data,[],'all','omitnan'),data_units);
                            fprintf('dates([1:3, end-2:end]):\n');
                            disp(dates([1:3, end-2:end]));
    
                            data_EM = cat(ndim+1, data_EM, data);  %/ concantenate along a new dimension.
                        end

                        %/ flag_skip_the_model == 0 means data reading has completed successfully
                        %/ Now storing data_EM to 'CMIPdata'
                        if flag_skip_the_model == 0

                            %/ Ensemble mean (based on realizations)
                            data_EM = mean(data_EM, ndim+1);
                            a = max(abs(data_EM), [], 'all', 'omitnan');
                            if a > 1e10   
                                error('*** DOUBLE-CHECK: max(abs(data_EM), [], ''all'', ''omitnan'') = %.3g %s ***\n', a,data_units);
                            else
                                fprintf('*** DOUBLE-CHECK: max(abs(data_EM), [], ''all'', ''omitnan'') = %.3g %s ***\n', a,data_units);
                            end
                            
                            %/ Seems no need now. 
                            % %/ [IMPORTANT] Finally, set data on land to nan if it is an ocean variable 
                            % %/             Why? Because some models output non-nan SST data over land!!
                            % if find(ismember(var, ocean_list))
                            %     fprintf('!!! Detected that ''%s'' is an ocean variable. Setting land values to NaN... !!!\n', var);
                            %     slct_reg = 'land';
                            %     [cond_land, ~, ~, ~, ~] = reg_extractor('lon', lon, 'lat', lat, 'slct_reg', slct_reg, 'saveload_cond_landocean', 0, 'data_folder', output_folder);
                            %     cond_land(isnan(cond_land)) = 0; %/ Convert nan to 0 since reg_extractor sets unassigned land/ocean grids to NaN. 
                            %     cond_land_2Dto1D = reshape(cond_land, [], 1);
                            %     ind_land = find(cond_land_2Dto1D);
                            % 
                            %     sz = size(data_EM);
                            %     if ndim == 3
                            %         data_EM = reshape(data_EM, sz(1)*sz(2), sz(3));
                            %         data_EM(ind_land,:) = nan;  
                            %         data_EM = reshape(data_EM, sz(1), sz(2), sz(3));
                            %     elseif ndim == 4
                            %         data_EM = reshape(data_EM, sz(1)*sz(2),  sz(3), sz(4));
                            %         data_EM(ind_land,:,:) = nan;  
                            %         data_EM = reshape(data_EM, sz(1), sz(2), sz(3), sz(4));
                            %     end
                            % end

                            %/ Store into CMIPdata
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep)) = squeeze(data_EM); %/ Finally, squeeze out any singlton dimension
                            CMIPdata.(model_strrep).(strcat(var_strrep,'_units')) = data_units;
                            CMIPdata.(model_strrep).('lon') = lon;
                            CMIPdata.(model_strrep).('lat') = lat;
                            CMIPdata.(model_strrep).(strcat(var,'_lon')) = lon;  %/ Since different var may have different griddings, even for the same model and exp.
                            CMIPdata.(model_strrep).(strcat(var,'_lat')) = lat;  %/ Since different var may have different griddings, even for the same model and exp.
                            if ndim == 4
                                CMIPdata.(model_strrep).('level') = level;
                            end
                            CMIPdata.(model_strrep).('ens') = ens_list;
                            CMIPdata.(model_strrep).(date_fld) = dates;  
                            
                            %/ Finally, read some global attributes
                            globalattr = {'institution', 'institution_id'};
                            for g = 1:length(globalattr)
                                CMIPdata.(model_strrep).(globalattr{g}) = ncreadatt(filename,"/",globalattr{g});
                            end
                        end
                    end

                    %/ Save as mat (if no errors)
                    if savemat && flag_skip_the_model == 0
                        disp(CMIPdata.(model_strrep));
                        data       = CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep));
                        data_units = CMIPdata.(model_strrep).(strcat(var_strrep,'_units'));
                        fprintf('Mean of data: %.3g %s\n', mean(data,'all','omitnan'),data_units);
                        fprintf('Min  of data: %.3g %s\n', min(data,[],'all','omitnan'),data_units);
                        fprintf('Max  of data: %.3g %s\n', max(data,[],'all','omitnan'),data_units);
                        fprintf('dates([1:3, end-2:end]):\n');
                        disp(CMIPdata.(model_strrep).(date_fld)([1:3, end-2:end]));

                        S = CMIPdata.(model_strrep);
                        tic;
                        disp(['Writing into ', output_filename, ' ...'])
                        save(output_filename, 'S','-v7.3');
                        clear S;
                        disp(strcat({'Time taken: '},num2str(toc),{' sec for writting.'}))
                    end
                end

                if load_or_not == 0 && flag_to_compute_MME == 0
                    CMIPdata.(model_strrep) = []; %/ Then clear up the field (spare memory)
                end
            end
            problem_files = unique(problem_files);  %/ remove redundant elements

            %/ Multi-model ensemble mean (MME) 
            if ~isempty(problem_files)
                if delete_prob_files
                    for i = 1:length(problem_files)
                        fprintf('*** Deleting problematic file: %s ***\n', problem_files{i});
                        delete(problem_files{i})
                    end
                end
            else
                if flag_to_compute_MME
                    if flag_skip_the_model
                        warning('Since there is at least one model''s variable being skipped due to data availability, no MME will be computed.');
                        continue;
                    end
                    
                    MME = 0; nonNaN_freq = 0;
                    for m = 1:length(model_list)
                        model        = model_list{m};
                        fprintf('Interpolating %s...\n', model)
                        model_strrep = strrep(model, '-', '_');  %/ Sadly, '-' is not allowed in the field name. Replacing it with '_'.
                        data         = CMIPdata.(model_strrep).(strcat(var_strrep,'_',exp_strrep));
                        data_units   = CMIPdata.(model_strrep).(strcat(var_strrep,'_units'));
                        a = max(abs(data), [], 'all', 'omitnan');
                        if a > 1e10  
                            error('*** DOUBLE-CHECK: max(abs(data), [], ''all'', ''omitnan'') = %.3g %s ***\n', a, data_units);
                        else
                            fprintf('*** DOUBLE-CHECK: max(abs(data), [], ''all'', ''omitnan'') = %.3g %s ***\n', a, data_units);
                        end

                        %/ Search for lon/lat exclusive for the variable (e.g. 'hus_lon', 'hus_lat')
                        %/ Otherwise, use ordinary 'lon' and 'lat' fields
                        if isfield(CMIPdata.(model_strrep), strcat(var,'_lon')) && ...
                           isfield(CMIPdata.(model_strrep), strcat(var,'_lat'))
                            lon_old = CMIPdata.(model_strrep).(strcat(var,'_lon'));
                            lat_old = CMIPdata.(model_strrep).(strcat(var,'_lat'));
                        else
                            lon_old = CMIPdata.(model_strrep).lon;
                            lat_old = CMIPdata.(model_strrep).lat;
                        end
                        lon_new = 0:MME_lon_res:360-MME_lon_res;               %/ exclude 360E
                        lat_new = -90+MME_lat_res:MME_lat_res:90-MME_lat_res;  %/ exclude two poles
    
                        is_global = 1;
                        lon_dim   = 1;
                        data(isinf(data)) = nan;  %/ Set inf as nan (Since variables like E/P may give Infs)
    
                        data_interp = my_interp('lon_old', lon_old, 'lat_old', lat_old, 'data', data,...
                                                'lon_new', lon_new, 'lat_new', lat_new, 'is_global', is_global, 'lon_dim', lon_dim,...
                                                'interp_method', interp_method, 'NumWorkers', NumWorkers);
                        size(data_interp)
                        if ~isempty(find(isnan(data_interp), 1)) 
                            warning('%s %s data contains %d NaN after interpolation! Replacing them with 0.', model, var, length(find(isnan(data_interp))));
                        end
                        a = max(abs(data_interp), [], 'all', 'omitnan');
                        if a > 1e10   
                            error('*** DOUBLE-CHECK: max(abs(data_interp), [], ''all'', ''omitnan'') = %.3g %s ***\n', a, data_units);
                        else
                            fprintf('*** DOUBLE-CHECK: max(abs(data_interp), [], ''all'', ''omitnan'') = %.3g %s ***\n', a, data_units);
                        end
    
                        logical_nonNaN = ~isnan(data_interp); %/ get a 3D logical for non-NaN data
                        data_interp(~logical_nonNaN) = 0;     
                        nonNaN_freq = nonNaN_freq + double(logical_nonNaN); %/ frequency (for taking average)
    
                        MME = MME + data_interp;  %/ reduction
                    end
                    if ~isempty(find(isnan(MME), 1)) 
                        error('MME still contains NaNs!'); 
                    end
                    MME = MME./nonNaN_freq; %/ Take average (equal weightings); 
                                            %/ It may then contain NaNs in places where all models have NaNs.
                    ndim = length(size(MME));
                    CMIPdata.(MME_name) = [];
                    CMIPdata.(MME_name).(strcat(var_strrep,'_',exp_strrep)) = MME;
                    CMIPdata.(MME_name).(strcat(var_strrep,'_units')) = CMIPdata.(model_strrep).(strcat(var_strrep,'_units'));
                    CMIPdata.(MME_name).('lon') = lon_new;
                    CMIPdata.(MME_name).('lat') = lat_new;
                    if ndim == 4
                        CMIPdata.(MME_name).('level') = CMIPdata.(model_strrep).('level');
                    end
                    CMIPdata.(MME_name).('ens_list') = slct_ens;
                    CMIPdata.(MME_name).('model_list') = model_list;
                    CMIPdata.(MME_name).(date_fld) = CMIPdata.(model_strrep).(date_fld);
        
                    % flds = {strcat(var_strrep,'_',exp_strrep), strcat(var_strrep,'_units'), 'lon', 'lat', 'ens_list', 'model_list'};
                    % if ndim == 4
                    %     flds{end+1} = 'level';
                    % end

                    %/ Save as mat
                    if savemat
                        S = CMIPdata.(MME_name);
                        tic;
                        disp(['Writing into ', MME_output_filename, ' ...'])
                        save(MME_output_filename, 'S','-v7.3');
                        clear S;
                        disp(strcat({'Time taken: '},num2str(toc),{' sec for writting.'}))
                    end
                    % disp(CMIPdata.(MME_name));
    
                    data = CMIPdata.(MME_name).(strcat(var_strrep,'_',exp_strrep));
                    data_units = CMIPdata.(MME_name).(strcat(var_strrep,'_units'));
                    fprintf('Mean of data: %.3g %s\n', mean(data,'all','omitnan'),data_units);
                    fprintf('Min  of data: %.3g %s\n', min(data,[],'all','omitnan'),data_units);
                    fprintf('Max  of data: %.3g %s\n', max(data,[],'all','omitnan'),data_units);
                    fprintf('dates([1:3, end-2:end]):\n');
                    disp(CMIPdata.(MME_name).(date_fld)([1:3, end-2:end]));
                end
            end
        end
    end

    if ~isempty(problem_files)
        fprintf('!!! read_CMIP is not completed successfully, please check the following problem files (and their siblings) !!!\n')
        disp(problem_files)
    else
        fprintf('!!! CONGRADULATIONS: read_CMIP is completed successfully !!!\n')
    end
%     disp(CMIPdata);
end
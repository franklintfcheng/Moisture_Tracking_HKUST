function [dataset, flag_exist, matfilename_full] = check_n_load_field(varargin)

    pnames       = {   'dataset', 'data_folder', 'var', 'model', 'exp', 'project_name', 'select_field',...
                       'noleap', 'middle_str', 'str_years', 'boxreg_lon', 'boxreg_lat', 'matfilename',...
                       'stlevel', 'edlevel', 'search_dataset', 'search_datafolder',...
                       'load_or_not',};
    
    dflts        = cell(length(pnames),1);
    
    [                   dataset,   data_folder,   var,   model,    exp,  project_name,  select_field,...
                        noleap,  middle_str,    str_years,  boxreg_lon,  boxreg_lat, matfilename,...
                        stlevel, edlevel, search_dataset,    search_datafolder,...
                        load_or_not] = ...
            internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
    %%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: Jan 10, 2024
    %/
    %/ Description: This function is to check if the request data field 
    %/              'select_field' of the variable 'var' exists in 
    %/              'dataset' or under 'data_folder'.
    %/
    %/      INPUT: 
    %/            'select_field': A string or cell
    %/                   'model': CMIP6 model name
    %/                     'exp': CMIP6 model experiment
    %/
    %/      OUTPUT:      
    %/                 'dataset': The dataset will contain the requested
    %/                            field if flag_exist == 1
    %/
    %/              'flag_exist': 0 or 1. Can be an array of the same
    %/                            size as 'select_field'
    %/
    %/        'matfilename_full': The corresponding file path
    %/=====================================================================
    
    if ischar(select_field)
        select_field = {select_field};
    end

    if ~isempty(model) && ~isempty(exp) 
        flag_cmip = 1;
    else
        flag_cmip = 0;
    end

    if ~isempty(stlevel) && ~isempty(edlevel)
        if stlevel > edlevel
            str_level_range = sprintf('_lv%d-%d', edlevel, stlevel);
        else
            str_level_range = sprintf('_lv%d-%d', stlevel, edlevel);
        end
    else
        str_level_range = '';
    end
    
    if isempty(search_dataset)
        search_dataset = 1;
    end
    if isempty(search_datafolder)
        search_datafolder = 1;
    end
    if isempty(load_or_not)
        load_or_not = 1;
    end
    flag_exist = zeros(length(select_field), 1);
    matfilename_full = cell(length(select_field), 1);

    %/ loop over select_field
    for k = 1:length(select_field)
        select_field_each = select_field{k};
        
        if noleap 
            if contains(select_field_each, 'daily')
                str_noleap = '_noleap'; 
            elseif contains(select_field_each, 'pentad')
                warning('%s data will not include the leap days. No need to set noleap = 1.', select_field_each)
                str_noleap = '';
            elseif contains(select_field_each, 'monthly')
                str_noleap = '';
                noleap = 0; 
                warning('For monthly data, resume ''noleap'' to 0.');
            else
                error('code not ready!');
            end
        else
            str_noleap = ''; 
        end

        if ~isempty(matfilename)
            str_matfilename = strcat('_', matfilename);
        else
            str_matfilename = [];
        end

        if isempty(boxreg_lon) && isempty(boxreg_lat)        %/ Set also boxreg_lon = [] and boxreg_lat = [] to load the raw cmip!
            str_domain = '_global'; 
        else
            str_domain = sprintf('_%d-%dE_%d-%dN', boxreg_lon(1),boxreg_lon(2),boxreg_lat(1),boxreg_lat(2)); %/ save processed 4D data -> time-savor!
        end
        
        %/ NOTE: the mat filename for the processed CMIP6 data (see read_CMIP.m) is different.
        if flag_cmip
            ens = read_CMIP_ens('project_name',project_name,'model_list',model,'var',var,'exp',exp,...
                                'ori_field',select_field_each,'remove_nest_cell',1); %/ Read the available ens 
            var_cmip6 = strcat(model,'_',exp,'_',var,'_',strjoin(ens,'_'));

            matfilename_full{k} = strcat(data_folder, var_cmip6, str_level_range, '_', select_field_each, middle_str, str_domain, '_', str_years, str_noleap, str_matfilename, '.mat');
        else
            matfilename_full{k} = strcat(data_folder, var,       str_level_range, '_', select_field_each, middle_str, str_domain, '_', str_years, str_noleap, str_matfilename, '.mat');
        end

        %/ First, check if the field has been loaded into 'dataset'
        dataset.placeholder = []; %/ make sure isfield(dataset, var) will work
        if search_dataset && isfield(dataset, var)
            if flag_cmip
                if isfield(dataset.(var), select_field_each) && isfield(dataset.(var), 'model') && isfield(dataset.(var), 'exp') && isfield(dataset.(var), 'ens')
                    if isequal(dataset.(var).model, model) && isequal(dataset.(var).exp, exp) && isequal(dataset.(var).ens, ens) 
                        flag_exist(k) = 1;
                        fprintf('!!! check_n_load_field: %s %s %s %s %s has already been loaded in ''dataset'' !!!\n', var, select_field_each, model, exp, strjoin(ens,'_') );
                    end
                end
            else
                %/ If matfilename is not empty, check if the data remark matches it.
                if ~isempty(matfilename) && isfield(dataset.(var), 'remark') 
                    if isequal(dataset.(var).remark, matfilename) && isfield(dataset.(var), select_field_each)
                        flag_exist(k) = 1;
                        fprintf('!!! check_n_load_field: %s %s has already been loaded in ''dataset'' !!!\n', var, select_field_each)
                    end
                else
                    if isfield(dataset.(var), select_field_each)
                        flag_exist(k) = 1;
                        fprintf('!!! check_n_load_field: %s %s has already been loaded in ''dataset'' !!!\n', var, select_field_each)
                    end
                end
            end
        end

        %/ Load the mat file if existed
        if flag_exist(k) == 0 && search_datafolder
            if isfile(matfilename_full{k}) 
                if load_or_not
                    fprintf('*** check_n_load_field: Data file found. Loading %s ***\n', matfilename_full{k})
                    load(matfilename_full{k}, 'S');  
                    % matfilename_full{k}
                    flds = fieldnames(S);
                    % disp(flds) %/ Checking

                    if ~isfield(dataset, var)
                        dataset.(var) = [];    %/ only when the field doesn't exist!
                    end
                    remaining_flds = flds;
                    if flag_cmip
                        if isequal(select_field_each, 'daily') || isequal(select_field_each, 'monthly')
                            if isempty(boxreg_lon) && isempty(boxreg_lat)
                                main_fld = strcat(var,'_',strrep(exp, '-', '_')); %/ Then it means we load the raw output of read_CMIP
                            else
                                main_fld = select_field_each;  %/ Otherwise, the data should have been restructured
                            end
                            dataset.(var).(select_field_each) = S.(main_fld); %/ Due to the different struct field format from the daily/monthly CMIP6 data processed by read_CMIP...
                            remaining_flds = setdiff(flds, main_fld);
                        end
                        dataset.(var).model  = model;
                        dataset.(var).exp    = exp;
                        dataset.(var).ens    = ens;
                    end
        
                    for f = 1:length(remaining_flds)
                        if contains(remaining_flds{f}, 'units')
                            dataset.(var).('units') = S.(remaining_flds{f});
                        else
                            dataset.(var).(remaining_flds{f}) = S.(remaining_flds{f});
                        end
                    end
                    dataset.(var).remark = matfilename; %/ To mark it, in case it has not been marked.(should be fine).
                    fprintf('*** Loading done ***\n')
                else
                    fprintf('*** check_n_load_field: Data file found, but skip Loading. ***\n')
                end
                flag_exist(k) = 1;  %/ Update it, since it was found from the post-processed mat file
            else
                warning('!!! check_n_load_field: Data file not found! %s !!!\n', matfilename_full{k})
            end
        end
    end
end

function var_list = cmip_availability(varargin)

    pnames = { 'var_list', 'model', 'exp', 'select_field'};
    
    dflts  = cell(length(pnames), 1);
      
    [           var_list,   model,      ~,   select_field] = ...
            internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
%%
    %/============================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Apr 19, 2024
    %/
    %/ DESCRIPTION:
    %/      This function was *hard-coded* to validate the cmip6 data
    %/  availability.
    %/
    %/=======================

    if ~iscell(var_list)
        var_list = {var_list};
    end
    
    for i = 1:length(var_list)
        if ~ischar(var_list{i})
            %/ Skip if the cell element contains no char (maybe number)
            continue;
        else
            if contains(select_field, 'daily')
                if ismember(var_list(i), {'wap'}) && ...
                    ismember(model, {'EC-Earth3-CC', 'EC-Earth3-Veg', 'EC-Earth3-Veg-LR', 'GISS-E2-1-G'})
                    var_list{i} = ''; %/ UPDATE
                elseif contains(var_list(i), {'tos'}) && ismember(model, {'GISS-E2-1-G', 'FGOALS-g3', 'MIROC-ES2H', 'MIROC-ES2L'})
                    var_list{i} = ''; %/ UPDATE
                end
            end
        end
    end
end
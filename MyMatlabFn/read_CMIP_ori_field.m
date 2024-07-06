function ori_field = read_CMIP_ori_field(varargin)

    pnames = {'model', 'var', 'output_field'};

    dflts =  cell(1, length(pnames));
    [          model,   var,   output_field] ...
                            = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 24 Mar 2024
    %/
    %/ NOTE: If the model has no ensemble for a certain variable, 
    %/       it is not a good practice to input an empty cell for slct_ens{m}.
    %/       Instead, you just shouldn't call this model in the first place.
    %/
    %/=====================================================================
    
    if isequal(output_field, 'daily') && isequal(model, 'GISS-E2-1-G') && isequal(var, 'pr') 
        ori_field = '3hr';   %/ no daily pr for GISS but 3hrly
    else
        ori_field = output_field; 
    end

end
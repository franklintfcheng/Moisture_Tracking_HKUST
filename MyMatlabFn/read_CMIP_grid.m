function grid = read_CMIP_grid(varargin)

    pnames = {'model', 'exp', 'var',  'ori_field'};

    dflts =  cell(1, length(pnames));
    [          model,   exp,   var,    ori_field ] ...
                            = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 29 Feb 2024
    %/
    %/ DESCRIPTION:: This function is to output the correct grid format (grid) 
    %/               string based on the input model, var and ori_field.
    %/
    %/=====================================================================
    
    if ismember(model, {'INM-CM4-8', 'INM-CM5-0'})
        grid = 'gr1';   %/ 'atmos data regridded from Cubed-sphere (c96) to 180,288; interpolation method: conserve_order1'
    
    elseif ismember(model, {'GFDL-ESM4', 'GFDL-CM4'})
        if isequal(var, 'tos')
            grid = 'gr';
        elseif isequal(ori_field, 'monthly')
            grid = 'gr1';
        elseif isequal(exp, 'ssp370')
            grid = 'gr1';
        else
            grid = 'gr2';   %/ 'atmos data regridded from Cubed-sphere (c96) to 180,288; interpolation method: conserve_order1'
        end
        
    elseif contains(model, 'EC-Earth3') || ...
           ismember(model, {'KACE-1-0-G', 'E3SM-2-0', 'IPSL-CM6A-LR', 'IPSL-CM6A-LR-INCA',...
                            'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'CNRM-ESM2-1'})

        if isequal(var, 'tos') && (contains(model, 'EC-Earth3') || ismember(model, {'IPSL-CM6A-LR', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'CNRM-ESM2-1'}))
            grid = 'gn';    
        else
            grid = 'gr';    %/ regridded to a CMIP6 standard 1x1 degree lonxlat grid from the native grid using an area-average preserving method.
        end
    else
        grid = 'gn';    %/ By default, may have exceptions
    end


end



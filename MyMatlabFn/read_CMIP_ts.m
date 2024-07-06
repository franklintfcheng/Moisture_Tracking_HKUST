function ts = read_CMIP_ts(varargin)

    pnames = {'model', 'var',  'ori_field'};

    dflts =  cell(1, length(pnames));
    [          model,    var,   ori_field ] ...
                            = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 29 Feb 2024
    %/
    %/ DESCRIPTION:: This function is to output the correct timescale (ts) 
    %/               string based on the input model, var and ori_field.
    %/
    %/=====================================================================
    
    if isequal(ori_field, 'monthly')
        if ismember(var, {'od550aer'})
            ts = 'AERmon';
        elseif ismember(var, {'tos'})
            ts = 'Omon';
        else
            ts = 'Amon';
        end

    elseif isequal(ori_field, 'daily')
        if ismember(var, {'tos'})
            ts = 'Oday';
        else
            ts = 'day';    %/ IMPORTANT: Do NOT mix 'day' with 'Eday' vars!!
                           %/  They are different/inconsistent and will mess up the analysis!
        end 

    elseif contains(ori_field, 'hr')
        ts = ori_field;

    else
        error('Invalid ori_field ''%s''!', ori_field);
    end


end
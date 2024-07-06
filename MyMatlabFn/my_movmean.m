function data_movmean = my_movmean(varargin)

    pnames = {'data',   'n',    'dim'};
    dflts =  cell(1, length(pnames));
    [          data,     n,      dim] = internal.stats.parseArgs(pnames, dflts, varargin{:});

%%
    %/==========================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 11 Dec 2023
    %/==========================================
    %/
    %/ This function differs from the built-in matlab function "movmean" by
    %/ the way that it handles NaNs. Here, we treat NaNs as boundaries and 
    %/ skip it when the center of the windows is a NaN.

    %/ Determine the indices of NaNs
    logical_nan = isnan(data);

    %/ Compute movmean by omitting nan first
    data_movmean = movmean(data, n, dim, 'omitnan');

    %/ Restore the smoothed values that were NaNs to NaNs
    data_movmean(logical_nan) = nan;
    
end
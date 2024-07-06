%%
function bnd = drawboundary(varargin)
    %/ create a set of valid parameters and their default value
    pnames = {'data2D'};
    dflts  = {[]};
    [data2D] = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

    [ind_row, ind_col] = find(~isnan(data2D));
    if isempty(ind_col)
        error('Please set the area outside your object to be NaN.');
    end

    disp('Drawing the boundary northeastward from the westernmost point ...')
    [ind_Wcorner, ~] = min(ind_col); %/ Assume data2D(lon, lat), then find the westernmost point.
    bnd = bwtraceboundary(~isnan(data2D), [ind_row(ind_Wcorner), ind_col(ind_Wcorner)], 'NE');

end
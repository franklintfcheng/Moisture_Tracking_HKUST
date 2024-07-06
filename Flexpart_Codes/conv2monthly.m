%%
function [X, X_dates] = conv2monthly(varargin)

pnames = {'data',   'dates',   'select_year', 'skip_conv', 'detrend_mode'};
dflts  = {[]        []          []                      0,              0};
[data, dates, select_year, skip_conv, detrend_mode] = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

time_dim = length(size(data)); %/ assume the last dimension is the time dimension.

%/ skip changing data's time scale. But detrend_mode will be still at work if toggled on
if skip_conv    
    X = data;
    X_dates = dates;
else
    if numel(num2str(dates(1))) == 8   
        CONV = 100; 
    elseif numel(num2str(dates(1))) == 6
        CONV = 1;
    else
        error('Check the format of the input date! The input data has to be daily (yyyymmdd)!')
    end

    %/ convert to monthly data
    X_dates = nan(length(select_year)*12, 1);

    if     time_dim == 2   X = nan(length(select_year)*12, 1);
    elseif time_dim == 3   X = nan(size(data,1), size(data,2), length(select_year)*12);  
    end

    tic
    for t = 1:length(select_year)
        for mth = 1:12

            X_dates(mth + (t-1)*12, 1) = select_year(t)*1e2 + mth;

            ind_date = find(floor(dates/CONV) == select_year(t)*1e2 + mth);
            if isempty(ind_date)   error('No date is retrieved!');  end

            if      time_dim == 2   X(mth+(t-1)*12, 1)    = mean(data(ind_date), 'omitnan');
            elseif  time_dim == 3   X(:, :, mth+(t-1)*12) = mean(data(:,:,ind_date), time_dim, 'omitnan');
            else                    error('The function only consider 1D, 2D or 3D data!');
            end
        end
    end
    toc;
end

if detrend_mode 
    [nlon, nlat, ntime] = size(X);

    ind = find(isnan(X));
    if ~isempty(ind)  
        warning('Detected %d NaNs in the data! Converting it into zeros...', length(ind)); 
        X(isnan(X)) = 0;
    end

    if time_dim == 2
        X_2D = X';                                                         %/ transpose to year x seasons
    %         size(X_2D)
        X_2D_detrend = detrend(X_2D);                                      %/ If x is a vector, then detrend subtracts the trend from the elements of x.
        X = X_2D_detrend';                                                 %/ transpose back

    elseif time_dim == 3

        X_2D = reshape(X, [], ntime);                                      %/ reshape   to grids x time
        X_2D = X_2D';                                                      %/ transpose to time x grids (for detrend to perform correctly)
    %         size(X_2D)
        X_2D_detrend = detrend(X_2D);                                      %/ If x is a matrix, then detrend operates on each column separately, subtracting each trend from the corresponding column.
        X_2D_detrend = X_2D_detrend';                                      %/ transpose back to grids x time
        X = reshape(X_2D_detrend, nlon, nlat, ntime);                      %/ reshape back to lon x lat x time
    end
end

end

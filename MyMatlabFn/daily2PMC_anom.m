%%
function [dailyAnom_PMC, PMC] = daily2PMC_anom(varargin)

    % create a set of valid parameters and their default value
    pnames = {'daily_data', 'dates'}; 
    dflts  = {      [],          []};
    [          daily_data,   dates] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %/ NOTE:
    %/       1. This function was called 'calc_PMC_anom'.
    %/          Have used 'Find files' from Matlab editor to confirm that no other .m files use this function 
    %/          -> safe to re-name.
    %/       2. Its sisters: e.g., daily2pentad.

    fprintf('*** Running daily2PMC_anom... ***\n')

    %/ check inputs
    if numel(num2str(dates(1))) ~= 8     error('Check if input date is in yyyymmdd!!');                           end
    if ~isempty(find(isnan(daily_data))) error('daily_data contains NaN!!');                                      end
    if length(size(daily_data)) ~= 3     error('The function handles 3D daily_data only. Modify it if needed.');  end
    
    date_mmdd_OneYr = date_array_gen('years', 1980, 'st_month', 1, 'st_day', 1, 'ed_month', 12, 'ed_day', 31, 'output_date_format', 'yyyymmdd');
    date_mmdd_OneYr = mod(date_mmdd_OneYr, 1e4);  %/ 366 days (we include the leap day)
    nday            = length(date_mmdd_OneYr);

    time_dim     = length(size(daily_data));
    [n1, n2, ~]  = size(daily_data);
    daily_clim   = nan(n1, n2, length(date_mmdd_OneYr));
    for i = 1:nday
        day_ind = find(date == date_mmdd_OneYr(i));
        daily_clim(:,:,i) = mean(daily_data(:,:,day_ind),time_dim);
    end

    % Pentad Moving Mean Clim (PMC)
    if date_mmdd_OneYr(1) == 101 && date_mmdd_OneYr(end) == 1231
        %/ Append 2 days to each of the two ends for 5-day moving mean -> circular moving mean
        a = cat(time_dim, daily_clim(:,:,end-1:end),...
                          daily_clim(:,:,:),...
                          daily_clim(:,:,1:2));
        PMC = movmean(a,5,time_dim); % movmean(X,movingdays,dim)
        PMC = PMC(:,:,3:end-2);
    else
        PMC = movmean(daily_clim,5,time_dim); % movmean(X,movingdays,dim)
    end
    
    % daily_data anom on PMC
    dailyAnom_PMC = nan(size(daily_data));
    for i = 1:nday
        day_ind = find(date == date_mmdd_OneYr(i));
        dailyAnom_PMC(:,:,day_ind) = daily_data(:,:,day_ind) - PMC(:,:,i);
    end

end
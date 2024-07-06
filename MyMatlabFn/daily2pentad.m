function [pentad_clim, pentad_allyr, pentad_anom, pentad_anom_std, date_yyyyptd] = daily2pentad(varargin)
    
    %/ create a set of valid parameters and their default value
    pnames = {'daily_data', 'dates'}; 
    dflts  = {     [],          []};
    %/ parse function arguments
    [daily_data, dates] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %/============================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Jan 10, 2023
    %/
    %/ NOTE:
    %/       1. Its sisters: e.g., daily2PMC_anom.
    %/============================================================================

    fprintf('*** Running daily2pentad... ***\n')
    
    %/ check inputs
    if numel(num2str(dates(1))) ~= 8     
        error('Check if input date is in yyyymmdd!!');  
    end
    if ~isempty(find(isnan(daily_data), 1)) 
        warning('daily_data contains NaN! Check if this is normal.');             
    end
    
    date_mmdd_OneYr = date_array_gen('year_list', 1979, 'st_month', 1, 'st_day', 1, 'ed_month', 12, 'ed_day', 31, 'output_date_format', 'yyyymmdd');
    date_mmdd_OneYr = mod(date_mmdd_OneYr, 1e4);  %/ 365 days (so we ignore the leap day)
    
    year         = int64(unique(floor(dates/1e4))); %/ make it in int64
    nyr          = length(year);
    nptd         = 73;
    
    %/ Since dates are in int64, we write int64([1:nptd] make sure the class is consistent
    % date_yyyyptd = reshape(repmat(year, 1, 73)', [], 1)*1e2 + reshape(repmat(int64(1:nptd), nyr, 1)', [], 1);
    date_yyyyptd = date_array_gen('year_list', year, 'output_date_format', 'yyyyptd');

    tic;
    dim_time     = length(size(daily_data));
    
    if dim_time == 2 
        [ni, ~]     = size(daily_data);   %/ always assume the last dim is time.
        pentad_clim = nan(ni, nptd);
        pentad_allyr = nan(ni, nptd * nyr);
        for ptd = 1:nptd
            ptd_dates = date_mmdd_OneYr((1:5)+5*(ptd-1));

            %/ pentad climatology
            ind                = findismember_loop(mod(dates,1e4), ptd_dates); %/ e.g., 210 days for a pentad for 42 years
            pentad_clim(:,ptd) = mean(daily_data(:,ind), dim_time);
       
            %/ pentad time series (full) - only 0.5s faster than for-loop.
            ptd_dates_ann    = reshape((year*1e4 + [date_mmdd_OneYr((1:5)+5*(ptd-1))]')', [], 1);
            ind              = findismember_loop(dates, ptd_dates_ann);
            daily_data_reshp = reshape(daily_data(:,ind), ni, 5, []);
            pentad_allyr(:,ptd + nptd * ((1:nyr)-1)) = squeeze(mean(daily_data_reshp, dim_time));
            
%             %/ pentad time series (full)
%             for y = 1:nyr
%                 ptd_dates_ann = year(y)*1e4 + date_mmdd_OneYr((1:5)+5*(ptd-1));
% 
%                 ind       = findismember_loop(dates, ptd_dates_ann);
%                 pentad_allyr(:, ptd + nptd * (y-1)) = mean(daily_data(:,ind), dim_time);
%             end
        end
        
        pentad_anom     = nan(ni, nptd * nyr);
        pentad_anom_std = nan(ni, nptd);
        for ptd = 1:nptd
            ind = ptd + nptd * ((1:nyr)- 1); %/ get the indices of a pentad p in all years.

            %/ pentad anomalies
            pentad_anom(:,ind) = pentad_allyr(:,ind) - pentad_clim(:,ptd);

            %/ (Interannual) S.D. of pentad anomalies
            pentad_anom_std(:,ptd) = std(pentad_anom(:,ind), 0, dim_time); 
        end
        
    elseif dim_time == 3
        
        [ni, nj, ~] = size(daily_data);
        
        pentad_clim = nan(ni, nj, nptd);
        pentad_allyr = nan(ni, nj, nptd * nyr);
        for ptd = 1:nptd
            ptd_dates = date_mmdd_OneYr((1:5)+5*(ptd-1));

            %/ pentad climatology
            ind                  = findismember_loop(mod(dates,1e4), ptd_dates); %/ e.g., 210 days for a pentad for 42 years
            pentad_clim(:,:,ptd) = mean(daily_data(:,:,ind), dim_time);

            %/ pentad time series (full) - only 0.5s faster than for-loop.
            ptd_dates_ann = reshape((year*1e4 + [date_mmdd_OneYr((1:5)+5*(ptd-1))]')', [], 1);
            ind           = findismember_loop(dates, ptd_dates_ann);
            daily_data_reshp = reshape(daily_data(:,:,ind), ni, nj, 5, []);
            pentad_allyr(:,:,ptd + nptd * ((1:nyr)-1)) = squeeze(mean(daily_data_reshp, dim_time));
        end
        
        pentad_anom     = nan(ni, nj, nptd * nyr);
        pentad_anom_std = nan(ni, nj, nptd);
        for ptd = 1:nptd
            ind = ptd + nptd * ((1:nyr)- 1); %/ get the indices of a pentad p in all years.

            %/ pentad anomalies
            pentad_anom(:,:,ind) = pentad_allyr(:,:,ind) - pentad_clim(:,:,ptd);

            %/ (Interannual) S.D. of pentad anomalies
            pentad_anom_std(:,:,ptd) = std(pentad_anom(:,:,ind), 0, dim_time); 
        end
        
    elseif dim_time == 4

        [ni, nj, nk, ~] = size(daily_data);
        
        pentad_clim  = nan(ni, nj, nk, nptd);
        pentad_allyr = nan(ni, nj, nk, nptd * nyr);
        for ptd = 1:nptd
            ptd_dates = date_mmdd_OneYr((1:5)+5*(ptd-1));
            
            %/ pentad climatology
            ind                     = findismember_loop(mod(dates,1e4), ptd_dates); %/ e.g., 210 days for a pentad for 42 years
            pentad_clim(:,:,:,ptd) = mean(daily_data(:,:,:,ind), dim_time);
            
            %/ pentad time series (full) - only 0.5s faster than for-loop.
            ptd_dates_ann = reshape((year*1e4 + [date_mmdd_OneYr((1:5)+5*(ptd-1))]')', [], 1);
            ind           = findismember_loop(dates, ptd_dates_ann);
            
            daily_data_reshp = reshape(daily_data(:,:,:,ind), ni, nj, nk, 5, []);
            pentad_allyr(:,:,:,ptd + nptd * ((1:nyr)-1)) = squeeze(mean(daily_data_reshp, dim_time));
        end
        
        pentad_anom     = nan(ni, nj, nk, nptd * nyr);
        pentad_anom_std = nan(ni, nj, nk, nptd);
        for ptd = 1:nptd
            ind = ptd + nptd * ((1:nyr)- 1); %/ get the indices of a pentad p in all years.

            %/ pentad anomalies
            pentad_anom(:,:,:,ind) = pentad_allyr(:,:,:,ind) - pentad_clim(:,:,:,ptd);

            %/ (Interannual) S.D. of pentad anomalies
            pentad_anom_std(:,:,:,ptd) = std(pentad_anom(:,:,:,ind), 0, dim_time); 
        end
    end
    toc;
    
    %/ double check
%     ptd_series = repmat([1:73], 1, NoOfYrs)';
%     test = nan(nlon, nlat, nptd);
%     for ptd = 1:nptd
%         
%         ind = find(ptd_series == ptd);
%         
%         test(:,:,ptd) = mean(pentad_data_full(:,:,ind),dim_time);
%     end
%  
%     a = test - pentad_data_clim;
%     max(abs(a), [], 'all')     %/ 5.6843e-13 --> so calculations are correct.
    
    
end
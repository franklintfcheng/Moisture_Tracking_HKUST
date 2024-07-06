%%
function X = conv2seasonal(varargin)

pnames = {'data',   'dates',   'select_year', 'detrend_mode', 'mean_or_sum', 'NumWorker'};
dflts  = {[]        []          []              0,                   'mean',          10};
[data, dates, select_year, detrend_mode, mean_or_sum, NumWorkers] = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

if ~isequal(mean_or_sum, 'mean') && ~isequal(mean_or_sum, 'sum')
    error('Wrong in put of ''mean_or_sum''!!');
end

isannual = 0;
if      numel(num2str(dates(1))) == 8   CONV = 100; 
elseif  numel(num2str(dates(1))) == 6   CONV = 1;   
elseif  numel(num2str(dates(1))) == 4   isannual = 1;
else    error('Check the format of the input date! It has to be either in yyyymmdd or yyyymm!') 
end
time_dim = length(size(data)); %/ assume the last dimension is the time dimension.

if     time_dim == 2   X = nan(5, length(select_year));                              
elseif time_dim == 3   X = nan(size(data,1), size(data,2), 5, length(select_year));  
end

flag = 0;
for t = 1:length(select_year) 
    for s = 1:5
        if isannual
            ind_date = find(dates == select_year(t)); %/ ANN
        else
            if     s == 1    ind_date = find(floor(dates/CONV) >= select_year(t)*100 + 3 & floor(dates/CONV) <= select_year(t)*100 + 5);
            elseif s == 2    ind_date = find(floor(dates/CONV) >= select_year(t)*100 + 6 & floor(dates/CONV) <= select_year(t)*100 + 8);
            elseif s == 3    ind_date = find(floor(dates/CONV) >= select_year(t)*100 + 9 & floor(dates/CONV) <= select_year(t)*100 + 11);
            elseif s == 4    
                             a = find(floor(dates/CONV) == select_year(t)*100 + 12);
                             b = find(floor(dates/CONV) >= (select_year(t)+1)*100 + 1 & floor(dates/CONV) <= (select_year(t)+1)*100 + 2);
                             if isempty(a) || isempty(b)
%                                  warning('Skip cos not enough data to compute DJF!')
                                 flag = 1;
                                 continue;                %/ skip when not enough data to compute DJF.
                             end   
                             ind_date = cat(1, a, b);
            elseif s == 5 
                ind_date = find(floor(dates/CONV) >= select_year(t)*100 + 1 & floor(dates/CONV) <= select_year(t)*100 + 12); %/ ANN
            end
        end
        
        if isempty(ind_date)   error('No date is retrieved!');  end

        if      time_dim == 2   
            if isequal(mean_or_sum, 'mean')
                X(s, t)       = mean(data(ind_date), 'omitnan');
            elseif isequal(mean_or_sum, 'sum')
                X(s, t)       = sum(data(ind_date), 'omitnan');
            end
            
        elseif  time_dim == 3   
            if isequal(mean_or_sum, 'mean')
                X(:, :, s, t) = mean(data(:,:,ind_date), time_dim, 'omitnan');
            elseif isequal(mean_or_sum, 'sum')
                X(:, :, s, t) = sum(data(:,:,ind_date), time_dim, 'omitnan');
            end
        else
            error('The function only consider 1D, 2D or 3D data!');
        end
    end
end

if detrend_mode 
    if time_dim == 2
        [nseason, ntime] = size(X);
    elseif time_dim == 3
        [nlon, nlat, nseason, ntime] = size(X);
    end
        
    ind = find(isnan(X));
    if ~isempty(ind)  error('Data contains %d NaNs!', length(ind)); end
    
    if time_dim == 2
        X_2D = X';                                                         %/ transpose to year x seasons
        for s = 1:nseason
            X_2D_bc = X_2D(:,s);
            
            %/ NOTE: DJF value in the last DJF could be unavailable (if so will then set to 0)
            %/       remove the element to make a 39-yr data before detrending.
            if flag == 1 && s == 4  
                X_2D_bc(end) = [];                                         %/ remove this last DJF value before detrending
            end
            
            X_2D_bc = detrend(X_2D_bc);                                    %/ If x is a vector, then detrend subtracts the trend from the elements of x.
            
            if flag == 1 && s == 4  
                X_2D_bc(end+1,:) = nan;                                    %/ add back last DJF value as nan to restore the data size
            end
            
            X_2D_detrend(:,s) = X_2D_bc;
        end
        X = X_2D_detrend';                                                 %/ transpose back
    
    elseif time_dim == 3
        for s = 1:nseason
            X_2D = reshape(X(:,:,s,:), [], ntime);                         %/ reshape   to grids x time
            X_2D = X_2D';                                                  %/ transpose to time x grids (for detrend to perform correctly)
            
            %/ NOTE: DJF value in the last DJF could be unavailable (if so will then set to 0)
            %/       remove the element to make a 39-yr data before detrending.
            if flag == 1 && s == 4  
                X_2D(end,:) = [];                                          %/ remove this last DJF value before detrending
            end
            
            X_2D_detrend = detrend(X_2D);                                  %/ If x is a matrix, then detrend operates on each column separately, subtracting each trend from the corresponding column.
            
            if flag == 1 && s == 4
                X_2D_detrend(end+1,:) = nan;                               %/ add back last DJF value as nan to restore the data size
            end
            
            X_2D_detrend = X_2D_detrend';                                  %/ transpose back to grids x time
            X(:,:,s,:) = reshape(X_2D_detrend, nlon, nlat, 1, ntime);      %/ reshape back to lon x lat x 1 x time
        end
    end
end

% if detrend_mode    %/ to handle omitnan, since detrend() in matlab2019 does not have nanflag argument!
%     
%     if isempty(gcp('nocreate')) && ~isempty(NumWorkers)
%         parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
%     end
% 
%     tic
%     N = length(X_2D_detrend_cell);
%     parfor i = 1:N
%     %     fprintf('i = %d/%d \n', i, N)
% 
%         X_2D_nonan_OneGrid = X_2D(:, i);                         %/ load the original data
%         ind_nan            = find(isnan(X_2D_nonan_OneGrid));    %/ find indices of NaNs
%         X_2D_nonan_OneGrid(ind_nan) = [];                        %/ remove NaNs
% 
%         if isempty(X_2D_nonan_OneGrid) || isempty(ind_nan) || length(X_2D_nonan_OneGrid) == 1     
%             continue;                                            %/ skip if all are NaNs, all are non-NaNs or only one value left.
%         end  
%         X_2D_nonan_OneGrid = detrend(X_2D_nonan_OneGrid);        %/ perform detrending
% 
%         ind_nonnan = setdiff(1:length(select_year), ind_nan);
%         a = nan(length(select_year), 1);
%         a(ind_nonnan) = X_2D_nonan_OneGrid;
% 
%         X_2D_detrend_cell{i} = a;
%     end
%     X_2D_detrend = [X_2D_detrend_cell{:}];                       %/ convert back to numerical matrix
%     toc
% 
%     X_2D_detrend = X_2D_detrend';                                                              %/ transpose back to grids x time
%     if time_dim == 2
%         X_detrend = X_2D_detrend;
% 
%     else if time_dim == 3
%         X_detrend = reshape(X_2D_detrend, size(X,1), size(X,2), size(X, 3), size(X, 4));       %/ reshape back to lon x lat x season x year
%     end
%     end
%     X = X_detrend;
% end

% tic
% [row, col] = find(isnan(X_2D));                                                            %/ year x grids
% uni_col    = unique(col)';                                                                  %/ since one column can have >1 nans, we avoid processing the same column repeatedly.
% find(diff(uni_col) < 0)
% for i = 1:length(uni_col)
%     fprintf('i = %d/%d \n', i, length(uni_col))
%     X_2D_nonan_OneGrid                            = X_2D(:, uni_col(i));
%     X_2D_nonan_OneGrid(isnan(X_2D_nonan_OneGrid)) = [];                                    %/ remove nan
% 
%     if isempty(X_2D_nonan_OneGrid) || length(X_2D_nonan_OneGrid) == 1     continue;   end  %/ skip if all are nans or only one value left.
% %         disp(X_2D_nonan_OneGrid)
%     X_2D_nonan_OneGrid = detrend(X_2D_nonan_OneGrid);                                      %/ perform detrending
% 
%     ind = find(col == uni_col(i));
%     rows_with_nan   = row(ind);
%     rows_with_nonan = setdiff(1:length(select_year), rows_with_nan);
%     X_2D_detrend(rows_with_nonan, uni_col(i)) = X_2D_nonan_OneGrid;
% 
% end
% toc

%     trend_array = (X - X_detrended);                                     %/ slope
%     trend = (trend_array(end) - trend_array(1))/length(year_list)*10*12; %/ trend (unit per decade)

end
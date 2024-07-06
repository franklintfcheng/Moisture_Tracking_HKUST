%%
function X_deseason = deseasonalize(X, X_dates)

    if      numel(num2str(X_dates(1))) == 8   mod_divider = 1e4; 
    elseif  numel(num2str(X_dates(1))) == 6   mod_divider = 1e2;   
    else    error('Check the format of the input date! It has to be either in yyyymmdd or yyyymm!') 
    end

    if any(isnan(X))     error('Input X contains NaN. Check if it is correct!'); end
    
    time_dim = length(size(X));
    if time_dim == 2 && size(X,1) ~= 1     X = X';     end                 %/ make it a row vector if not so
    
    
    %/ A simple deseasonality method by removing the daily/monthly climatology.
    b = mod(X_dates, mod_divider);
    b_unique = unique(b);

    tic;
    X_deseason = nan(size(X));
    for i = 1:length(b_unique)
    
        ind  = find(b == b_unique(i));
        if time_dim == 2
            clim            = mean(X(ind), time_dim, 'omitnan');
            X_deseason(ind) = X(ind) - clim;
            
        elseif time_dim == 3
            clim                = mean(X(:,:,ind), time_dim, 'omitnan');
            X_deseason(:,:,ind) = X(:,:,ind) - clim;
            
        end
    end
    fprintf('*** Time cost in deseasonalizing data: %.2f sec ***\n', toc)
end
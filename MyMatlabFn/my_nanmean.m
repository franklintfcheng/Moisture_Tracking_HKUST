%%

function [X_nanmean] = my_nanmean(X, dim)
    
    %/ Inf (due to 1/0) will cause incorrect results of nanmean()! Set Inf to NaN;
    if ~isempty(isinf(X))
        X(isinf(X)) = NaN;   
%         warning('Inf is found, replaced with nan.');
    end
    
    X_nanmean = mean(X, dim, 'omitnan');
    
end




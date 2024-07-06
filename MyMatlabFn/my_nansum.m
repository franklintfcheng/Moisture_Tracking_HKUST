%%

function [X_nansum] = my_nansum(X, dim)
    
    %/ Inf (due to 1/0) will cause incorrect results of nanmean()! Set Inf to NaN;
    if ~isempty(isinf(X))
        X(isinf(X)) = NaN;   
%         warning('Inf is found, replaced with nan.');
    end
    
    X_nansum = nansum(X, dim);
    
end
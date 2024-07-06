
function dX = centdiff(X, h) % X has to be 1D array

    n = length(X);
%   i.e., elementwise computation on dX(k) = (X(k+1) - X(k-1))/(2*h);
    dX = (X(3:end) - X(1:end-2))/(2*h);
end
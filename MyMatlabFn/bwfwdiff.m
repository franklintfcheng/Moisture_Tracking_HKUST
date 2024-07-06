
function dX = bwfwdiff(X, h) % X has to be 1D array

%   i.e., elementwise computation on dX(k) = (X(k+1) - X(k))/h;
    dX = (X(2:end) - X(1:end-1))/(h);
end
%%
function D = EulerDis(x1, x2)
 

    sq_dev    = squeeze((x1 - x2).^2);                      %/ square of deviation
    sq_dev_1D = reshape(sq_dev, [], 1);                     %/ reshape to 1D before summing.
    D         = sqrt(sum(sq_dev_1D, 'omitnan'));            %/ summing up and take sqrt. (by omitting nans)
    

end
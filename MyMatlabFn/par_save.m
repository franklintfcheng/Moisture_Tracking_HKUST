
function par_save(fname, x) 
    warning('If the variable to save is not named by x, stop using this function!!!');
    save(fname, 'x', '-v7.3')
end
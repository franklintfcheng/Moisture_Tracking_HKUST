function theta = compute_PT(varargin)

    pnames = {'T', 'p_dim',  'p',  'p_unit', 'debug'};

    dflts =  cell(1, length(pnames));
    [          T,   p_dim,     p,    p_unit,  debug ] ...
                            = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 10 Apr 2024
    %/
    %/ Description: This function outputs potential temperature based on
    %/              the given 3D/4D temperature data
    %/
    %/      'p_dim': the dimension index of pressure level in T 
    %/               (set it [] if there is no such dimension, then p must be a const only)
    %/=====================================================================

    p = double(p);  %/ make sure it's double, not integer

    nlevel = length(p);
    if length(size(T)) >= 4
        if isempty(p_dim)
            p_dim = 3;
            warning('Detected T is a 4D or higher dimensional data. Since p_dim is not given, assuming ''p_dim == 3''.');
        end
    end

    %/ If p_dim is given, that indicates T has a dimension of plevel,
    %/ check consistency then
    if ~isempty(p_dim)
        T_nlevel = size(T, p_dim);
        if ~isequal(T_nlevel, nlevel)
            error('The p-level dimension of ''T'' is not consistent with the length of ''p''!');
        end
        p = reshape(p, [ones(1,p_dim-1),nlevel]); %/ 1 x nlevel -> 1 x 1 x nlevel (if p_dim == 3)
    else
        if nlevel ~= 1
            error('p should be a constant!')
        end
    end
    
    if debug
        disp('Size of p (after reshape):')
        size(p)
        disp('Size of T:')
        size(T)
    end

    if isequal(p_unit, 'hPa')
        p_CONV = 1;
    elseif isequal(p_unit, 'Pa')
        p_CONV = 100;
    else
        error('Specific ''p_unit''! Either ''Pa'' or ''hPa''.');
    end
    cp = 1004;         %/ Specific heat at constant pressure [J kg-1 K-1]
    R  = 287;          %/ Gas constant for dry air [J kg-1 K-1]
    ps = 1000*p_CONV;  %/ Pressure at reference level (generally taken as 1000 hPa) 
    
    theta = T.*(ps./p).^(R/cp); %/ Use elementwise multiplication since p can be an array
    
    if ~isempty(find(isnan(theta), 1))
        warning('compute_PT: theta contains %d NaNs!', length(find(isnan(theta))));
    end

    %/ Checking
    fprintf('compute_PT: min(theta)  = %.2f K \n', min(theta, [], 'all', 'omitnan'));
    fprintf('compute_PT: mean(theta) = %.2f K \n', mean(theta,    'all', 'omitnan'));
    fprintf('compute_PT: max(theta)  = %.2f K \n', max(theta, [], 'all', 'omitnan'));
 
end
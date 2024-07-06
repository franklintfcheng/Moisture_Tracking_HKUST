%%
function MSE = MoistStaticEnergy(varargin)
    % ----------------------------
    %Author: Fandy Cheng
    %Date of creation: 9 May 2022

    %/ Compute the Moist Static Energy (MSE or m).

    %----- About the input------%
    %/ Temperature (K)
    %/ Geopotential Height Z (m) 
    %/ Specific humidity q (kg/kg)
    
    pnames = {'T', 'Z', 'q'};
    dflts  = cell(1, length(pnames));
    [T, Z, q] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    fprintf('*** Running MoistStaticEnergy... ***\n')
    
    if any(T < 100)                      error('Check the input T, it should be in K!');    end
    if mean(Z(:,:,end,1), 'all') > 1e5   warning('Z appears too large. Check if you have divided it by 9.81 from geopotential.'); end
    
    ind = find(q < 0);
    if ~isempty(ind)
        warning('Specific humidity has %d negative values! Auto correct them to zeros...\n', length(ind))
        q(q < 0) = 0;
    end
    
    cp = 1004;  %/ J kg-1 K-1
    g  = 9.81;  %/ m s-2
    Lc = 2.5e6; %/ J kg-1 = Latent heat of condensation at 0 degC (See Maloney 2009; Holton 4th Ed. P. 291, 492).
    
    MSE = cp*T + g*Z + Lc*q;  %/ J kg-1
    
end







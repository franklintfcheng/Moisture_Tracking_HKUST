%%
function P = z2p_hydrostatic(varargin)
    
    pnames       = {   'z',    'unit'};
    dflts        = {    [],       'm'};
    [                   z,      unit] = ...
            internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
    
%     fprintf('*** Running p2z_hydrostatic... ***\n')
    if isequal(unit, 'm')
        conv = 1e-3;
    elseif isequal(unit, 'km')
        conv = 1;
    end

    H  = 8.5;                    %/ in km. 
                                 %/ NOTE: Typical scale height of Earth's atmosphere (mean T = 290K) = 8.5 km (by NASA);
                                 %/       The value was said to be 7.5 km by Wallace and Hobbs (2006, Ch.3),
                                 %/       but it gives 813 hPa when z = 1550 m;
                                 %/       Setting 8.5 km of scale height gives 833 hPa for z = 1550 m. (more reasonable).
    P0 = 1000;                   %/ in hPa.  (Suggest not changing this)
    P  = P0*exp(-z*conv/H);      %/ in hPa
end
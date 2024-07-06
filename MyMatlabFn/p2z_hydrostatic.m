function z = p2z_hydrostatic(varargin)
    
    pnames       = {   'P',    'unit'};
    dflts        = {    [],     'hPa'};
    [                   P,      unit] = ...
            internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
%%
    %/ NOTE: Typical scale height of Earth's atmosphere (mean T = 290K) = 8.5 km (by NASA);
    %/       The value was said to be 7.5 km by Wallace and Hobbs (2006, Ch.3),
    %/       but it gives 813 hPa when z = 1550 m;
    %/       Setting 8.5 km of scale height gives 833 hPa for z = 1550 m. (more reasonable).

    if isequal(unit, 'hPa')
        conv = 1;
    elseif isequal(unit, 'Pa')
        conv = 100;
    end
    H  = 8.5;                    %/ in km                      
    P0 = 1000*conv;              %/ in hPa
    z  = -H * log (P/P0)*1000;   %/ in m
    
end
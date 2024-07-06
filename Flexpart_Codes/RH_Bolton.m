function [RH, p] = RH_Bolton(z, q, T) 
    
    %/ Author: Fandy Cheng
    %/ Date of last update: 4 Jan 2023

    %/==== References ====
    %/ Geoffrey K. Vallis (2019). Essentials of Atmospheric and Oceanic
    %/      Dynamics. Chapter 11.
    %/ Bolton, D. 1980. The computation of equivalent potential temperature. Mon. Wea. Rev.. 108. 1046?1053.
    %/ https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
    %/ https://glossary.ametsoc.org/wiki/Clausius-clapeyron_equation
    
    %/==== Input variables ====
    %/ z = altitude (m)
    %/ q = specific humdity (kg/kg)
    %/ T = temperature (K)
    
    R = 287;          % specific gas constant for dry air (J K**-1 kg**-1)
    g = 9.81;         % (m s**-2)
    H = R*T/g;        % (m)    %/ CAVEAT: using T from parcel to compute H may be physically incorrect!
    p0 = 1013.25;     % (hPa)
    
    %/ Hydrostatic balance (Assumption only!)
    p = p0*exp(-z./H); % (hPa)
    
    md = 28.97;            % molecular weight of dry air     (g/mol)
    mv = 18.0153;          % molecular weight of water vapor (g/mol)
    eps = round(mv/md, 3); % = 0.622
    
    %/ Water vapor pressure (e)
    e = q.*p./(eps + (1-eps)*q); % See Vallis (2019)
    
    %/ Saturated water vapor pressure (es)
    T0 = 273.16;       % (K)
    e0 = 6.112;        % (hPa)
    T = T - T0;        % convert to (deg C)
    
    %/ CAVEAT: This estimation of es is most emprically accurate for -30C <= T <= 35C 
    %/         See Eq.10 in Bolton (1980).
    es = e0*exp(17.67*T./(T+243.5));        % (hPa). T is in (deg C). 
    
    %/ Compute RH
    RH = e./es*100;
    
    
    
%------------ old codes -------------%
%     Rv = 462;         % specific gas constant for water vapor (J K**-1 kg**-1)
%     L  = 2.44e6;      % latent heat of vapor condensation (J kg**-1)
    
    %/ Saturated water vapor pressure
%     es = e0*exp(L/Rv*(1/T0 - 1./T));  % this eq. is not too accurate since L = L(T) (can vary by ~10%)

end
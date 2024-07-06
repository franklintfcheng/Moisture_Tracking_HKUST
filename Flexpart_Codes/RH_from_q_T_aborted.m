function [RH, p] = RH_from_q_T(z, q, T)    

    R = 287;          % J K**-1 kg**-1
    g = 9.81;         % m s**-2
    H = R*T/g;        % m
    p0 = 1013.25*100; % Pa
    
    %/ hydrostatic balance
    p = p0*exp(-z./H);

    T0 = 273.16;
    RH = 0.263.*p.*q.*exp((17.67.*(T-T0))./(T-29.65)).^(-1); % RH in %

    %/ https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
    % Variables used:
    % q specific humidity or the mass mixing ratio of water vapor to total air (kg/kg)
    % mv specific mass of water vapor (kg)
    % mvs specific mass of water vapor at equilibrium (kg)
    % md specific mass of dry air (kg)
    % w mass mixing ratio of water vapor to dry air (dimensionless)
    % ws mass mixing ratio of water vapor to dry air at equilibrium (dimensionless)
    % es(T) saturation vapor pressure (Pa)
    % es0 saturation vapor pressure at T0 (Pa)
    % Rd specific gas constant for dry air (J kg?1 K?1)
    % Rv specific gas constant for water vapor (J kg?1 K?1)
    % p pressure (Pa)
    % T temperature (K)
    % T0 reference temperature (typically 273.16 K) (K)
end
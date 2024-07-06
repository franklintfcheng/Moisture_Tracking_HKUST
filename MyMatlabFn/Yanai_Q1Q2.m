%%
function [Q1, Q2] = Yanai_Q1Q2(varargin)
    % ----------------------------
    %Author: Fandy Cheng
    %Date of creation: 17 May 2022

    %/ Compute Q1 and Q2 (Yanai et al. 1973; Johnson and Ciesielski (2000) 

    %----- About the input------%
    %/ T    = Temperature       (K)
    %/ q    = Specific humidity (kg/kg)
    %/ U, V = horizontal winds  (m/s)
    %/ W    = pressure velocity (Pa/s)
    %/ P    = pressure levels   (hPa)
    
    pnames = {'T', 'q', 'U', 'V', 'W', 'lon', 'lat', 'P'};
    dflts  = cell(1, length(pnames));
    
    [          T,   q,   U,   V,   W,   lon,   lat,   P] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    
    fprintf('*** Running Yanai_Q1Q2... ***\n')
    
    %/ reshape input pressure (P) level to 1 x 1 x nplev
    reshp_P_size      = ones(1, length(size(q))-1);
    reshp_P_size(end) = length(P);
    P                 = reshape(P, reshp_P_size);     %/ 1 x 1 x nplev
    
    if isempty(P)                        error('This function only handles pressure level data. Input the corresp. pressure plz.');  end
    if any(P > 1e4)                      error('Check the input P, it should be in hPa!');  end
    if any(T < 100)                      error('Check the input T, it should be in K!');    end
    
    ind = find(q < 0);
    if ~isempty(ind)
        warning('Specific humidity has %d negative values! Auto correct them to zeros...\n', length(ind))
        q(q < 0) = 0;
    end
    
    cp = 1004;   %/ [J kg-1 K-1]
    g  = 9.81;   %/ [m s-2]
    Lc = 2.5e6;  %/ [J kg-1],       Latent heat of condensation at 0 degC (See Maloney 2009; Holton 4th Ed. P. 291, 492).
    R  = 287;    %/ [J kg-1 K-1],   Specific gas constant for dry air.
    r  = 6371e3; %/ [m],            Earth's radius.
    P0 = 1000;   %/ hPa
    
    dt  = 24*3600;                          %/ NOTE: Assume daily time step (dt = 24*3600 s)
    hx  = diff(lon(1:2))*pi/180;            %/ rmb to change degree to radian!!!!
    hy  = diff(lat(1:2))*pi/180;            %/ rmb to change degree to radian!!!! NOTE: hy can be -ve. <- correct.
    
    %/ Derivatives of T
    [~, ~, ~, dT_dt] = gradient(T, dt);                      %/ K s-1
    
    [~, dT_dx, ~, ~] = gradient(T, hx);    
    unit_conv        = (1./(r*cosd(lat)))';                  %/ 1 x lat
    dT_dx            = unit_conv.*dT_dx;                     %/ K m-1.   div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))
    
    [dT_dy, ~, ~, ~] = gradient(T, hy);         
    unit_conv        = 1./r;                     
    dT_dy            = unit_conv.*dT_dy;                     %/ K m-1.   div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))
    
    theta                = T.*(P0./P).^(R/cp);               %/ K, potential temp.
    [~, ~, dtheta_dA, ~] = gradient(theta);                  %/ Chain rule to do derivate with *non-uniform* spacings.
    dP_dA                = gradient(squeeze(P) * 100);       %/ get non-uniform gradient of pressure level first.
    dP_dA                = reshape(dP_dA, reshp_P_size);     %/ 1 x 1 x plev
    dtheta_dP            = dtheta_dA./dP_dA;                 %/ elementwise division.  J kg-1 Pa-1
    
    %/ Derivatives of q
    [~, ~, ~, dq_dt] = gradient(q, dt);                      %/ kg/kg s-1
    
    [~, dq_dx, ~, ~] = gradient(q, hx);    
    unit_conv        = (1./(r*cosd(lat)))';                  %/ 1 x lat
    dq_dx            = unit_conv.*dq_dx;                     %/ kg/kg  m-1.   div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))

    [dq_dy, ~, ~, ~] = gradient(q, hy);         
    unit_conv        = 1./r;                     
    dq_dy            = unit_conv.*dq_dy;                     %/ kg/kg  m-1.   div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))
    
    [~, ~, dq_dA, ~] = gradient(q);                          %/ Chain rule to do derivate with *non-uniform* spacings.
    dq_dP            = dq_dA./dP_dA;                         %/ elementwise division.  J kg-1 Pa-1
    
    %/ Compute Q1
    Q1 =  cp * ( dT_dt + U.*dT_dx + V.*dT_dy + (P./P0).^(R/cp).*W.*dtheta_dP );  %/ J kg-1 s-1 
    
    %/ Compute Q2
    Q2 = -Lc * ( dq_dt + U.*dq_dx + V.*dq_dy + W.*dq_dP );                       %/ J kg-1 s-1 
    

    
    %==========================================================================================
    %/ Following Luo and Yanai (1984), you can plot the vertical profile of Q1/cp [K day-1].
    %==========================================================================================
    
    
%     clc
%     mean(dq_dP, 'all')
%     squeeze(dtheta_dA(100, 40, :, 1))
%     squeeze(Q1(100, 40, :, 1))
%     squeeze(Q2(100, 40, :, 1))
%     squeeze(q(100, 40, :, 1)) 
    
    
    
end







function output = calc_IVTdiv(varargin)
    
    % create a set of valid parameters and their default value
    pnames = {'uIVT', 'vIVT', 'lon', 'lat',    'output_unit'};
    dflts  = {    [],     [],    [],    [],         'mm/day'};

    %/ parse function arguments
    [         uIVT,    vIVT,   lon,   lat,       output_unit] ...
           = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %==========================================================================
    %/      Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 17 May 2024
    %/
    %/ Description: This function is designd to calculate vertically
    %/              integrated moisture divergence (IVTdiv)
    %==========================================================================
    
    if isequal(uIVT, vIVT)
        error('Did you just input uIVT as vIVT or vice versa?');
    end
    
    %/ compute divergence of IVT
    %/ [FX,FY,FZ,...,FN] = gradient(F) returns the N components of the numerical gradient of F, where F is an array with N dimensions.
    %/ The first output FX is always the gradient along the *2nd* dimension of F, going across columns.
    %/ The second output FY is always the gradient along the *1st* dimension of F, going across rows.
    %/ For the third output FZ and the outputs that follow, the Nth output is the gradient along the Nth dimension of F.
    
    if isequal(output_unit, 'mm/day')
        rho_w = 997;                      %/ kg m**-3
        unitconv = 1/rho_w*1000*24*3600;  %/ since rho_w is not perfectly 1000 kg m**-3.
    else
        error('code not set yet!');
    end

    %/ NOTE: the calculation of IVT has already incorporated rho_w and g. -> THat's why [IVT] = kg/m/s. 
       
    fprintf('*** Running calc_IVTdiv... ***\n')

    h_lambda = diff(lon(1:2))*pi/180;   %/ rmb to change degree to radian!!!!
    [~, du_dlambda, ~] = gradient(uIVT, h_lambda);

    h_theta = diff(lat(1:2))*pi/180;   %/ rmb to change degree to radian!!!!
    [dv_dtheta, ~, ~] = gradient(vIVT, h_theta);

    [~, lat_2D] = meshgrid(lon, lat);
    lat_2D = lat_2D';

    r     = 6371e3;   %/ in m    
    output = 1./(r*cosd(lat_2D)).*(du_dlambda + dv_dtheta.*cosd(lat_2D))*unitconv;  %/ == du/dx + dv/dy in cartesian coor.
    
end
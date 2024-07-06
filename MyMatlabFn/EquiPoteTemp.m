%%
function EPT = EquiPoteTemp(varargin)
    % ----------------------------
    %Author: Fandy Cheng
    %Date of creation: 27 Feb 2020

    %/ Approximate equivalent potential temperature 

    %----- About the input------%
    %/ Pressure P (hPa) 
    %/ Specific humidity q (kg/kg)
    %/ Temperature (K)

    %/ mixing ratio will be converted to (g/kg).

    pnames = {'P', 'q', 'T'};
    dflts  = {[], [] []};
    [P, q, T] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    fprintf('*** Running EquiPoteTemp... ***\n')
    
    %/ reshape input pressure (P) level to 1 x 1 x nplev
    reshp_size      = ones(1, length(size(q))-1);
    reshp_size(end) = length(P);
    P               = reshape(P, reshp_size);       
    
    if isempty(P)      error('This function only handles pressure level data. Input the corresp. pressure plz.');  end
    if any(P > 1e4)    error('Check the input P, it should be in hPa!');  end
    if any(T < 100)    error('Check the input T, it should be in K!');    end
    
    ind = find(q < 0);
    if ~isempty(ind)
        warning('Specific humidity has %d negative values! Auto correct them to zeros...\n', length(ind))
        q(q < 0) = 0;
    end
    
    %/ See EAOD Vallis (2019)
    r = (q./(1 - q))*1000;                          %/ mixing ratio (convert into g/kg)
    e = (P.*r)./(622 + r);                          %/ vapor pressure (hPa)
    
    %/ See Bolton (1980)
    TL = 2840./(3.5*log(T) - log(e) - 4.805) + 55;  %/ Eq. (21). The Parcel's temperature at lifting condensation level
    EPT = T.*(1000./P).^(0.2854*(1 - r*0.28e-3))...
                 .* exp((3.376./TL - 0.00254).*r.*(1 + r*0.81e-3)); %/ Eq. (43). equivalent potential temperature (K)

end







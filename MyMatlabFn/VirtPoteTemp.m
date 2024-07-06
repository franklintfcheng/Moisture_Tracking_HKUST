%%
function VPT = VirtPoteTemp(varargin)

    % ----------------------------
    %Author: Fandy Cheng
    %Date of creation: 9 Nov 2021

    %/ Approximate virtual potential temperature 
    %/ REF BOOK: Mesoscale Meteorology in Mid-latitudes 
    %               by Markowski and Richardson 2010.

    %----- About the input------%
    %/ Pressure P (hPa) 
    %/ Specific humidity q (kg/kg)
    %/ Temperature (K)
    
    pnames = {'T', 'q', 'P', 'sp', 'a', 'b'};
    dflts  = cell(length(pnames), 1);
    [T, q, P, sp, a, b] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %/ check input and compute pressure
    if isempty(P) && isempty(sp) && isempty(a) && isempty(b)
        error('No input for pressure!')
    elseif isempty(P)
        P = (a + b*sp)/100;   %/ Compute the air pressure of the model level given sfc pressure. Convert Pa into hPa.
    end
    
    if any(P > 1e4)          error('Check the input P, it should be in hPa!'); end
    if any(T < 100)          error('Check the input T, it should be in K!');   end

    ind = find(q < 0);
    if ~isempty(ind)
        warning('Specific humidity has %d negative values! Auto correct them to zeros...\n',length(ind))
        q(q < 0) = 0;
    end

    P0 = 1000;  %/ reference pressure (hPa)
    Rd = 287;   %/ gas constant for dry air (J kg-1 K-1)
    cp = 1004;   %/ specific heat capacity of air at constant pressure (J kg-1 K-1)
    rv = q./(1-q);
    
    PT  = T.*(P0./P).^(Rd/cp);  %/ potential temp
    VPT = PT.*(1 + 0.61*rv); %/ virtual potential temp 

end
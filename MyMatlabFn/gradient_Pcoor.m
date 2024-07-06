%%
function dAdp = gradient_Pcoor(varargin)

pnames = {'data', 'level', 'level_unit', 'lev_dim', 'lnP_coor'};
dflts  = cell(length(pnames),1);
[           data,  level,   level_unit,   lev_dim,   lnP_coor] = internal.stats.parseArgs(pnames, dflts, varargin{:});

%/ NOTE: this function is designed for a non-uniform spacing of pressure
%         level when computing vertical gradient in pressure coordinate.

if isempty(level_unit)   
    error('You must specify the unit of the input pressure level!');
else
    if isequal(level_unit, 'hPa')
        unitconv = 100;               %/ convert from hPa to Pa
    elseif isequal(level_unit, 'Pa')
        unitconv = 1;
    else
        error('P-levels are only in hPa or Pa. Modify the function code if needed.')
    end
end

if ~isvector(level)     error('level should be a vector only!'); end
if size(level,1) == 1   level = level';                          end %/ make it a column vector

[nlon,nlat,nlev,ntime] = size(data);
level = level*unitconv;
% dP    = diff(level);   %/ it can be -ve, if level is not in ascending order
lnP   = log(level);
% dlnP  = diff(lnP);
% dAdp  = nan(size(data));

%/ Q: How to do centered difference for uneven spacing of p-levels???
if lev_dim == 3
    if lnP_coor 
        %/ A smarter way to do centered difference on a non-uniform gridding
        level_nD         = reshape(level, 1, 1, nlev);            %/ [nlev, 1] -> [1 1 nlev]
        lnP_nD           = reshape(lnP, 1, 1, nlev);              %/ [nlev, 1] -> [1 1 nlev]
        lnP_nD           = repmat(lnP_nD, nlon, nlat, 1, ntime);  %/ [nlev, 1] -> [nlon nlat nlev ntime]
        [~, ~, dlnP]     = gradient(lnP_nD);
        [~, ~, dA]       = gradient(data);
        dAdp             = 1./level_nD.*dA./dlnP;
        
%         for p = 1:length(level)   %/ forward difference, and lastly backward difference
%             if p == length(level)
%                 dAdp(:,:,p) = 1/level(p)*(data(:,:,p) - data(:,:,p-1))/(lnP(p)-lnP(p-1));
%             else
%                 dAdp(:,:,p) = 1/level(p)*(data(:,:,p+1) - data(:,:,p))/dlnP(p);
%             end
%         end
    else
        %/ A smarter way to do centered difference on a non-uniform gridding
        level_nD         = reshape(level, 1, 1, nlev);              %/ [nlev, 1] -> [1 1 nlev]
        level_nD         = repmat(level_nD, nlon, nlat, 1, ntime);  %/ [nlev, 1] -> [nlon nlat nlev ntime]
        [~, ~, dP]  = gradient(level_nD);
        [~, ~, dA]       = gradient(data);
        dAdp             = dA./dP;
        
%         dAdp_old  = nan(size(data));
%         for p = 1:length(level)   %/ forward difference, and lastly backward difference
%             if p == length(level)
%                 dAdp_old(:,:,p) = (data(:,:,p) - data(:,:,p-1))/(level(p)-level(p-1));
%             else
%                 dAdp_old(:,:,p) = (data(:,:,p+1) - data(:,:,p))/dP(p);
%             end
%         end
%         mean(abs((dAdp - dAdp_old))./dAdp_old*100, 'all')
    end
else
    error('The function only works for lev_dim = 3. Please add the correspond code for the input lev_dim %d\n', lev_dim);
end

























end
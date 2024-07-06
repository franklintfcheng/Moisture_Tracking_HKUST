%%
function X_detrend = detrend_3D(varargin)

pnames = {'X', 'NumWorkers'};
dflts  = cell(length(pnames), 1);
[X, NumWorkers]= internal.stats.parseArgs(pnames, dflts, varargin{:});
    
if length(size(X)) ~= 3   error('This function only works for 3D matrix, with time in the 3rd dim.');  end

if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
    parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
end

%/ assume X is a 3D matrix with dimensions (lon, lat, time)
[nlon, nlat, ~] = size(X);

X_detrend = nan(size(X));
parfor ii = 1:nlon
    for jj = 1:nlat
        
        X_detrend(ii,jj,:) = detrend(squeeze(X(ii,jj,:)));
    end
end


end
function output = gavg(varargin)

pnames = {'data_3D', 'lat'}; 
dflts  = repmat([], length(pnames), 1);
[data_3D, lat] = internal.stats.parseArgs(pnames, dflts, varargin{:});

wgt = cosd(lat);                                %/ Create area weights first because the grid area decreases with latitude
wgt = wgt/sum(wgt);                             %/ Normalize wgt so that its sum is 1 
wgt = reshape(wgt, 1, [], 1);                   %/ Reshape to perform 3D element-wise mult.

 %/ we sum along lat dim since it has been weighted by wgt which has sum = 1.
output = squeeze(sum(mean(wgt.*data_3D, 1, 'omitnan'), 2, 'omitnan'));
        
end
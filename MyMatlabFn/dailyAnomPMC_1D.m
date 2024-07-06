%%

function dailyAnom_PMC = dailyAnomPMC_1D(varargin)

% create a set of valid parameters and their default values
pnames = {'dailydata', 'date_mmdd_AllYr', 'date_mmdd_OneYr'};

dflts  = {[], [], []};

%# parse function arguments
[dailydata, date_mmdd_AllYr, date_mmdd_OneYr] = internal.stats.parseArgs(pnames, dflts, varargin{:});

%====================================
[n1, n2] = size(dailydata);
if n1 == 1
    dailydata = dailydata'; %/ just to ensure time is in the 1st dim.
end

nday = length(date_mmdd_OneYr);
for i = 1:nday
    cday_ind = find(date_mmdd_AllYr == date_mmdd_OneYr(i));
    daily_clim(i, 1) = mean(dailydata(cday_ind, 1));

end
%disp(daily_clim(1:10))

% Pentad Moving Mean Clim (PMC)
if date_mmdd_OneYr(1) == 101 && date_mmdd_OneYr(end) == 1231
    %then we can append 2 days to each of the two ends for 5-day moving mean
    a = cat(1, daily_clim(:,:,end-1:end),...
                      daily_clim(:,:,:),...
                      daily_clim(:,:,1:2));
    PMC = movmean(a,5); % movmean(X,movingdays,dim)
    PMC = PMC(:,:,3:end-2);

else
    PMC = movmean(daily_clim,5); % movmean(X,movingdays,dim)
end
% disp(PMC(2))


% daily anom on PMC
dailyAnom_PMC = nan(size(dailydata));
for i = 1:nday
    cday_ind = find(date_mmdd_AllYr == date_mmdd_OneYr(i));
    dailyAnom_PMC(cday_ind, 1) = dailydata(cday_ind, 1) - PMC(i, 1);
end



end
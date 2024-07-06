function events = CI2events(varargin)

%/ create a set of valid parameters and their default value
pnames = {'CI_dates', 'CI_data', 'thres', 'n'};  
dflts  = cell(1, length(pnames));

[CI_dates, CI_data, thres, n] ...
               = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
    %%
    %/ NOTE:
    %/  thres = threshold value to define +ve and -ve events
    %/  n = min. consecutive months
    
    fprintf('*** Constructing climate events with threshold = %.4f, minimum consecutive months = %d ... ***\n', thres, n);
    
    M_hankel = hankel(CI_data(1:n), CI_data(n:end));
    events = zeros(size(CI_data));                  %/ 0 = neutral, 1 = El, -1 = La
    events(:,1) = CI_dates;
    
    %/ +ve event
    pve_ind = find(all(M_hankel >= thres, 1) == 1);  
    for i = 1:length(pve_ind)
        events(pve_ind(i):pve_ind(i) + n - 1, 2) = 1; %/ NOTE: whatever a column of 5 months met the criteria, we mark *all the 5 months* instead of one.
    end
    
    %/ -ve event
    nve_ind = find(all(M_hankel <= -thres, 1) == 1);  
    for i = 1:length(nve_ind)
        events(nve_ind(i):nve_ind(i) + n - 1, 2) = -1;
    end

end

%%
function [X_diff, X_diff_sig, str_drywet] = compute_drywet_anom(varargin)
    
    % create a set of valid parameters and their default value
    pnames = {'drywet_mode', 'alpha', 'X_daily', 'X_dates', 'dry_days', 'wet_days', 'normal_days'};
    dflts  = cell(1, length(pnames));
    [          drywet_mode,    alpha,  X_daily,   X_dates,   dry_days,   wet_days,   normal_days] ...
                        = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    %/ NOTE:
    %/ X_dates, dry_days, wet_days, normal_days in yyyymmdd.
    
    
    
    X_diff = []; X_diff_sig = []; str_drywet = [];
    if     drywet_mode == 1           str_drywet = '_Dry-Norm';  
    elseif drywet_mode == 2           str_drywet = '_Wet-Norm';  
    elseif drywet_mode == 3           str_drywet = '_Dry-Wet';      
    elseif drywet_mode == 4           str_drywet = '_JJA-DJF';          end    %/ sources during dry days are more interesting.

    %/ Will not proceed if X_daily is empty. Return str_drywet.
    if isempty(X_daily)               warning('Empty X_daily. Skip computing drywet anomalies.');   return;  end
    
    %/ Only handle 3D data
    if length(size(X_daily)) ~= 3     error('this function only handles 3D X_daily !');  end
         
    ind_X_dry_days       = findismember_loop(X_dates, dry_days);
    ind_X_wet_days       = findismember_loop(X_dates, wet_days);
    ind_X_normal_days    = findismember_loop(X_dates, normal_days);
    
    if isempty(ind_X_dry_days)                                                error('ind_dry_days is empty!!');        end
    if isempty(ind_X_wet_days)                                                error('ind_wet_days is empty!!');        end
    if isempty(ind_X_normal_days) && (drywet_mode == 1 || drywet_mode == 2)   error('ind_normal_days is empty!!');     end
 
    fprintf('*** Compute %s anomalies... ***\n', strrep(str_drywet, '_', ''));
    
    if drywet_mode == 1   %/ Dry minus Normal
        dataX = X_daily(:,:,ind_X_dry_days);
        dataY = X_daily(:,:,ind_X_normal_days);
        
    elseif drywet_mode == 2   %/ Wet minus Normal
        dataX = X_daily(:,:,ind_X_wet_days);
        dataY = X_daily(:,:,ind_X_normal_days);
    
    elseif drywet_mode == 3 || drywet_mode == 4   %/ 3]: Dry minus Wet. 4]: JJA minus DJF.
        dataX = X_daily(:,:,ind_X_dry_days);
        dataY = X_daily(:,:,ind_X_wet_days);

    end
    
    [X_diff_sig, X_diff] = ttest2_sig_fn(dataX, dataY, alpha, 3); %/ two-sample t-test
    
end
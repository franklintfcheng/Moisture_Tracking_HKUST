function leadlag_data = hovmoller2leadlagcorr(varargin)
    
    pnames = {'hovdata', 'hori', 'hori_refpt', 'along_timelag_or_hori', 'sliding_range', 'zm_or_mm', 'timelag', 'timelag_full',...
              };
    
    dflts  = cell(length(pnames), 1);
    
    [           hovdata,   hori,   hori_refpt,   along_timelag_or_hori,   sliding_range,  zm_or_mm,   timelag,   timelag_full,...
               ] = ...
                                            internal.stats.parseArgs(pnames, dflts, varargin{:});
%%

    if length(size(hovdata)) ~= 2
        error('The input hovdata must be 2D (lag-lon or lag-lat hovdata)!');
    end
    hori_res     = diff(abs(hori(1:2)));
    nhori        = size(hovdata, 1);  %/ Assume in (hori, timelag)
    [B, i0]      = min(abs(hori-hori_refpt));
    t0           = find(timelag_full == 0); %/ the index of lag 0 in timelag_bc
    
    if isempty(t0) || B > hori_res
        error('Check inputs of timelag and contf_hori!');
    end

    leadlag_data = nan(size(hovdata));   %/ hori, time
    if along_timelag_or_hori == 1
        if isempty(sliding_range)
            sliding_range = length(timelag);
        end
        sliding_grids = sliding_range;  %/ The sliding window for lead lag correation

        for lag = timelag
            t   = find(timelag_full == lag);
            tw  = (t-(sliding_grids-1)/2):(t+(sliding_grids-1)/2);
            tw0 = (t0-(sliding_grids-1)/2):(t0+(sliding_grids-1)/2); %/ window centered at lag 0
            for i = 1:nhori  %/ Lead-lag correlation (w.r.t. (lon0, lag0))
                X = squeeze(hovdata(i0, tw0))';
                Y = squeeze(hovdata(i,  tw))';
                leadlag_data(i,t) = corr(X, Y);
            end
        end
        %/ Restore the inqueried range of the lead-lag dim
        ind = findismember_loop(timelag_full, timelag);
        leadlag_data = leadlag_data(:,ind);
    
    elseif along_timelag_or_hori == 2   %/ ** Recommended **
        if isempty(sliding_range)
            if zm_or_mm == 1
                sliding_range = 40;   %/ in terms of lat deg 
            elseif zm_or_mm == 2
                sliding_range = 160;  %/ in terms of lon deg 
            end
        end
        sliding_grids = floor(sliding_range/hori_res);  %/ total number of grids for 80E (varies across different resolutions)
                        
        if mod(sliding_grids, 2) == 0
            iw0 = (i0-(sliding_grids/2)):(i0+(sliding_grids/2));
        else
            iw0 = (i0-((sliding_grids-1)/2)):(i0+((sliding_grids-1)/2));
        end
        
        if ~isempty(find(iw0 < 1, 1))
            error('The sliding_range at the reference point (hori_refpt) exceeds the horizontal dim!')
        end
        
        for lag = timelag
            t = find(timelag == lag);
            for i = 1:nhori   %/ Lead-lag correlation (w.r.t. (lon0, lag0))
                if mod(sliding_grids, 2) == 0
                    iw = (i-(sliding_grids/2)):(i+(sliding_grids/2));
                else
                    iw = (i-((sliding_grids-1)/2)):(i+((sliding_grids-1)/2));
                end
                ind_out = (iw < 1 | iw > nhori);
                iw(ind_out) = [];
                X = squeeze(hovdata(iw0, t0));
                Y = squeeze(hovdata(iw,   t));
                X(ind_out,:) = [];
                leadlag_data(i,t) = corr(X, Y, 'Rows', 'complete');
            end
        end
    else
        error('Invalid input of ''along_timelag_or_hori''!');
    end
end
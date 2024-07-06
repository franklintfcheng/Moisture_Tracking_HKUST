%%
function [related_TC] = IBTrACS_related_TCdays(varargin)

% create a set of valid parameters and their default values
pnames = {'allbtlons', 'allbtlats', 'allbtdates', 'allbtwind', 'inputdates', 'sinkreglon', 'sinkreglat', 'effradius_indeg'};

dflts  = {[] [] [] [] [] [] []};

%# parse function arguments
[allbtlons, allbtlats, allbtdates, allbtwind, inputdates, sinkreglon, sinkreglat, effradius_indeg] = internal.stats.parseArgs(pnames, dflts, varargin{:});

missing_TCdays = 0;
related_TCdates = [];
related_TCcount = [];
related_TCtrackid = {};
related_TCwind = {};
for t = 1:length(inputdates)
%     disp(t);
    [ind_row, ind_col] = find(ismember(allbtdates, inputdates(t)));
    Tracks = unique(ind_col);
    
    Ncount = 0;
    Ntracks = {};
    Nwind = {};
    for i = 1:length(Tracks)

        Tracks_timesteps = ind_row(find(ind_col == Tracks(i)));

        lons = allbtlons(Tracks_timesteps, Tracks(i)); %/ retrieve timesteps of a TC track concurrent with eventdates
        lats = allbtlats(Tracks_timesteps, Tracks(i));
        wind = allbtwind(Tracks_timesteps, Tracks(i));
        
        nonnan_ind = find(~isnan(lons)); %/ extract non-NaN timesteps
        lons = lons(nonnan_ind);
        lats = lats(nonnan_ind);
        wind = wind(nonnan_ind);
%         disp(wind)
        if isempty(lons)
            missing_TCdays = missing_TCdays + 1;
            continue %/skip if the selected best track dataset have no data on this day.
        end
        
        for pt = 1:length(lons)
            [arclen, az] = distance('rh', sinkreglat, sinkreglon, lats(pt), lons(pt)); %rumbline distance in degree
%                     distkm = deg2km(arclen);
            if arclen < effradius_indeg %/ in degree
                Ncount = Ncount + 1;
                Ntracks = [Ntracks; {Tracks(i)}];
                Nwind = [Nwind; {wind(pt)}];
%                 if  %/ if more than one track within the effective range on the day
%                     related_TCcount(end) = Ncount; %/then modify the last (most recent) element.
%                     related_TCtrackid(end) = [related_TCtrackid(end), {Tracks(i)}];
                    
%                     related_TCtrackid(end) = {cat(1, related_TCtrackid{end}, {Tracks(i)})}; %/ append cell array (e.g. 1x1 cell to 2x1 cell
%                     related_TCtrackid = cat(1, related_TCtrackid, {Tracks(i)});
%                 end
                break;
            end
        end
    end

%     disp(Ntracks)
%     disp(related_TCtrackid)
    related_TCdates = cat(1, related_TCdates, inputdates(t));
    related_TCcount = cat(1, related_TCcount, Ncount);
    related_TCtrackid = [related_TCtrackid; {[Ntracks{:}]}]; %/ {[A{:}]} to remove unnecessary inner brackets
    related_TCwind = [related_TCwind; {[Nwind{:}]}]; %/ {[A{:}]} to remove unnecessary inner brackets

end
% B = unique(related_TCdates);
% Ncount = histc(related_TCdates, B);
% related_TC = [B, Ncount];
% disp(size(related_TCdates))
% disp(size(related_TCcount))
% disp(size(related_TCtrackid))
% disp(related_TCtrackid)
related_TC = [num2cell(related_TCdates), num2cell(related_TCcount), related_TCtrackid, related_TCwind];

% ind = find(ismember(inputdates, [related_TC{:,1}]));
% related_TC_full = cell(length(inputdates), 3);
% related_TC_full(:,1) = num2cell(inputdates);
% related_TC_full(ind,2) = related_TC(:,2);
% related_TC_full(ind,3) = related_TC(:,3);
% disp(related_TC_full)
% 
% ind = find(isempty(related_TC_full(:,2)));
% related_TC_full(ind,2) = {0};
% disp(related_TC_full)

% fprintf('Missing TC days = %d', missing_TCdays);
end
%%
function [lead, season] = hs_lead_season_set(varargin)

pnames = {'indexname',   'hotspot'};
dflts  = {[], []};
[indexname, hotspot] = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

%/ lead: 0 to 3
%/ season: 1. MAM, 2. JJA, 3. SON, 4. DJF

lead = []; season = [];         
if ismember(indexname, {'NINO34'})
    
    if ismember(hotspot, {'Amazon'})                lead = 0; season = 1;   end
    if ismember(hotspot, {'Orinoco'})               lead = 1; season = 1;   end
    if ismember(hotspot, {'La Plata'})              lead = 0; season = 1;   end
    if ismember(hotspot, {'NE. South America'})     lead = 0; season = 1;   end
    if ismember(hotspot, {'Shebelle'})              lead = 1; season = 3;   end
    if ismember(hotspot, {'Lake Turkana-Swash'})    lead = 1; season = 3;   end
    if ismember(hotspot, {'Ganges'})                lead = 3; season = 2;   end
    if ismember(hotspot, {'New Guinea'})            lead = 0; season = 3;   end
    if ismember(hotspot, {'Pearl'})                 lead = 0; season = 3;   end
end
    
if ismember(indexname, {'IOD'})
    if ismember(hotspot, {'Congo'})                 lead = 2; season = 3;   end
    if ismember(hotspot, {'Orinoco'})               lead = 1; season = 4;   end
    if ismember(hotspot, {'Shebelle'})              lead = 0; season = 3;   end
    if ismember(hotspot, {'Lake Turkana-Swash'})    lead = 0; season = 2;   end
    if ismember(hotspot, {'Ganges'})                lead = 0; season = 3;   end
    if ismember(hotspot, {'New Guinea'})            lead = 0; season = 3;   end
    if ismember(hotspot, {'Irrawaddy-Salween'})     lead = 0; season = 3;   end
end

if ismember(indexname, {'AO'})
    if ismember(hotspot, {'Nile'})                  lead = 3; season = 4;   end
    if ismember(hotspot, {'NE. South America'})     lead = 2; season = 3;   end
    if ismember(hotspot, {'Shebelle'})              lead = 2; season = 3;   end
    if ismember(hotspot, {'Lake Turkana-Swash'})    lead = 2; season = 3;   end
    if ismember(hotspot, {'New Guinea'})            lead = 2; season = 3;   end
    if ismember(hotspot, {'Mekong'})                lead = 1; season = 2;   end
    if ismember(hotspot, {'Irrawaddy-Salween'})     lead = 0; season = 4;   end
end

if ismember(indexname, {'WP'})
    if ismember(hotspot, {'Yangtze'})               lead = 3; season = 4;   end
    if ismember(hotspot, {'Pearl'})                 lead = 0; season = 1;   end
end

if ismember(indexname, {'PDO'})
%     if ismember(hotspot, {'Amazon'})                lead = 3; season = 1;   end
%     if ismember(hotspot, {'Orinoco'})               lead = 2; season = 1;   end
%     if ismember(hotspot, {'NE. South America'})     lead = 0; season = 2;   end
    if ismember(hotspot, {'Yangtze'})               lead = 0; season = 5;   end
    if ismember(hotspot, {'New Guinea'})            lead = 1; season = 5;   end
end

if ismember(indexname, {'PNA'})
    if ismember(hotspot, {'Orinoco'})               lead = 0; season = 4;   end
    if ismember(hotspot, {'Ganges'})                lead = 2; season = 2;   end
    if ismember(hotspot, {'Mekong'})                lead = 2; season = 2;   end
    if ismember(hotspot, {'Pearl'})                 lead = 3; season = 2;   end
end

if ismember(indexname, {'NAO'})
    if ismember(hotspot, {'Nile'})                  lead = 3; season = 4;   end
end

end

















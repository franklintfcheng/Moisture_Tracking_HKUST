%%
function W = omega2w(varargin)
% ----------------------------
%Author: Fandy Cheng
%Date of creation: 25 Feb 2020

%/ Convert pressure velocity (omega) to vertical velocity /%
%/ Knowledege of the layer's temperature is required /%

pnames = {'level', 'omega', 'T'};
dflts  = {[] [] []};
[level, omega, T] = internal.stats.parseArgs(pnames, dflts, varargin{:});

W = nan(size(omega));
for i = 1:length(level)
    R = 287;
    g = 9.81;
    P = level(i)*100; %hPa -> Pa

    %MAKE SURE T is in Kelvin.
    W(:,:,i,:) = -1* omega(:,:,i,:) * R .* T(:,:,i,:)/P/g; % w = -omega*RT/(Pg)

end
end



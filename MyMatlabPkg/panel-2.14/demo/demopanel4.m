
% Panels can be any size.
%
% (a) Create an asymmetrical grid of panels.
% (b) Create another.
% (c) Use select('all') to load them all with axes
% (d) Get handles to all the axes and modify them.
close all
clc
figure
p = panel();
% p.pack(2, 2);
% p(1).repack(0.05);
% p(2,1).pack(3,3);
% p(2,2).repack(0.05);
% p(2,2).pack(3,1);


p.pack(2, 2);
p(1).repack(1e-10); %/ placeholder to avoid over-cutting at the top.
p(2,1).pack(panel_row,panel_col);
p(2,2).repack(0.05);
p(2,2).pack(3,1);
        
        
p.select('all');
p.identify();
%% Panels can be any size.
%
% (a) Create an asymmetrical grid of panels.
% (b) Create another.
% (c) Use select('all') to load them all with axes
% (d) Get handles to all the axes and modify them.
close all
% clf
clc
fig = figure;
p = panel(fig);
p.pack(2, 2);
p(1).repack(1e-10); %/ placeholder to avoid over-cutting at the top.
p(2,1).pack(3,3);
p(2,2).repack(0.5);
% p.pack({15 [] 15}, {70 30});
% p(2,1).pack({90 10}, {[]});
% p(2,2).pack({50 50}, {[]}); % create subpanels 
% 
% p(2,2,1,1).pack({[]}, {80 20});
% p(2,2,2,1).pack({[]}, {80 20});

% p(6).repack(0.05);
% p.de.margin = 5;
p.select('all');
p.identify();
% p2 = panel(u2,'add');
% p2.pack(3,1);
% p2.de.margin = 5;
% p2.select('all');

%%
% p(1, 1, 1, 1).margin = [0 0 0 0];% main figure

% p(1, 1, 1, 2).pack({50 50}, {[]}); % create subpanels 
% p(1, 1, 1, 2).margin = [0 0 0 0]; % colorbar, %left bot right top
% p(1, 1, 1, 3).margin = [0 0 0 0]; % placeholder
% p(1, 1, 2).select();
% p(1, 1, 1, 1).select();
% p(1, 1, 1, 1, 2, 1).marginright = 10;
% p(1, 2, 1, 2).marginright = 10;
% 
% p(1, 2).pack({[]}, {[] 5 5});
% p(1, 2, 1, 1).marginright = 3;
% p(1, 2, 1, 2).marginleft = 0;
% p(1, 2, 1, 3).marginleft = 0;
% p(1, 2, 1, 2).marginright = 10;

p.select('all');
h_axes = p.de.axis;
p.identify();
% so then we might want to set something on them.
% set(h_axes, 'color', [0 0 0]);
% (a)

% create a 2x2 grid in gcf with different fractionally-sized
% rows and columns. a row or column sized as "[]" will
% stretch to fill the remaining unassigned space.
p = panel();
p.pack({1/3 []}, {1/3 []});




%% (b)

% pack a 2x3 grid into p(2, 2). note that we can pack by
% percentage as well as by fraction - the interpretation is
% just based on the size of the numbers we pass in (1 to
% 100 for percentage, or 0 to 1 for fraction).
p(2, 2).pack({30 70}, {20 20 []});



%% (c)

% use select('all') to quickly show the layout you've achieved.
% this commits all uncommitted panels as axis panels, so
% they can't be parents anymore (i.e. they can't have more
% children pack()ed into them).
%
% this is no use at all once you've got organised - look at
% the first three demos, which don't use it - but it may help
% you to see what you're doing as you're starting out.
p.select('all');



%% (d)

% whilst we're here, we can get all the axes within a
% particular panel like this. there are three "groups"
% associated with a panel: (fa)mily, (de)scendants, and
% (ch)ildren. see "help panel/descendants", for instance, to
% see who's in them. they're each useful in different
% circumstances. here, we use (de)scendants.
h_axes = p.de.axis;

% so then we might want to set something on them.
set(h_axes, 'color', [0 0 0]);

% yeah, real gothic.


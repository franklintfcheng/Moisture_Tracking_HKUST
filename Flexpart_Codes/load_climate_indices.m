%% Load climate Indices
clc; close all;

year_list = 1971:2010;

ClimIndices = [];
datafolder = '/disk/r128/tfchengac/fandyr128SOM/Climate Indices/';

indexname =  {'MJO',           'GMT',     'SAM',   'IOD_HadISST',  'EP_ENSO', 'CP_ENSO', 'Mixed_ENSO',    'ONI', 'NINO34_HadISST',  'NINO3_HadISST',  'NINO4_HadISST',    'AMO',    'PDO',     'NP',    'PNA',     'WP',     'AO',    'NAO',  'SolarFlux'};
dataformat = {'%f',             '%f',      '%f',            '%f',       '%f',      '%f',         '%f',     '%f',             '%f',             '%f',             '%f',     '%f',     '%f',     '%f',     '%f',     '%f',     '%f',    '%f,',     '%f'};
type =       {'MJO_format', 'Stndrd',  'Stndrd',        'Stndrd',   'Stndrd',  'Stndrd',     'Stndrd', 'Stndrd',         'Stndrd',         'Stndrd',         'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'TwoCol', 'Stndrd', 'OneCol',   'OneCol'};
skipheadings=[   2,                0,        0,                0,          0,         0,            0,        0,                0,                0,                0,        0,        0,        0,        0,        9,        0,        0,        0];

% indexname =  {'SAM',    'ONI',    'NINO34', 'NINO4',  'AMO',    'PDO',    'NP',     'PNA',    'WP',     'AO',     'NAO',    'SolarFlux',  'EuraSnowCover', 'NAGLSnowCover', 'NHSnowCover'};
% dataformat = {'%f',     '%f',      '%f',     '%f',     '%f',     '%f',     '%f',     '%f',     '%f',     '%f',     '%f,',    '%f',         '%f,',           '%f,',           '%f,'};
% type =       {'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'TwoCol', 'Stndrd', 'OneCol', 'OneCol',     'OneCol',        'OneCol',        'OneCol'};
% skipheadings=[   0,        0,        0,        0,        0,        0,        0,        0,        9,        0,        0,        0,            4,               4,               4];
for i = 1:length(indexname)
    a = read_ClimIndex('type', type{i}, 'datafolder', datafolder, 'skipheadings', skipheadings(i), 'indexname', indexname{i},...
                       'select_year', year_list, 'dataformat', dataformat{i});
    fld = fieldnames(a);
    for f= 1:length(fld)
        ClimIndices.(fld{f}) = a.(fld{f});
    end
end

%/ EP ENSO (ONI)
%  Warm (red) and cold (blue) periods based on a threshold of +/- 0.5oC for the Oceanic Nino Index (ONI)
%  for at least 5 consecutive overlapping seasons
%  https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php

% NOTE: my result is slightly different from that in the website which
%       is based on a ONI index rounded to 1 digit only.

ClimIndices.ENSO_event = CI2events('CI_dates', ClimIndices.ONI(:,1), 'CI_data', ClimIndices.ONI(:,2), 'thres', 0.5, 'n', 5);
% ClimIndices.EP_ENSO_event = CI2events('CI_dates', ClimIndices.EP_ENSO(:,1), 'CI_data', ClimIndices.EP_ENSO(:,2), 'thres', 0.7, 'n', 5); %/ my guess
% ClimIndices.CP_ENSO_event = CI2events('CI_dates', ClimIndices.CP_ENSO(:,1), 'CI_data', ClimIndices.CP_ENSO(:,2), 'thres', 0.7, 'n', 5); %/ my guess

%/ IOD
%  For monitoring the IOD, Australian climatologists consider sustained values above +0.4 °C 
%  as typical of a positive IOD, and values below −0.4 °C as typical of a negative IOD.
%  http://www.bom.gov.au/climate/enso/indices/about.shtml

%/ NOTE: Since no statement on how 'sustain' it should be. Assume n = 5 months.

ClimIndices.IOD_event = CI2events('CI_dates', ClimIndices.IOD_HadISST(:,1), 'CI_data', ClimIndices.IOD_HadISST(:,2), 'thres', 0.4, 'n', 5);

ClimIndices.AMO_event = CI2events('CI_dates', ClimIndices.AMO(:,1), 'CI_data', ClimIndices.AMO(:,2), 'thres', 0.1, 'n', 5); %/ here 0.1 refers to 0.5sd

%/ MJO events 
%/ Amp > 1 and persistent for >= 3 days (my choice)
MJO_dates = ClimIndices.MJO(:,1);
MJO_amp   = ClimIndices.MJO(:,2);
MJO_phase = ClimIndices.MJO(:,3);

thres = 1;
% n = 1;     
n = 3;
M_hankel_amp   = hankel(MJO_amp(1:n),   MJO_amp(n:end));
M_hankel_phase = hankel(MJO_phase(1:n), MJO_phase(n:end));

MJO_mature_phase = zeros(size(MJO_phase));
for k = 1:8
    %/ find each moving n length of data with the same phase and strong enough amp
    mjo_phase_ind = find(all(M_hankel_amp > thres, 1) == 1   &...
                         all(M_hankel_phase == k, 1) == 1);  
                     
    for i = 1:length(mjo_phase_ind)
        MJO_mature_phase(mjo_phase_ind(i):(mjo_phase_ind(i) + n - 1), 1) = k; %/ NOTE: whatever a column of n data met the criteria, we mark *all the n data* instead of one.                                                       
    end
end
ClimIndices.MJO_event(:,1) = MJO_dates;
ClimIndices.MJO_event(:,2) = MJO_mature_phase;

% for k = 1:8
%     disp(length(find( MJO_mature_phase  == k)))
% end

%/ Zhang et al 2019's method to classify EP, CP and Mixed ENSO (Assume N3 and N4 have the same date)
N3_dates = ClimIndices.NINO3_HadISST(:,1);
N4_dates = ClimIndices.NINO4_HadISST(:,1);
if ~isequal(N3_dates, N4_dates)  error('Make sure N3 and N4 have the same date!!'); end

N3 = ClimIndices.NINO3_HadISST(:,2);
N4 = ClimIndices.NINO4_HadISST(:,2);

%/ 3-month running mean
N3 = movmean(N3, 3);
N4 = movmean(N4, 3);

ClimIndices.UCEI(:,1) = N3_dates;
ClimIndices.UCEI(:,2) = sqrt(2*(N3.^2 + N4.^2));                           %/ r
ClimIndices.UCEI(:,3) = atand((N3-N4)./(N3+N4));                           %/ theta, atand in degree

ind_La = find(N3+N4 < 0);
ClimIndices.UCEI(ind_La,3) = ClimIndices.UCEI(ind_La,3) - 180;
                              
complex_ENSO_event      = zeros(length(N3_dates), 2);                      
complex_ENSO_event(:,1) = N3_dates;

r        = ClimIndices.UCEI(:,2);
thres    = 0.5;   n = 5;                                                   %/ Condition: r > 0.5 for at least 5 months
M_hankel = hankel(r(1:n), r(n:end));
enso_ind = find(all(M_hankel > thres, 1) == 1);  
for i = 1:length(enso_ind)
    complex_ENSO_event(enso_ind(i):(enso_ind(i) + n - 1), 2) = 999;        %/ NOTE: whatever a column of 5 months met the criteria, we mark *all the 5 months* instead of one.
                                                                           %/       mark as 999 since we will further classify the type based on theta
end
theta_range =     [  15   90    1;        %/ EP    El (1)
                    -15   15    2;        %/ Mixed El (2)
                    -90  -15    3;        %/ CP    El (3)
                   -165  -90   -1;        %/ EP    La (-1)
                   -195  -165  -2;        %/ Mixed La (-2)
                   -270  -195  -3;];      %/ CP    La (-3)
                   
theta = ClimIndices.UCEI(:,3);
for i = 1:length(theta_range)
    
    ind_phase = find(complex_ENSO_event(:,2) == 999   &...
                     theta > theta_range(i,1) &...
                     theta < theta_range(i,2));
                 
    complex_ENSO_event(ind_phase,2) = theta_range(i,3);
end

ClimIndices.complex_ENSO_event = complex_ENSO_event;

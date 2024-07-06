function ClimateIndices = load_all_Clim_Events(varargin)
    
    % create a set of valid parameters and their default values
    pnames = {'datafolder', 'indexname', 'get_ENSO_events', 'get_complex_ENSO_events', 'get_IOD_events', 'get_AMO_events', 'get_MJO_phases', 'save_csv', 'fields_to_save', 'csv_path'};
    dflts  = cell(1, length(pnames));
    [          datafolder,   indexname,   get_ENSO_events,   get_complex_ENSO_events,   get_IOD_events,   get_AMO_events,    get_MJO_phases,   save_csv,   fields_to_save,   csv_path] ...
        = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    %%
    %/============================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Apr 24, 2024
    %/============================================================================
    
    %/ MAIN
    if isempty(datafolder)  
        datafolder = '/disk/r128/tfchengac/fandyr128SOM/Climate Indices/';    
    end
    
    %/ Database
    indexname_db    = {       'MJO',    'GMT',    'SAM', 'IOD_HadISST', 'EP_ENSO', 'CP_ENSO', 'Mixed_ENSO',    'ONI', 'NINO34_ersstv5', 'NINO3_ersstv5', 'NINO4_ersstv5', 'NINO34_HadISST', 'NINO3_HadISST', 'NINO4_HadISST',    'AMO',    'PDO',     'NP',    'PNA',     'WP',     'AO',    'NAO', 'SolarFlux'};
    dataformat_db   = {        '%f',     '%f',     '%f',          '%f',      '%f',      '%f',         '%f',     '%f',             '%f',            '%f',            '%f',             '%f',            '%f',            '%f',     '%f',     '%f',     '%f',     '%f',     '%f',     '%f',    '%f,',        '%f'};
    type_db         = {'MJO_format', 'Stndrd', 'Stndrd',      'Stndrd',  'Stndrd',  'Stndrd',     'Stndrd', 'Stndrd',         'Stndrd',        'Stndrd',        'Stndrd',         'Stndrd',        'Stndrd',        'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'Stndrd', 'TwoCol', 'Stndrd', 'OneCol',    'OneCol'};
    skipheadings_db = [           2,        0,        0,             0,         0,         0,            0,        0,                1,               1,               1,                0,               0,               0,        0,        0,        0,        0,        9,        0,        0,           0];
    missing_val_db  = [         nan,      nan,      nan,           nan,       nan,       nan,          nan,      nan,           -99.99,          -99.99,          -99.99,              nan,             nan,             nan,      nan,      nan,      nan,      nan,      nan,      nan,      nan,         nan];

    if isempty(indexname)
        indexname = indexname_db;
    end
    if get_ENSO_events
        ENSO_index = {'ONI'};
        indexname = [indexname, ENSO_index]; 
    end
    if get_complex_ENSO_events
        indexname = [indexname, {'NINO3_HadISST', 'NINO4_HadISST'}]; 
    end
    if get_IOD_events
        IOD_index = {'IOD_HadISST'};
        indexname = [indexname, IOD_index]; 
    end
    if get_AMO_events
        AMO_index = {'AMO'};
        indexname = [indexname, AMO_index]; 
    end
    if get_MJO_phases
        MJO_index = {'MJO'};
        indexname = [indexname, MJO_index]; 
    end

    ClimateIndices = [];
    ind_indexname  = findismember_loop(indexname_db, indexname);
    for k = ind_indexname
        X = read_ClimateIndex('type', type_db{k}, 'datafolder', datafolder, 'skipheadings', skipheadings_db(k), 'indexname', indexname_db{k},...
                              'dataformat', dataformat_db{k});
        fld = fieldnames(X);
        for f = 1:length(fld)
            a = X.(fld{f});
            if ~isnan(missing_val_db(k))
                a(a == missing_val_db(k)) = nan;
            end
            ClimateIndices.(fld{f}) = a;
        end
    end
    
    if get_ENSO_events
        %  Warm (red) and cold (blue) periods based on a threshold of +/- 0.5oC for the Oceanic Nino Index (ONI)
        %  for at least 5 consecutive overlapping seasons
        %  https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php
        
        % NOTE: my result is slightly different from that in the website, which is based on a ONI index rounded up to 1 digit.
    
        ClimateIndices.ENSO_event = CI2events('CI_dates', ClimateIndices.(ENSO_index{:})(:,1), 'CI_data', ClimateIndices.(ENSO_index{:})(:,2), 'thres', 0.5, 'n', 5);
        % ClimIndices.EP_ENSO_event = CI2events('CI_dates', ClimIndices.EP_ENSO(:,1), 'CI_data', ClimIndices.EP_ENSO(:,2), 'thres', 0.7, 'n', 5); %/ my choice
        % ClimIndices.CP_ENSO_event = CI2events('CI_dates', ClimIndices.CP_ENSO(:,1), 'CI_data', ClimIndices.CP_ENSO(:,2), 'thres', 0.7, 'n', 5); %/ my choice
    end

    if get_IOD_events
        %  For monitoring the IOD, Australian climatologists consider sustained values above +0.4 °C 
        %  as typical of a positive IOD, and values below −0.4 °C as typical of a negative IOD.
        %  http://www.bom.gov.au/climate/enso/indices/about.html
        
        %/ NOTE: Since no refs on how 'persistent' it should be. Assume n = 5 months.
        ClimateIndices.IOD_event = CI2events('CI_dates', ClimateIndices.(IOD_index{:})(:,1), 'CI_data', ClimateIndices.(IOD_index{:})(:,2), 'thres', 0.4, 'n', 5);
    end

    if get_AMO_events
        ClimateIndices.AMO_event = CI2events('CI_dates', ClimateIndices.(AMO_index{:})(:,1), 'CI_data', ClimateIndices.(AMO_index{:})(:,2), 'thres', 0.1, 'n', 5); %/ here 0.1 refers to 0.5 sd
    end

    if get_MJO_phases
        %/ Amp > 1 and persistent for >= 3 days (my choice)
        MJO_dates = ClimateIndices.(MJO_index{:})(:,1);
        MJO_amp   = ClimateIndices.(MJO_index{:})(:,2);
        MJO_phase = ClimateIndices.(MJO_index{:})(:,3);
        thres     = 1;
        n         = 3;
        M_hankel_amp   = hankel(MJO_amp(1:n),   MJO_amp(n:end));
        M_hankel_phase = hankel(MJO_phase(1:n), MJO_phase(n:end));
        
        MJO_mature_phase = zeros(size(MJO_phase));
        for k = 1:8
            %/ find each moving n length of data with the same phase and strong enough amp
            mjo_phase_ind = find(all(M_hankel_amp > thres, 1) &...
                                 all(M_hankel_phase == k,  1));        
            for i = 1:length(mjo_phase_ind)
                MJO_mature_phase(mjo_phase_ind(i):(mjo_phase_ind(i) + n - 1), 1) = k; %/ NOTE: whatever a column of n data met the criteria, we mark *all the n data* instead of one.                                                       
            end
        end
        ClimateIndices.MJO_event(:,1) = MJO_dates;
        ClimateIndices.MJO_event(:,2) = MJO_mature_phase;
        % for k = 1:8
        %     disp(length(find( MJO_mature_phase  == k)))
        % end
    end

    if get_complex_ENSO_events
        %====== EP, mixed, CP ENSO ======%
        %/ Zhang et al.'s (2019) method to classify EP, CP and Mixed ENSO (Assume N3 and N4 have the same date)
        N3_dates = ClimateIndices.NINO3_HadISST(:,1);
        N4_dates = ClimateIndices.NINO4_HadISST(:,1);
        if ~isequal(N3_dates, N4_dates)  
            error('Make sure N3 and N4 have the same date!!'); 
        end
    
        N3 = ClimateIndices.NINO3_HadISST(:,2);
        N4 = ClimateIndices.NINO4_HadISST(:,2);
    
        %/ 3-month running mean
        N3 = movmean(N3, 3);
        N4 = movmean(N4, 3);
    
        ClimateIndices.UCEI(:,1) = N3_dates;
        ClimateIndices.UCEI(:,2) = sqrt(2*(N3.^2 + N4.^2));                       %/ r
        ClimateIndices.UCEI(:,3) = atand((N3-N4)./(N3+N4));                       %/ theta, atand in degree
    
        ind_La = find(N3+N4 < 0);
        ClimateIndices.UCEI(ind_La,3) = ClimateIndices.UCEI(ind_La,3) - 180;
    
        complex_ENSO_event      = zeros(length(N3_dates), 2);                      
        complex_ENSO_event(:,1) = N3_dates;
        
        r        = ClimateIndices.UCEI(:,2);
        thres    = 0.5;   n = 5;                                               %/ Condition: r > 0.5 for at least 5 months
        M_hankel = hankel(r(1:n), r(n:end));
        enso_ind = find(all(M_hankel > thres, 1));  
        for i = 1:length(enso_ind)
            complex_ENSO_event(enso_ind(i):(enso_ind(i) + n - 1), 2) = 999;    %/ NOTE: whatever a column of 5 months met the criteria, we mark *all the 5 months* instead of one.
                                                                               %/ mark as 999 since we will further classify the type based on theta
        end
        theta_range =     [  15   90    1;        %/ EP    El (1)
                            -15   15    2;        %/ Mixed El (2)
                            -90  -15    3;        %/ CP    El (3)
                           -165  -90   -1;        %/ EP    La (-1)
                           -195  -165  -2;        %/ Mixed La (-2)
                           -270  -195  -3;];      %/ CP    La (-3)
                       
        theta = ClimateIndices.UCEI(:,3);
        for i = 1:length(theta_range)
            ind_phase = complex_ENSO_event(:,2) == 999   &...
                             theta > theta_range(i,1)    &...
                             theta < theta_range(i,2);
    
            complex_ENSO_event(ind_phase,2) = theta_range(i,3);
        end
        ClimateIndices.complex_ENSO_event = complex_ENSO_event;
        %     theta_range_name = {'EP El'; ...
        %                         'Mixed El';  ...
        %                         'CP El'; ...
        %                         'EP La'; ...
        %                         'Mixed La'; ...
        %                         'CP La'};        
            
        %     complex_ENSO_event{ind_phase,3} = theta_range_name{i};
    end

    %/ Save into csv (for other ppl to use)
    if save_csv
       for i = 1:length(fields_to_save)
           S = ClimateIndices.(fields_to_save{i});
           csv_filename = strcat(csv_path, fields_to_save{i}, '.csv');
           writematrix(S, csv_filename);
           fprintf('*** The climate index/event is saved into %s! ***\n', csv_filename);
       end
    end
end



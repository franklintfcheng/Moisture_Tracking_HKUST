function [dqc, RHc, RHc_map_lon, RHc_map_lat, str_RHc_dqc] = read_RHc_dqc_scheme(varargin)
    
    switch nargin
        case 1   %/ given 1 inputs
            RHc_dqc_scheme = varargin{1};
            forcing        = 'CERA';         %/ Assume CERA-20C forcing by default
            output_res     = 1;              %/ Use 1x1deg RHc map by default
        case 2   %/ given 2 inputs
            RHc_dqc_scheme = varargin{1};
            forcing        = varargin{2};
            output_res     = 1;
        case 3 
            RHc_dqc_scheme = varargin{1};
            forcing        = varargin{2};
            output_res     = varargin{3};
    otherwise
        error('Unexpected inputs')
    end
    % fprintf('*** Running read_RHc_dqc_scheme... ***\n')
    fprintf('*** Use RHc Scheme %d... ***\n', RHc_dqc_scheme)

    %/ For ERA5 RHc map
    if isequal(forcing, 'EA')  
        if RHc_dqc_scheme == 0  
            dqc = 0.1;                                      %/ critical dq for prcp to occur [g/kg]
            RHc = 80;                                       %/ critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = strcat('RHc', num2str(RHc), '_dqc', num2str(dqc)); %/ This is the old format.
            RHc_map_lon = [];
            RHc_map_lat = [];

        elseif RHc_dqc_scheme == 1  %/ Control setting
            dqc = 0.1;                                      %/ critical dq for prcp to occur [g/kg]
            RHc = 85;                                       %/ critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = strcat('RHc', num2str(RHc), '_dqc', num2str(dqc)); %/ This is the old format.
            RHc_map_lon = [];
            RHc_map_lat = [];

        elseif RHc_dqc_scheme == 2   
            dqc = 0.1;                                      %/ critical dq for prcp to occur [g/kg]
            RHc = 90;                                       %/ critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = strcat('RHc', num2str(RHc), '_dqc', num2str(dqc)); %/ This is the old format.
            RHc_map_lon = [];
            RHc_map_lat = [];

        elseif RHc_dqc_scheme == 3  
            dqc = 0.1;                                      %/ critical dq for prcp to occur [g/kg]
            RHc = 95;                                       %/ critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = strcat('RHc', num2str(RHc), '_dqc', num2str(dqc)); %/ This is the old format.
            RHc_map_lon = [];
            RHc_map_lat = [];

        else
            error('Not prescribed setting for ''RHc_dqc_scheme == %d''!', RHc_dqc_scheme) 
        end
        
    %/ For CERA-20C RHc map
    elseif isequal(forcing, 'CERA') && output_res == 1  

        if RHc_dqc_scheme == 1  %/ Control setting
            dqc = 0.1;                                      %/ critical dq for prcp to occur [g/kg]
            RHc = 85;                                       %/ critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = strcat('RHc', num2str(RHc), '_dqc', num2str(dqc)); %/ This is the old format.
            RHc_map_lon = [];
            RHc_map_lat = [];

        elseif RHc_dqc_scheme == 2
            dqc = 0.1;                                      %/ critical dq for prcp to occur [g/kg]
            %/ Load RHc_map
            RHc_levels   = [80, 85, 90];  
            thres_upbnd  = [-30,   30,   inf];
            thres_lowbnd = [-inf, -30,    30];  
            var          = 'P_LA_P_prct_diff';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map_%s_%s_thres_%s_1971-2010.mat',...
                                           strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'))];
            disp(RHc_map_matfilename)                         
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme == 3
            %/ Load RHc_map
            dqc = 0.1;                                      %/ critical dq for prcp to occur [g/kg]
            RHc_levels   = [75, 85, 95];  
            thres_upbnd  = [-1,   1,   inf];
            thres_lowbnd = [-inf, -1,    1];  
            var          = 'P_LA_P_diff';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map_%s_%s_thres_%s_1971-2010.mat',...
                                           strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'))];
            disp(RHc_map_matfilename)                         
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme == 4
            dqc = 0.1;                                      %/ critical dq for prcp to occur [g/kg]
            %/ Load RHc_map
            RHc_levels  = [70, 80, 85, 90, 100];
            thres_upbnd  = [-200,  -30,  30, 200, inf];
            thres_lowbnd = [-inf, -200, -30,  30, 200]; 
            var          = 'P_LA_P_prct_diff';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map_%s_%s_thres_%s_1971-2010.mat',...
                                           strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'))];
            disp(RHc_map_matfilename)                         
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme == 5
            dqc = 0.15;                                      %/ critical dq for prcp to occur [g/kg]
            RHc = 85;                                        %/ critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = strcat('RHc', num2str(RHc), '_dqc', num2str(dqc)); %/ don't change this format!
            RHc_map_lon = [];
            RHc_map_lat = [];

        elseif RHc_dqc_scheme == 6 
            dqc = 0.2;                                     %/ critical dq for prcp to occur [g/kg]
            RHc_levels    = [75, 80, 85, 90, 100];   %/ RHc for underestimation, normal, overestimation.
            thres_upbnd   = [ -40,   40,  inf];     %/ For prct diff (upper bound)
            thres_lowbnd  = [-inf,  -40,   40];     %/ For prct diff (lower bound)
            thres2_upbnd  = [  -2,    2,  inf];     %/ For abs  diff (upper bound)
            thres2_lowbnd = [-inf,   -2,    2];     %/ For abs  diff (lower bound)
            var           = 'P_LA_P_prct_diff';
            var2          = 'P_LA_P_diff';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map_%s_%s_thres_%s_%s_thres2_%s_1971-2010.mat',...
                                           strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'),...
                                           var2, strjoin(string(thres2_upbnd),'_'))];
            disp(RHc_map_matfilename)                         
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme == 7
            dqc = 0.2;                                     %/ critical dq for prcp to occur [g/kg]
            RHc_levels    = [75, 80, 85, 90, 95, 100];   %/ RHc for underestimation, normal, overestimation.
            thres_upbnd   = [ -80,   80,  inf];     %/ For prct diff (upper bound)
            thres_lowbnd  = [-inf,  -80,   80];     %/ For prct diff (lower bound)
            thres2_upbnd  = [  -2,    2,  inf];     %/ For abs  diff (upper bound)
            thres2_lowbnd = [-inf,   -2,    2];     %/ For abs  diff (lower bound)
            var           = 'P_LA_P_prct_diff';
            var2          = 'P_LA_P_diff';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map_%s_%s_thres_%s_%s_thres2_%s_1971-2010.mat',...
                                           strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'),...
                                           var2, strjoin(string(thres2_upbnd),'_'))];
            disp(RHc_map_matfilename)                         
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme == 8
            dqc = 0.2;                              %/ critical dq for prcp to occur [g/kg]
            thres_upbnd   = [ -80,   80,  inf];     %/ For prct diff (upper bound)
            thres_lowbnd  = [-inf,  -80,   80];     %/ For prct diff (lower bound)
            thres2_upbnd  = [  -2,    2,  inf];     %/ For abs  diff (upper bound)
            thres2_lowbnd = [-inf,   -2,    2];     %/ For abs  diff (lower bound)
            thres3_topo   = [1500];

            RHc_levels_2D = [75 80 85;
                             80 85 95;
                             85 90 95;];

            RHc_levels = unique(RHc_levels_2D)';   %/ RHc for underestimation, normal, overestimation.
            RHc_levels = [RHc_levels, 120];        %/ for areas with altitude > 1500 m.

            var           = 'P_LA_P_prct_diff';
            var2          = 'P_LA_P_diff';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map_%s_%s_thres_%s_%s_thres2_%s_topo%d_1971-2010.mat',...
                                           strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'),...
                                           var2, strjoin(string(thres2_upbnd),'_'), thres3_topo)];
            disp(RHc_map_matfilename)                         
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme == 9
            %/ NOTE: Scheme 9 is idnetical to Scheme 8, but I have
            %/       corrected the calculation of RH using both
            %/       particle's z-pos and topo height.
            dqc = 0.2;                              %/ critical dq for prcp to occur [g/kg]
            thres_upbnd   = [ -80,   80,  inf];     %/ For prct diff (upper bound)
            thres_lowbnd  = [-inf,  -80,   80];     %/ For prct diff (lower bound)
            thres2_upbnd  = [  -2,    2,  inf];     %/ For abs  diff (upper bound)
            thres2_lowbnd = [-inf,   -2,    2];     %/ For abs  diff (lower bound)
            thres3_topo   = [1500];

            RHc_levels_2D = [75 80 85;
                             80 85 95;
                             85 90 95;];

            RHc_levels = unique(RHc_levels_2D)';   %/ RHc for underestimation, normal, overestimation.
            RHc_levels = [RHc_levels, 120];        %/ for areas with altitude > 1500 m.

            var           = 'P_LA_P_prct_diff';
            var2          = 'P_LA_P_diff';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map_%s_%s_thres_%s_%s_thres2_%s_topo%d_1971-2010.mat',...
                                           strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'),...
                                           var2, strjoin(string(thres2_upbnd),'_'), thres3_topo)];
            disp(RHc_map_matfilename)                         
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme == 10
            dqc = 0.1;                                      %/ critical dq for prcp to occur [g/kg]
            RHc = 85;                                       %/ critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);
            RHc_map_lon = [];
            RHc_map_lat = [];

        elseif RHc_dqc_scheme == 11
            %/ NOTE: Scheme 9 is identical to Scheme 8, but I have
            %/       corrected the calculation of RH using both
            %/       particle's z-pos and topo height.
            dqc = 0.1;                              %/ critical dq for prcp to occur [g/kg]
            thres_upbnd   = [ -20,   20,  inf];     %/ For prct diff (upper bound)
            thres_lowbnd  = [-inf,  -20,   20];     %/ For prct diff (lower bound)
            thres2_upbnd  = [  -1,    1,  inf];     %/ For abs  diff (upper bound)
            thres2_lowbnd = [-inf,   -1,    1];     %/ For abs  diff (lower bound)
        %         thres3_topo   = [1500];

            %/              ---> i (cols): prct diff
            RHc_levels_2D = [77.5  80  85;       % | 
                             82.5  85  87.5;     % | j (rows): abs diff
                             85    90  92.5;];   % v

            RHc_levels = unique(RHc_levels_2D)';   %/ RHc for underestimation, normal, overestimation.

            var                 = 'P_LA_P_prct_diff';
            var2                = 'P_LA_P_diff';
            str_baseline        = '_BaselineScheme10';
            str_baseline_years  = '2010-2010';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map%s_%s_%s_thres_%s_%s_thres2_%s_%s.mat',...
                                           str_baseline, strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'),...
                                           var2, strjoin(string(thres2_upbnd),'_'), str_baseline_years)];
            disp(RHc_map_matfilename)                         
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme == 12
            dqc           = 0.1;                    %/ critical dq for prcp to occur [g/kg]
            thres_upbnd   = [ -20,   20,  inf];     %/ For prct diff (upper bound)
            thres_lowbnd  = [-inf,  -20,   20];     %/ For prct diff (lower bound)
            thres2_upbnd  = [  -1,    1,  inf];     %/ For abs  diff (upper bound)
            thres2_lowbnd = [-inf,   -1,    1];     %/ For abs  diff (lower bound)
            thres3_topo   = [3000];

            %/              ---> i (cols): prct diff
            RHc_levels_2D = [80    80  85;     % | 
                             77.5  85  90;     % | j (rows): abs diff
                             85    90  90;];   % v

            RHc_levels = unique(RHc_levels_2D)';   %/ RHc for underestimation, normal, overestimation.
            RHc_TP = 75;        %/ a lower RHc for TP
            RHc_levels = unique([RHc_levels, RHc_TP]);

            var           = 'P_LA_P_prct_diff';
            var2          = 'P_LA_P_diff';
            str_baseline        = '_BaselineScheme10';
            str_baseline_years  = '2010-2010';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map%s_%s_%s_thres_%s_%s_thres2_%s_topo%d_%s.mat',...
                                           str_baseline, strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'),...
                                           var2, strjoin(string(thres2_upbnd),'_'), thres3_topo, str_baseline_years)];
            disp(RHc_map_matfilename)                         
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme == 13
            dqc           = 0.1;
            thres_upbnd   = [ -20,   20,  inf];     %/ For prct diff (upper bound)
            thres_lowbnd  = [-inf,  -20,   20];     %/ For prct diff (lower bound)
            thres2_upbnd  = [  -1,    1,  inf];     %/ For abs  diff (upper bound)
            thres2_lowbnd = [-inf,   -1,    1];     %/ For abs  diff (lower bound)

            %/              ---> i (cols): prct diff
            RHc_levels_2D = [80    80  85;     % | 
                             80    85  90;     % | j (rows): abs diff
                             85    90  90;];   % v

            RHc_levels = unique(RHc_levels_2D)';   %/ RHc for underestimation, normal, overestimation.

            %/ RHc map for Pan-TP region (>1500m)
            thres3_TP_upbnd  = [  -2, -1.5,   -1, inf]; %/ For abs  diff (upper bound)
            thres3_TP_lowbnd = [-inf,   -2, -1.5,  -1]; %/ For abs  diff (lower bound)
            thres3_topo      = [1500];
            RHc_TP           = [65, 70, 75, 85];        %/ a lower RHc for TP
            RHc_levels       = unique([RHc_levels, RHc_TP]);

            var           = 'P_LA_P_prct_diff';
            var2          = 'P_LA_P_diff';
            str_baseline        = '_BaselineScheme10';
            str_baseline_years  = '2010-2010';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map%s_%s_%s_thres_%s_%s_thres2_%s_topo%d%s_%s.mat',...
                                           str_baseline, strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'),...
                                           var2, strjoin(string(thres2_upbnd),'_'), thres3_topo, '_pan_TP', str_baseline_years)];
            disp(RHc_map_matfilename)                         
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme == 14
            dqc           = 0.1;
            thres_upbnd   = [ -20,   20,  inf];     %/ For prct diff (upper bound)
            thres_lowbnd  = [-inf,  -20,   20];     %/ For prct diff (lower bound)
            thres2_upbnd  = [  -1,    1,  inf];     %/ For abs  diff (upper bound)
            thres2_lowbnd = [-inf,   -1,    1];     %/ For abs  diff (lower bound)

            %/              ---> i (cols): prct diff
            RHc_levels_2D = [80    80  85;     % | 
                             80    85  90;     % | j (rows): abs diff
                             85    90  90;];   % v

            RHc_levels = unique(RHc_levels_2D)';   %/ RHc for underestimation, normal, overestimation.

            %/ RHc map for Pan-TP region (>1500m)
            thres3_TP_upbnd  = [ -60,  -40,  -20, inf]; %/ For prct  diff (upper bound)
            thres3_TP_lowbnd = [-inf,  -60,  -40, -20]; %/ For prct  diff (lower bound)
            thres4_TP_upbnd  = [-2.5,  -2, -1.5,   -1, inf]; %/ For abs  diff (upper bound)
            thres4_TP_lowbnd = [-inf, -2.5,   -2, -1.5,  -1]; %/ For abs  diff (lower bound)
            thres3_topo      = [1500];

            %/              ---> i (cols): prct diff
            RHc_TP        = [
                             65 70 75 80;     % | 
                             70 75 80 85;     % | 
                             70 75 80 85;     % | j (rows): abs diff
                             70 75 80 85;     % |
                             70 75 80 85;];   % v

            RHc_levels    = unique([RHc_levels, unique(RHc_TP)']);
            var           = 'P_LA_P_prct_diff';
            var2          = 'P_LA_P_diff';
            str_baseline        = '_BaselineScheme10';
            str_baseline_years  = '2010-2010';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map%s_%s_%s_thres_%s_%s_thres2_%s_topo%d_pan_TP_scheme%d_%s.mat',...
                                           str_baseline, strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'),...
                                           var2, strjoin(string(thres2_upbnd),'_'), thres3_topo, RHc_dqc_scheme, str_baseline_years)];
            disp(RHc_map_matfilename)
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme == 15
            %/ NOTE: Scheme 15 uses an averaged topo map consistent
            %/       with that used in 'show_land_or_ocean'.
            %/       Also, the scheme is again fine-tuned.
            dqc           = 0.1;
            thres_upbnd   = [ -20,   20,  inf];     %/ For prct diff (upper bound)
            thres_lowbnd  = [-inf,  -20,   20];     %/ For prct diff (lower bound)
            thres2_upbnd  = [  -1,    1,  inf];     %/ For abs  diff (upper bound)
            thres2_lowbnd = [-inf,   -1,    1];     %/ For abs  diff (lower bound)

            %/              ---> i (cols): prct diff
            RHc_levels_2D = [80    80  85;     % | 
                             80    85  90;     % | j (rows): abs diff
                             85    90  90;];   % v

            RHc_levels = unique(RHc_levels_2D)';   %/ RHc for underestimation, normal, overestimation.

            %/ RHc map for Pan-TP region (>1500m)
            thres3_TP_upbnd  = [ -60,  -40,  -20, inf]; %/ For prct  diff (upper bound)
            thres3_TP_lowbnd = [-inf,  -60,  -40, -20]; %/ For prct  diff (lower bound)
            thres4_TP_upbnd  = [-2.5,  -2, -1.5,   -1, inf]; %/ For abs  diff (upper bound)
            thres4_TP_lowbnd = [-inf, -2.5,   -2, -1.5,  -1]; %/ For abs  diff (lower bound)
            thres3_topo      = [1500];

            %/              ---> i (cols): prct diff
            RHc_TP        = [
                             65 70 75 80;     % | 
                             70 75 80 85;     % | 
                             70 75 80 85;     % | j (rows): abs diff
                             70 75 80 85;     % |
                             79 81 83 85;];   % v

            RHc_levels    = unique([RHc_levels, unique(RHc_TP)']);
            var           = 'P_LA_P_prct_diff';
            var2          = 'P_LA_P_diff';
            str_baseline        = '_BaselineScheme10';
            str_baseline_years  = '2010-2010';
            RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map%s_%s_%s_thres_%s_%s_thres2_%s_topo%d_pan_TP_scheme%d_%s.mat',...
                                           str_baseline, strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'),...
                                           var2, strjoin(string(thres2_upbnd),'_'), thres3_topo, RHc_dqc_scheme, str_baseline_years)];
            disp(RHc_map_matfilename)
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);

        elseif RHc_dqc_scheme >= 16
            %/ NOTE: Scheme 16 uses a more accurate outline of the Pan-TP
            %/       together with 'show_land_or_ocean_hydrosheds'.
            %/       Also, the scheme is again fine-tuned.
            %/
            %/       Scheme 17 is the same as Scheme 16, but now we fine tune
            %/       the RHc over TP based on the known bias from Scheme 16

            dqc           = 0.1;
            thres_upbnd   = [ -20,   20,  inf];     %/ For prct diff (upper bound)
            thres_lowbnd  = [-inf,  -20,   20];     %/ For prct diff (lower bound)
            thres2_upbnd  = [  -1,    1,  inf];     %/ For abs  diff (upper bound)
            thres2_lowbnd = [-inf,   -1,    1];     %/ For abs  diff (lower bound)

            %/              ---> i (cols): prct diff
            RHc_levels_2D = [80    80  85;     % | 
                             80    85  90;     % | j (rows): abs diff
                             85    90  90;];   % v

            RHc_levels = unique(RHc_levels_2D)';   %/ RHc for underestimation, normal, overestimation.

            %/ RHc map for Pan-TP region (>1500m)
            thres3_TP_upbnd  = [ -60,  -40,  -20, inf]; %/ For prct  diff (upper bound)
            thres3_TP_lowbnd = [-inf,  -60,  -40, -20]; %/ For prct  diff (lower bound)
            thres4_TP_upbnd  = [-2.5,  -2, -1.5,   -1, inf]; %/ For abs  diff (upper bound)
            thres4_TP_lowbnd = [-inf, -2.5,   -2, -1.5,  -1]; %/ For abs  diff (lower bound)
            thres3_topo      = [1500];

            %/              ---> i (cols): prct diff
            RHc_TP        = [
                             65 70 75 80;     % | 
                             70 75 80 85;     % | 
                             70 75 80 85;     % | j (rows): abs diff
                             70 75 80 85;     % |
                             75 81 83 85;];   % v

            RHc_levels    = unique([RHc_levels, unique(RHc_TP)']);
            var           = 'P_LA_P_prct_diff';
            var2          = 'P_LA_P_diff';
            str_baseline        = '_BaselineScheme10';
            str_baseline_years  = '2010-2010';
            if RHc_dqc_scheme >= 20
                RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                   sprintf('RHc_map%s_%s_thres_%s_%s_thres2_%s_topo%d_pan_TP_scheme%d_TPtuned_seasonal_%s_global.mat',...
                                           str_baseline, var, strjoin(string(thres_upbnd),'_'),...
                                           var2, strjoin(string(thres2_upbnd),'_'), thres3_topo, RHc_dqc_scheme, str_baseline_years)];

            else
                RHc_map_matfilename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/',...
                                       sprintf('RHc_map%s_%s_%s_thres_%s_%s_thres2_%s_topo%d_pan_TP_scheme%d_%s.mat',...
                                               str_baseline, strjoin(string(RHc_levels), '_'), var, strjoin(string(thres_upbnd),'_'),...
                                               var2, strjoin(string(thres2_upbnd),'_'), thres3_topo, RHc_dqc_scheme, str_baseline_years)];
            end
            disp(RHc_map_matfilename)
            load(RHc_map_matfilename, 'RHc_map', 'RHc_map_lon', 'RHc_map_lat');
            RHc = RHc_map;                                  %/ Patially varying critical RH to indicate prcp occurrence [%]
            str_RHc_dqc = sprintf('RHc_dqc_scheme%d', RHc_dqc_scheme);
        else
            error('Not prescribed setting for ''RHc_dqc_scheme == %d''!', RHc_dqc_scheme) 
        end
    else
        error('code not ready for resolution (res) = %.2f!\n', output_res)
    end
    fprintf('Max RHc = %.2f\n', max(RHc, [], 'all'));
    fprintf('Min RHc = %.2f\n', min(RHc, [], 'all'));
end

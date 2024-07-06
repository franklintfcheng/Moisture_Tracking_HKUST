%% working code
addpath(genpath('/home/tfchengac/MyMatlabPkg'));
addpath(genpath('/home/tfchengac/MyMatlabFn'));
addpath(genpath('/home/tfchengac/Flexpart_Codes'));
expmnt = 'domfill_CERA_MPI';  optimal_rr = []; optimal_prcnt = []; fwd = [];
%>>>>>>>>>>>>>>>
%/ WaterSip setting
% year_list = 1971:2010;
year_list = 2010;  %/ for r37

% trajdirect = 'bwd';     from_hs = 0;              %/ backward (bwd) or forward (fwd)
% maxtraj_day = 20;                                 %/ max traj days
% optimal_rr  = 0.9;                                %/ if not [], then track recursively until rr >= optimal value

trajdirect = 'fwd';     from_hs = 1;                %/ backward (bwd) or forward (fwd)
maxtraj_day = 15;    %/15 days seems to be enough   %/ max traj days
optimal_prcnt  = 0.1;                               %/ if not [], then forward track the contribution until f < optimal value

% traj_rm_jump = 1; rm_when_dz = 5000;                %/ toggle on to slice traj of a particle with a *big* jump
traj_rm_jump = 0; rm_when_dz = [];                  %/ toggle on to slice traj of a particle with a *big* jump
dqc = 0.1;                                          %/ critical dq for prcp to occur [g/kg]
RHc = 85;                                           %/ critical RH for prcp to occur [%]
BLH_factor = 1;                                     %/ multiplier on BLH to account for small-scale variability
% str_remark = 'RH2_';                                %/ REMEMBER to change that in "slct_traj_parfor" function !!!
% str_remark = 'RH2_correct_';                        %/ REMEMBER to change that in "slct_traj_parfor" function !!!
str_remark = 'testing_';                              %/ REMEMBER to change that in "slct_traj_parfor" function !!!

%>>>>>>>>>>>>>>
overwrite_traj               = 1; %/ [0]: skip if file exists  [1]: overwrite file in any case
overwrite_watersip_products  = 2; %/ [0]: skip if file exists  [1]: overwrite file in any case  [2]: just to skip 
process_jan_only             = 0; %/ for testing use
matlab_id = 1;                    % <--- change this for different matlab processes
job_divi = 1;
NumWorkers = 20;  %20;                              %/ optimal no. of worker

for year = year_list

    %/ Main folder to save processed data (traj + WSV)
    if ismember(trajdirect, {'bwd'})
        if year >= 1991
            masterfolder = '/disk/r034/tfchengac/FLEXPART/';  %/ occupies 10 TB
        else
            masterfolder = '/disk/r059/tfchengac/FLEXPART/';  %/ occupies 10 TB
        end
    else if ismember(trajdirect, {'fwd'})
        masterfolder = '/disk/r034/tfchengac/FLEXPART/';   
    end
    end
    
    %/ FLEXPART output dir
    if year >= 2001 && year <= 2010
        flexpart_output_dir = strcat('/disk/r037/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
    else if year >= 1991 && year <= 2000
        flexpart_output_dir = strcat('/disk/r034/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
    else if year >= 1981 && year <= 1990
        flexpart_output_dir = strcat('/disk/r012/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
    else if year >= 1971 && year <= 1980
        flexpart_output_dir = strcat('/disk/r014/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
    end
    end
    end
    end

    %==== Always read the header file ====%
    nest = 0;
    readp = 1; %/ read release points (0/1)
    calcarea = 1;

    [header, fail] = flex_header(flexpart_output_dir, nest, readp, calcarea); %/ Postprocessing tool from flexpart
    header.dates_str = string(importdata([flexpart_output_dir 'dates'])); %/ read dates file in string

    header.dates_dt = datetime(header.dates_str,'InputFormat','yyyyMMddHHmmss', 'Format', 'yyyy-MM-dd HH:mm:ss');

    ind = find(ismember(string(header.dates), string(header.ireleasestart)));
    header.ireleasestart_dt = header.dates_dt(ind, :);

    ind = find(ismember(string(header.dates), string(header.ireleaseend)));
    header.ireleaseend_dt = header.dates_dt(ind, :);
    

    %==== Derive 10-day trajectories from domain-fill and the uptake map ====%
    dt_slct = 3;                                                           %/ traj time interval (in hr; user-defined)
    dt  = 3;                                                               %/ FLEXPART time interval (in hr)
    maxtraj_ts = 24*maxtraj_day/dt_slct + 1;                               %/ traj time steps, +1 to include the very last time step
    res = 1;
    lon = 0:res:359;
    lat = -90:res:90;
    area = calc_grid_area_header(header);                                  %/ my function to calculating gird area.
    
    if ismember(trajdirect, {'bwd'})
        trajtime_loop_direct = 0:-1:-(maxtraj_ts-1);                       %/ backward traj time loop
        trajtime             = trajtime_loop_direct*dt_slct;

    else if ismember(trajdirect, {'fwd'})
        trajtime_loop_direct = -1:1:maxtraj_ts-1;                          %/ forward traj time loop (starting from -1)
        trajtime             = trajtime_loop_direct(1:end)*dt_slct;
        
    end
    end

    if ismember(trajdirect, {'bwd'})
        if process_jan_only                                                %/ in the test mode, process only 2010 Jan data.
            slct_date_dt = header.dates_dt(1:496);
        else
            slct_date_dt = header.dates_dt(1:end);
        end
        dates_yyyymmdd = yyyymmdd(slct_date_dt);
        endofbuffer = max(find(dates_yyyymmdd == (year-1)*1e4 + 1231));    %/ get the time index of the end of the buffer period (i.e. yyyy-12-31 21:00)
        
        %/ split the job
        t_split = floor((length(slct_date_dt) - endofbuffer)/job_divi);
    
        if matlab_id == job_divi
            t_list = (1 + endofbuffer + t_split*(matlab_id-1)):length(slct_date_dt); %/ to workaround the uneven job split.
        else
            t_list = (1:t_split) + endofbuffer + t_split*(matlab_id-1);
        end
    end
    
    if ismember(trajdirect, {'fwd'})
        if process_jan_only                                                %/ in the test mode, process only 2010 Jan data.
            a = yyyymmdd(header.dates_dt);
            ind_st = min(find(a == year*1e4 + 101)) - 1;                              %/ from one timestep before Jan 1.
            ind_ed = max(find(a == year*1e4 + 131));                                  %/ up to Jan 31.
            slct_date_dt = header.dates_dt(ind_st:ind_ed+max(trajtime_loop_direct));  %/ extend the dates for maxtraj_day.
        else
            a = yyyymmdd(header.dates_dt);
            if year == 1971
                ind_st = min(find(a == year*1e4 + 101)) - 1;                          %/ from one timestep before Jan 1, 1971
            else
                ind_st = min(find(a == (year-1)*1e4 + 1216)) - 1;                     %/ from one timestep before Dec 16 in year(-1)
            end
            slct_date_dt = header.dates_dt(ind_st:end);                              
        end
        
        dates_yyyymmdd = yyyymmdd(slct_date_dt);
        if process_jan_only
            endofbuffer = length(slct_date_dt) - max(find(dates_yyyymmdd == year*1e4 + 131)); %/ the buffer is from Feb 1 to Feb xx (max. traj day)
        else
            endofbuffer = max(trajtime_loop_direct);                       %/ to track till 16 Dec (or till 15 Dec in 2010).
        end
        
        %/ split the job
        t_split = floor((length(slct_date_dt) - endofbuffer)/job_divi);
        
        if matlab_id == job_divi
            t_list = (1 + t_split*(matlab_id-1)):length(slct_date_dt)-endofbuffer; 
        else
            t_list = (1:t_split) + t_split*(matlab_id-1);
        end
        
        if matlab_id == 1     t_list(1) = [];    end   %/ since it starts from 1231 21:00, we skip it and start from the 2nd time step.
    end
    
    disp(slct_date_dt(t_list,:))

    prev_t = []; numpart_esti = 'unknown'; 
    if ismember(trajdirect, {'fwd'})
        for t = t_list(2)
            if t == length(slct_date_dt)
                fprintf('*** NOTE: The last time step of the experiment is skipped cos the output vars contain NaNs only. ***\n');
                continue;
            end

            domfill = []; %/ release memory
            fprintf('\n*** overwrite_traj   = %d *** \n',   overwrite_traj);
            fprintf('*** overwrite_watersip_products = %d *** \n', overwrite_watersip_products);
            if process_jan_only 
                fprintf('*** WARNING: test mode is on. Process only Jan data !! *** \n'); 
            end

            src_z   = [];  src_lon = [];  src_lat = [];
            if ~isempty(optimal_rr) || ~isempty(optimal_prcnt)  str_optimal      = '_optimal';                       else  str_optimal      = [];  end
            if traj_rm_jump == 0                                str_traj_rm_jump = 'intact_';                        else  str_traj_rm_jump = [];  end
            if BLH_factor ~= 1                                  str_BLH_factor   = sprintf('%.1fBLH_', BLH_factor);  else  str_BLH_factor   = [];  end
            if isempty(src_z)                                   str_src_z        = 'nozbound_';                      else  str_src_z        = [];  end

            %>>>>>>>>> Put a for-loop over hotspots here? <<<<<<<<<
            if from_hs && ~exist('AYR_hotspot', 'var')
                AYR_hs_data_filename = string(strcat('/disk/r059/tfchengac/FLEXPART/', expmnt, '/prcssd_data_4plotting/',...
                                              'AYR_hs_p95_area_sum_BL_Pm_RHc85_dqc0.1_bwd_20d_optimal_3h_glbland_nozbound_intact_RH2_1971-2010.mat'));
                fprintf('*** The queried data are found. *** \n*** Loading: %s *** \n', AYR_hs_data_filename)
                load(AYR_hs_data_filename);   %/ output is AYR_hotspot
            end

            %===== Set filepaths =====%
            str_domain = sprintf('AYR_hs%d', length(AYR_hotspot));

            traj_suffix = strcat('RHc', num2str(RHc),'_dqc', num2str(dqc), '_', trajdirect,'_',...
                            num2str(maxtraj_day), 'd', str_optimal, '_', num2str(dt_slct), 'h_', str_domain,'_', str_src_z, str_BLH_factor, str_remark,...
                            datestr(slct_date_dt(t,:), 'yyyymmddHHMM'), '.mat');

            suffix = strcat('RHc', num2str(RHc),'_dqc', num2str(dqc), '_', trajdirect,'_',...
                            num2str(maxtraj_day), 'd', str_optimal, '_', num2str(dt_slct), 'h_', str_domain,'_', str_traj_rm_jump, str_src_z, str_BLH_factor, str_remark,...
                            datestr(slct_date_dt(t,:), 'yyyymmddHHMM'), '.mat');

            %/ dir to store processed data
            if yyyymmdd(slct_date_dt(t,:)) < year*1e4     %/ then modify the output dir for year (-1)
                prcssd_dir = strcat(masterfolder, expmnt, '/prcssd_data_', num2str(year-1),'/');
            else
                prcssd_dir = strcat(masterfolder, expmnt, '/prcssd_data_', num2str(year),'/');
            end
            mkdir(char(prcssd_dir));
            
            traj_fullpath             = strcat(prcssd_dir,  'traj_',         traj_suffix);

            %/ for fwd mode
            contr_map_fullpath        = strcat(prcssd_dir,  'contr_map_',         suffix);
            optimal_trajtime_fullpath = strcat(prcssd_dir,  'optimal_trajtime_',  suffix);

            %===== Skip if all data exist and do not overwrite them =======%
            if isfile(traj_fullpath) && isfile(contr_map_fullpath) && isfile(optimal_trajtime_fullpath) &&...
                overwrite_traj == 0 && overwrite_watersip_products == 0

                fprintf('!!! [t = %d/%d] All traj, uptake and watersip files exist, and no overwriting is needed. So skip this loop. !!! \n', t, t_list(end))
                continue;
            end

            %===== Check if traj data exists, load from it if overwrite = 0 =======%
            traj_exist = 0;
            if overwrite_traj == 0
                if isfile(traj_fullpath)
                    fprintf('!!! [t = %d/%d] traj & EminusP_LA data are found: %s. Loading traj data... !!! \n', t, t_list(end), traj_fullpath)
                    load(traj_fullpath);
                    traj_exist = 1;
                end
            end

            if traj_exist == 0
                tic
                fprintf('*** t = %d, prev_t = %d ***\n', t, prev_t);
                if (isempty(prev_t)) || ( t - prev_t ~= 1)                 %/ only two scenarios when direct loading is needed.
                    plist = struct('partoutput', cell(length(trajtime_loop_direct), 1));
                end
                
                fwd = {};
                trajtime_loop_direct_init = trajtime_loop_direct(1:2);
                for cnt = 1:length(trajtime_loop_direct_init)              %/ run the first 2 timesteps
                    j = trajtime_loop_direct_init(cnt);

                    date_flag = datestr(slct_date_dt(t+j,:), 'yyyymmddHHMMSS');
                    
                    if (isempty(prev_t)) || ( t - prev_t ~= 1)             %/ load the traj if it's the first time of the year, OR if the difference in two timesteps ~= 1 (e.g. across two years, or some timesteps were skipped)
                        if ismember(numpart_esti, 'unknown')
                            numpart_esti = readpart10(flexpart_output_dir, date_flag, numpart_esti); %/ If numpart is unknown, input 'unknown' to estimate it.
                        end

                        plist(cnt).partoutput = readpart10(flexpart_output_dir, date_flag, numpart_esti);           
                        str_load_or_retrieve = 'Direct Loading:';
                    elseif cnt == 1
                        [plist(1:end-1).partoutput] = plist(2:end).partoutput;   %/ shift the previous plist one timestep forward.
                        str_load_or_retrieve = 'Retrieve From Memory:';
                    end

                    fprintf('*** [t = %d/%d] [%.1f days %s from %s]: %s %d trajs on %s (%d/%d) are loaded. *** \n',...
                                t, t_list(end), maxtraj_day, trajdirect, slct_date_dt(t,:), str_load_or_retrieve, numpart_esti, slct_date_dt(t+j,:), cnt, length(trajtime_loop_direct))
                            
                    %/ loop over hotspots
                    for top = 1:length(AYR_hotspot)
                        fprintf('*** AYR hotspot #%d/%d ... ***\n', top, length(AYR_hotspot))
                        bndry_data = AYR_hotspot(top).bndry;   %<<-- important!

                        if cnt == 1
                            %/ preallocate struct field for EACH hotspot!
                            domfill_parfor = struct('slct_id_list', cell(1,1),...
                                                     'traj',cell(1,length(trajtime_loop_direct)),...
                                                     'q',   cell(1,length(trajtime_loop_direct)),...
                                                     'BLH', cell(1,length(trajtime_loop_direct)),...
                                                     'dates_dt', cell(1,length(trajtime_loop_direct)));
                        else
                            domfill_parfor = fwd{top};                     %/retrieve the struct of the hotspot in the cell
                        end

                        %======================== allocation ======================
                        %/ store q and id one timestep before the start of forward traj.
                        id = plist(cnt).partoutput.npoint;
                        domfill_parfor(cnt).traj = plist(cnt).partoutput.xyz;
                        domfill_parfor(cnt).q    = plist(cnt).partoutput.vars(:,1) * 1000; %/ from kg/kg to g/kg.
                        domfill_parfor(cnt).BLH  = plist(cnt).partoutput.vars(:,2);
                        domfill_parfor(cnt).T    = plist(cnt).partoutput.vars(:,3);
                        domfill_parfor(cnt).dates_dt = slct_date_dt(t+j);

                        if cnt == 2
                            if ~isempty(bndry_data(:,1) < 0)               %/ convert to [0, 360) if not so.
                                ind = find(bndry_data(:,1) < 0);
                                bndry_data(ind,1) = bndry_data(ind,1) + 360;  
                            end 

                            traj_x = plist(cnt).partoutput.xyz(:,1);
                            traj_y = plist(cnt).partoutput.xyz(:,2);
                            if ~isempty(traj_x < 0)                        %/ convert to [0, 360) if not so.
                                ind = find(traj_x < 0);
                                traj_x(ind) = traj_x(ind) + 360;  
                            end

                            %/ First, select the traj within the hotspot's boundary
                            disp('*** Matching trajs within the hotspot boundary (The length of bndry_data determines the speed! Not the # of particles!)... ***')
                            [in, ~] = inpoly2([traj_x, traj_y], bndry_data);    %/ inpoly2 is 600xx faster than inpolygon!!
%                             tic;   [in,~] = inpolygon(traj_x, traj_y, bndry_data(:,1), bndry_data(:,2));  toc
                            ind = find(in == 1);

                            domfill_parfor(1).slct_id_list = id(ind);
                            domfill_parfor(1).mass     = plist(cnt).partoutput.xmass(ind,1); 

                            for k = 1:2
                                domfill_parfor(k).traj   = domfill_parfor(k).traj(ind,:);
                                domfill_parfor(k).q      = domfill_parfor(k).q(ind);
                                domfill_parfor(k).BLH    = domfill_parfor(k).BLH(ind);
                                domfill_parfor(k).T      = domfill_parfor(k).T(ind);
                            end

                            [domfill_parfor, ~] = slct_traj_parfor('trajdirect', trajdirect, 'lon', lon, 'lat', lat, 'domfill', domfill_parfor,...
                                                                   'RHc', RHc, 'dqc', dqc, 'BLH_factor', BLH_factor, 'area', area); 
                            ind = []; %/ release memory
                        end

                        fwd{top} = domfill_parfor;                         %/update the struct in the cell
                    end
%                     partoutput = []; %/ release memory 
                end
                toc

                %/ [IMPORTANT] replicate the initial domfill_parfor for each of the rest trajtime!!
                for k = cnt+1:length(trajtime_loop_direct)
                    fwd(k,:) = fwd(1,:);
                end

                %======== Parfor loop over the rest traj time !!! ============%
                if isempty(gcp('nocreate')) && ~isempty(NumWorkers)        %/ if set worker number
                    parpool('local', NumWorkers)                               %/ use process-based parpool (threads-based parpool is only available after R2020a :((
                end
%%
                tic
                for k = cnt+1:length(trajtime_loop_direct)
%                 parfor k = cnt+1:length(trajtime_loop_direct)              %/ parfor for the rest (starting from cnt+1)
%                     plist_bc = plist;
                    
                    j = trajtime_loop_direct(k);

                    date_flag = datestr(slct_date_dt(t+j,:), 'yyyymmddHHMMSS');
                    
                    if ( (isempty(prev_t)) || ( t - prev_t ~= 1) ) || ...
                       ( (~isempty(prev_t)) && ( t - prev_t == 1) && k == length(trajtime_loop_direct) )          %/ only two scenarios when direct loading is needed.
                   
                        plist(k).partoutput = readpart10(flexpart_output_dir, date_flag, numpart_esti);   
                        str_load_or_retrieve = 'Direct Loading:';
                    else
                        str_load_or_retrieve = 'Retrieve From Memory:';
                    end
                    partoutput = plist(k).partoutput;    %/ then we load the correpsonding traj at timestep k
                    
                    fprintf('*** [t = %d/%d] [%.1f days %s from %s]: %s %d trajs on %s (%d/%d) are loaded. *** \n',...
                                t, t_list(length(t_list)), maxtraj_day, trajdirect, slct_date_dt(t,:), str_load_or_retrieve, numpart_esti, slct_date_dt(t+j,:), k, length(trajtime_loop_direct))

                    %======================== allocation ======================
                    a = fwd(k, :);      %/ [IMPORTANT] use a temporary variable to avoid classification problem in an inner loop of a parfor!
                    
                    for top = 1:length(AYR_hotspot)
                        
                        domfill_parfor = a{top};
                        slct_id_list = domfill_parfor(1).slct_id_list;

                        ind = findismember_loop([partoutput.npoint], slct_id_list);  %/ get indices based on id list (no auto sorting!!)

                        domfill_parfor(k).traj    = partoutput.xyz(ind, :);  
                        domfill_parfor(k).q       = partoutput.vars(ind,1) * 1000; %/ from kg/kg to g/kg.
                        domfill_parfor(k).BLH     = partoutput.vars(ind,2);
                        domfill_parfor(k).T       = partoutput.vars(ind,3);
                        domfill_parfor(k).dates_dt = slct_date_dt(t+j);

                        a{top} = domfill_parfor; %/update the struct in the cell
                    end
                    
                    partoutput = [];                                           %/ release memory 
                    fwd(k, :) = a; 
                end
                toc
%%
                %/ resemble the processed struct data into one, for each hotspot
                tic
                fprintf('*** resemble struct data (domfill_parfor) ... ***\n')
                fwd_final = {};
                for top = 1:length(AYR_hotspot)
                    for k = cnt+1:length(trajtime_loop_direct) 
                        domfill_parfor   = fwd{1, top};
                        domfill_parfor_k = fwd{k, top};

                        fld = fieldnames(domfill_parfor);
                        for i = 1:length(fld)
                            domfill_parfor(k).(fld{i}) = domfill_parfor_k(k).(fld{i});
                        end

                        fwd{1, top} = domfill_parfor;                          %/ update it 
                    end

                    domfill_parfor   = fwd{1, top};                            %/ now use this complete domfill_parfor.

                    domfill.traj     = cat(3, domfill_parfor(:).traj);         %/ change nonscalar domfill_parfor to scalar domfill (compatible for other codes) 
                    domfill.q        = cat(2, domfill_parfor(:).q);            %/ and concatenate the numbered fields 
                    domfill.BLH      = cat(2, domfill_parfor(:).BLH);
                    domfill.dates_dt = cat(2, domfill_parfor(:).dates_dt);
                    domfill.T        = cat(2, domfill_parfor(:).T);
%                     domfill.T        = domfill_parfor(1).T;
                    domfill.mass     = domfill_parfor(1).mass;

        %             domfill_parfor = [];   %/ release memory

                    domfill.trajdirect   = trajdirect;
                    domfill.src_lon      = src_lon;
                    domfill.src_lat      = src_lat;
                    domfill.src_z        = src_z;
                    domfill.dt_slct      = dt_slct;
                    domfill.maxtraj_ts   = maxtraj_ts;
                    domfill.maxtraj_day  = maxtraj_day;
                    domfill.trajtime     = trajtime;

                    fwd_final{top} = domfill;
                end
                toc

                %===== Save traj data =====%
                fprintf('*** [t = %d/%d] Saving traj into: %s ... *** \n', t, t_list(end), traj_fullpath)
                save(traj_fullpath, 'fwd_final', '-v7.3');                 %/ save/overwrite
                
                prev_t = t;                                                %/ update prev_t
            end

            %===== WaterSip algorithm =====%
            if overwrite_watersip_products == 2
                fprintf('!!! [t = %d/%d] uptake and watersip data are skipped as per user''s request. !!!', t, t_list(end))
            else
                
                %===== Check if the uptake_map data exists, load from it if overwrite = 0 =======%
                contr_map_exist = 0; watersip_fwd_exist = 0; optimal_trajtime_exist = 0; 
                if overwrite_watersip_products == 0
                    %/ fwd mode
                    if isfile(contr_map_fullpath)         contr_map_exist = 1;          end
                    if isfile(optimal_trajtime_fullpath)  optimal_trajtime_exist = 1;   end
                end

                %/ Either overwrite is required or any watersip files not exists, we will run the algorithm.
                if overwrite_watersip_products == 1 || (contr_map_exist == 0 || watersip_fwd_exist == 0 || optimal_trajtime_exist == 0)
                    tic
                    contr_map = {}; mean_optimal_trajtime_map = {};
                    for top = 1:length(AYR_hotspot)
                        fprintf('*** AYR hotspot #%d/%d ... ***\n', top, length(AYR_hotspot))
                        domfill = fwd_final{top};
                        
                        fprintf('*** [t = %d/%d] Start WaterSip (Forward) algorithm... *** \n', t, t_list(end))
                        [contr_map{top}, watersip_fwd_stacktable, mean_optimal_trajtime_map{top}] ...
                                = WaterSip_fwd('domfill', domfill, 'lon', lon, 'lat', lat, 'dqc', dqc, 'RHc', RHc, 'BLH_factor', BLH_factor, 'area', area,...
                                               'NumWorkers', NumWorkers, 'optimal_prcnt', optimal_prcnt);
                    end
                    toc
                    
                    if contr_map_exist == 0             save(contr_map_fullpath,           'contr_map',                 '-v7.3'); end
                    if optimal_trajtime_exist == 0      save(optimal_trajtime_fullpath,    'mean_optimal_trajtime_map', '-v7.3'); end
                    fprintf('*** [t = %d/%d] WaterSip products saved under %s *** \n', t, t_list(end), prcssd_dir)
                end
            end
        end
    end

    if ismember(trajdirect, {'bwd'}) 
        for t = t_list
            if t == length(slct_date_dt)
                fprintf('*** NOTE: The last time step of the experiment is skipped cos the output vars contain NaNs only. ***\n');
                continue;
            end

            domfill = []; %/ release memory
            fprintf('\n*** overwrite_traj   = %d *** \n',   overwrite_traj);
            fprintf('*** overwrite_watersip_products = %d *** \n', overwrite_watersip_products);
            if process_jan_only 
                fprintf('*** WARNING: test mode is on. Process only Jan data !! *** \n'); 
            end

            src_z   = [];  src_lon = [];  src_lat = [];
            if ~isempty(optimal_rr) || ~isempty(optimal_prcnt)  str_optimal      = '_optimal';                       else  str_optimal      = [];  end
            if traj_rm_jump == 0                                str_traj_rm_jump = 'intact_';                        else  str_traj_rm_jump = [];  end
            if BLH_factor ~= 1                                  str_BLH_factor   = sprintf('%.1fBLH_', BLH_factor);  else  str_BLH_factor   = [];  end
            if isempty(src_z)                                   str_src_z        = 'nozbound_';                      else  str_src_z        = [];  end

            %>>>>>>>>> Put a for-loop over hotspots here? <<<<<<<<<
            domain = 3;
            if domain == 1
                %/ Regional domain
                src_lon = [85  90];                
                src_lat = [10  15];
                src_z   = [0  9000]; %/ 9000 --> 300 hPa
                str_domain = strcat('lon', join(string(src_lon), '_'), ...
                                    '_lat', join(string(src_lat), '_'),...
                                    '_z',   join(string(src_z),   '_'));
            else if domain == 2
                %/ Global domain
                src_lon = [0  360];             
                src_lat = [-90  90];
                src_z   = [0  9000]; %/ 9000 --> 300 hPa
                str_domain = strcat('lon', join(string(src_lon), '_'), ...
                                    '_lat', join(string(src_lat), '_'),...
                                    '_z',   join(string(src_z),   '_'));
            else if domain == 3
                %/ global land region
                src_lon = 'glbland';
                src_lat = 'glbland';
                src_z   = []; 
                str_domain = 'glbland';
            end
            end
            end

            traj_suffix = strcat('RHc', num2str(RHc),'_dqc', num2str(dqc), '_', trajdirect,'_',...
                            num2str(maxtraj_day), 'd', str_optimal, '_', num2str(dt_slct), 'h_', str_domain,'_', str_src_z, str_BLH_factor, str_remark,...
                            datestr(slct_date_dt(t,:), 'yyyymmddHHMM'), '.mat');

            suffix = strcat('RHc', num2str(RHc),'_dqc', num2str(dqc), '_', trajdirect,'_',...
                            num2str(maxtraj_day), 'd', str_optimal, '_', num2str(dt_slct), 'h_', str_domain,'_', str_traj_rm_jump, str_src_z, str_BLH_factor, str_remark,...
                            datestr(slct_date_dt(t,:), 'yyyymmddHHMM'), '.mat');

            prcssd_dir = strcat(masterfolder, expmnt, '/prcssd_data_', num2str(year),'/');
                        
            traj_fullpath             = strcat(prcssd_dir,  'traj_',         traj_suffix);
            %/ for bwd mode
            EminusP_LA_fullpath       = strcat(prcssd_dir,  'EminusP_LA_',        suffix);
            uptake_fullpath           = strcat(prcssd_dir,  'uptake_',            suffix);
            Pm_fullpath               = strcat(prcssd_dir,  'Pm_',                suffix);
            BL_uptake_fullpath        = strcat(prcssd_dir,  'BL_uptake_',         suffix);
            BL_Pm_fullpath            = strcat(prcssd_dir,  'BL_Pm_',             suffix);
            watersip_rr_fullpath      = strcat(prcssd_dir,  'watersip_',          suffix);
            optimal_trajtime_fullpath = strcat(prcssd_dir,  'optimal_trajtime_',  suffix);
            P_LA_fullpath             = strcat(prcssd_dir,  'P_LA_',              suffix);

            %/ for fwd mode
            contr_map_fullpath        = strcat(prcssd_dir,  'contr_map_',         suffix);
            watersip_fwd_fullpath     = strcat(prcssd_dir,  'watersip_fwd_',      suffix);

            %===== Skip if all data exist and do not overwrite them =======%
            if ismember(trajdirect, {'bwd'}) 
                if isfile(traj_fullpath) && isfile(EminusP_LA_fullpath) && isfile(uptake_fullpath) && isfile(Pm_fullpath)...
                   && isfile(BL_uptake_fullpath) && isfile(BL_Pm_fullpath)...
                   && isfile(watersip_rr_fullpath) && isfile(optimal_trajtime_fullpath) && isfile(P_LA_fullpath) ...
                   && overwrite_traj == 0 && overwrite_watersip_products == 0

                    fprintf('!!! [t = %d/%d] All traj, uptake and watersip files exist, and no overwriting is needed. So skip this loop. !!! \n', t, t_list(end))
                    continue;
                end
            end

            if ismember(trajdirect, {'fwd'}) 
                if isfile(traj_fullpath) && isfile(contr_map_fullpath) && isfile(watersip_fwd_fullpath) && isfile(optimal_trajtime_fullpath) &&...
                    overwrite_traj == 0 && overwrite_watersip_products == 0

                    fprintf('!!! [t = %d/%d] All traj, uptake and watersip files exist, and no overwriting is needed. So skip this loop. !!! \n', t, t_list(end))
                    continue;
                end
            end

            %===== Check if traj data exists, load from it if overwrite = 0 =======%
            traj_exist = 0;
            if overwrite_traj == 0
                if (isfile(traj_fullpath) && isfile(EminusP_LA_fullpath) && ismember(trajdirect, {'bwd'})) || ...
                   (isfile(traj_fullpath) && ismember(trajdirect, {'fwd'}))

                    fprintf('!!! [t = %d/%d] traj & EminusP_LA data are found: %s. Loading traj data... !!! \n', t, t_list(end), traj_fullpath)
                    load(traj_fullpath);
                    traj_exist = 1;
                end
            end

            %===== Derive Traj =====%
            if traj_exist == 0

                partoutput = []; %/ release memory
                numpart_esti = 'unknown';

                %/ preallocate struct field for parfor loop later.
                domfill_parfor = struct('slct_id_list', cell(1,1),...
                                         'traj',cell(1,length(trajtime_loop_direct)),...
                                         'q',   cell(1,length(trajtime_loop_direct)),...
                                         'BLH', cell(1,length(trajtime_loop_direct)),...
                                         'dates_dt', cell(1,length(trajtime_loop_direct)));

                trajtime_loop_direct_init = trajtime_loop_direct(1:2);
                for cnt = 1:length(trajtime_loop_direct_init)                  %/ run the first 2 timesteps
                    j = trajtime_loop_direct_init(cnt);

                    fprintf('*** Loading traj ... ***\n')
                    disp(numpart_esti)
                    date_flag = datestr(slct_date_dt(t+j,:), 'yyyymmddHHMMSS');
                    if ismember(numpart_esti, 'unknown')
                        numpart_esti = readpart10(flexpart_output_dir, date_flag, numpart_esti); %/ If numpart is unknown, input 'unknown' to estimate it.
                    end

                    partoutput = readpart10(flexpart_output_dir, date_flag, numpart_esti);           
                    fprintf('*** [t = %d/%d] [%.1f days %s from %s]: %d trajs on %s (%d/%d) are loaded. *** \n',...
                            t, t_list(end), maxtraj_day, trajdirect, slct_date_dt(t,:), numpart_esti, slct_date_dt(t+j,:), cnt, length(trajtime_loop_direct))

                    %======================== allocation ======================
                    if ismember(trajdirect, {'fwd'}) %/ i.e., store q and id one timestep before the start of forward traj.
                        id = partoutput.npoint;
                        domfill_parfor(cnt).traj = partoutput.xyz;
                        domfill_parfor(cnt).q    = partoutput.vars(:,1) * 1000; %/ from kg/kg to g/kg.
                        domfill_parfor(cnt).BLH  = partoutput.vars(:,2);
                        domfill_parfor(cnt).T    = partoutput.vars(:,3);
                        domfill_parfor(cnt).dates_dt = slct_date_dt(t+j);

                        if cnt == 2
                            if from_hs
                                if ~isempty(bndry_data(:,1) < 0)               %/ convert to [0, 360) if not so.
                                    ind = find(bndry_data(:,1) < 0);
                                    bndry_data(ind,1) = bndry_data(ind,1) + 360;  
                                end 

                                traj_x = partoutput.xyz(:,1);
                                traj_y = partoutput.xyz(:,2);
                                if ~isempty(traj_x < 0)                        %/ convert to [0, 360) if not so.
                                    ind = find(traj_x < 0);
                                    traj_x(ind) = traj_x(ind) + 360;  
                                end

                                %/ First, select the traj within the hotspot's boundary
                                disp('*** Matching trajs within the hotspot boundary ... ***')
                                tic;   [in, ~] = inpoly2([traj_x, traj_y], bndry_data);  toc     %/ inpoly2 is 600xx faster than inpolygon!!
    %                             tic;   [in,~] = inpolygon(traj_x, traj_y, bndry_data(:,1), bndry_data(:,2));  toc
                                ind = find(in == 1);

                            else
                                if domain == 3
                                    ind = which_on_land('pos_lon_array', partoutput.xyz(:,1), 'pos_lat_array', partoutput.xyz(:,2),...
                                                        'pos_hgt_array', partoutput.xyz(:,3), 'hgt_range', src_z);

                                else if domain == 1 || domain == 2
                                    ind = find(partoutput.xyz(:,1) >= src_lon(1) & partoutput.xyz(:,1) <= src_lon(2) & ...
                                               partoutput.xyz(:,2) >= src_lat(1) & partoutput.xyz(:,2) <= src_lat(2) & ...
                                               partoutput.xyz(:,3) >= src_z(1)   & partoutput.xyz(:,3) <= src_z(2));
                                end
                                end
                            end
                            domfill_parfor(1).slct_id_list = id(ind);
                            domfill_parfor(1).mass     = partoutput.xmass(ind,1); 

                            for k = 1:2
                                domfill_parfor(k).traj   = domfill_parfor(k).traj(ind,:);
                                domfill_parfor(k).q      = domfill_parfor(k).q(ind);
                                domfill_parfor(k).BLH    = domfill_parfor(k).BLH(ind);
                                domfill_parfor(k).T      = domfill_parfor(k).T(ind);
                            end

                            [domfill_parfor, ~] = slct_traj_parfor('trajdirect', trajdirect, 'lon', lon, 'lat', lat, 'domfill', domfill_parfor,...
                                                                   'RHc', RHc, 'dqc', dqc, 'BLH_factor', BLH_factor, 'area', area); 
                            ind = []; %/ release memory
                        end
                    end

                    if ismember(trajdirect, {'bwd'}) 
                        if cnt == 1
                            %/ Selection 1: handle only those particles from the land grids
                            if domain == 3

                                ind = which_on_land('pos_lon_array', partoutput.xyz(:,1), 'pos_lat_array', partoutput.xyz(:,2),...
                                                    'pos_hgt_array', partoutput.xyz(:,3), 'hgt_range', src_z);

                            else if domain == 1 || domain == 2
                                ind = find(partoutput.xyz(:,1) >= src_lon(1) & partoutput.xyz(:,1) <= src_lon(2) & ...
                                           partoutput.xyz(:,2) >= src_lat(1) & partoutput.xyz(:,2) <= src_lat(2) & ...
                                           partoutput.xyz(:,3) >= src_z(1)   & partoutput.xyz(:,3) <= src_z(2));
                            end
                            end
                            domfill_parfor(1).slct_id_list = partoutput.npoint(ind);

                            %/ Preallocation (will help the codes run faster)
                            domfill_parfor(cnt).traj  = partoutput.xyz(ind, :);
                            domfill_parfor(cnt).q     = partoutput.vars(ind, 1) * 1000;
                            domfill_parfor(cnt).BLH   = partoutput.vars(ind, 2);
                            domfill_parfor(cnt).T     = partoutput.vars(ind, 3);
                            domfill_parfor(1).mass    = partoutput.xmass(ind,1); 
                            domfill_parfor(cnt).dates_dt = slct_date_dt(t+j);

                        else if cnt == 2 
                            ind = findismember_loop(partoutput.npoint, domfill_parfor(1).slct_id_list); %/ get indices based on id list (no auto sorting)
                            domfill_parfor(cnt).traj  = partoutput.xyz(ind, :);
                            domfill_parfor(cnt).q     = partoutput.vars(ind, 1)*1000;
                            domfill_parfor(cnt).BLH   = partoutput.vars(ind, 2);
                            domfill_parfor(cnt).T     = partoutput.vars(ind, 3);
                            domfill_parfor(cnt).dates_dt = slct_date_dt(t+j);
                            %/ Selection 2: only handle trajs ended up with RH > RHc and dq > dqc
                            %/    --> update domfill.slct_id_list (speedup by > 1 s for each loop!)

                            [domfill_parfor, EminusP_LA_map] = slct_traj_parfor('trajdirect', trajdirect, 'lon', lon, 'lat', lat, 'domfill', domfill_parfor,...
                                                              'RHc', RHc, 'dqc', dqc, 'BLH_factor', BLH_factor, 'area', area); 

                            %===== Save EminusP_LA data =====%
                            fprintf('*** Saving EminusP_LA into: %s ... *** \n', EminusP_LA_fullpath)
                            save(EminusP_LA_fullpath, 'EminusP_LA_map', '-v7.3');             %/ save/overwrite               
                        end
                        end
                    end

                    partoutput = []; %/ release memory 
                end

                %======== Parfor loop over the rest traj time !!! ============%
                if isempty(gcp('nocreate')) && ~isempty(NumWorkers)            %/ if set worker number
                    parpool('local', NumWorkers)                               %/ use process-based parpool (threads-based parpool is only available after R2020a :((
                end

                slct_id_list = domfill_parfor(1).slct_id_list;
                tic
    %             for k = cnt+1:length(trajtime_loop_direct)                   %/ testing
                parfor k = cnt+1:length(trajtime_loop_direct)                  %/ parfor for the rest (starting from cnt+1)
                    j = trajtime_loop_direct(k);

                    fprintf('*** Loading traj ... ***\n')
                    date_flag = datestr(slct_date_dt(t+j,:), 'yyyymmddHHMMSS');

                    partoutput = readpart10(flexpart_output_dir, date_flag, numpart_esti);            
                    fprintf('*** [t = %d/%d] [%.1f days %s from %s]: %d trajs on %s (%d/%d) are loaded. *** \n',...
                            t, t_list(end), maxtraj_day, trajdirect, slct_date_dt(t,:), numpart_esti, slct_date_dt(t+j,:), k, length(trajtime_loop_direct))

                    %======================== allocation ======================
                    ind = findismember_loop(partoutput.npoint, slct_id_list);  %/ get indices based on id list (no auto sorting)

                    domfill_parfor(k).traj    = partoutput.xyz(ind, :);  
                    domfill_parfor(k).q       = partoutput.vars(ind,1) * 1000; %/ from kg/kg to g/kg.
                    domfill_parfor(k).BLH     = partoutput.vars(ind,2);

                    domfill_parfor(k).dates_dt = slct_date_dt(t+j);
                    partoutput = [];                                           %/ release memory 
                end
                toc

                domfill.traj     = cat(3, domfill_parfor(:).traj);             %/ change nonscalar domfill_parfor to scalar domfill (compatible for other codes) 
                domfill.q        = cat(2, domfill_parfor(:).q);                %/ and concatenate the numbered fields 
                domfill.BLH      = cat(2, domfill_parfor(:).BLH);
                domfill.dates_dt = cat(2, domfill_parfor(:).dates_dt);
                domfill.T        = domfill_parfor(1).T;
                domfill.mass     = domfill_parfor(1).mass;

                domfill_parfor = [];   %/ release memory

                domfill.trajdirect      = trajdirect;
                domfill.src_lon      = src_lon;
                domfill.src_lat      = src_lat;
                domfill.src_z        = src_z;
                domfill.dt_slct      = dt_slct;
                domfill.maxtraj_ts   = maxtraj_ts;
                domfill.maxtraj_day  = maxtraj_day;
                domfill.trajtime     = trajtime;

                %===== Save traj data =====%
                fprintf('*** [t = %d/%d] Saving traj into: %s ... *** \n', t, t_list(end), traj_fullpath)
                save(traj_fullpath, 'domfill', '-v7.3');             %/ save/overwrite
            end

            %===== WaterSip algorithm =====%
            if overwrite_watersip_products == 2
                fprintf('!!! [t = %d/%d] uptake and watersip data are skipped as per user''s request. !!!', t, t_list(end))
            else
                %===== Check if the uptake_map data exists, load from it if overwrite = 0 =======%
                uptake_exist = 0; Pm_exist = 0; BL_uptake_exist = 0; BL_Pm_exist = 0; watersip_exist = 0; optimal_trajtime_exist = 0; P_LA_exist = 0;
                contr_map_exist = 0; watersip_fwd_exist = 0;
                if overwrite_watersip_products == 0
                    %/ bwd mode
                    if isfile(uptake_fullpath)            uptake_exist = 1;             end
                    if isfile(Pm_fullpath)                Pm_exist = 1;                 end
                    if isfile(BL_uptake_fullpath)         BL_uptake_exist = 1;          end
                    if isfile(BL_Pm_fullpath)             BL_Pm_exist = 1;              end
                    if isfile(watersip_rr_fullpath)       watersip_exist = 1;           end
                    if isfile(optimal_trajtime_fullpath)  optimal_trajtime_exist = 1;   end
                    if isfile(P_LA_fullpath)              P_LA_exist = 1;               end

                    %/ fwd mode
                    if isfile(contr_map_fullpath)         contr_map_exist = 1;          end
                    if isfile(watersip_fwd_fullpath)      watersip_fwd_exist = 1;       end
                end

                %/ Either overwrite is required or any watersip files not exists, we will run the algorithm.
                tic
                if ismember(trajdirect, {'bwd'}) 
                    if overwrite_watersip_products || ...
                       (uptake_exist == 0 || Pm_exist == 0 || BL_uptake_exist == 0 || BL_Pm_exist == 0 || watersip_exist == 0 ...
                        || optimal_trajtime_exist == 0 || P_LA_exist == 0)

                        fprintf('*** [t = %d/%d] Start WaterSip algorithm... *** \n', t, t_list(end))
                        [uptake_map, Pm_map, BL_uptake_map, BL_Pm_map, watersip_rr, watersip_stacktable, mean_optimal_trajtime_map, P_LA_map] ...
                                = WaterSip('domfill', domfill, 'lon', lon, 'lat', lat, 'dqc', dqc, 'BLH_factor', BLH_factor, 'area', area,...
                                           'traj_rm_jump', traj_rm_jump, 'rm_when_dz', rm_when_dz, 'NumWorkers', NumWorkers, 'optimal_rr', optimal_rr);

                        if uptake_exist == 0             save(uptake_fullpath,             'uptake_map',                '-v7.3'); end
                        if Pm_exist == 0                 save(Pm_fullpath,                 'Pm_map',                    '-v7.3'); end
                        if BL_uptake_exist == 0          save(BL_uptake_fullpath,          'BL_uptake_map',             '-v7.3'); end
                        if BL_Pm_exist == 0              save(BL_Pm_fullpath,              'BL_Pm_map',                 '-v7.3'); end
                        if watersip_exist == 0           save(watersip_rr_fullpath,        'watersip_rr',               '-v7.3'); end
                        if optimal_trajtime_exist == 0   save(optimal_trajtime_fullpath,   'mean_optimal_trajtime_map', '-v7.3'); end
                        if P_LA_exist == 0               save(P_LA_fullpath,               'P_LA_map',                  '-v7.3'); end
                        fprintf('*** [t = %d/%d] WaterSip products saved under %s *** \n', t, t_list(end), prcssd_dir)
                    end
                end

                if ismember(trajdirect, {'fwd'}) 
                    if overwrite_watersip_products || (contr_map_exist == 0 || watersip_fwd_exist == 0 || optimal_trajtime_exist == 0)

                        fprintf('*** [t = %d/%d] Start WaterSip (Forward) algorithm... *** \n', t, t_list(end))
                        [contr_map, watersip_fwd, watersip_fwd_stacktable, mean_optimal_trajtime_map] ...
                                = WaterSip_fwd('domfill', domfill, 'lon', lon, 'lat', lat, 'dqc', dqc, 'BLH_factor', BLH_factor, 'area', area,...
                                           'traj_rm_jump', traj_rm_jump, 'rm_when_dz', rm_when_dz, 'NumWorkers', NumWorkers, 'optimal_prcnt', optimal_prcnt);

                        if contr_map_exist == 0             save(contr_map_fullpath,           'contr_map',                 '-v7.3'); end
                        if watersip_fwd_exist == 0          save(watersip_fwd_fullpath,        'watersip_fwd',              '-v7.3'); end
                        if optimal_trajtime_exist == 0      save(optimal_trajtime_fullpath,    'mean_optimal_trajtime_map', '-v7.3'); end
                        fprintf('*** [t = %d/%d] WaterSip products saved under %s *** \n', t, t_list(end), prcssd_dir)
                    end
                end
                toc
            end
            disp('stop here')
        end
    end
end
fprintf('\n !!!!!!!!!!!!!!!!!!!!!!')
fprintf('\n !!! Job completed. !!!\n')
fprintf(' !!!!!!!!!!!!!!!!!!!!!!\n')

%% Distribution of recylcing ratio 
% close all
% figname = '/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/Plottings/fixed_vs_optimal';
% plot_hist_cdf('x', 0:0.05:1, 'y', watersip_rr(:,5),  'y2', watersip_rr_optimal(:,5),...
%               'yname', '10-d', 'y2name', 'optimal till 20-d', 'y_lim', [0, 12000], 'savepath', figname);

%%
% a_cell = cell(1000, 1);
% a = ones(1, 360).*ones(181,1);
% 
% for i = 1:length(a_cell)
%     disp(i)
%     a_cell{i} = a;
% end
% b = cell2mat(a_cell);
% b = sum(b,1);
% 
% 
% b = reshape(b, size(a, 1), [])';
% 
% %%
% b = sum(cat(3, a_cell{:}),3);

%% compare 6h and 3h accuracy
% close all;
% masterfolder = '/disk/r059/tfchengac/FLEXPART/';
% 
% x = 0:0.05:1;
% x_label = 'Recycling ratio';
% 
% rr_3h = cellfun(@(x) sum(x(1,5:6)), watersip, 'UniformOutput', false);
% rr_3h = [rr_3h{:}]';
% 
% titlename = 'hist_rr_3h';
% savepath = strcat(masterfolder, expmnt, '/Plottings/', strrep(titlename, ' ', '_'));
% plot_hist_cdf('x', x, 'y', rr_3h, 'x_label', x_label, 'titlename', titlename, 'savepath', savepath)
% 
% 
% rr_6h = cellfun(@(x) sum(x(1,5:6)), watersip_6h, 'UniformOutput', false);
% rr_6h = [rr_6h{:}]';
% 
% titlename = 'hist_rr_6h';
% savepath = strcat(masterfolder, expmnt, '/Plottings/', strrep(titlename, ' ', '_'));
% plot_hist_cdf('x', x, 'y', rr_6h, 'x_label', x_label, 'titlename', titlename, 'savepath', savepath)

%% Rename files starting with xxx
% cd('/disk/r034/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_2010/')
% files = dir('EminusP_LA_RHc80_dqc0.2_bwd_20d_optimal_3h_glbland_*.mat'); % Get all files starting with xxx
% % Loop through each
% for id = 1:length(files)
% %     % Get the file name (minus the extension)
%     
%     [~, f] = fileparts(files(id).name);
%     
%     filename_parts = strsplit(f, '_');
%     
%     filename_parts_new = {filename_parts{1:9}, 'intact', filename_parts{10:end}};
%     
%     f_new = strjoin(filename_parts_new, '_');
%     disp(f_new)
% 
%     movefile(files(id).name, sprintf('%s.mat', f_new));
% 
% end






function [uptake_map, Pm_map, BL_uptake_map, BL_Pm_map, watersip_rr, watersip_stacktable,...
          mean_optimal_trajtime_map, mean_CWRT_map, P_LA_map, rr_L_tot_map, rr_NLL_tot_map, rr_NLO_tot_map, RH2_map,...
          LA_3D, Pm_AR_map] = WaterSip(varargin)

    %/ create a set of valid parameters and their default value
    pnames = {       'domfill',           'lon',             'lat',             'dqc',    'BLH_factor',       'area',...
                   'cond_land',    'cond_ocean',    'traj_rm_jump',      'rm_when_dz',    'NumWorkers', 'optimal_rr',...
               'output_rr_map',  'output_LA_3D', 'basin_lon_range', 'basin_lat_range', 'basin_z_range', 'basin_name',...
                          'AR'};  
           
    dflts  = {              [],              [],                [],                [],               1,           [],...
                            [],              [],                 0,                [],              [],           [],...
                             0,               0,                [],                [],              [],           [],...
                            []};
    
    [                  domfill,             lon,               lat,               dqc,      BLH_factor,         area,...
                     cond_land,      cond_ocean,                 ~,                 ~,      NumWorkers,   optimal_rr,...
                 output_rr_map,    output_LA_3D,   basin_lon_range,   basin_lat_range,   basin_z_range,   basin_name,...
                            AR] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
%%
    %/=====================================================================
    %/      Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 4 Jul 2024
    %/
    %/  
    %/ 4 Jul 2024: Allow to compute AR-related moisture source 
    %/             contribution (Pm_AR_map)
    %/
    %/
    %/       INPUT:
    %/           traj_rm_jump (deprecated):  [0 or 1] toggle on to slice traj of a particle with a *big* jump
    %/
    %/           rm_when_dz (deprecated):    dz threshold to remove traj with *big* jump
    %/
    %/           output_rr_map: To compute recycling ratios (Local, non-local
    %/                          land, non-local oceans) and output the final ratios on map.
    %/
    %/           output_LA_3D:  Output 3D traj freq, contribution, water loss
    %/                          and the BLH.
    %/
    %/    Reference: Sodemann et al. (2008)
    %/=====================================================================
    
    %/ OUTPUT: 
    %/      no need to process uptake_map, BL_uptake_map, BL_Pm_map anymore
    
    if output_LA_3D && (isempty(basin_lon_range) || isempty(basin_lat_range) || isempty(basin_z_range))
        error('Input basin_lon_range, basin_lat_range and basin_z_range if toggling ''output_LA_3D'' on!!');   
    end
    
    if optimal_rr > 1 || optimal_rr < 0
        error('optimal_rr must be from 0 to 1.');
    end
    
    % fprintf('*** Running WaterSip algorithm for the input trajs... *** \n');
    
    %/ Create reduction variables (we can just set 0, it speeds up and save memory!!!)
    %/ NOTE: whether or not they will be used, have to assign them to return the value.

    %/ Make sure consistent dimensions of the output even when no trajs
    uptake_map                = zeros(length(lon), length(lat));      
    Pm_map                    = zeros(length(lon), length(lat));
    Pm_AR_map                 = zeros(length(lon), length(lat));
    BL_uptake_map             = zeros(length(lon), length(lat));   
    BL_Pm_map                 = zeros(length(lon), length(lat));   
    mean_optimal_trajtime_map = zeros(length(lon), length(lat));
    mean_CWRT_map             = zeros(length(lon), length(lat));
    traj_startpt_count_map    = zeros(length(lon), length(lat));
    P_LA_map                  = zeros(length(lon), length(lat));
    RH2_map                   = zeros(length(lon), length(lat));
    rr_L_tot_map              = zeros(length(lon), length(lat));
    rr_NLL_tot_map            = zeros(length(lon), length(lat));
    rr_NLO_tot_map            = zeros(length(lon), length(lat));
    freq_3D                   = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
    rr_3D                     = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
    contr_3D                  = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
    loss_3D                   = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
    u_3D                      = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
    v_3D                      = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
    w_3D                      = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
    BLH_2D                    = zeros(length(basin_lon_range), length(basin_lat_range));
    watersip_rr               = cell(size(domfill.traj, 1), 1);
    watersip_stacktable       = []; %/ only save this for debugging.
    
    kkn      = size(domfill.traj, 1);
    dates_dt = domfill.dates_dt';    %/ for computing AR contr.
    dates    = datetime2int(dates_dt, 'yyyyMMddHHmm');
    trajtime = domfill.trajtime';
    dt       = abs(trajtime(1)-trajtime(2));
    
    if kkn > 1
        x   = squeeze(domfill.traj(:,1,:));     %/ from ntraj x pos x trajtime -> ntraj x trajtime
        y   = squeeze(domfill.traj(:,2,:)); 
        z   = squeeze(domfill.traj(:,3,:)); 
    else
        x   = squeeze(domfill.traj(:,1,:))';    %/ transpose to avoid dimension order being mistakenedly switched by squeeze() when kkn == 1.
        y   = squeeze(domfill.traj(:,2,:))'; 
        z   = squeeze(domfill.traj(:,3,:))'; 
    end
    %/ IMPORTANT: Always convert into [0 360) for the convenience of the following codes 
    x(x < 0) = x(x < 0) + 360;
    
    q    = domfill.q(:,:);                    %/ g/kg
    dq   = -1*diff(q, 1, 2);                  %/ +ve: water uptake. -ve: water loss.
    BLH  = domfill.BLH(:,:);
    topo = domfill.topo(:,:);
    mass = domfill.mass;
    RH2  = domfill.RH2;

    % ============================================================
    % %% test
    % dx = nan(kkn,1); dy = nan(kkn,1); dz = nan(kkn,1); 
    % for kk = 1:kkn
    %     mean_x12_kk   = lon_movpairmean_360(x(kk,:), 300); %/ my function to handle the mean across the lon boundary 
    %     mean_y12_kk   = movmean(y(kk,:),   2, 'Endpoints','discard');
    %     mean_z12_kk   = movmean(z(kk,:),   2, 'Endpoints','discard');
    % %     mean_BLH12_all = movmean(BLH, 2, 'Endpoints','discard');
    % 
    %     dx(kk) = mean_x12_kk(1) - x(kk,1); %/ in deg
    %     dy(kk) = mean_y12_kk(1) - y(kk,1); %/ in deg
    %     dz(kk) = mean_z12_kk(1) - z(kk,1); %/ in m
    % end
    % %%
    % clc;
    % fprintf('Mean dx: %.2f deg\n', mean(dx)) 
    % fprintf('    %.2f   %.2f   %.2f   %.2f   %.2f\n',...
    %          quantile(dx, 0.1), quantile(dx, 0.25), quantile(dx, 0.5), quantile(dx, 0.75), quantile(dx, 0.9)) 
    % 
    % fprintf('Mean dy: %.2f deg\n', mean(dy)) 
    % fprintf('    %.2f   %.2f   %.2f   %.2f   %.2f\n',...
    %          quantile(dy, 0.1), quantile(dy, 0.25), quantile(dy, 0.5), quantile(dy, 0.75), quantile(dy, 0.9)) 
    %      
    % fprintf('Mean dz: %.2f m\n', mean(dz)) 
    % fprintf('    %.2f   %.2f   %.2f   %.2f   %.2f\n',...
    %          quantile(dz, 0.1), quantile(dz, 0.25), quantile(dz, 0.5), quantile(dz, 0.75), quantile(dz, 0.9)) 
    %      
    % disp('==============================')
    % fprintf('Mean abs. dx: %.2f deg\n', mean(abs(dx))) 
    % fprintf('    %.2f   %.2f   %.2f   %.2f   %.2f\n',...
    %          quantile(abs(dx), 0.1), quantile(abs(dx), 0.25), quantile(abs(dx), 0.5), quantile(abs(dx), 0.75), quantile(abs(dx), 0.9)) 
    % 
    % fprintf('Mean abs. dy: %.2f deg\n', mean(abs(dy))) 
    % fprintf('    %.2f   %.2f   %.2f   %.2f   %.2f\n',...
    %          quantile(abs(dy), 0.1), quantile(abs(dy), 0.25), quantile(abs(dy), 0.5), quantile(abs(dy), 0.75), quantile(abs(dy), 0.9)) 
    % 
    % fprintf('Mean abs. dz: %.2f m\n', mean(abs(dz))) 
    % fprintf('    %.2f   %.2f   %.2f   %.2f   %.2f\n',...
    %          quantile(abs(dz), 0.1), quantile(abs(dz), 0.25), quantile(abs(dz), 0.5), quantile(abs(dz), 0.75), quantile(abs(dz), 0.9)) 
    
    %/ Printed Results
    % Mean dx: -0.13 deg
    %     -1.20   -0.71   -0.23   0.03   0.28
    % Mean dy: -0.15 deg
    %     -0.73   -0.28   -0.06   0.09   0.30
    % Mean dz: -189.78 m
    %     -461.94   -255.44   -119.30   -30.56   28.79
    % ==============================
    % Mean abs. dx: 0.74 deg
    %     0.05   0.13   0.35   0.77   1.23
    % Mean abs. dy: 0.30 deg
    %     0.03   0.07   0.18   0.39   0.75
    % Mean abs. dz: 213.46 m
    %     19.30   53.33   129.30   262.28   465.37
    % ============================================================
    flag_shutdown = 0;
    if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if parpool is not open
        parpool('Threads', NumWorkers) %/ much faster
        flag_shutdown = 1;
    end

    tic;
    parfor kk = 1:kkn 
    % for kk = 1:kkn      %/ NOTE: although for-loop uses multi-thread, but parfor is still faster!
        % fprintf('kk = %d/%d \n', kk, kkn);
    
        %====== Set broadcast vars (using 'parallel.pool.Constant' does not speedup!) ======%
        trajtime_bc = trajtime;
        x_bc        = x(kk,:)';
        y_bc        = y(kk,:)';
        z_bc        = z(kk,:)';
        q_bc        = q(kk,:)';
        dq_bc       = dq(kk,:)';
        BLH_bc      = BLH(kk,:)';
        topo_bc     = topo(kk,:)';
        mass_bc     = mass(kk);
        RH2_bc      = RH2(kk);

        %====== Testing my result using the example in Sodemann 2008) ========%
        % id = 'test';
        % trajtime = [0:-6:-54]';
        % x = repmat(100, length(trajtime), 1); 
        % y = repmat(20, length(trajtime), 1); 
        % z = repmat(500, length(trajtime), 1); %/ z = 1000m
        % z(4) = 1000;
        % q = [2.1, 2.6, 2.6, 2.6, 2.3, 2.3, 2.5, 1, 1, 0.2]'; 
        % dq = -1*diff(q);   %/ +ve: water uptake. -ve: water loss.
        % BLH = repmat(1200, length(trajtime), 1); 
        % BLH(4) = 100; %/ particle above BL at t = -18 hr
        % RH = 85; %/ q needs to be in kg/kg
        % 
        % %/ retain only the traj *after* the big jump 
        % if traj_rm_jump
        %     dz_ipart = diff(z_bc);
        %     first_change_pt = min(find(dz_ipart > rm_when_dz));
        % 
        %     if first_change_pt == 1
        %         disp('*** The sliced trajectory only has one time step. Skip it. ***')
        %         continue;  
        %     end
        % 
        %     if ~isempty(first_change_pt) 
        %         trajtime_bc = trajtime_bc(1:first_change_pt);
        %         x_bc = x_bc(1:first_change_pt);
        %         y_bc = y_bc(1:first_change_pt);
        %         z_bc = z_bc(1:first_change_pt);
        %         q_bc = q_bc(1:first_change_pt);
        %         dq_bc = -1*diff(q_bc);
        %         BLH_bc = BLH_bc(1:first_change_pt)';
        % 
        %         NoOf_sliced_traj = NoOf_sliced_traj + 1;
        %     else
        %         NoOf_nonsliced_traj = NoOf_nonsliced_traj + 1;
        %     end
        % end
    
        %=======================================================================%
        [ind_lon_xy,  ind_lat_xy]  = point2map('pt_x', x_bc, 'pt_y', y_bc, 'lon', lon, 'lat', lat);
        ind_lon_des = ind_lon_xy(1);
        ind_lat_des = ind_lat_xy(1);
        
        %/ Input the pre-processed cond_land & cond_ocean [speed-up]
        flag_land = show_land_or_ocean_hydrosheds('pos_lon_array', x_bc, 'pos_lat_array', y_bc,...
                                       'lon_grids', lon, 'lat_grids', lat, 'land_or_ocean', 1, 'output_as_logical', 1,...
                                       'cond_land', cond_land, 'cond_ocean', cond_ocean);
        
        MaxTrajTime = length(trajtime);         
        if ~isempty(optimal_rr)                 %/ [optimal_mode is on]: recursively run WaterSip over different lengths of trajtime.
            recursive_trajtime = 3:MaxTrajTime; 
            if isempty(recursive_trajtime)
                error('MaxTrajTime is too short! \n');
            end
        else
            recursive_trajtime = MaxTrajTime;   %/ [optimal_mode is off]: run WaterSip at a fixed maxtrajtime.
        end
    
        %/ IMPORTANT: The current recursive_trajtime loop is optimized (?). 
        %/            Putting all one m_loop is found to be even SLOWER.
    
        %/ Since the recursive loop efficiently scans over the needed indices,
        %/ here just preallocate arrays with the largest length. 10% speedup!
        N_max  = recursive_trajtime(end);
        loss   = zeros(N_max,1);  dq_BL  = nan(N_max,1);    dq_FT  = nan(N_max,1);  
        f      = nan(N_max,1);    e      = nan(N_max,1);  
        f_tot  = zeros(N_max,1);  e_tot  = zeros(N_max,1);  rr_tot = zeros(N_max,1);
    
        if output_rr_map  %/ do not initilize useless vars --> speedup!
            dq_L     = nan(N_max,1);   dq_NLL     = nan(N_max,1);   dq_NLO     = nan(N_max,1);
            rr_L     = nan(N_max,1);   rr_NLL     = nan(N_max,1);   rr_NLO     = nan(N_max,1);
            rr_L_tot = zeros(N_max,1); rr_NLL_tot = zeros(N_max,1); rr_NLO_tot = zeros(N_max,1); 
            rr_tot2  = zeros(N_max,1);
        end
    
        % time_elapse_block1 = 0;
        % time_elapse_block2 = 0;
        % time_elapse_block3 = 0;
        for N = recursive_trajtime
            q_nearest_uptake = []; 
            for j = N-1:-1:1 %/ loop from the 2nd-to-last to the 1st point (since no dq for the last index)
                %/ WATER LOSS
                if dq_bc(j) < 0         %/ we consider all the loss (regardless of the magnitude!) -> as did in Cheng and Lu 2022
                                        %/ Why? Because recycling ratio (f, e) depends also on current q. The ratio
                                        %/ will be excergerate when there has been a continuous, weak loss of moisture.
                                        %/ To workaround this, we just deplete the previous source whenever there is
                                        %/ a drop in the parcel's moisture. 
                    loss(j) = dq_bc(j); %/ record any water loss. 
    
                    %/ update dq_BL and dq_FT if ...
                    %/      1. There is water uptake before this loss,
                    %/      2. And it is not the water loss at the starting timepoint (mind this!)
                    if ~isempty(q_nearest_uptake) && j ~= 1
    %                     tic;
                        for m = N-1:-1:j 
                            %/ all previous uptakes are discounted proportionally with water loss
                            %/ NOTE: recycling ratios will stay the same (see the proof in PDF)
                            dq_BL(m)  = dq_BL(m)  + loss(j)*f(m); 
                            dq_FT(m)  = dq_FT(m)  + loss(j)*e(m);
    
                            if output_rr_map
                                dq_L(m)   = dq_L(m)   + loss(j)*rr_L(m); 
                                dq_NLL(m) = dq_NLL(m) + loss(j)*rr_NLL(m);
                                dq_NLO(m) = dq_NLO(m) + loss(j)*rr_NLO(m);
                            end
    
                            %/ As we skipped those small |dq| > dqc, the previous water uptake (dq_BL, dq_FT) can be depleted to -ve.
                            %/ In this case, reset it to zero.
                            %/ So, you sometimes may see zero dq_BL and dq_FT, but non-zero f and e.
                            if dq_BL(m)  < 0                      dq_BL(m)  = 0;  end
                            if dq_FT(m)  < 0                      dq_FT(m)  = 0;  end
    
                            if output_rr_map
                                if dq_L(m)   < 0                  dq_L(m)   = 0;  end
                                if dq_NLL(m) < 0                  dq_NLL(m) = 0;  end
                                if dq_NLO(m) < 0                  dq_NLO(m) = 0;  end
                            end
                        end
                    end
    
                %/ WATER UPTAKE
                elseif dq_bc(j) > dqc            %/ NOTE: Trivial dq (< dqc) is not considered but adds up into q value.
                                                 %/       --> may slightly overestimate the sources' contribution.
                    q_nearest_uptake = q_bc(j);  %/ update the current q value
    
                    %/ classify the recycling source by BLH
                    if z_bc(j) <= BLH_factor*BLH_bc(j)
                        dq_BL(j) = dq_bc(j);
                    else
                        dq_FT(j) = dq_bc(j);
                    end
    
                    if output_rr_map
                        %/ classify the recycling source by local (L), non-local land (NLL) and non-local oceans (NLO)
                        if ind_lon_xy(j) == ind_lon_des && ind_lat_xy(j) == ind_lat_des  %/ still in local grid?
                            dq_L(j) = dq_bc(j);
                        elseif flag_land(j)                    %/ over land grid?
                            dq_NLL(j) = dq_bc(j);
                        else                                   %/ else then should be over oceans
                            dq_NLO(j) = dq_bc(j);
                        end
                    end
    
                    for m = N-1:-1:j %/ update the importance of previous water uptakes till time j.
                        f(m)      = dq_BL(m) /q_nearest_uptake;
                        e(m)      = dq_FT(m) /q_nearest_uptake;
    
                        if output_rr_map
                            rr_L(m)   = dq_L(m)  /q_nearest_uptake;
                            rr_NLL(m) = dq_NLL(m)/q_nearest_uptake;
                            rr_NLO(m) = dq_NLO(m)/q_nearest_uptake;
                        end
                    end
                end
    
                % tic;
                if j == 1  %/ IMPORTANT: This is a speed-critical block. To SAVE TIME, sum up ratios ONLY at j = 1.
                    f_tot(j)      = sum(f, 'omitnan'); 
                    e_tot(j)      = sum(e, 'omitnan'); 
                    rr_tot(j)     = f_tot(j) + e_tot(j);
    
                    if output_rr_map
                        rr_L_tot(j)   = sum(rr_L, 'omitnan'); 
                        rr_NLL_tot(j) = sum(rr_NLL, 'omitnan'); 
                        rr_NLO_tot(j) = sum(rr_NLO, 'omitnan');  
                        rr_tot2(j)    = rr_L_tot(j) + rr_NLL_tot(j) + rr_NLO_tot(j);
                    end
                end
                % time_elapse_block3 = time_elapse_block3 + toc;
            end
    
            %/ Optimal tracking
            if ~isempty(optimal_rr) && rr_tot(1) > optimal_rr
                loss   = loss(1:N);     dq_BL  = dq_BL(1:N);    dq_FT  = dq_FT(1:N);
                f      = f(1:N);        e      = e(1:N);        
                f_tot  = f_tot(1:N);    e_tot  = e_tot(1:N);    rr_tot = rr_tot(1:N);
                if output_rr_map
                    dq_L       = dq_L(1:N);       dq_NLL     = dq_NLL(1:N);      dq_NLO     = dq_NLO(1:N);
                    rr_L       = rr_L(1:N);       rr_NLL     = rr_NLL(1:N);      rr_NLO     = rr_NLO(1:N);
                    rr_L_tot   = rr_L_tot(1:N);   rr_NLL_tot = rr_NLL_tot(1:N);  rr_NLO_tot = rr_NLO_tot(1:N);
                    rr_tot2    = rr_tot2(1:N);
                end
                break;
            end
            % time_elapse_block2 = time_elapse_block2 + toc;
        end
        % fprintf('\ntime_elapse_block1: %.5f sec (x60,000)\n', time_elapse_block1)
        % fprintf('time_elapse_block2: %.5f sec (x60,000)\n', time_elapse_block2)
        % fprintf('time_elapse_block3: %.5f sec (x60,000)\n', time_elapse_block3)
    
        if ismember(kk, 1:10) 
            T_vars = {'Time','x','y', 'z', 'BLH', 'topo', 'q', 'dq0', 'loss', 'dq_BL', 'dq_FT', 'f', 'e', 'f_tot', 'e_tot', 'rr_tot'};
            watersip_table = nan(N, length(T_vars)); %/ recording water uptake from the corrected traj
            watersip_table(:, 1)  = trajtime_bc(1:N);
            watersip_table(:, 2)  = x_bc(1:N); 
            watersip_table(:, 3)  = y_bc(1:N);
            watersip_table(:, 4)  = z_bc(1:N);
            watersip_table(:, 5)  = BLH_bc(1:N);
            watersip_table(:, 6)  = topo_bc(1:N);
            watersip_table(:, 7)  = q_bc(1:N);
            watersip_table(1:end-1, 8)  = dq_bc(1:N-1);
            watersip_table(:, 9)  = loss;
            watersip_table(:, 10) = dq_BL;
            watersip_table(:, 11) = dq_FT;
            watersip_table(:, 12) = f;
            watersip_table(:, 13) = e;
            watersip_table(:, 14) = f_tot;
            watersip_table(:, 15) = e_tot;
            watersip_table(:, 16) = rr_tot;
    
            if output_rr_map
                T_vars = [T_vars, {'rr_tot2', 'dq_L', 'dq_NLL', 'dq_NLO', 'rr_L', 'rr_NLL', 'rr_NLO', 'rr_L_tot', 'rr_NLL_tot', 'rr_NLO_tot'}];
                watersip_table(:, 17) = rr_tot2;
                watersip_table(:, 18) = dq_L;
                watersip_table(:, 19) = dq_NLL;
                watersip_table(:, 20) = dq_NLO;
                watersip_table(:, 21) = rr_L;
                watersip_table(:, 22) = rr_NLL;
                watersip_table(:, 23) = rr_NLO;
                watersip_table(:, 24) = rr_L_tot;
                watersip_table(:, 25) = rr_NLL_tot;
                watersip_table(:, 26) = rr_NLO_tot;
            end
    
            watersip_table      = array2table(watersip_table, 'VariableNames', T_vars);
            watersip_stacktable = [watersip_stacktable, {watersip_table}];                          %/ just for checking the initial water loss.
            % watersip_stacktable{kk} = watersip_table;                          %/ just for checking the initial water loss.
        end
    
        %/ Compute contribution-weighted residence time (CWRT) (Laderach and Sodemann 2016)
        f_bc = f; f_bc(isnan(f_bc)) = 0;  %/ nan to zeros
        e_bc = e; e_bc(isnan(e_bc)) = 0;  %/ nan to zeros
        fe = f_bc + e_bc;
        ind_nonnan = find(~isnan(fe));
        CWRT = sum(ind_nonnan.*fe)/sum(fe); %/ CAVEAT: CWRT can be NaN or Inf!
    
        %/ Store the final recycling ratios
        if output_rr_map
            watersip_rr{kk} = [N, CWRT, q_bc(1), loss(1), f_tot(1), e_tot(1), rr_tot(1), rr_L_tot(1), rr_NLL_tot(1), rr_NLO_tot(1)];   %/ (index = 1 as WaterSip goes from earliest to latest)
        else
            watersip_rr{kk} = [N, CWRT, q_bc(1), loss(1), f_tot(1), e_tot(1), rr_tot(1)];
        end
    
        %>>>>>>>>>>>>>>>>> optimal trajtime map (use the true destination pt) <<<<<<<<<<<<<<<<<<<%
        mean_optimal_trajtime_map_kk = zeros(length(lon),length(lat));
        mean_CWRT_map_kk             = zeros(length(lon),length(lat));
        traj_startpt_count_map_kk    = zeros(length(lon),length(lat));
    
        if ~isnan(CWRT) && ~isinf(CWRT)
            mean_optimal_trajtime_map_kk(ind_lon_des, ind_lat_des) = N;
            mean_CWRT_map_kk(ind_lon_des, ind_lat_des)             = CWRT;
            traj_startpt_count_map_kk(ind_lon_des, ind_lat_des)    = 1;
    
            mean_optimal_trajtime_map = mean_optimal_trajtime_map  + mean_optimal_trajtime_map_kk;  %/ Reduction in parfor is allowed! (save memory!)
            mean_CWRT_map             = mean_CWRT_map              + mean_CWRT_map_kk;  %/ Reduction in parfor is allowed! (save memory!)
            traj_startpt_count_map    = traj_startpt_count_map     + traj_startpt_count_map_kk;
        end

        %>>>>>>>>>>>>>>>>>>> P_LA map <<<<<<<<<<<<<<<<<<<%
        P_LA_map_kk = zeros(length(lon),length(lat));
        P_LA = mass_bc*-1*loss(1)/1000;                                 %/ Using Stohl and James's method (actually this is E-P)
        P_LA_map_kk(ind_lon_des, ind_lat_des) = P_LA;                
        P_LA_map = P_LA_map + P_LA_map_kk;                              %/ Reduction in parfor is allowed! (save memory!)
    
        %>>>>>>>>>>>>>>>>>>> RH2 map <<<<<<<<<<<<<<<<<<<%
        RH2_map_kk = zeros(length(lon),length(lat));
        RH2_map_kk(ind_lon_des, ind_lat_des) = RH2_bc;
        RH2_map = RH2_map + RH2_map_kk;                                 %/ Reduction in parfor is allowed! (save memory!)
        
        
        %/ combine f and e tgt (count all their contributions)
        rr = sum([f, e], 2, 'omitnan');  %/ recycling ratio along the trajectory

        %/ Double-check
        if rr_tot(1) > 1
            error('total recycling ratio > 1!');
        elseif rr_tot(1) < 0
            error('total recycling ratio < 0!');
        end

        if any(rr > 1)
            error('rr contains values > 1!');
        elseif any(rr < 0)
            error('rr contains values < 0!');
        end

        %>>>>>>>>>>>>>>>>>>> Moisture source (Pm_map) <<<<<<<<<<<<<<<<<<<%
        Pm_map_kk    = zeros(length(lon),length(lat));
        Pm_AR_map_kk = zeros(length(lon),length(lat));
        % uptake_map_kk = zeros(length(lon),length(lat));
        
        ind = find(rr ~= 0);    %/ rr is in (ntrajtime, 1)
        for i = 1:length(ind)
            contr  = -1 * rr(ind(i)) * loss(1); %/ consider water uptake at different location along the trajectory (g/kg)
            Pm     = mass_bc*contr/1000;        %/ g/kg -> kg/kg

            Pm_map_kk(ind_lon_xy(ind(i)), ind_lat_xy(ind(i))) = Pm_map_kk(ind_lon_xy(ind(i)), ind_lat_xy(ind(i))) + Pm;
            % uptake_map_kk(ind_lon_meanxy(ind(i)), ind_lat_meanxy(ind(i))) = uptake_map_kk(ind_lon_meanxy(ind(i)), ind_lat_meanxy(ind(i))) + src_contr;

            %/ Locate source contributions from AR if a logical matrix is given
            if ~isempty(AR)
                ind_t     = findismember_loop(AR.dates, dates(ind(i)));
                if isempty(ind_t)
                    error('Missing date %d in AR.dates!', dates(ind(i)));
                end
                AR_logical_t = squeeze(AR.logical(:,:,ind_t));
                
                [ind_AR_lon_xy,  ind_AR_lat_xy] = point2map('pt_x', x_bc(ind(i)), 'pt_y', y_bc(ind(i)), 'lon', AR.lon, 'lat', AR.lat);
                
                if AR_logical_t(ind_AR_lon_xy, ind_AR_lat_xy) == 1  %/ If the particle is at the grid cell where AR is identified
                    % disp('Identified AR contribution!')
                    Pm_AR_map_kk(ind_lon_xy(ind(i)), ind_lat_xy(ind(i))) = Pm_AR_map_kk(ind_lon_xy(ind(i)), ind_lat_xy(ind(i))) + Pm;
                end
            end
        end

        %/ Double-check
        eps = 1e-10;
        if abs(sum(Pm_map_kk,'all')/P_LA - rr_tot(1)) > eps
            sum(Pm_map_kk,'all')/P_LA
            rr_tot(1)
            error('The ratio of P_LA to the summed sources does not equal to total recycling ratio!');
        end
        
        Pm_map     = Pm_map    + Pm_map_kk; 
        Pm_AR_map  = Pm_AR_map + Pm_AR_map_kk; 
        % uptake_map = uptake_map + uptake_map_kk;     %/ Reduction in parfor is allowed! (save memory!)
        
        % plot_contfmap('contf_data', Pm_AR_map, 'contf_lon', lon, 'contf_lat', lat);

        %>>>>>>>>>>>>>>>>>>> rr_L_tot, rr_NLL_tot, rr_NLO_tot map (use the true destination pt) <<<<<<<<<<<<<<<<<<<%
        if output_rr_map
            rr_L_tot_map_kk    = zeros(length(lon),length(lat));
            rr_NLL_tot_map_kk  = zeros(length(lon),length(lat));
            rr_NLO_tot_map_kk  = zeros(length(lon),length(lat));
    
            rr_L_tot_map_kk(ind_lon_des, ind_lat_des) = rr_L_tot(1);
            rr_L_tot_map = rr_L_tot_map + rr_L_tot_map_kk;                              %/ Reduction in parfor is allowed! (save memory!)
    
            rr_NLL_tot_map_kk(ind_lon_des, ind_lat_des) = rr_NLL_tot(1);
            rr_NLL_tot_map = rr_NLL_tot_map + rr_NLL_tot_map_kk;                              %/ Reduction in parfor is allowed! (save memory!)
    
            rr_NLO_tot_map_kk(ind_lon_des, ind_lat_des) = rr_NLO_tot(1);
            rr_NLO_tot_map = rr_NLO_tot_map + rr_NLO_tot_map_kk;                              %/ Reduction in parfor is allowed! (save memory!)
        end
    
        %>>>>>>>>>>>> output 3D matrics in lagrangian (speed-critical step) >>>>>>>>%
        if output_LA_3D
            %/       - Retrieve the cell index in the given basin_lon_range, basin_lat_range and basin_z_range.
            [~, ind_box_lon] = find_nearest_val(basin_lon_range, x_bc(1:N), 'lon');
            [~, ind_box_lat] = find_nearest_val(basin_lat_range, y_bc(1:N));
            [~,   ind_box_z] = find_nearest_val(basin_z_range,   z_bc(1:N)+topo_bc(1:N));         
    
            %/ combine f and e tgt (count all their contributions)
            % uptake = sum([dq_BL, dq_FT], 2, 'omitnan');   %/ here 'uptake' means how much water the parcel gains, but NOT the contribution (it is rr!).

            %/ loop only those in the 3D box.
            I = strfind((~isnan(ind_box_lon) & ~isnan(ind_box_lat) & ~isnan(ind_box_z))', 1); %/ faster than find()!
    
            freq_3D_kk   = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
            rr_3D_kk     = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
            contr_3D_kk  = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
            loss_3D_kk   = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
            u_3D_kk      = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
            v_3D_kk      = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
            w_3D_kk      = zeros(length(basin_lon_range), length(basin_lat_range), length(basin_z_range));
            BLH_2D_kk    = zeros(length(basin_lon_range), length(basin_lat_range));
    
            %/ let nans to 0 in loss prior to reduction operation.
            loss_bc = loss;
            loss_bc(isnan(loss_bc)) = 0;
    
            %/ Compute Traj velocities
            %/ using the actual xyz, then the velocity will correspond to mean_xyz.
            %/ since it's bwd now, we let the ending point to be t + 1.
            %/ Note that uvw has one less element as other data, but the last
            %/ row is always NaN, so doesn't matter.
            [u, v, w] = traj_velocity('lon1', x_bc(2:N),   'lat1', y_bc(2:N),   'z1', z_bc(2:N),...
                                      'lon2', x_bc(1:N-1), 'lat2', y_bc(1:N-1), 'z2', z_bc(1:N-1), 'dt', dt); %/ (u,v) in m/s, w in cm/s
    
            for i = 1:length(I)
                freq_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                            freq_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))   + 1;
    
                rr_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                            rr_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))     + rr(I(i));
    
                contr_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                            contr_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))  + -1 * rr(I(i)) * loss_bc(1);
    
                loss_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                            loss_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))   + loss_bc(I(i));       
    
                %/ Since no velocity can be computed for the last index. 
                %/ This should not affect the result much.
                if i ~= length(I)
                    u_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                                u_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))   + u(I(i));   
    
                    v_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                                v_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))   + v(I(i));   
    
                    w_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                                w_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))   + w(I(i));   
                end
                BLH_2D_kk(ind_box_lon(I(i)), ind_box_lat(I(i))) = ...
                            BLH_2D_kk(ind_box_lon(I(i)), ind_box_lat(I(i))) + BLH_bc(I(i));
            end
            % size(freq_3D_kk)
            % size(freq_3D)
            freq_3D   = freq_3D   + freq_3D_kk;   %/ frequency of the traj position (3D)
            rr_3D     = rr_3D     + rr_3D_kk;
            contr_3D  = contr_3D  + contr_3D_kk;
            loss_3D   = loss_3D   + loss_3D_kk;
            u_3D      = u_3D      + u_3D_kk;
            v_3D      = v_3D      + v_3D_kk;
            w_3D      = w_3D      + w_3D_kk;
            BLH_2D    = BLH_2D    + BLH_2D_kk;
        end
    
        %>>>>>>>>>>>> update the water uptake map (from BL only) <<<<<<<<<<<<<%
        BL_uptake_map_kk = zeros(length(lon),length(lat));
        BL_Pm_map_kk     = zeros(length(lon),length(lat));
    
        ind = find(~isnan(f));
        for i = 1:length(ind)
            uptake = -1 * f(ind(i)) * loss(1); %/ consider water uptake at different location along the trajectory (g/kg)
            Pm     = mass_bc*uptake/1000;
    
            BL_uptake_map_kk(ind_lon_xy(ind(i)), ind_lat_xy(ind(i))) = BL_uptake_map_kk(ind_lon_xy(ind(i)), ind_lat_xy(ind(i))) + uptake;
            BL_Pm_map_kk    (ind_lon_xy(ind(i)), ind_lat_xy(ind(i))) = BL_Pm_map_kk    (ind_lon_xy(ind(i)), ind_lat_xy(ind(i))) + Pm;
        end
    
        %/ Reduction in parfor is allowed! (save memory!)
        BL_uptake_map = BL_uptake_map + BL_uptake_map_kk;
        BL_Pm_map     = BL_Pm_map     + BL_Pm_map_kk;
    end
    fprintf('WaterSip Loop Time: %.2f s\n', toc);
    
    if kkn ~= 0
        % size(Pm_map)
        % size(area)
        Pm_map     = Pm_map(:,2:end-1)./area;         %/ remove two poles (90N and 90S) as they have very small area.
        Pm_AR_map  = Pm_AR_map(:,2:end-1)./area;      %/ remove two poles (90N and 90S) as they have very small area.
        P_LA_map   = P_LA_map(:,2:end-1)./area;       %/ remove two poles (90N and 90S) as they have very small area.
        BL_Pm_map  = BL_Pm_map(:,2:end-1)./area;      %/ remove two poles (90N and 90S) as they have very small area.
        
        if output_rr_map
            rr_L_tot_map   = rr_L_tot_map./traj_startpt_count_map;    %/ it will contain nans, but retrieve_WSV function can deal with it.
            rr_NLL_tot_map = rr_NLL_tot_map./traj_startpt_count_map;  %/ it will contain nans, but retrieve_WSV function can deal with it.
            rr_NLO_tot_map = rr_NLO_tot_map./traj_startpt_count_map;  %/ it will contain nans, but retrieve_WSV function can deal with it.
            
            T_vars2     = {'Ntrajtime', 'Wgted_Resi_Time', 'q', 'loss', 'f_tot', 'e_tot', 'rr_tot', 'rr_L_tot', 'rr_NLL_tot', 'rr_NLO_tot'};
        else
            T_vars2     = {'Ntrajtime', 'Wgted_Resi_Time', 'q', 'loss', 'f_tot', 'e_tot', 'rr_tot'};
        end
        
        watersip_rr = cat(1,watersip_rr{:});
        watersip_rr = array2table(watersip_rr, 'VariableNames', T_vars2);
    
        RH2_map                   = RH2_map./traj_startpt_count_map;                   %/ it will contain nans, but retrieve_WSV function can deal with it.
        mean_optimal_trajtime_map = mean_optimal_trajtime_map./traj_startpt_count_map; %/ it will contain nans, but retrieve_WSV function can deal with it.
        mean_CWRT_map             = mean_CWRT_map./traj_startpt_count_map;             %/ it will contain nans, but retrieve_WSV function can deal with it.
    else
        %/ If no traj (kkn == 0), we need to make zero matrices to ensure consistent
        %/ dimension of the output!
        Pm_map    = zeros(length(lon),length(lat)-2);
        Pm_AR_map = zeros(length(lon),length(lat)-2);
        P_LA_map  = zeros(length(lon),length(lat)-2);
        BL_Pm_map = zeros(length(lon),length(lat)-2);
        
        if output_rr_map
            rr_L_tot_map = zeros(length(lon),length(lat));
            rr_NLL_tot_map = zeros(length(lon),length(lat));
            rr_NLO_tot_map = zeros(length(lon),length(lat));
            T_vars2     = {'Ntrajtime', 'Wgted_Resi_Time', 'q', 'loss', 'f_tot', 'e_tot', 'rr_tot', 'rr_L_tot', 'rr_NLL_tot', 'rr_NLO_tot'};
        else
            T_vars2     = {'Ntrajtime', 'Wgted_Resi_Time', 'q', 'loss', 'f_tot', 'e_tot', 'rr_tot'};
        end
    end
    
    %/ Double-check the global moisture attribution (relative to P_LA)
    RR = sum(Pm_map, 'all')/sum(P_LA_map, 'all');
    if RR > 1
        error('Global moisture attribution exceeds 1 (100%)!');
    elseif RR < 0
        error('Global moisture attribution is negative!');
    end
    
    LA_3D = [];
    if output_LA_3D
        rr_3D     = rr_3D./freq_3D;                 %/ 3hrly traj-mean (will contain nans)
        contr_3D  = contr_3D./freq_3D;              %/ 3hrly traj-mean (will contain nans)
        loss_3D   = loss_3D./freq_3D;               %/ 3hrly traj-mean (will contain nans)
        BLH_2D    = BLH_2D./sum(freq_3D,3);         %/ 3hrly traj-mean (will contain nans)
        LA_3D.('freq_3D')         = freq_3D;
        LA_3D.('rr_3D')           = rr_3D;
        LA_3D.('contr_3D')        = contr_3D;
        LA_3D.('loss_3D')         = loss_3D;
        LA_3D.('u_3D')            = u_3D;
        LA_3D.('v_3D')            = v_3D;
        LA_3D.('w_3D')            = w_3D;
        LA_3D.('BLH_2D')          = BLH_2D;
        LA_3D.('basin_lon_range') = basin_lon_range;
        LA_3D.('basin_lat_range') = basin_lat_range;
        LA_3D.('basin_z_range')   = basin_z_range;
        LA_3D.('basin_name')      = basin_name;
    end
    
    % fprintf('Aggregation step: %.5f sec\n', toc);
    % if traj_rm_jump
    %     disp(['NoOf_sliced_traj = ',    num2str(NoOf_sliced_traj)])
    %     disp(['NoOf_nonsliced_traj = ', num2str(NoOf_nonsliced_traj)])
    % end
    
    %/ Shut down the thread-based parpool initiated within the function
    %/ to avoid conflicting with the need of process-based parpool outside
    if flag_shutdown
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
    disp('!!! WaterSip (backward) is completed !!!')  

end

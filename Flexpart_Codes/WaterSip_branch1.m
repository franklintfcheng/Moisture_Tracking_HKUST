%%
function [uptake_map, Pm_map, BL_uptake_map, BL_Pm_map, watersip_rr, watersip_stacktable,...
          mean_optimal_trajtime_map, P_LA_map, rr_L_tot_map, rr_NLL_tot_map, rr_NLO_tot_map, LA_3D...
                        ] = WaterSip_branch1(varargin)

%/ create a set of valid parameters and their default value
pnames = {       'domfill', 'lon', 'lat', 'dqc', 'BLH_factor', 'area', 'traj_rm_jump', 'rm_when_dz', 'NumWorkers', 'optimal_rr', 'land_grids',...
           'output_rr_map', 'output_LA_3D', 'hs_lon_range', 'hs_lat_range', 'hs_z_range', 'hs_name'};  
       
dflts  = {              [],    [],    [],    [],            1,     [],             0,            [],          [],           [],           [],...
                         0,              0,             [],             [],           [],        []};

[domfill, lon, lat, dqc, BLH_factor, area, traj_rm_jump, rm_when_dz, NumWorkers, optimal_rr, land_grids,...
    output_rr_map, output_LA_3D, hs_lon_range, hs_lat_range, hs_z_range, hs_name] ...
               = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

%/ INPUT:
%/      traj_rm_jump:       [0 or 1] toggle on to slice traj of a particle with a *big* jump
%/
%/      rm_when_dz:         dz threshold to remove traj with *big* jump
%/
%/      output_rr_map:      To compute recycling ratios (Local, non-local
%/                          land, non-local oceans) and output the final ratios on map.
%/
%/      output_LA_3D:    
%/

%/ OUTPUT: 
%/      no need to process uptake_map, BL_uptake_map, BL_Pm_map anymore

if output_LA_3D && (isempty(hs_lon_range) || isempty(hs_lat_range) || isempty(hs_z_range))
    error('Input hs_lon_range, hs_lat_range and hs_z_range if toggling ''output_LA_3D'' on!!');   
end

fprintf('*** Running WaterSip algorithm for these trajs... *** \n');

%/ create reduction variables (we can just set 0, it speeds up and save memory!!!)
%/ NOTE: whether or not they will be used, have to assign them to return the value.
uptake_map                = 0;      
Pm_map                    = 0;      
BL_uptake_map             = 0;   
BL_Pm_map                 = 0;   
mean_optimal_trajtime_map = 0;
traj_startpt_count_map    = 0;
P_LA_map                  = 0;
rr_L_tot_map              = 0;
rr_NLL_tot_map            = 0;
rr_NLO_tot_map            = 0;
freq_3D                   = 0;
rr_3D                     = 0;
contr_3D                  = 0;
loss_3D                   = 0;
u_3D                      = 0;
v_3D                      = 0;
w_3D                      = 0;
BLH_2D                    = 0;


watersip_rr         = cell(size(domfill.traj, 1), 1);
watersip_stacktable = []; %/ only save this for debugging.

kkn      = size(domfill.traj, 1);
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
q    = domfill.q(:,:);                    %/ g/kg
dq   = -1*diff(q, 1, 2);                  %/ +ve: water uptake. -ve: water loss.
BLH  = domfill.BLH(:,:);
mass = domfill.mass;


% IMPORTANT: Convert to [0 360) for the convenience of the following codes if the lon of traj is in [-179 180].
if ~isempty(find(x < 0))
   x(x < 0) = x(x < 0) + 360;
end

if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
    parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
end
% fprintf('Initiation step: %.5f sec\n', toc);

% tic;
parfor kk = 1:kkn
% for kk = 327
%     fprintf('kk = %d/%d \n', kk, kkn);

    %====== Set broadcast vars (using 'parallel.pool.Constant' does not speedup!) ======%
    trajtime_bc = trajtime;
    x_bc        = x(kk,:)';
    y_bc        = y(kk,:)';
    z_bc        = z(kk,:)';
    q_bc        = q(kk,:)';
    dq_bc       = dq(kk,:)';
    BLH_bc      = BLH(kk,:)';
    mass_bc     = mass(kk);

    %====== Testing my result using the example in Sodemann 2008) ========%
%     id = 'test';
%     trajtime = [0:-6:-54]';
%     x = repmat(100, length(trajtime), 1); 
%     y = repmat(20, length(trajtime), 1); 
%     z = repmat(500, length(trajtime), 1); %/ z = 1000m
%     z(4) = 1000;
%     q = [2.1, 2.6, 2.6, 2.6, 2.3, 2.3, 2.5, 1, 1, 0.2]'; 
%     dq = -1*diff(q);   %/ +ve: water uptake. -ve: water loss.
%     BLH = repmat(1200, length(trajtime), 1); 
%     BLH(4) = 100; %/ particle above BL at t = -18 hr
%     RH = 85; %/ q needs to be in kg/kg

    %/ retain only the traj *after* the big jump 
%     if traj_rm_jump
%         dz_ipart = diff(z_bc);
%         first_change_pt = min(find(dz_ipart > rm_when_dz));
% 
%         if first_change_pt == 1
%             disp('*** The sliced trajectory only has one time step. Skip it. ***')
%             continue;  
%         end
% 
%         if ~isempty(first_change_pt) 
%             trajtime_bc = trajtime_bc(1:first_change_pt);
%             x_bc = x_bc(1:first_change_pt);
%             y_bc = y_bc(1:first_change_pt);
%             z_bc = z_bc(1:first_change_pt);
%             q_bc = q_bc(1:first_change_pt);
%             dq_bc = -1*diff(q_bc);
%             BLH_bc = BLH_bc(1:first_change_pt)';
% 
%             NoOf_sliced_traj = NoOf_sliced_traj + 1;
%         else
%             NoOf_nonsliced_traj = NoOf_nonsliced_traj + 1;
%         end
%     end
    
    mean_x   = lon_movpairmean_360(x_bc, 300); %/ my function to handle the mean across the lon boundary 
    mean_y   = movmean(y_bc,   2, 'Endpoints','discard');
    mean_z   = movmean(z_bc,   2, 'Endpoints','discard');
    mean_BLH = movmean(BLH_bc, 2, 'Endpoints','discard');

    
    %/ retrieve the grid cell index at the destination point.
    [ind_lon_des, ind_lat_des] = point2map('pt_x', x_bc(1), 'pt_y', y_bc(1), 'lon', lon, 'lat', lat);  

    %/ retrieve the grid cell index along trajectory (based on mean_x, mean_y). Detect if over land.
    %/ NOTE: no difference in speed compared to 'find_nearest_val'
    [ind_lon_meanxy, ind_lat_meanxy] = point2map('pt_x', mean_x, 'pt_y', mean_y, 'lon', lon, 'lat', lat);
    
    flag_land = which_on_land('pos_lon_array', mean_x, 'pos_lat_array', mean_y, 'full_logic_array', 1, 'land_grids', land_grids); %< time-consuming fn!!

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
    loss   = nan(N_max,1);    dq_BL  = nan(N_max,1);    dq_FT  = nan(N_max,1);  
    f      = nan(N_max,1);    e      = nan(N_max,1);  
    f_tot  = zeros(N_max,1);  e_tot  = zeros(N_max,1);  rr_tot = zeros(N_max,1);

    if output_rr_map  %/ do not initilize useless vars --> speedup!
        dq_L     = nan(N_max,1);   dq_NLL     = nan(N_max,1);   dq_NLO     = nan(N_max,1);
        rr_L     = nan(N_max,1);   rr_NLL     = nan(N_max,1);   rr_NLO     = nan(N_max,1);
        rr_L_tot = zeros(N_max,1); rr_NLL_tot = zeros(N_max,1); rr_NLO_tot = zeros(N_max,1); 
        rr_tot2  = zeros(N_max,1);
    end
    
%     time_elapse_block1 = 0;
%     time_elapse_block2 = 0;
%     time_elapse_block3 = 0;
    for N = recursive_trajtime
        q_nearest_uptake = []; 
%         tic;
        for j = N-1:-1:1 %/ loop from the 2nd-to-last to the 1st point (since no dq for the last index)
            %/ WATER LOSS
            if dq_bc(j) < 0         %/ we consider all the loss (regardless of the magnitude!)
                loss(j) = dq_bc(j); %/ record any sig. water loss.

                %/ update dq_BL and dq_FT if ...
                %/      1. there is water uptake before this loss
                %/      2. it is not the water loss at the starting timepoint (mind this!)
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
                        
                        %/ Since we skipped those small |dq| > dqc, dq_BL sometimes can be -ve.
                        %/ In this case, reset it to zero.
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
                if mean_z(j) <= BLH_factor*mean_BLH(j)
                    dq_BL(j) = dq_bc(j);
                else
                    dq_FT(j) = dq_bc(j);
                end
                
                if output_rr_map
                    %/ classify the recycling source by local (L), non-local land (NLL) and non-local oceans (NLO)
                    if ind_lon_meanxy(j) == ind_lon_des && ind_lat_meanxy(j) == ind_lat_des  %/ still in local grid?
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

%             tic;
            if j == 1  %/ IMPORTANT: This is a speed-critical block. To SAVE TIME, sum up ratios ONLY at j = 1.
                f_tot(j)      = nansum(f); 
                e_tot(j)      = nansum(e); 
                rr_tot(j)     = f_tot(j) + e_tot(j);

                if output_rr_map
                    rr_L_tot(j)   = nansum(rr_L); 
                    rr_NLL_tot(j) = nansum(rr_NLL); 
                    rr_NLO_tot(j) = nansum(rr_NLO); 
                    rr_tot2(j)    = rr_L_tot(j) + rr_NLL_tot(j) + rr_NLO_tot(j);
                end
            end
%             time_elapse_block3 = time_elapse_block3 + toc;
        end
        
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
%         time_elapse_block2 = time_elapse_block2 + toc;
    end
%     fprintf('\ntime_elapse_block1: %.5f sec (x60,000)\n', time_elapse_block1)
%     fprintf('time_elapse_block2: %.5f sec (x60,000)\n', time_elapse_block2)
%     fprintf('time_elapse_block3: %.5f sec (x60,000)\n', time_elapse_block3)
%     toc
    
    if ismember(kk, [1:1000:kkn]) 
        T_vars = {'Time','mean_x','mean_y', 'mean_z','mean_BLH', 'q', 'dq0', 'loss', 'dq_BL', 'dq_FT', 'f', 'e', 'f_tot', 'e_tot', 'rr_tot'};
        watersip_table = nan(N, length(T_vars)); %/ recording water uptake from the corrected traj
        watersip_table(:,1) = trajtime_bc(1:N);
        watersip_table(1:end-1,2) = mean_x(1:N-1);  %/ since mean_x is based on x at two time steps, so its length is N-1.
        watersip_table(1:end-1,3) = mean_y(1:N-1);
        watersip_table(1:end-1,4) = mean_z(1:N-1);
        watersip_table(1:end-1,5) = mean_BLH(1:N-1);
        watersip_table(:, 6)  = q_bc(1:N);
        watersip_table(1:end-1, 7) = dq_bc(1:N-1);
        watersip_table(:, 8)  = loss;
        watersip_table(:, 9)  = dq_BL;
        watersip_table(:, 10) = dq_FT;
        watersip_table(:, 11) = f;
        watersip_table(:, 12) = e;
        watersip_table(:, 13) = f_tot;
        watersip_table(:, 14) = e_tot;
        watersip_table(:, 15) = rr_tot;

        if output_rr_map
            T_vars = [T_vars, {'rr_tot2', 'dq_L', 'dq_NLL', 'dq_NLO', 'rr_L', 'rr_NLL', 'rr_NLO', 'rr_L_tot', 'rr_NLL_tot', 'rr_NLO_tot'}];
            watersip_table(:, 16) = rr_tot2;
            watersip_table(:, 17) = dq_L;
            watersip_table(:, 18) = dq_NLL;
            watersip_table(:, 19) = dq_NLO;
            watersip_table(:, 20) = rr_L;
            watersip_table(:, 21) = rr_NLL;
            watersip_table(:, 22) = rr_NLO;
            watersip_table(:, 23) = rr_L_tot;
            watersip_table(:, 24) = rr_NLL_tot;
            watersip_table(:, 25) = rr_NLO_tot;
        end
        
        watersip_table      = array2table(watersip_table, 'VariableNames', T_vars);
        watersip_stacktable = [watersip_stacktable, {watersip_table}];                          %/ just for checking the initial water loss.
%         watersip_stacktable{kk} = watersip_table;                          %/ just for checking the initial water loss.
    end

%     watersip_rr{kk} = [N, q_bc(1), loss(1), f_tot(1), e_tot(1), rr_tot(1)];   %/ final recycling ratios (index = 1 as WaterSip goes from earliest to latest)
    
    if output_rr_map
        watersip_rr{kk} = [N, q_bc(1), loss(1), f_tot(1), e_tot(1), rr_tot(1), rr_L_tot(1), rr_NLL_tot(1), rr_NLO_tot(1)];   %/ final recycling ratios (index = 1 as WaterSip goes from earliest to latest)
    else
        watersip_rr{kk} = [N, q_bc(1), loss(1), f_tot(1), e_tot(1), rr_tot(1)];
    end
    
    %>>>>>>>>>>>>>>>>> optimal trajtime map (use the true destination pt) <<<<<<<<<<<<<<<<<<<%
    avg_optimal_trajtime_map_kk = zeros(length(lon),length(lat));
    traj_startpt_count_map_kk   = zeros(length(lon),length(lat));

    avg_optimal_trajtime_map_kk(ind_lon_des, ind_lat_des) = N;
    traj_startpt_count_map_kk(ind_lon_des, ind_lat_des)   = 1;
    
    mean_optimal_trajtime_map = mean_optimal_trajtime_map  + avg_optimal_trajtime_map_kk;  %/ Reduction in parfor is allowed! (save memory!)
    traj_startpt_count_map    = traj_startpt_count_map     + traj_startpt_count_map_kk;

    %>>>>>>>>>>>>>>>>>>> P_LA map <<<<<<<<<<<<<<<<<<<%
    P_LA_map_kk  = zeros(length(lon),length(lat));
    P_LA = mass_bc*-1*loss(1)/1000;                                %/ Using Stohl and James's method (actually this is E-P)
    P_LA_map_kk(ind_lon_meanxy(1), ind_lat_meanxy(1)) = P_LA;
    P_LA_map = P_LA_map + P_LA_map_kk;                              %/ Reduction in parfor is allowed! (save memory!)
        
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
    
    %>>>>>>>>>>>>>>>>>>> update the water uptake map <<<<<<<<<<<<<<<<<<<%
    %/ combine f and e tgt (count all their contributions)
    rr = nansum([f, e], 2);

%     uptake_map_kk           = zeros(length(lon),length(lat));
    Pm_map_kk               = zeros(length(lon),length(lat));
    ind = find(rr ~= 0);
    for i = 1:length(ind)
        contr = -1 * rr(ind(i)) * loss(1); %/ consider water uptake at different location along the trajectory (g/kg)
        Pm     = mass_bc*contr/1000;
        Pm_map_kk    (ind_lon_meanxy(ind(i)), ind_lat_meanxy(ind(i))) = Pm_map_kk(ind_lon_meanxy(ind(i)), ind_lat_meanxy(ind(i))) + Pm;
        
%         uptake_map_kk(ind_lon_meanxy(ind(i)), ind_lat_meanxy(ind(i))) = uptake_map_kk(ind_lon_meanxy(ind(i)), ind_lat_meanxy(ind(i))) + src_contr;
    end
    Pm_map     = Pm_map     + Pm_map_kk; 
%     uptake_map = uptake_map + uptake_map_kk;     %/ Reduction in parfor is allowed! (save memory!)

    %>>>>>>>>>>>> output 3D matrics in lagrangian (speed-critical step) >>>>>>>>%
%     tic
    if output_LA_3D
        
        %/ NOTE: - For now, we consider ONLY the optimally-tracked traj (1:N-1), NOT the full traj!
        %/       - Retrieve the cell index in the given hs_lon_range, hs_lat_range and hs_z_range.
        [~, ind_box_lon] = find_nearest_val(hs_lon_range, mean_x(1:N-1), 'lon');
        [~, ind_box_lat] = find_nearest_val(hs_lat_range, mean_y(1:N-1));
        [~,   ind_box_z] = find_nearest_val(hs_z_range,   mean_z(1:N-1));         
        
        %/ combine f and e tgt (count all their contributions)
%         uptake = nansum([dq_BL, dq_FT], 2);   %/ here 'uptake' means how much water the parcel gains, but NOT the contribution (it is rr!).
        
        %/ loop only those in the 3D box.
%         I = find(~isnan(ind_box_lon) & ~isnan(ind_box_lat) & ~isnan(ind_box_z));
        I = strfind([~isnan(ind_box_lon) & ~isnan(ind_box_lat) & ~isnan(ind_box_z)]', 1); %/ faster than find()!
        
        freq_3D_kk   = zeros(length(hs_lon_range), length(hs_lat_range), length(hs_z_range));
        rr_3D_kk     = zeros(length(hs_lon_range), length(hs_lat_range), length(hs_z_range));
        contr_3D_kk  = zeros(length(hs_lon_range), length(hs_lat_range), length(hs_z_range));
        loss_3D_kk   = zeros(length(hs_lon_range), length(hs_lat_range), length(hs_z_range));
        u_3D_kk      = zeros(length(hs_lon_range), length(hs_lat_range), length(hs_z_range));
        v_3D_kk      = zeros(length(hs_lon_range), length(hs_lat_range), length(hs_z_range));
        w_3D_kk      = zeros(length(hs_lon_range), length(hs_lat_range), length(hs_z_range));
        BLH_2D_kk    = zeros(length(hs_lon_range), length(hs_lat_range));

        %/ let nans to 0 in loss prior to reduction operation.
        loss_bc = loss;
        loss_bc(isnan(loss_bc)) = 0;

        %/ Compute Traj velocities
        %/ using the actual xyz, then the velocity will correspond to mean_xyz.
        %/ since it's bwd now, we let the ending point to be t + 1.
        %/ Note that uvw has one less element as other data, but the last
        %/ row is always NaN, so doesn't matter.
        [u, v, w] = traj_velocity('lon1', x_bc(2:N),   'lat1', y_bc(2:N),   'z1', z_bc(2:N),...
                                  'lon2', x_bc(1:N-1), 'lat2', y_bc(1:N-1), 'z2', z_bc(1:N-1), 'dt', dt);
                              
        for i = 1:length(I)
            freq_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                        freq_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))   + 1;
          
            rr_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                        rr_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))     + rr(I(i));
                    
            contr_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                        contr_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))  + -1 * rr(I(i)) * loss_bc(1);

            loss_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                        loss_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))   + loss_bc(I(i));       
            
            u_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                        u_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))   + u(I(i));   
            
            v_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                        v_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))   + v(I(i));   
                    
            w_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i))) = ...
                        w_3D_kk(ind_box_lon(I(i)), ind_box_lat(I(i)), ind_box_z(I(i)))   + w(I(i));   
                    
            BLH_2D_kk(ind_box_lon(I(i)), ind_box_lat(I(i))) = ...
                        BLH_2D_kk(ind_box_lon(I(i)), ind_box_lat(I(i))) + mean_BLH(I(i));
        end

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
%     BL_uptake_map_kk = zeros(length(lon),length(lat));
%     BL_Pm_map_kk     = zeros(length(lon),length(lat));
%     
%     ind = find(~isnan(f));
%     for i = 1:length(ind)
%         uptake = -1 * f(ind(i)) * loss(1); %/ consider water uptake at different location along the trajectory (g/kg)
%         Pm     = mass_bc*uptake/1000;
%         
% %         [ind_lon, ind_lat] = point2map('pt_x', mean_x(ind(i)), 'pt_y', mean_y(ind(i)), 'lon', lon, 'lat', lat);
% % 
% %         BL_uptake_map_kk(ind_lon, ind_lat) = BL_uptake_map_kk(ind_lon, ind_lat) + uptake;
% %         BL_Pm_map_kk(ind_lon, ind_lat)     = BL_Pm_map_kk(ind_lon, ind_lat)     + Pm;

%         BL_uptake_map_kk(ind_lon_meanxy(ind(i)), ind_lat_meanxy(ind(i))) = BL_uptake_map_kk(ind_lon_meanxy(ind(i)), ind_lat_meanxy(ind(i))) + uptake;
%         BL_Pm_map_kk    (ind_lon_meanxy(ind(i)), ind_lat_meanxy(ind(i))) = BL_Pm_map_kk    (ind_lon_meanxy(ind(i)), ind_lat_meanxy(ind(i))) + Pm;

%     end
    
    %/ Reduction in parfor is allowed! (save memory!)
%     BL_uptake_map = BL_uptake_map + BL_uptake_map_kk;
%     BL_Pm_map     = BL_Pm_map     + BL_Pm_map_kk;

end
% toc;

% tic;
if kkn ~= 0
    Pm_map     = Pm_map(:,2:end-1)./area;         %/ remove two poles (90N and 90S) as they have no area.
    P_LA_map   = P_LA_map(:,2:end-1)./area;       %/ remove two poles (90N and 90S) as they have no area.
    % BL_Pm_map  = BL_Pm_map(:,2:end-1)./area;      %/ remove two poles (90N and 90S) as they have no area.

    if output_rr_map
        rr_L_tot_map   = rr_L_tot_map./traj_startpt_count_map;    %/ NOTE: it will contain nans, but Retrieve_WSV function can deal with it.
        rr_NLL_tot_map = rr_NLL_tot_map./traj_startpt_count_map;  %/ NOTE: it will contain nans, but Retrieve_WSV function can deal with it.
        rr_NLO_tot_map = rr_NLO_tot_map./traj_startpt_count_map;  %/ NOTE: it will contain nans, but Retrieve_WSV function can deal with it.
        T_vars2     = {'Ntrajtime', 'q', 'loss', 'f_tot', 'e_tot', 'rr_tot', 'rr_L_tot', 'rr_NLL_tot', 'rr_NLO_tot'};
    else
        T_vars2     = {'Ntrajtime', 'q', 'loss', 'f_tot', 'e_tot', 'rr_tot'};
    end

    watersip_rr = cat(1,watersip_rr{:});
    watersip_rr = array2table(watersip_rr, 'VariableNames', T_vars2);

    if ~isempty(optimal_rr)
        mean_optimal_trajtime_map = mean_optimal_trajtime_map./traj_startpt_count_map; %/ it will contain nans, but Retrieve_WSV function can deal with it.
    end
end


LA_3D = [];
if output_LA_3D
    rr_3D     = rr_3D./freq_3D;                 %/ 3hrly traj-mean (will contain nans)
    contr_3D  = contr_3D./freq_3D;              %/ 3hrly traj-mean (will contain nans)
    loss_3D   = loss_3D./freq_3D;               %/ 3hrly traj-mean (will contain nans)
    BLH_2D    = BLH_2D./sum(freq_3D,3);         %/ 3hrly traj-mean (will contain nans)
    
    LA_3D.('freq_3D')      = freq_3D;
    LA_3D.('rr_3D')        = rr_3D;
    LA_3D.('contr_3D')     = contr_3D;
    LA_3D.('loss_3D')      = loss_3D;
    LA_3D.('u_3D')         = u_3D;
    LA_3D.('v_3D')         = v_3D;
    LA_3D.('w_3D')         = w_3D;
    LA_3D.('BLH_2D')       = BLH_2D;
    LA_3D.('hs_lon_range') = hs_lon_range;
    LA_3D.('hs_lat_range') = hs_lat_range;
    LA_3D.('hs_z_range')   = hs_z_range;
    LA_3D.('hs_name')      = hs_name;
end


% fprintf('Aggregation step: %.5f sec\n', toc);

% if traj_rm_jump
%     disp(['NoOf_sliced_traj = ',    num2str(NoOf_sliced_traj)])
%     disp(['NoOf_nonsliced_traj = ', num2str(NoOf_nonsliced_traj)])
% end

disp('!!! Finished WaterSip. !!!')  


end

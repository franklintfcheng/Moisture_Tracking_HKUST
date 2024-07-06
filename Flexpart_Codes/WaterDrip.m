function [Cf_map, watersip_fwd_stacktable, mean_optimal_trajtime_map] = WaterDrip(varargin)

    %/ create a set of valid parameters and their default value
    pnames = {'domfill', 'lon', 'lat', 'dqc', 'RHc', 'BLH_factor', 'area', 'NumWorkers', 'optimal_rr'};  
    dflts  = {[],           [],    [],   [],     [],            1,     [],           [],           []};
    
    [domfill, lon, lat, dqc, RHc, BLH_factor, area, NumWorkers, optimal_rr] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
    
    %%
    %/=====================================================================
    %/      Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 29 Jan 2024
    %/
    %/   Reference: Cheng and Lu (2023)
    %/
    %/=====================================================================

    if optimal_rr > 1 || optimal_rr < 0
        error('optimal_rr must be from 0 to 1.');
    end

    trajtime                  = domfill.trajtime';
    Cf_map                    = zeros(length(lon),length(lat),length(trajtime)-1);      %/ no need to create cell array since they are reduction variables (save memory!!!)
    traj_startpt_count_map    = zeros(length(lon),length(lat));
    watersip_fwd_stacktable   = cell(1, 1); %/ only save this for debugging.
    % watersip_fwd              = cell(size(domfill.traj, 1), 1);
    
    kkn      = size(domfill.traj, 1);
    if kkn == 0
        warning('No trajs input. Skip WaterDrip.')
        mean_optimal_trajtime_map = nan(length(lon),length(lat));  %/ then set the mean trajtime map to be all nans.
        return;
    else
        mean_optimal_trajtime_map = zeros(length(lon),length(lat));  %/ if not, then set to be all zeros for computation.
    end
    
    if kkn > 1
        x   = squeeze(domfill.traj(:,1,:));     %/ ntraj * pos * trajtime -> ntraj * trajtime
        y   = squeeze(domfill.traj(:,2,:)); 
        z   = squeeze(domfill.traj(:,3,:)); 
    else
        x   = squeeze(domfill.traj(:,1,:))';    %/ transpose to avoid dimension order being mistakenedly switched by squeeze() when kkn == 1.
        y   = squeeze(domfill.traj(:,2,:))'; 
        z   = squeeze(domfill.traj(:,3,:))'; 
    end
    % x   = squeeze(domfill.traj(:,1,:));
    % y   = squeeze(domfill.traj(:,2,:));
    % z   = squeeze(domfill.traj(:,3,:)); 
    q     = domfill.q(:,:); 
    dq    = diff(q, 1, 2);                  %/ +ve: water uptake. -ve: water loss.
    BLH   = domfill.BLH(:,:);
    T     = domfill.T(:,:);
    topo  = domfill.topo(:,:);
    mass  = domfill.mass;
    
    % IMPORTANT: Convert to [0 360) for the convenience of the following codes if the lon of traj is in [-179 180].
    if ~isempty(find(x < 0, 1))
       x(x < 0) = x(x < 0) + 360;
    end
    
    flag_shutdown = 0;
    if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if parpool is not open
        parpool('Threads', NumWorkers) %/ much faster
        flag_shutdown = 1;
    end

    % fprintf('*** Running WaterDrip (forward) algorithm for the input trajs... *** \n');
    cnt_earlystop = 0;
    ind_lon_list = [];
    ind_lat_list = [];
    ind_trajtime_list = [];
    Em_list = [];
    parfor kk = 1:kkn 
    % for kk = 1:kkn     %/ for-loop is already the fastest, for it incorporates multi-thread computing
    %     fprintf('kk = %d/%d \n', kk, kkn);  
    
        %====== Set broadcast vars ======%
        trajtime_bc = trajtime;
        x_bc    = x(kk,:)';
        y_bc    = y(kk,:)';
        z_bc    = z(kk,:)';
        q_bc    = q(kk,:)';
        dq_bc   = dq(kk,:)';
        BLH_bc  = BLH(kk,:)';
        T_bc    = T(kk,:)';
        topo_bc = topo(kk,:)';
        
        %/ Identify water loss along the forward trajectory 
        N = length(trajtime_bc); %/ if got sliced, the max track time is shorter; otherwise it should be 241 time steps for 30 days
        loss  = nan(N-1,1);  contr = nan(N-1,1);  res_contr = nan(N-1,1);  prcnt_res_contr = nan(N-1,1); f = nan(N-1,1);   
        stoptime = N-1; %/ assume no early stop.
        for j = 1:N-1                                                           %/ loop from the start
    
            %/ RH at t-1 == RH at j
            [RH2, ~] = RH_Bolton_v2(z_bc(j)+topo_bc(j), q_bc(j)/1000, T_bc(j)); %/ Great water loss occurs with a high RH at t = -1, but extremely low RH at t = 0.. q needs to be in kg/kg
            
    %         fprintf('j = %d, dq = %.4f,  RH = %.4f \n', j, dq_bc(j), RH2)
            if j == 1 && dq_bc(j) > dqc && z_bc(j) <= BLH_factor*BLH_bc(j)                                            
                res_contr(j)       = dq_bc(j);                                
                res_contr_latest   = res_contr(j);                              %/ initialize latest_res_contr
               
            elseif dq_bc(j) < -dqc && RH2 > RHc                                 %/ we consider only the loss < -dqc with RH(t-1) > RHc
                loss(j)            = dq_bc(j);                                  %/ record the sig. moisture loss.
                contr(j)           = f(j-1)*loss(j)*-1;                         %/ use the previous f to compute contr.
                res_contr_latest   = res_contr_latest - contr(j);
                res_contr(j)       = res_contr_latest;                          %/ update the latest residual contribution
                
            elseif dq_bc(j) < 0                                                 %/ consider any 'minor/unsaturated' loss that causes res_contr to drop correspondingly. 
                loss(j)            = dq_bc(j);                                  %/ record the sig. moisture loss.
                loss_of_src        = f(j-1)*loss(j)*-1;                         %/ use the previous f to compute loss of the source's moisture.
                res_contr_latest   = res_contr_latest - loss_of_src;
                res_contr(j)       = res_contr_latest;                          %/ update the latest residual contribution
            end
    
            %/ always update res_contr, prcnt_res_contr and f for next water loss event.
            current_q          = q_bc(j+1);
            res_contr(j)       = res_contr_latest;
            prcnt_res_contr(j) = res_contr_latest/res_contr(1);
            f(j)               = res_contr_latest/current_q;                   
            
            if ~isempty(optimal_rr) && prcnt_res_contr(j) < (1 - optimal_rr)
    %             fprintf('!!! [kk = %d/%d] Forward tracking stops as the remaining ratio of contriution = %.4f < %.4f !!! \n', kk, kkn, prcnt_res_contr(j), optimal_rr);
                cnt_earlystop = cnt_earlystop + 1;                             %/ record # of traj being stopped.
                stoptime = j - 1;                                              %/ record the time step at the stopper.
                break;
            end
        end
        
        if kk <= 10
            T_vars = {'Time','x','y', 'z','BLH', 'q', 'dq0', 'loss', 'contr', 'res_contr', 'prcnt_res_contr', 'f'};
            watersip_table = nan(N, length(T_vars)); %/ recording water uptake from the corrected traj
            watersip_table(:,1)       = trajtime_bc;
            watersip_table(:,2)   = x_bc;  
            watersip_table(:,3)   = y_bc;
            watersip_table(:,4)   = z_bc;
            watersip_table(:,5)   = BLH_bc;
            watersip_table(:, 6)      = q_bc;
            watersip_table(2:end, 7)  = dq_bc;
            watersip_table(2:end, 8)  = loss;
            watersip_table(2:end, 9)  = contr;
            watersip_table(2:end, 10) = res_contr;
            watersip_table(2:end, 11) = prcnt_res_contr;
            watersip_table(2:end, 12) = f;
            watersip_table = array2table(watersip_table, 'VariableNames', T_vars);
            watersip_fwd_stacktable{kk} = watersip_table;
        end
    %     watersip_fwd{kk} = [q_bc(1), loss(1), contr(1), res_contr(1), f(1)];
       
        %>>>>>>>>>>>>>>>>> optimal trajtime map  <<<<<<<<<<<<<<<<<<<%
        avg_optimal_trajtime_map_kk = zeros(length(lon),length(lat));
        traj_startpt_count_map_kk   = zeros(length(lon),length(lat));
    
        [ind_lon, ind_lat] = point2map('pt_x', x_bc(2), 'pt_y', y_bc(2), 'lon', lon, 'lat', lat);   %/ x_bc(2) and y_bc(2) is the start pt in the source region.
        
        avg_optimal_trajtime_map_kk(ind_lon, ind_lat) = stoptime;
        traj_startpt_count_map_kk(ind_lon, ind_lat)   = 1;
        
        mean_optimal_trajtime_map = mean_optimal_trajtime_map  + avg_optimal_trajtime_map_kk;  %/ Reduction in parfor is allowed! (save memory!)
        traj_startpt_count_map    = traj_startpt_count_map     + traj_startpt_count_map_kk;
        
        % %>>>>>>>>>>>>>>>>> Forward Contr. (Cf) map (lon, lat, trajtime)  <<<<<<<<<<<<<<<<<<<%
        % % Cf_map_kk = zeros(length(lon),length(lat));
        % Cf_map_kk = zeros(length(lon),length(lat), N-1);
        % ind       = find(~isnan(contr));
        % for i = 1:length(ind)
        %     Em = mass(kk)*contr(ind(i))/1000;                             %/ Using Stohl and James's method (actually this is E-P)
        % 
        %     [ind_lon, ind_lat] = point2map('pt_x', x_bc(ind(i)), 'pt_y', y_bc(ind(i)), 'lon', lon, 'lat', lat);
        % 
        %     % Cf_map_kk(ind_lon, ind_lat) = Cf_map_kk(ind_lon, ind_lat) + Em;
        %     Cf_map_kk(ind_lon, ind_lat, ind(i)) = Cf_map_kk(ind_lon, ind_lat, ind(i)) + Em;
        % end
        % %/ Reduction 
        % Cf_map = Cf_map + Cf_map_kk;  


        %>>>>>>>>>>>>>>>>> Forward Contr. (Cf) map (lon, lat, trajtime)  <<<<<<<<<<<<<<<<<<<%
        ind = find(~isnan(contr));
        for i = 1:length(ind)
            Em = mass(kk)*contr(ind(i))/1000;                             %/ Using Stohl and James's method (actually this is E-P)
            
            [ind_lon, ind_lat] = point2map('pt_x', x_bc(ind(i)), 'pt_y', y_bc(ind(i)), 'lon', lon, 'lat', lat);
        
            %/ Reduction (appending), the ordering doesn't matter. 
            %/ Avoid making a huge 3D matrix in a parfor 
            %/    -> extremely slow and heavy
            ind_lon_list      = [ind_lon_list; ind_lon];
            ind_lat_list      = [ind_lat_list; ind_lat];
            ind_trajtime_list = [ind_trajtime_list; ind(i)];
            Em_list           = [Em_list; Em];
            % Cf_map_kk(ind_lon, ind_lat) = Cf_map_kk(ind_lon, ind_lat) + Em;
            % Cf_map(ind_lon, ind_lat, ind(i)) = Cf_map(ind_lon, ind_lat, ind(i)) + Em;
        end
        %/ Reduction 
        % Cf_map = Cf_map + Cf_map_kk;  
    end
    
    %/ Can't use parfor, but a simple for-loop is fast enough
    L = length(ind_lon_list);
    for i = 1:L  
        X = ind_lon_list(i);
        Y = ind_lat_list(i);
        TJ = ind_trajtime_list(i);
        Em = Em_list(i);
        Cf_map(X,Y,TJ) = Cf_map(X,Y,TJ) + Em;
    end

    if ~isempty(optimal_rr)
        mean_optimal_trajtime_map = mean_optimal_trajtime_map./traj_startpt_count_map;  %/ 0/0 = nan
        mean_optimal_trajtime_map(isnan(mean_optimal_trajtime_map)) = 0;                %/ set nan to zeros.
    end
    
    Cf_map     = Cf_map(:,2:end-1,:)./area;         %/ convert into kg/m2 per 3h, remove two poles (90N and 90S) as they have no area.
    % watersip_fwd = cat(1,watersip_fwd{:});

    %/ Shut down the thread-based parpool initiated within the function
    %/ to avoid conflicting with the need of process-based parpool outside
    if flag_shutdown
        poolobj = gcp('nocreate');
        delete(poolobj);
    end

    disp('!!! WaterDrip (forward) is completed !!!')  
    fprintf('*** %d/%d of fwd tracking stopped when %.3G%% of initial water gain has been precipitated *** \n', cnt_earlystop, kkn, optimal_rr*100)

end

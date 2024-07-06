%%
function domfill = slct_traj_parfor(varargin)

    %/ create a set of valid parameters and their default value
    pnames = {'domfill', 'trajdirect', 'RHc', 'dqc', 'str_RHc_dqc', 'str_remark', 'RHc_map_lon', 'RHc_map_lat', 'slct_one_date_dt', 'BLH_factor'}; 
    dflts  = cell(length(pnames),1);
    [domfill, trajdirect, RHc, dqc, str_RHc_dqc, str_remark, RHc_map_lon, RHc_map_lat, slct_one_date_dt, BLH_factor] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %%
    %/ dqc                min dq for prcp to occur [g/kg]
    %/ RHc                min RH for prcp to occur [%] It can be a single value or a matrix.
    %/ RHc_map_lon        useful only when RHc is a matrix.
    %/ RHc_map_lat        useful only when RHc is a matrix.

    %=========================================================================%
    %===== Pre-select trajs that satisfy RHc and dqc (which imply prcp) ======%
    %=========================================================================%
    % EminusP_LA_map = zeros(length(lon),length(lat));

    if isempty(domfill(1).slct_id_list)
        warning('slct_traj_parfor: No traj id is input. Skipping...')
        return
    end
    if isempty(BLH_factor)   BLH_factor = 1;  end

    %=========================================================================%
    %/ NOTE: The particle's z-pos is measured by *height above terrain*
    %/       You'll need to add the topo height to its z-pos to restore
    %/       the actual altitude!
    %=========================================================================%
    traj0  = cat(3, domfill(1:2).traj); 
    x0     = squeeze(traj0(:,1,:));
    y0     = squeeze(traj0(:,2,:));
    z0     = squeeze(traj0(:,3,:));   
    % mean_x = squeeze(mean(traj0(:,1,:),3));
    % mean_y = squeeze(mean(traj0(:,2,:),3));
    % mean_z = mean(z0, 2);

    q0       = cat(2, domfill(1:2).q);
    T0       = cat(2, domfill(1:2).T); 
    BLH0     = cat(2, domfill(1:2).BLH);
    % mean_BLH = mean(BLH0, 2);

    topo0    = squeeze(cat(3, domfill(1:2).topo)); 
    % mass     = domfill(1).mass;
    % a = [z0(:,2), topo0(:,2), z0(:,2)-topo0(:,2)]; %/ checking

    if ismember(trajdirect, {'bwd'})
    %     [RH2_old, ~] = RH_Bolton(z0(:,2), q0(:,2)/1000, T0(:,2));                  %/ Great water loss occurs with a high RH at t = -1, but extremely low RH at t = 0.. q needs to be in kg/kg
    %     [RH2_old, ~] = RH_Bolton_v2(z0(:,2), q0(:,2)/1000, T0(:,2));                  %/ Great water loss occurs with a high RH at t = -1, but extremely low RH at t = 0.. q needs to be in kg/kg
        
        if isempty(str_remark)
            [RH2, ~] = RH_Bolton_v2(z0(:,1)+topo0(:,1), q0(:,1)/1000, T0(:,1));  
        elseif contains(str_remark, 'RH2')
            [RH2, ~] = RH_Bolton_v2(z0(:,2)+topo0(:,2), q0(:,2)/1000, T0(:,2));                  %/ Great water loss occurs with a high RH at t = -1, but extremely low RH at t = 0.. q needs to be in kg/kg
        else
            error('Invalid ''str_remark == %s''!', str_remark);
        end
        dq0 = -1*diff(q0, 1, 2);                                               %/ +ve: water uptake. -ve: water loss.
    
    %     %===================================================================%
        if length(RHc) > 1   %/ RHc is a matrix (RHc_map)
            if size(RHc, 3) == 4  %/ seasonal RHc map
                if any(ismember([3:5], slct_one_date_dt.Month)) %/ MAM
                    RHc_bc = RHc(:,:,1);
                elseif any(ismember([6:8], slct_one_date_dt.Month)) %/ JJA
                    RHc_bc = RHc(:,:,2);
                elseif any(ismember([9:11], slct_one_date_dt.Month)) %/ SON
                    RHc_bc = RHc(:,:,3);
                elseif any(ismember([12,1,2], slct_one_date_dt.Month)) %/ DJF
                    RHc_bc = RHc(:,:,4);
                else
                    error('check your code!')
                end
            else
                RHc_bc = RHc;
            end
            
            %/ Retrieve the grid cell index at the destination point.
            [ind_lon_des, ind_lat_des] = point2map('pt_x', x0(:,1), 'pt_y', y0(:,1), 'lon', RHc_map_lon, 'lat', RHc_map_lat);  
    
            %/ This simple loop is already fast
            cond_RH_RHc = zeros(length(ind_lon_des), 1);
            for i = 1:length(ind_lon_des)
                if isnan(ind_lon_des(i)) || isnan(ind_lat_des(i))
                    %/ There are nans in ind_lon_des / ind_lat_des
                    %/ Since the particle may reach the north / south pole when lat is only from -89 to 89. 
                    continue;
                end
                cond_RH_RHc(i) = (RH2(i) > RHc_bc(ind_lon_des(i), ind_lat_des(i)));
                
                %/ Double-check
    %             fprintf('i = %d, RH2(i)=%.2f%%, RHc = %.2f%%, cond = %d \n', i, RH2(i), RHc_bc(ind_lon_des(i), ind_lat_des(i)), cond_RH_RHc(i));
            end
            traj_prcp_list = find(cond_RH_RHc == 1 & dq0 < -1*dqc);  
    
            %/ Double-check
    %         for j = 1:100 %length(traj_prcp_list)
    %             fprintf('i = %d, RH2(i)=%.2f%%, dq0 = %.2f \n', traj_prcp_list(j), RH2(traj_prcp_list(j)), dq0(traj_prcp_list(j)));
    %         end
    
        else  %/ RHc is a single value
            traj_prcp_list = find(RH2 > RHc & dq0 < -1*dqc);   
        end
    
        for k = 1:2
            domfill(k).q    = domfill(k).q   (traj_prcp_list);
            domfill(k).traj = domfill(k).traj(traj_prcp_list,:);
            domfill(k).BLH  = domfill(k).BLH (traj_prcp_list);
            domfill(k).T    = domfill(k).T   (traj_prcp_list);
            domfill(k).topo = domfill(k).topo (traj_prcp_list);
        end
        %/ update other fields
    %     domfill(1).T            = domfill(1).T           (traj_prcp_list);   
        domfill(1).slct_id_list = domfill(1).slct_id_list(traj_prcp_list);
        domfill(1).mass         = domfill(1).mass        (traj_prcp_list);
        domfill(1).RH2          = RH2(traj_prcp_list);  %/ New! RH2 when dq < -dqc.
        
        fprintf('*** [%s mode] Detected %d/%d trajs starting w/ %s *** \n',...
                trajdirect, length(traj_prcp_list), length(RH2), str_RHc_dqc)
    end

    if ismember(trajdirect, {'fwd'})
        dq0 = diff(q0, 1, 2);                                                  %/ +ve: water uptake. -ve: water loss.
    %     traj_src_list = find(dq0 > dqc & mean_z <= BLH_factor*mean_BLH);       %/ select trajs that intially gain moisture (> dqc) within BL.
        traj_src_list = find(dq0 > dqc & z0(:,1) <= BLH_factor*BLH0(:,1));       %/ select trajs that intially gain moisture (> dqc) within BL.

        for k = 1:2
            domfill(k).q    = domfill(k).q   (traj_src_list);
            domfill(k).traj = domfill(k).traj(traj_src_list,:);
            domfill(k).BLH  = domfill(k).BLH (traj_src_list);
            domfill(k).T    = domfill(k).T   (traj_src_list);
            domfill(k).topo = domfill(k).topo(traj_src_list);
        end
        %/ update other fields
    %     domfill(1).T            = domfill(1).T           (traj_src_list);  
    %     domfill(2).T            = [];                                          %/ no use.
        domfill(1).slct_id_list = domfill(1).slct_id_list(traj_src_list);
        domfill(1).mass         = domfill(1).mass        (traj_src_list);

        fprintf('*** [%s mode] Detected %d/%d trajs starting w/ dq > %.2f g/kg and z <= %.2f*BLH *** \n',...
                trajdirect, length(traj_src_list), length(dq0), dqc, BLH_factor)
    end

end
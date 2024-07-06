function domfill = slct_traj(varargin)

%/ create a set of valid parameters and their default value
pnames = {'domfill', 'trajdirect', 'RHc', 'dqc'}; 
dflts  = {[], [], [], []};
%/ parse function arguments
[domfill, trajdirect, RHc, dqc ] = internal.stats.parseArgs(pnames, dflts, varargin{:});

%/ dqc                min dq for prcp to occur [g/kg]
%/ RHc                min RH for prcp to occur [%]

%=========================================================================%
%===== Pre-select trajs that satisfy RHc and dqc (which imply prcp) ======%
%=========================================================================%
z0 = squeeze(domfill.traj(:, 3, 1)); 
q0 = domfill.q(:, 1:2);    %/ convert from g/kg to kg/kg for calculating RH
T0 = domfill.T(:, 1); 

if ismember(trajdirect, {'bwd'})
    [RH, ~] = RH_Bolton(z0, q0(:,1)/1000, T0); %/ q needs to be in kg/kg
    
    dq0 = -1*diff(q0, 1, 2);                         %/ +ve: water uptake. -ve: water loss.
    traj_prcp_list = find(RH > RHc & dq0 < -1*dqc);  %/ Assume prcp occurs when RH > 80%, dq0 > dqc
    domfill.q   = domfill.q(traj_prcp_list,:);  
    
    fprintf('*** [%s mode] Detected %d/%d trajs starting w/ RH > %.2f%% and dq < -%.2f g/kg *** \n',...
            trajdirect, length(traj_prcp_list), length(RH), RHc, dqc)
        
else if ismember(trajdirect, {'fwd'})
    
    dq0 = diff(q0, 1, 2);                %/ +ve: water uptake. -ve: water loss.
    traj_prcp_list = find(dq0 > dqc);    %/ select trajs that gain moisture (> dqc) at the initial grid compared to its earlier state.
    domfill.dq0 = dq0(traj_prcp_list);  %/ store this moisture gain for further analysis.
    domfill.q   = domfill.q(traj_prcp_list,2);  %/ discard the q before start. 
    
    fprintf('*** [%s mode] Detected %d/%d trajs starting w/ dq > %.2f g/kg *** \n',...
            trajdirect, length(traj_prcp_list), length(dq0),dqc)
end
end

%/ Output domfill with the selected trajectories 
domfill.traj = domfill.traj(traj_prcp_list,:,:);
domfill.BLH  = domfill.BLH (traj_prcp_list,:);
domfill.T    = domfill.T   (traj_prcp_list,:);   
    
%/ update id list and the corresp. masses
domfill.slct_id_list = domfill.slct_id_list(traj_prcp_list);
domfill.mass         = domfill.mass(traj_prcp_list);
end
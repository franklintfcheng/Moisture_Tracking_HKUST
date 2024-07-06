%% Read WaterSip Parameterization

function [str_RHc_dqc, ldirect, str_remark, optimal_tracking, traj_rm_jump, BLH_factor, str_src_z,...
          str_years, WSV_lon, WSV_lat, maxtraj_day, str_BLH_factor, str_optimal, str_traj_rm_jump, dt_slct, data_folder,...
          plotting_folder, WSV_dir, str_expmntinfo, basin_catalog, str_domain, ...
          ldirect_fwd, str_remark_fwd, maxtraj_day_fwd, WSV_dir_fwd, str_expmntinfo_fwd,...
          global_basin_2D, global_basin_bndry_list, global_basin_name_list, global_basin_centr_list,...
          cond_TP, TP_bndry, TP_centr, dataset] = param_watersip(varargin)

    pnames = {'expmnt', 'output_res', 'dt', 'year_list', 'optimal_rr', 'RHc_dqc_scheme', 'from_basin', 'dataset', 'masterfolder'};  
    dflts  = cell(1, length(pnames));
    [          expmnt,   output_res,   dt,   year_list,   optimal_rr,   RHc_dqc_scheme,   from_basin,   dataset,   masterfolder] = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
    %%
    
    %/ Derive forwad trajectories from domain-fill and the uptake map
    fprintf('*** NOTE: For the derivation of trajs and prcp footprint maps, see ''moisture_tracking.m''. ***\n')

    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    [~, ~, ~, ~, str_RHc_dqc] = read_RHc_dqc_scheme(RHc_dqc_scheme);
    ldirect = 'bwd'; str_remark = '_RH2'; optimal_tracking = 1; traj_rm_jump = 0; BLH_factor = 1; str_src_z = '_nozbound';  % rm_when_dz = []; 

    [str_years, WSV_lon, WSV_lat, maxtraj_day, str_BLH_factor, str_optimal, str_traj_rm_jump, dt_slct, data_folder,...
     plotting_folder, WSV_dir, str_expmntinfo, basin_catalog, str_domain] = ...
        set_flexpart_fpath('expmnt', expmnt, 'output_res', output_res, 'dt', dt, 'year_list', year_list, 'ldirect', ldirect, 'str_remark', str_remark, 'str_RHc_dqc', str_RHc_dqc,...
                           'optimal_tracking', optimal_tracking, 'optimal_rr', optimal_rr, 'traj_rm_jump', traj_rm_jump, 'BLH_factor', BLH_factor, 'str_src_z', str_src_z,...
                           'masterfolder', masterfolder, 'from_basin', from_basin);

    %/ By default, also load str_expmntinfo_fwd (can later plot bwd and fwd products tgt)
    %/ NOTE: '_RH2_correct' is just to indicate rectified code. Don't mind.
    if RHc_dqc_scheme == 1
        ldirect_fwd = 'fwd'; str_remark_fwd = '_RH2_correct'; optimal_tracking = 1; traj_rm_jump = 0;  BLH_factor = 1; str_src_z = '_nozbound'; % rm_when_dz = [];
        [~, ~, ~, maxtraj_day_fwd, ~, ~, ~, ~, ~, ~, WSV_dir_fwd, str_expmntinfo_fwd, ~, str_domain] = ...
            set_flexpart_fpath('expmnt', expmnt, 'output_res', output_res, 'dt', dt, 'year_list', year_list, 'ldirect', ldirect_fwd, 'str_remark', str_remark_fwd, 'str_RHc_dqc', str_RHc_dqc,...
                               'optimal_tracking', optimal_tracking, 'optimal_rr', optimal_rr, 'traj_rm_jump', traj_rm_jump, 'BLH_factor', BLH_factor, 'str_src_z', str_src_z,...
                               'masterfolder', masterfolder, 'from_basin', from_basin);
    else
        ldirect_fwd = []; 
        str_remark_fwd = []; 
        maxtraj_day_fwd     = [];  %/ return empty (for now)
        WSV_dir_fwd         = [];
        str_expmntinfo_fwd  = [];
%         str_domain          = [];
    end

    %/ Load the global basins first (save time for computing SR network if needed)
    recompute = 0;
    slct_reg  = 'hydrosheds_TPbasins_oceans';

    [global_basin_2D, global_basin_bndry_list, global_basin_name_list, ~, global_basin_centr_list] = reg_extractor('slct_reg', slct_reg, 'lon', WSV_lon,...
                                               'lat', WSV_lat(2:end-1), 'data_folder', data_folder, 'savemat', 1, 'recompute', recompute, 'verbose', 1);

    %/ Load the TP basin boundary
    slct_reg = 'TP';
    [cond_TP, TP_bndry, ~, ~, TP_centr] = reg_extractor('slct_reg', slct_reg, 'lon', WSV_lon, 'lat', WSV_lat(2:end-1),...
                                                'data_folder', data_folder, 'savemat', 1, 'recompute', recompute, 'verbose', 1);
    cond_TP(isnan(cond_TP)) = 0;
    cond_TP = logical(cond_TP);

    %/ Load Land Cover Type from CERA-20C
    if contains(expmnt, '_CERA_')
        LC_name = {'HiVeg', 'LoVeg'};
        for k = 1:length(LC_name)
            if isfield(dataset, LC_name{k})       %/ skip it if it has been loaded.
                disp(['!!! ', LC_name{k}, ' has already been loaded. Skip data retrieval. !!!'])
            else
                select_field = 'daily';
                domain_str = '_global';

                loadfile = strcat(data_folder, LC_name{k}, '_', select_field, domain_str, '.mat');
                disp(['Loading ', loadfile, ' ...']);
                load(loadfile, 'S');

                fld = fieldnames(S);
                for f = 1:length(fld)
                    dataset.(LC_name{k}).(fld{f}) = S.(fld{f});
                end
                clear S;
            end
        end
    end
    %/ HiVeg and LoVeg can overlap in some grids.
    % dataset.('VegType').daily = dataset.('HiVeg').daily + dataset.('LoVeg').daily;
    % unique(reshape(dataset.('VegType').daily,[], 1))
    fprintf('=== NOTE: RHc_dqc_scheme = %d, from_basin = %d. Period: %d-%d. ===\n', RHc_dqc_scheme, from_basin, year_list(1), year_list(end))

end
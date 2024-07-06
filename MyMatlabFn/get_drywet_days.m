%%
function [dataset, dry_days, wet_days, normal_days] = get_drywet_days(varargin)
    
    % create a set of valid parameters and their default value
    pnames = {'year_list', 'bndry_data', 'bndry_name', 'dataset', 'prcp_datasource', 'set_prcntl',  'mth_list',  'plot_hist', 'savefig',  'plotting_folder'};
    dflts  = {         [],           [],           [],        [],            'CERA',   [0.1, 0.9],           0,            0,         0,                 []};
    [          year_list,   bndry_data,   bndry_name,   dataset,   prcp_datasource,    set_prcntl,    mth_list,    plot_hist,   savefig,    plotting_folder] ...
                        = internal.stats.parseArgs(pnames, dflts, varargin{:});

    str_mth = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'MAM', 'JJA', 'SON', 'DJF'};

    %/ load prcp first
    if isfield(dataset, 'P')                                   %/ skip it if it has been loaded.
        disp(['!!! P has already been loaded. Skip data retrieval. !!!'])
    else
        if isequal(prcp_datasource, 'CERA')
            loadfile = '/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/P_daily_global_1971-2010.mat';
            disp(['''get_drywet_days'': Loading ', loadfile, ' ...']);
            load(loadfile, 'S');

            fld = fieldnames(S);
            for f = 1:length(fld)
                dataset.('P').(fld{f}) = S.(fld{f});
            end
        end
    end
    

    %/ find data and area within the boundaries
    [map_data_in, area_in, ~] = which_in_boundary('bndry_data_inpoly2', bndry_data, 'map_data', dataset.('P').daily,...
                                                  'lon', dataset.('P').lon, 'lat', dataset.('P').lat, 'output_1D_or_2D', '1D');

    for mth = mth_list

        if mth >= 13 && mth <= 16
            slct_dates = date_array_gen('years', year_list, 'season', mth-12,...   %/ omit 20101231 if using season function.
                                       'output_date_format', 'yyyymmdd');
            str_season = strcat('_', str_mth{mth});

        elseif mth >= 1 && mth <= 12
            slct_dates = date_array_gen('years', year_list, 'season', mth,...   %/ omit 20101231 if using season function.
                                       'output_date_format', 'yyyymmdd');
            str_season = strcat('_', str_mth{mth});

        else
            slct_dates = dataset.('P').date_yyyymmdd_AllYr;
            str_season = [];
        end  

        ind_dates       = findismember_loop(dataset.('P').date_yyyymmdd_AllYr, slct_dates);
        P_subset_dates  = dataset.('P').date_yyyymmdd_AllYr(ind_dates);
        P_AWA           = [sum(map_data_in(:,ind_dates).*area_in, 1)/sum(area_in)]';  %/ AWA = area-weighted average. mm/day

        %/ Plot histogram
        if plot_hist
            close all;
            titlename = strrep(strcat(bndry_name, str_season), '_', ' ');
            if savefig      
                if isempty(plotting_folder)   error('empty plotting_folder!!');  end
                savepath = strrep(strcat(plotting_folder, 'hist_P_AWA', '_', titlename), ' ', '_');
            else
                savepath = [];                                                                  
            end

            plot_hist_cdf('y', P_AWA, 'y_lim', [], 'set_prcntl', set_prcntl, 'basic_stat', 1,...
                          'zoomin_mode', 0, 'draw_cdf', 0, 'titlename', titlename, 'savepath', savepath)
            pause(0.5)
        end

        %/ derive dry and wet days
        P_lower_prcntl = quantile(P_AWA, set_prcntl(1));
        P_upper_prcntl = quantile(P_AWA, set_prcntl(2));


        ind_dry_days = find(P_AWA <= P_lower_prcntl);
        ind_wet_days = find(P_AWA >= P_upper_prcntl);
        ind_normal_days = setdiff(1:length(P_AWA), [ind_dry_days; ind_wet_days])';

        %/ obtain the exact dates in yyyymmdd (to minimize mistakes)
        dry_days    = P_subset_dates(ind_dry_days);
        wet_days    = P_subset_dates(ind_wet_days);
        normal_days = P_subset_dates(ind_normal_days);

        %/ no need to savemat. data is not large.
%         if savemat
%             str_set_prcntl = sprintf('%dp-%dp', set_prcntl(1)*100, set_prcntl(2)*100);
%             matfilename = strcat(data_folder, 'P_AWA_drywetdays_', str_set_prcntl, '_', bndry_name, str_season, '.mat');
%             fprintf('*** Saving dry wet normal days into %s ***\n', matfilename);
%             save(matfilename, 'dry_days', 'wet_days', 'normal_days');
%         end

    end
end
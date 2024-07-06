%%
function [dataset] = calc_WAF(varargin)

    %/ create a set of valid parameters and their default value
    pnames = {'dataset', 'plev', 'perturb_fld', 'str_years', 'prcssd_data_folder', 'savemat'}; 
    dflts  = cell(length(pnames), 1);
    [           dataset,   plev,   perturb_fld,   str_years,   prcssd_data_folder,  savemat] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    
    %%
    %/ Input: 
    %/  perturbation_field:    'pentad_anom', 'dailyAnom_PMC' 
    
    fprintf('*** Running calc_WAF...***\n')
    dataset.(sprintf('uWAF%d', plev))      = [];   %/ remove leftover WAF fields from last run.
    dataset.(sprintf('vWAF%d', plev))      = [];   %/ remove leftover WAF fields from last run.
        
    for k = 1:length(perturb_fld)
        perturb_fld_bc = perturb_fld{k};
        
        windspharm_vars    = {sprintf('psi%d',        plev);...
                              sprintf('d_psi%d_dx',   plev);...
                              sprintf('d_psi%d_dy',   plev);...
                              sprintf('d_psi%d_dxdy', plev);...
                              sprintf('d_psi%d_dxdx', plev);...
                              sprintf('d_psi%d_dydy', plev)};

        %/ Read psia and its derivatives (calc_derivatives_psia)
        for j = 1:length(windspharm_vars)
            fprintf('*** Loading %s field of %s ...***\n', perturb_fld_bc, windspharm_vars{j})

            load_ncfilename = strcat(prcssd_data_folder, windspharm_vars{j}, '_', perturb_fld_bc, '_', str_years, '.nc');
    %         ncdisp(load_ncfilename)

            dataset_temp.(windspharm_vars{j})                         = [];
            dataset_temp.(windspharm_vars{j}).(perturb_fld_bc) = double(ncread(load_ncfilename, windspharm_vars{j}));
            dataset_temp.(windspharm_vars{j}).lon                     = ncread(load_ncfilename, 'lon');
            dataset_temp.(windspharm_vars{j}).lat                     = ncread(load_ncfilename, 'lat');
            dataset_temp.(windspharm_vars{j}).date     = ncread(load_ncfilename, 'time');

            date_AllYr = dataset_temp.(windspharm_vars{j}).date;

            %/ Basic check
            if ~isempty(find(isnan(dataset_temp.(windspharm_vars{j}).(perturb_fld_bc))))  error('%s data contains nan!!', windspharm_vars{j});   end
        end

        %/ get the corresponding mean field of U,V
        uv_vars = {sprintf('u%d', plev); sprintf('v%d', plev)};
        if isequal(perturb_fld_bc, 'dailyAnom_PMC')

            %/ treat the final WAF as daily, not daily anom. (see if this is ok)
            WAF_field     = 'daily';
            uv_clim_field = 'PMC';
            for j = 1:length(uv_vars)
                fprintf('*** Loading %s field of %s ...***\n', uv_clim_field, uv_vars{j})

                load_ncfilename = strcat(prcssd_data_folder, uv_vars{j}, '_', uv_clim_field, '_', str_years, '.nc');
                dataset_temp.(uv_vars{j})                  = [];
                dataset_temp.(uv_vars{j}).(perturb_fld_bc)    = double(ncread(load_ncfilename, uv_vars{j}));
                dataset_temp.(uv_vars{j}).lon              = ncread(load_ncfilename, 'lon');
                dataset_temp.(uv_vars{j}).lat              = ncread(load_ncfilename, 'lat');
                dataset_temp.(uv_vars{j}).date             = ncread(load_ncfilename, 'time');

                %/ get PMC
                [~, dataset_temp.(uv_vars{j}).(uv_clim_field)] = daily2PMC_anom('daily_data', dataset_temp.(uv_vars{j}).daily,...
                                                                              'date',      dataset_temp.(uv_vars{j}).date);
            end

        elseif ismember(perturb_fld_bc, {'pentad_anom', 'pentad_anom_AC', 'pentad_anom_CISO'})

            %/ treat the final WAF as pentad, not pentad anom. (see if this is ok)
            if isequal(perturb_fld_bc, 'pentad_anom')
                WAF_field      = 'pentad';
                WAF_clim_field = 'pentad_clim';
                uv_clim_field  = 'pentad_clim';

            elseif isequal(perturb_fld_bc, 'pentad_anom_AC')
                WAF_field      = 'pentad_AC';
                WAF_clim_field = 'pentad_clim_AC';
                uv_clim_field  = 'pentad_clim_AC';

            elseif isequal(perturb_fld_bc, 'pentad_anom_CISO')
                WAF_field      = 'pentad_CISO';
                WAF_clim_field = 'pentad_clim_CISO';
                uv_clim_field  = 'pentad_clim_CISO';
            end

            %/ load clim. U,V
            for j = 1:length(uv_vars)
                fprintf('*** Loading %s field of %s ...***\n', uv_clim_field, uv_vars{j})

                load_ncfilename = strcat(prcssd_data_folder, uv_vars{j}, '_', uv_clim_field, '_', str_years, '.nc');
                dataset_temp.(uv_vars{j})                 = [];
                dataset_temp.(uv_vars{j}).(uv_clim_field) = double(ncread(load_ncfilename, uv_vars{j}));
                dataset_temp.(uv_vars{j}).lon             = ncread(load_ncfilename, 'lon');
                dataset_temp.(uv_vars{j}).lat             = ncread(load_ncfilename, 'lat');
                dataset_temp.(uv_vars{j}).date_yyyyptd    = ncread(load_ncfilename, 'time');
            end
        else
            error('code not set yet!');
        end

        dataset.(sprintf('uWAF%d', plev)).date = date_AllYr;
        dataset.(sprintf('vWAF%d', plev)).date = date_AllYr;

        %/ Calculate WAF
        if isequal(perturb_fld_bc, 'dailyAnom_PMC')
            mod_number = 1e4;  %/ 8-digit yyyymmdd for daily
        else
            mod_number = 1e2;  %/ 6-digit yyyypt for pentad
        end
        date_OneYr = unique(mod(date_AllYr, mod_number)); %/ yyyypt
        ndate = length(date_OneYr);

        dataset.(sprintf('uWAF%d', plev)).lon = dataset_temp.(sprintf('psi%d', plev)).lon;
        dataset.(sprintf('uWAF%d', plev)).lat = dataset_temp.(sprintf('psi%d', plev)).lat;

        dataset.(sprintf('vWAF%d', plev)).lon = dataset_temp.(sprintf('psi%d', plev)).lon;
        dataset.(sprintf('vWAF%d', plev)).lat = dataset_temp.(sprintf('psi%d', plev)).lat;

        dataset.(sprintf('uWAF%d', plev)).(WAF_field) = nan(size(dataset_temp.(sprintf('psi%d', plev)).(perturb_fld_bc)));
        dataset.(sprintf('vWAF%d', plev)).(WAF_field) = nan(size(dataset_temp.(sprintf('psi%d', plev)).(perturb_fld_bc)));

        dataset.(sprintf('uWAF%d', plev)).(WAF_clim_field) = nan(size(dataset_temp.(uv_vars{j}).(uv_clim_field)));
        dataset.(sprintf('vWAF%d', plev)).(WAF_clim_field) = nan(size(dataset_temp.(uv_vars{j}).(uv_clim_field)));

        for d = 1:ndate
            ind_date        = find(mod(date_AllYr, mod_number) == date_OneYr(d));
            psia            = dataset_temp.(sprintf('psi%d',        plev)).(perturb_fld_bc)(:,:,ind_date);
            d_psia_dx       = dataset_temp.(sprintf('d_psi%d_dx',   plev)).(perturb_fld_bc)(:,:,ind_date);
            d_psia_dy       = dataset_temp.(sprintf('d_psi%d_dy',   plev)).(perturb_fld_bc)(:,:,ind_date);
            d_psia_dxdy     = dataset_temp.(sprintf('d_psi%d_dxdy', plev)).(perturb_fld_bc)(:,:,ind_date);
            d_psia_dxdx     = dataset_temp.(sprintf('d_psi%d_dxdx', plev)).(perturb_fld_bc)(:,:,ind_date);
            d_psia_dydy     = dataset_temp.(sprintf('d_psi%d_dydy', plev)).(perturb_fld_bc)(:,:,ind_date);
            u_mean          = dataset_temp.(sprintf('u%d',          plev)).(uv_clim_field)(:,:,d);
            v_mean          = dataset_temp.(sprintf('v%d',          plev)).(uv_clim_field)(:,:,d);
            uv_mean_speed   = sqrt(u_mean.^2 + v_mean.^2);

            %/ compute WAF 
            dataset.(sprintf('uWAF%d', plev)).(WAF_field)(:,:,ind_date) = 1./(2*uv_mean_speed).*(u_mean .* (d_psia_dx.^2           - psia .* d_psia_dxdx) ...
                                                             + v_mean .* (d_psia_dx .* d_psia_dy - psia .* d_psia_dxdy));

            dataset.(sprintf('vWAF%d', plev)).(WAF_field)(:,:,ind_date) = 1./(2*uv_mean_speed).*(u_mean .* (d_psia_dx .* d_psia_dy - psia .* d_psia_dxdy) ...
                                                             + v_mean .* (d_psia_dy.^2           - psia .* d_psia_dydy));

            %/ compute climatology of WAF                                             
            dataset.(sprintf('uWAF%d', plev)).(WAF_clim_field)(:,:,d) = mean(dataset.(sprintf('uWAF%d', plev)).(WAF_field)(:,:,ind_date), 3);
            dataset.(sprintf('vWAF%d', plev)).(WAF_clim_field)(:,:,d) = mean(dataset.(sprintf('vWAF%d', plev)).(WAF_field)(:,:,ind_date), 3);

        end

        %/ Final check
        if ~isempty(find(isnan(dataset.(sprintf('uWAF%d', plev)).(WAF_field))))         error('data contains nan!!');  end
        if ~isempty(find(isnan(dataset.(sprintf('vWAF%d', plev)).(WAF_field))))         error('data contains nan!!');  end
        if ~isempty(find(isnan(dataset.(sprintf('uWAF%d', plev)).(WAF_clim_field))))    error('data contains nan!!');  end
        if ~isempty(find(isnan(dataset.(sprintf('vWAF%d', plev)).(WAF_clim_field))))    error('data contains nan!!');  end
    end
    
    if savemat
        str_WAF = strcat({'uWAF','vWAF'}, num2str(plev));
        for i = 1:length(str_WAF)
            %/ Remove pentad, keep pentad_clim.
            S = dataset.(str_WAF{i});  
            filename = strcat(prcssd_data_folder, str_WAF{i}, '_pentad_', str_years, '.mat');
            disp(['Writing into ', filename, ' ...'])
            save(filename, 'S','-v7.3');
        end
    end
end
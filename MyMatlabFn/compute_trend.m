function [trend_data, trend_data_sig, ins_or_acc, pval_data, intercept_data, CI95_data, Rsquared_data, Rsquared_adj_data, RMSE_data] = compute_trend(varargin)

    pnames = {'data_daily', 'data_yearly', 'ins_or_acc', 'lon', 'lat', 'data_dates', 'data_name', 'year_list', 'mth',  'trend_test', 'trend_alpha',  'save_trend', 'recompute_trend', 'data_folder', 'NumWorkers'}; 
    dflts  = {          [],            [],           [],    [],    [],           [],          [],          [],     0,            [],          0.05,             0,                 0,            [],           60 };
    [           data_daily,   data_yearly,   ins_or_acc,   lon,   lat,   data_dates,   data_name,   year_list,   mth,    trend_test,   trend_alpha,    save_trend,   recompute_trend,   data_folder,   NumWorkers  ] ...
               = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    
    %/===== Description =====%
    %/  This function uses simple linear regression to estimate the long-term trends.
    %/  'mth' = 0:      compute trends in annual mean time series.
    %/  'mth' = 1:12:   compute trends in monthly mean time series.
    %/  'mth' = 13:16:  compute trends in seasonal mean (MAM/JJA/SON/DJF) time series.
    
    %/  Input Data:   'data_daily' (2D or 3D), will be converted into yearly
    %/                data based on 'mth'
    %/
    %/                'data_yearly' (2D or 3D), should has the same length
    %/                with length(year_list)
    %/
    %/                'lon' and 'lat' are optional, only for saving.
    %/  Output Data:  'trend_data' [unit per yr]
    
    if isempty(data_daily) && isempty(data_yearly)
        error('Neither ''data_daily'' nor ''data_yearly'' is given!');
    end
    
    %/ Set year_list to a row vector
    if size(year_list, 2) == 1   
        year_list = year_list';  
    end  

    %/ Set the strings
    str_mth = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec',....
               'MAM', 'JJA', 'SON', 'DJF', 'AMJJAS', 'ONDJFM', 'JFD', 'nonJJA', 'MJJASO', 'NDJFMA', 'MJJAS'};
%     str_mth = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'MAM', 'JJA', 'SON', 'DJF'};
    seq = sequence('data', year_list);
    str_years = '';
    for i = 1:length(seq)
        a = strjoin({num2str(seq{i}(1)), num2str(seq{i}(end))}, '-');
        str_years = strjoin({str_years, a}, '_');
    end
    if isempty(trend_test)
        str_trend_test = '_t-test';
    elseif isequal(trend_test, 'MK')
        str_trend_test = '_MK_test';
    else
        error('Invalid input of ''trend_test''!');
    end
    time_dim = length(size(data_daily));
    
    %/ If yearly data is not given
    if isempty(data_yearly) 
        %/ Convert data_daily to data_yearly (based on 'mth' and 'ins_or_acc')
        %/ NOTE: daily2any also works for monthly data
        [data_yearly, ~, ~] = daily2any('data_daily', data_daily, 'dates',  data_dates, 'mth', mth, 'ins_or_acc', ins_or_acc);
        
        if ~isequal(size(data_yearly, time_dim), length(year_list))
            size(data_yearly, time_dim)
            length(year_list)
            error('Time dim (assumed to be at the last dim) is not consistent with length(year_list)!');
        end
    end
    
    if time_dim == 2 %/ Assume in (vars, time)
        fprintf('============================================================\n')
        fprintf('======= Start computing the trend/slope (%s)... ======\n', strrep(str_trend_test, '_', ' '))
        fprintf('============================================================\n')
        [nvar, ~]         = size(data_yearly);
        trend_data        = nan(nvar,1);
        trend_data_sig    = nan(nvar,1);
        intercept_data    = nan(nvar,1);
        pval_data         = nan(nvar,1);
        CI95_data         = nan(nvar,1);
        Rsquared_data     = nan(nvar,1);
        Rsquared_adj_data = nan(nvar,1);
        RMSE_data         = nan(nvar,1);
        for i = 1:nvar
            model     = fitlm(year_list', squeeze(data_yearly(i,:)), 'y ~ x1 + 1');
            % disp(model)
            beta         = table2array(model.Coefficients(2,1));   %/ [unit per yr]
            se           = table2array(model.Coefficients(2,2));   %/ standard error of the regression coef (beta)
            CI95         = 1.96*se;                                %/ 95% Confidence Interval of the regression coef (beta)
            Rsquared     = model.Rsquared.Ordinary;                %/ For prediction, adjusted Rsquared may be a better option; otherwise, use the oridinary one
            Rsquared_adj = model.Rsquared.Adjusted;
            RMSE         = model.RMSE;
            if isempty(trend_test)
                beta_pval = table2array(model.Coefficients(2,4));
            elseif isequal(trend_test, 'MK')  %/ Mann-Kendall test
                %/ CAVEAT:  if data are not temporally evenly spaced, Sen's slope becomes
                %/          inaccurate (a future TODO).
                if all(isnan(data_yearly(i,:)))
                    continue;
                else
                    datain = [year_list', data_yearly(i,:)']; %/ must be in (nObs, 2)
                    [~, ~, ~, MK_pval] = ktaub(datain, trend_alpha, 0);
                    beta_pval = MK_pval;
                end
            else
                error('Invalid input of ''trend_test''!');
            end
            intercept = table2array(model.Coefficients(1,1));

            trend_data(i) = beta;
            if beta_pval < trend_alpha  %/ Store sig trend.
                trend_data_sig(i) = beta;  
            end
            pval_data(i)         = round(beta_pval, 3);
            intercept_data(i)    = intercept;
            CI95_data(i)         = CI95;
            Rsquared_data(i)     = Rsquared;
            Rsquared_adj_data(i) = Rsquared_adj;
            RMSE_data(i)         = RMSE;
        end
        
    elseif time_dim == 3 %/ Assume in (lon, lat, time)  
        if mth == 0
            str_mth_bc = [];
        else
            str_mth_bc = strcat('_', str_mth{mth});
        end 
        trend_data_filename = char(strcat(data_folder, data_name, '_trend_', ins_or_acc, str_trend_test, '_a', num2str(trend_alpha), str_years, str_mth_bc, '.mat'));
        [nlon, nlat, ~]     = size(data_yearly);

        tic;
        if isfile(trend_data_filename) && recompute_trend == 0
            fprintf('*** The queried trend data are found. *** \n*** Loading: %s *** \n', trend_data_filename)
            load(trend_data_filename, 'trend_data', 'trend_data_sig');    
        else
            disp(trend_data_filename)
            fprintf('=============================================\n')
            fprintf('======= Start computing trends (%s)... ======\n', strrep(str_trend_test, '_', ' '))
            fprintf('=============================================\n')

            if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
                 parpool('Threads', NumWorkers) %/ faster
            end

            %/ linear regression on years (to get the slope -> trends)
            trend_data_cell     = cell(nlon, 1);
            trend_data_sig_cell = cell(nlon, 1);
            intercept_data      = [];   %/ no need for 3D data
            CI95_data           = [];   %/ no need for 3D data
            Rsquared_data       = [];   %/ no need for 3D data
            Rsquared_adj_data   = [];   %/ no need for 3D data
            RMSE_data           = [];   %/ no need for 3D data
            parfor i = 1:nlon
                fprintf('*** i = %d/%d ***\n', i, nlon);
                temp_1D     = nan(1,nlat);
                temp_1D_sig = nan(1,nlat);
                for j = 1:nlat
                    model     = fitlm(year_list', squeeze(data_yearly(i,j,:)), 'y ~ x1 + 1');
                    beta      = table2array(model.Coefficients(2,1));   %/ [unit per yr]
                    
                    if isempty(trend_test)
                        beta_pval = table2array(model.Coefficients(2,4));
                    elseif isequal(trend_test, 'MK')  %/ Mann-Kendall test
                        %/ CAVEAT:  if data are not temporally evenly spaced, Sen's slope becomes
                        %/          inaccurate (a future TODO).
                        if all(isnan(data_yearly(i,j,:)))
                            continue;
                        else
                            datain = [year_list', squeeze(data_yearly(i,j,:))]; %/ must be in (nObs, 2)
%                             size(datain)
                            [~, ~, ~, MK_pval] = ktaub(datain, trend_alpha, 0);
                            beta_pval = MK_pval;
                        end
                    else
                        error('Invalid input of ''trend_test''!');
                    end
%                     intercept = table2array(model.Coefficients(1,1));
                    
                    temp_1D(j) = beta;
                    if beta_pval < trend_alpha  %/ Store sig trend.
                        temp_1D_sig(j) = beta;  
                    end
                end
                trend_data_cell{i}     = temp_1D;
                trend_data_sig_cell{i} = temp_1D_sig;
            end
            trend_data     = cat(1, trend_data_cell{:});
            trend_data_sig = cat(1, trend_data_sig_cell{:});

            if save_trend 
                fprintf('*** Saving trend data: %s *** \n', trend_data_filename)
                save(trend_data_filename, 'trend_data', 'trend_data_sig', 'lon', 'lat', '-v7.3');
            end
        end
        toc;
    else
        error('Input data is not a 3D data!');   
    end
    
    fprintf('======= Finished computing the trend/slope (%s) ======\n', strrep(str_trend_test, '_', ' '))
    fprintf('============================================================\n')

    %/ Close parpool when not using
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj)
    end


end
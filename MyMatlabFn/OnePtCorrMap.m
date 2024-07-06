%%

function [CorrMap, CorrMap_sig, BetaMap, BetaMap_sig, BetaMap_SE_sig, Stepwise_X, Stepwise_R2] = OnePtCorrMap(varargin)

pnames = {'CorrdataX', 'CorrdataY', 'CorrType', 'alpha', 'corr_or_reg', 'onept_or_local', 'stepwise_mode', 'local_partcorr', 'beta_in_prec', 'NumWorkers'};
dflts  = {         [],          [],         {},    0.05,             1,                1,               0,                0,              0,           20};
[           CorrdataX,   CorrdataY,   CorrType,   alpha,   corr_or_reg,   onept_or_local,   stepwise_mode,   local_partcorr,   beta_in_prec,   NumWorkers] = ...
        internal.stats.parseArgs(pnames, dflts, varargin{:});

%/ Matrix Initialization
CorrMap = [];  CorrMap_sig = [];
BetaMap = [];  BetaMap_sig = []; BetaMap_SE_sig = [];

if local_partcorr
    if length(size(CorrdataX)) == 4 || length(size(CorrdataX)) == 3
        [nlon, nlat, ntime, nvar] = size(CorrdataX);  %/ if CorrdataX is 3D, nvar = 1; 4D CorrdataX is for computing partial corr.
        CorrMap     = squeeze(nan(nlon,nlat,nvar));   %/ could be 3D or 4D.            4D CorrdataX is for computing partial corr.
        CorrMap_sig = squeeze(nan(nlon,nlat,nvar));

    elseif length(size(CorrdataX)) == 2
        [~, nvar] = size(CorrdataY);
        CorrMap     = nan;
        CorrMap_sig = nan;
    else
        error('Check the size of CorrdataX')
    end
    Stepwise_X = [];
    Stepwise_R2 = [];
else
    %/ Here, CorrdataY --> predictor(s)
    if length(size(CorrdataY)) == 2
        [~, nvar] = size(CorrdataY);
        
    elseif length(size(CorrdataY)) == 3
        [~, ~, nvar] = size(CorrdataY);
    end
    
    if length(size(CorrdataX)) == 3
        [nlon, nlat, ntime] = size(CorrdataX);  %/ if CorrdataX is 3D, nvar = 1; 4D CorrdataX is for computing partial corr.
        BetaMap        = nan(nlon,nlat,nvar);
        BetaMap_sig    = nan(nlon,nlat,nvar);
        BetaMap_SE_sig = nan(nlon,nlat,nvar);
        Stepwise_X     = cell(nlon,nlat);
        Stepwise_R2    = nan(nlon,nlat);

    elseif length(size(CorrdataX)) == 2
        BetaMap        = nan(1,nvar);
        BetaMap_sig    = nan(1,nvar);
        BetaMap_SE_sig = nan(1,nvar);
        Stepwise_X     = {};
        Stepwise_R2    = nan;

    else
        error('Check the size of CorrdataX')
    end
    CorrMap        = nan(nlon,nlat);
    CorrMap_sig    = nan(nlon,nlat);
end

%================== correlation (One-Point) ==================%
if corr_or_reg == 1 && onept_or_local == 1
    
    if length(size(CorrdataY)) ~= 2   error('Check the size of CorrdataY! it has to be a vector!');  end
    if size(CorrdataY, 1) == 1        CorrdataY = CorrdataY';           end         %/ transpose it if not a column vector.

    if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
        parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
    end

    parfor ii = 1:nlon
%     for ii = 1:nlon
        for jj= 1:nlat
%             fprintf('ii = %d, jj = %d\n', ii, jj);

            [rho,pval]= corr(CorrdataY, squeeze(CorrdataX(ii,jj,:)),...
                            'Tail','both','Type',CorrType, 'Rows','complete'); %/ 'Rows','complete' to ignore NaNs.
            
            CorrMap(ii,jj) = rho;
            if pval < alpha
                CorrMap_sig(ii,jj) = rho;
            end
        end
    end
    
%================== correlation (local), partial or ordinary ==================%
elseif corr_or_reg == 1 && onept_or_local == 2
 
    %--------------------
    %/ no need to reshape. Simple and fast.
    if local_partcorr                                                      %/ parital correlation (local)

        if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
            parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
        end

        fprintf('** Performing local partial correlation ... ***\n');
        %/ Assume CorrdataX is a 4D data in (lon, lat, time, vars)
        parfor ii = 1:nlon
            CorrdataX_bc = CorrdataX;
            CorrdataY_bc = CorrdataY;
            for jj= 1:nlat
                
                if all(isnan(squeeze(CorrdataY_bc(ii,jj,:)))) continue;   end   %/ skip if all are nans.
                for kk = 1:nvar

                    %/ partialcorr(x,y,z): linear partial correlation coefficients between pairs of variables in x and y, controlling for the variables in z.
                    if nvar == 1   
                        [rho, pval] = partialcorr(squeeze(CorrdataY_bc(ii,jj,:)), squeeze(CorrdataX_bc(ii,jj,:,kk)), zeros(ntime,1),...
                                                  'Tail','both','Type',CorrType, 'Rows','complete'); %/ should be equivalent to corr(A(:,i), Z) if only one variable
                    else
                        remaining_kk = find(~ismember(1:nvar, kk)); 
                        [rho, pval] = partialcorr(squeeze(CorrdataY_bc(ii,jj,:)), squeeze(CorrdataX_bc(ii,jj,:,kk)), squeeze(CorrdataX_bc(ii,jj,:,remaining_kk)),...
                                                  'Tail','both','Type',CorrType, 'Rows','complete');
                    end

                    CorrMap(ii,jj,kk) = rho;
                    if pval < alpha
                        CorrMap_sig(ii,jj,kk) = rho;
                    end
                end
            end
        end
    else
        fprintf('** Performing local correlation ... ***\n');
        if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
            parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
        end

        parfor ii = 1:nlon
            for jj= 1:nlat
                if all(isnan(squeeze(CorrdataY(ii,jj,:)))) continue;   end   %/ skip if all are nans.
                
                [rho,pval]= corr(squeeze(CorrdataY(ii,jj,:)), squeeze(CorrdataX(ii,jj,:)),...
                                'Tail','both','Type',CorrType, 'Rows','complete'); %/ 'Rows','complete' to ignore NaNs.

                CorrMap(ii,jj) = rho;
                if pval < alpha
                    CorrMap_sig(ii,jj) = rho;
                end
            end
        end
    end
    
%================== linear regression (One-Point) ==================%
elseif corr_or_reg == 2

%     tic
%     warning('!!! Performing CorrdataX = beta * CorrdataY + C... (see if this is what you want) !!!')
    if length(size(CorrdataX)) == 2   %/ if so, then onept or local are just the same.
        
        if size(CorrdataX, 1) == 1        CorrdataX = CorrdataX';           end         %/ transpose it if not a column vector.
        if size(CorrdataY, 1) == 1        CorrdataY = CorrdataY';           end         %/ transpose it if not a column vector.
        
        %/ NOTE: model = fitlm(X, y, 'y ~ x1')
        if stepwise_mode
            model             = stepwiselm(CorrdataY, CorrdataX, 'constant','Upper','linear', ...
                                           'Criterion', 'bic', 'PEnter', 0, 'PRemove', 0.01, 'Verbose',0);  %/ 'constant','Upper','linear' to avoid including interaction terms in model
%             disp(model)
            slcted_Xs = model.PredictorNames;
            ori_Xs    = split(sprintf('x%d/', 1:nvar), '/')';
            ori_Xs    = ori_Xs(1:end-1);
            ind_which_X = findismember_loop(ori_Xs, slcted_Xs);

            %/ Only record Xs and R-square if Xs is not empty.
            if ~isempty(slcted_Xs)   
                Stepwise_X  = {ind_which_X};                   %/ since the size of Xs can be different among grids, make it a cell
                Stepwise_R2 = model.Rsquared.Ordinary;      
            end
        else 
            model       = fitlm(CorrdataY, CorrdataX,  'linear');          %/ mdl = fitlm(X,y)
%             disp(model)
            ind_which_X = 1:nvar;
        end
        beta      = nan(1,nvar);
        beta_SE   = nan(1,nvar);
        beta_pval = nan(1,nvar);
        
        intercept              = table2array(model.Coefficients(1,1));
        beta(ind_which_X)      = table2array(model.Coefficients(2:end,1));  %/ get all beta  
        beta_SE(ind_which_X)   = table2array(model.Coefficients(2:end,2));  %/ get all the standard errors of beta  
        beta_pval(ind_which_X) = table2array(model.Coefficients(2:end,4));  %/ get all the pvalue of beta  
        
        for k = 1:nvar
            BetaMap(k) = beta(k);
            if  beta_pval(k) < alpha    
                BetaMap_sig(k)    = beta(k);  
                BetaMap_SE_sig(k) = beta_SE(k);  
            end
        end
    else
        if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
            parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
        end
        parfor ii = 1:nlon
            for jj= 1:nlat
                CorrdataX_bc = squeeze(CorrdataX(ii,jj,:));
                
                if onept_or_local == 1
                    CorrdataY_bc = CorrdataY;                       %/ onept -> time nvar
                elseif onept_or_local == 2
                    CorrdataY_bc = squeeze(CorrdataY(ii,jj,:,:));   %/ local -> lon lat time nvar
                end
               
                if all(isnan(CorrdataX_bc)) || all(CorrdataX_bc == 0)    continue;    end   %/ skip if all nans or zeros, otherwise BIC will remove all predictors!

                %/ NOTE: model = fitlm(X,y, 'y ~ x1')
                if stepwise_mode
                    model             = stepwiselm(CorrdataY_bc, CorrdataX_bc, 'constant','Upper','linear', ...
                                                   'Criterion', 'bic', 'PEnter', 0, 'PRemove', 0.01, 'Verbose',0); %/ 'constant','Upper','linear' to avoid including interaction terms in model
%                     disp(model)                                            %/ display anova table
                    slcted_Xs = model.PredictorNames';
                    ori_Xs    = split(sprintf('x%d/', 1:nvar), '/')';
                    ori_Xs    = ori_Xs(1:end-1);
                    ind_which_X = findismember_loop(ori_Xs, slcted_Xs);
                    
                    %/ Only record Xs and R-square if Xs is not empty.
                    if ~isempty(slcted_Xs)   
                        Stepwise_X(ii,jj) = {ind_which_X};                   %/ since the size of Xs can be different among grids, make it a cell
                        Stepwise_R2(ii,jj) = model.Rsquared.Ordinary;      
                    end
                else 
                    model     = fitlm(CorrdataY_bc, CorrdataX_bc,  'linear');    
                    ind_which_X = 1:nvar;
                end
                
                beta      = nan(1,nvar);
                beta_SE   = nan(1,nvar);
                beta_pval = nan(1,nvar);

                intercept              = table2array(model.Coefficients(1,1));
                beta(ind_which_X)      = table2array(model.Coefficients(2:end,1));  %/ get all the beta coefs 
                beta_SE(ind_which_X)   = table2array(model.Coefficients(2:end,2));  %/ get all the beta coefs 
                beta_pval(ind_which_X) = table2array(model.Coefficients(2:end,4));  %/ get all the pvalue of beta coefs 

                for k = 1:nvar
                    BetaMap(ii,jj,k) = beta(k);           
                    if  beta_pval(k) < alpha    
                        BetaMap_sig(ii,jj,k)    = beta(k);
                        BetaMap_SE_sig(ii,jj,k) = beta_SE(k);
                    end
                end
                
            end
        end
    end
    BetaMap        = squeeze(BetaMap);          %/ in case nvar == 1
    BetaMap_sig    = squeeze(BetaMap_sig);      %/ in case nvar == 1
    BetaMap_SE_sig = squeeze(BetaMap_SE_sig);   %/ in case nvar == 1
%     fprintf('!!! Time taken for local linear regression: %.2f s !!!\n', toc)
end

end


%                 if beta_in_prec                 %/ beta in terms of % change
%                     BetaMap_sig(ii,jj) = beta/(min(squeeze(CorrdataY(ii,jj,:)))*beta + intercept)*100;    %/ does not make any sense. The min value of x can be far from the fitted line!
%                 end



%%
function [CorrMap, CorrMap_sig] = LocalCorr(varargin)
    
    pnames = {'CorrdataX', 'CorrdataY', 'CorrType', 'alpha', 'NumWorkers'};
    dflts  = {         [],          [],  'Pearson',    0.05,           20};
    [           CorrdataX,   CorrdataY,   CorrType,   alpha,   NumWorkers] = ...
            internal.stats.parseArgs(pnames, dflts, varargin{:});
    
        
    %/ NOTE: This function borrows some lines from 'OnePtCorrMap.m'
    %/       It is simpler if you just want to do local correlation for
    %/       whatever matrix it is (e.g. geomap, hovmoller)
    
                        
    fprintf('** Performing local correlation ... ***\n');
%                     if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
%                         parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
%                     end
    if isempty(alpha)    error('do Not input an empty alpha!!!');     end

    time_dim = length(size(CorrdataX));  %/ always assume the last dim is time.
    
    if time_dim == 3
        [n1, n2, ~] = size(CorrdataX);
        CorrMap     = squeeze(nan(n1, n2)); 
        CorrMap_sig = squeeze(nan(n1, n2));
        
        for ii = 1:n1
            for jj= 1:n2
                if all(isnan(squeeze(CorrdataY(ii,jj,:))))   continue;   end   %/ skip if all are nans.

                [rho,pval]= corr(squeeze(CorrdataY(ii,jj,:)), squeeze(CorrdataX(ii,jj,:)),...
                                'Tail','both','Type', CorrType, 'Rows','complete'); %/ 'Rows','complete' to ignore NaNs.

                CorrMap(ii,jj) = rho;
                if pval < alpha
                    CorrMap_sig(ii,jj) = rho;
                end
            end
        end
    else
        error('code not yet set.');
    end        
end
%%
function A_sig = AS_ttest_wrapper(varargin)

    %/ create a set of valid parameters and their default value
    pnames = {'A',     'S',  'X',  'dates',  'alpha'};
    dflts  = { [],      [],   [],       [],    0.05};
              [ A,       S,    X,    dates,    alpha] = ...
                    internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
    
    %------------------------------
    if isempty(alpha)    error('do Not input an empty alpha!!!');     end
    
    dim_time = length(size(X));
    if dim_time ~= 2 && dim_time ~= 3    error('code not yet for >3 dim data');  end
    
    %/ Retrieve years from the date format (e.g., yyyymmdd -> yyyy)
    L = numel(num2str(dates(1)));
    if L-4 >= 0
        conv = 1./10^(L-4);
        years = unique(floor(dates.*conv));
    else
        error('Invalid format of dates (e.g., %d)!', dates(1));
    end
    nyr   = length(years);
    
    %/ if X is not empty, always use X to compute A and S.
    if ~isempty(X)
        if ~isempty(A) || ~isempty(S)
            error('do not input A and S when inputting X at the same time! Will always use X to compute A and S.')
        end
        
        if dim_time == 2
            [ni, ntime] = size(X);
            nptd        = ntime/nyr;
            S           = nan(ni, nptd);      %/ ni x 73
            A           = nan(ni, nptd);      %/ ni x 73
            
        elseif dim_time == 3
            [ni, nj, ntime] = size(X);
            nptd            = ntime/nyr;
            S               = nan(ni, nj, nptd);  %/ ni x nj x 73
            A               = nan(ni, nj, nptd);  %/ ni x nj x 73
        end
        
        for ptd = 1:nptd
            ind = find(mod(date_yyyyptd, 1e2) == ptd);
            
            if dim_time == 2
                S(:,ptd)   = std (X(:,ind),   0, dim_time);  %/ std(X) = std(X - X_mean). No difference.
                A(:,ptd)   = mean(X(:,ind),      dim_time);  %/ if X = pentad_CISO, A will = pentad_clim_CISO. (have shown already)
            elseif dim_time == 3
                S(:,:,ptd) = std (X(:,:,ind), 0, dim_time);
                A(:,:,ptd) = mean(X(:,:,ind),    dim_time);  
            end
        end
%         fprintf('*** Max. interannual SD of %s: %.2f ***\n', select_field, max(S, [], 'all'))

        if any(isnan(S))  warning('S contains nan! Check if this is normal.');  end
    end
    
    if isempty(A)  error('empty A!'); end
    if isempty(S)  error('empty S!'); end
    
    if ~isequal(size(A), size(S))   
        size(A)
        size(S)
        error('dimension of A and S inconsistent!');  
    end
    
    %/ perform A/S t-test on zm mm data. NaN -> Insig. value
    A_sig = AS_ttest('A', A, 'S', S, 'n', nyr, 'alpha', alpha);

%---- backup code ----%
%     %/ create a set of valid parameters and their default value
%     pnames = {'A',     'S', 'pentad_allyr',  'date_yyyyptd', 'n_list',  'alpha'};
%     dflts  = { [],      [],             [],              [],       [],     0.05};
%               [ A,       S,   pentad_allyr,    date_yyyyptd,   n_list,    alpha] = ...
%                     internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
%                 
%     %------------------------------
%     if isempty(alpha)    error('do Not input an empty alpha!!!');     end
%     
%     dim_time = length(size(pentad_allyr));
%     if dim_time ~= 2 && dim_time ~= 3    error('code not yet for >3 dim data');  end
%     
%     years = unique(floor(date_yyyyptd/1e2));
%     nyr   = length(years);
%     
%     if isempty(S)
%         if ~isempty(n_list)
%             pentad_allyr_bc = nan(size(pentad_allyr));
%             
%              %/ Compute Harmonics by *EACH* year, and then append them tgt.
%             for t = 1:nyr
%                 ind = find(floor(date_yyyyptd/1e2) == years(t));
%                 if dim_time == 2
%                     pentad_allyr_bc(:,ind)   = my_fourier('f', pentad_allyr(:,ind),   'n_list', n_list); 
%                 elseif dim_time == 3
%                     pentad_allyr_bc(:,:,ind) = my_fourier('f', pentad_allyr(:,:,ind), 'n_list', n_list);  
%                 end
%             end
%         else
%             %/ If n_list is empty, directly use the input pentad_allyr to compute S.
%             pentad_allyr_bc = pentad_allyr;
%         end
%         
%         if any(isnan(pentad_allyr_bc))  warning('pentad_allyr_bc contains nan! Check if this is normal.');  end
%         
%         %/ Obtain year-to-year variation
%         if dim_time == 2
%             [ni, ntime] = size(pentad_allyr_bc);
%             nptd        = ntime/nyr;
%             S           = nan(ni, nptd);  %/ ni x 73
%         elseif dim_time == 3
%             [ni, nj, ntime] = size(pentad_allyr_bc);
%             nptd        = ntime/nyr;
%             S           = nan(ni, nj, nptd);  %/ ni x 73
%         end
%         
%         for ptd = 1:nptd
%             ind = find(mod(date_yyyyptd, 1e2) == ptd);
% 
%             if dim_time == 2
%                 S(:,ptd)   = std(pentad_allyr_bc(:,ind),   0, dim_time);  %/ std(X) = std(X - X_mean). No difference.
%             elseif dim_time == 3
%                 S(:,:,ptd) = std(pentad_allyr_bc(:,:,ind), 0, dim_time);
%             end
%         end
%         if any(isnan(S))  warning('S contains nan! Check if this is normal.');  end
%     end
%     
%     
%     if ~isequal(size(A), size(S))   
%         size(A)
%         size(S)
%         error('dimension of A and S inconsistent!');  
%     end
%     
%     %/ perform A/S t-test on zm mm data. NaN -> Insig. value
%     A_sig = AS_ttest('A', A, 'S', S, 'n', nyr, 'alpha', alpha);
    
end
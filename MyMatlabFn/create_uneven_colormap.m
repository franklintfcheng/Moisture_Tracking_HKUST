%%
function [contf_levels_new, colmap_new] = create_uneven_colormap(varargin)

    % create a set of valid parameters and their default value
    pnames = {'contf_levels_uneven', 'colmap_uneven'};

    dflts  = cell(length(pnames), 1);

    [ contf_levels_uneven, colmap_uneven] = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/ testing
%     contf_levels = [0, 0.01, 0.05, 0.1, 0.5, 1:4];
%     colmap = brewermap(length(contf_levels)-1, '*RdYlBu');

    if length(unique(contf_levels_uneven)) ~= length(contf_levels_uneven)
        error('There are redundant values in ''unique(contf_levels_uneven)''!');
    end

    %/ get the minimum interval of the levels
    intvl = min(abs(diff(contf_levels_uneven)));
    contf_levels_new = min(contf_levels_uneven):intvl:max(contf_levels_uneven);
    colmap_new       = nan(length(contf_levels_new)-1, 3);
    
    tol = 1e-8;  
    for i = 1:length(contf_levels_uneven)-1
        %/ We have to use %/ a tolerance value (not too small!) 
        %/ to avoid bug due to imprecision when comparing two floating values!!!
        ind1 = find(abs(contf_levels_new - contf_levels_uneven(i)) < tol); 
        ind2 = find(abs(contf_levels_new - contf_levels_uneven(i+1)) < tol);
        ind = ind1:(ind2-1);
        
%         ind = find(contf_levels_new >= contf_levels_uneven(i) & ...
%                    contf_levels_new < contf_levels_uneven(i+1));  %/ this is not robust due to floating number issue.    
%         disp(i)
%         disp(ind)
        colmap_new(ind, :) = repmat(colmap_uneven(i,:), length(ind), 1);
    end

    if ~isempty(find(isnan(colmap_new)))
        error('colmap_new contains NaN! Check your function!');
    end
    
end
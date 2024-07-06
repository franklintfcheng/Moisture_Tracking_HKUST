function [data_diff_sig, data_diff] = ttest2_sig_fn(data1, data2, alpha, dim)

    %/=====================================================================
    %/ Author: Franklin Cheng
    %/ Last update: 1 Mar 2024
    %/
    %/ Two-sample t test on meteorological variable 
    %/
    %/   Output: difference in the temporal means of dataX and dataY 
    %/           with sig values shown only ====%
    %/
    %/           Raw and anomaly field are both ok.
    %/=====================================================================

    fprintf('*** Running ttest2_sig_fn... ***\n')
    
    %/ Check if the gridding is consistent
    sz1 = size(data1);
    sz2 = size(data2);
    if ~isequal(sz1(1:2), sz2(1:2))
        error('Inconsistent gridding between data1 (%d x %d)and data2 (%d x %d)!', sz1(1), sz1(2), sz2(1), sz2(2));
    end

    %/ Inf -> NaN if exists (before ttest2)
    if ~isempty(find(isinf(data1), 1))
        data1(isinf(data1)) = NaN;   
        warning('Inf is found in data1, replaced with nan.');
    end
    if ~isempty(find(isinf(data2), 1))
        data2(isinf(data2)) = NaN;   
        warning('Inf is found in data2, replaced with nan.');
    end
    
    %/ X, Y are from distr with unknown and unequal variance.
    ttest_sig = ttest2(data1, data2, 'Alpha',alpha,'Dim',dim, 'Tail','both', 'Vartype','unequal'); %/ output is a logical matrix
    
    cond = (ttest_sig==1);
    data_diff = mean(data1, dim, 'omitnan') - mean(data2, dim, 'omitnan');
    
    data_diff_sig = data_diff;
    data_diff_sig(~cond) = nan;
    
%     [row_sig, col_sig] = find(ttest_sig==1);
% 
%     data_diff     = nanmean(dataX, dim) - nanmean(dataY, dim);
%     data_diff_sig = nan(size(data_diff));
%     
%     for j = 1:length(row_sig)
%         data_diff_sig(row_sig(j),col_sig(j)) = data_diff(row_sig(j),col_sig(j)); 
%     end
%     fprintf('*** ttest2_sig_fn costs %.2f s ***\n', toc);
end
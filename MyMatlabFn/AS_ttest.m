%%
function A_sig = AS_ttest(varargin)
    
    %/ create a set of valid parameters and their default value
    pnames = {'A', 'S', 'n', 'alpha'}; 
    dflts  = { [],  [],  [],    0.05};
    %/ parse function arguments
    [A, S, n, alpha] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    %/ Author: Fandy Cheng
    %/
    %/ NOTE: This t-test was designed by Wang and Xu 1997 to test the
    %/       significance of CISO component.
    %/       E.g.,  p    = at a given pentad
    %/              A(p) = CISO component
    %/              S(p) = Interannual S.D. of pentad anomalies (Since mu = 0 for pentad anomalies.)
    %/              H0: A(p) = 0
    %/              HA: A(p) ~= 0
    
    fprintf('*** Performing AS t-test (Wang and Xu 1997)... ***\n');
    if isempty(alpha)    error('do Not input an empty alpha!!!');     end
    
    t = A./(S/sqrt(n));
    p = (1-tcdf(abs(t),n))*2; %/ two-tailed p-value; Here take absolute value of t since t-distribution is symmetric anyway.

    cond = (p > alpha);
    A_sig = A;
    A_sig(cond) = nan;
    
% https://www.mathworks.com/matlabcentral/answers/163578-obtaining-the-p-value-from-the-t-stat
%     tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
%     tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;              % 1-tailed t-distribution
%     
%     p1tail = 1 - tdist1T(t, v)
%     p2tail = 1 - tdist2T(t, v)
    
    %/ ttest
%     if length(h0_mu) ~= 1
%         if size(anomaly_data, 1) ~= size(h0_mu, 1) || size(anomaly_data, 2) ~= size(h0_mu, 2)
%             error('Dimension sizes of the data and h0_mu are different!');
%         end
%         
%         %/ have to do ttest grid by grid.
%         tic
%         ttest_sig = zeros(size(h0_mu));
%         for i = 1:size(h0_mu, 1)
%         for j = 1:size(h0_mu, 2)
%             ttest_sig(i,j) = ttest(anomaly_data(i,j,:), h0_mu(i,j),'Alpha',alpha,'Dim', dim);
%         end
%         end
%         toc
%     else
%         ttest_sig = ttest(anomaly_data, h0_mu,'Alpha',alpha,'Dim',dim);
%     end
%     
%     anomaly_data_mean   = nanmean(anomaly_data, dim) - h0_mu;  %/ get the difference from the ho mu.
%     nsize               = size(anomaly_data_mean);
%     anomaly_data_sig    = nan(nsize(1), nsize(2));  
%     
%     [row_sig, col_sig]  = find(ttest_sig==1);
%     for j = 1:length(row_sig)
%         anomaly_data_sig(row_sig(j),col_sig(j)) = anomaly_data_mean(row_sig(j),col_sig(j)); 
%     end

end
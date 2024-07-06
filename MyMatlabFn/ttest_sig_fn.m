function [data_mean_sig, data_mean, ttest_ci, ttest_p] = ttest_sig_fn(varargin)
    
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 27 Feb 2024
    %/
    %/ Description: One-sample t-test on anomalies of meteorological variable.
    %/ Return:      Time mean of anomalous dataX with sig values shown only.
    %/=====================================================================

    switch nargin
        case 3   %/ given 3 inputs
            data    = varargin{1};
            alpha   = varargin{2};
            dim     = varargin{3};
            h0_mu   = 0;           
            mode    = [];
            
        case 4   %/ given 4 inputs
            data    = varargin{1};
            alpha   = varargin{2};
            dim     = varargin{3};
            h0_mu   = varargin{4}; 
            mode    = [];
            
        case 5   %/ given 5 inputs
            data    = varargin{1};
            alpha   = varargin{2};
            dim     = varargin{3};
            h0_mu   = varargin{4}; 
            mode    = varargin{5}; %/ 'nanincluded' -> set nans to h0_mu before ttest and restore them afterwards
            
        otherwise
            error('Unexpected inputs')
    end
    %%
    data_mean_sig = []; data_mean = []; ttest_ci = [];
    if isempty(data)
        return;
    end
    fprintf('*** Running ttest_sig_fn... ***\n')
    
    if ~isempty(mode) 
        if isequal(mode, 'nanincluded')
            fprintf('*** ttest_sig_fn: ''nanincluded'' is on -> set nans to h0_mu before ttest and restore them afterwards ***\n')
            cond_nan_3D       = isnan(data);
            data(cond_nan_3D) = h0_mu;  %/ all nans are then set to h0_mu. 
        else
            error('Invalid input of the 5th input argument ''mode''!');
        end
    end
    
    %/ Perform ttest.
    if length(h0_mu) ~= 1    %/ if h0_mu is a matrix
        if size(data, 1) ~= size(h0_mu, 1) || size(data, 2) ~= size(h0_mu, 2)
            error('Dimension sizes of the data and h0_mu are different!');
        end
        
        if dim > 3   
            error('[ttest_sig_fn]: The function only works for 3D, 2D or 1D data!');   
        end
        
        %/ Have to do ttest grid by grid.
        ttest_sig = zeros(size(h0_mu));
        ttest_p   = zeros(size(h0_mu));
        ttest_ci  = zeros(size(h0_mu));
        for i = 1:size(h0_mu, 1)
            for j = 1:size(h0_mu, 2)
                [ttest_sig(i,j), p, ci, ~] = ttest(data(i,j,:), h0_mu(i,j),'Alpha',alpha,'Dim', dim); %/ Note: ttest() does ignore NaN and inf values.
                ci = ci - mean(data(i,j,:), 'omitnan');  %/ remove the data mean to return the true confidence interval (plus and minus)
                ttest_ci(i,j) = ci(2);        %/ take the +ve interval (should be of the same magnitude as the -ve one)
                ttest_p(i,j)  = p;
            end
        end
    else %/ then h0_mu is a const (e.g., 0)
        % size(data)
        [ttest_sig, p, ci, ~] = ttest(data, h0_mu,'Alpha',alpha,'Dim',dim);  %/ Note: ttest() does ignore NaN and inf values.
        ci = ci - mean(data, dim, 'omitnan');  %/ remove the data mean to return the true confidence interval (plus and minus)

        ttest_p = p;
        if dim == 1
            ttest_ci = ci(2);
        elseif dim == 2
            ttest_ci = ci(:,2);
        elseif dim == 3
            ttest_ci = ci(:,:,2);        %/ take the +ve interval (should be of the same magnitude as the -ve one)
        end
    end
%     size(ttest_sig)
    
    %/ Restore back to nans after ttest if mode = 'nanincluded'.
    if isequal(mode, 'nanincluded') 
        data(cond_nan_3D) = nan;   
    end
    
    %/ NEW CODE
    data_mean = mean(data, dim, 'omitnan');

    data_mean_sig = data_mean;
    data_mean_sig(ttest_sig == 0) = nan;
    

    %/ OLD CODE
    % if h0_mu ~= 0 
    %     anomaly_data_mean   = mean(data, dim, 'omitnan') - h0_mu;  %/ get the difference from h0_mu. (-> anomaly)
    % else
    %     anomaly_data_mean   = mean(data, dim, 'omitnan');          %/ then get the mean (ignore the naming here).
    % end    
    % anomaly_data_sig    = anomaly_data_mean;  
    % anomaly_data_sig(ttest_sig == 0) = nan;

end
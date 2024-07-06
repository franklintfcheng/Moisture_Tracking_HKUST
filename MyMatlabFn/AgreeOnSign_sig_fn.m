function [data_mean_diff_sig, data_mean_diff, data_diff_pve_sum, data_diff_nve_sum] = AgreeOnSign_sig_fn(varargin)
    
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 5 Apr 2024
    %/
    %/ Description: Check agreement on signs
    %/
    %/ Return:      'data_mean_diff_sig': Significant Mean difference of data1 minus data2
    %/                  'data_mean_diff': Mean difference of data1 minus data2
    %/                            'pval': (N-1)-dimensional pvalue matrix
    %/=====================================================================

    switch nargin
        case 3   %/ given 3 inputs
            data_diff = varargin{1};
            thres     = varargin{2};
            dim       = varargin{3};  %/ along which dimension to perform KS-test
        otherwise
            error('Unexpected number of inputs. Input data1, data2, alpha, dim.')
    end
    %%
    data_mean_diff_sig = []; data_mean_diff = []; data_diff_sign_sum = [];
    if isempty(data_diff)
        return;
    end
    fprintf('*** Running AgreeOnSign_sig_fn... ***\n')

    %/ If the input data is a column vector, convert it into a row vector
    if iscolumn(data_diff)
        data_diff = data_diff';  %/ Convert a column vector into a row vector
        dim = 2;  %/ Correct the dim
    end
    sz = size(data_diff);
    SampleSize = sz(dim);

    data_mean_diff = mean(data_diff, dim, 'omitnan');
    % data_mean_diff_sign = sign(data_mean_diff);

    data_diff_pve = (data_diff > 0);
    data_diff_nve = (data_diff < 0);

    data_diff_pve_sum = sum(data_diff_pve, dim, 'omitnan');
    data_diff_nve_sum = sum(data_diff_nve, dim, 'omitnan');
    
    cond_pass_pve = (data_diff_pve_sum >= SampleSize*thres) & (data_mean_diff > 0);
    cond_pass_vve = (data_diff_nve_sum >= SampleSize*thres) & (data_mean_diff < 0);
    cond_pass = (cond_pass_pve | cond_pass_vve);
    
    data_mean_diff_sig = data_mean_diff;
    data_mean_diff_sig(~cond_pass) = nan;
end
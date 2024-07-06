function [data_mean_diff_sig, data_mean_diff, pval] = kstest2_sig_fn(varargin)
    
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 5 Apr 2024
    %/
    %/ Description: Two-sample KS-test on N-dimensional arrays
    %/
    %/ Return:      'data_mean_diff_sig': Significant Mean difference of data1 minus data2
    %/                  'data_mean_diff': Mean difference of data1 minus data2
    %/                            'pval': (N-1)-dimensional pvalue matrix
    %/=====================================================================

    switch nargin
        case 4   %/ given 4 inputs
            data1   = varargin{1};
            data2   = varargin{2};
            alpha   = varargin{3};
            dim     = varargin{4};  %/ along which dimension to perform KS-test
        otherwise
            error('Unexpected number of inputs. Input data1, data2, alpha, dim.')
    end
    %%
    data_mean_diff_sig = []; data_mean_diff = []; pval = [];
    if isempty(data1) || isempty(data2)
        return;
    end
    fprintf('*** Running kstest2_sig_fn... ***\n')
    sz1 = size(data1);
    sz2 = size(data2);

    %/ If the input data is a column vector, convert it into a row vector
    if isvector(data1) && isvector(data2)
        [~, pval] = kstest2(data1',data2');
        % pval
    elseif numel(sz1) == dim && numel(sz2) == dim
        if ~isequal(sz1(1:end-1), sz2(1:end-1))
            disp(sz1(1:end-1))
            disp(sz2(1:end-1))
            error('The gridding of data1 appears different from that of data2!');
        end

        %/ NOTE: KS-test allows different sample sizes of data1 and data2
        data1_reshp = reshape(data1, prod(sz1(1:end-1)), sz1(end)); %/ e.g., data in (lon, lat, time) -> (grids, time)
        data2_reshp = reshape(data2, prod(sz2(1:end-1)), sz2(end)); %/ e.g., data in (lon, lat, time) -> (grids, time)

        pval_reshp = nan(size(data1_reshp, 1), 1);
        for i = 1:size(data1_reshp, 1)

            %/ kstest2 will fail if the data contains all NaNs
            if ~all(isnan(data1_reshp(i,:))) && ~all(isnan(data2_reshp(i,:)))
                [~, pval_reshp(i)] = kstest2(data1_reshp(i,:),data2_reshp(i,:));
            end
        end

        %/ Restore the dimension
        pval = reshape(pval_reshp, sz1(1:end-1));
    else
        error('code not set for ''dim'' that is not the last dim of input data!');
    end
    
    cond = (pval <= alpha);
    data_mean_diff = mean(data1, dim, 'omitnan') - mean(data2, dim, 'omitnan');  %/ by default, data2 - data1
    data_mean_diff_sig = data_mean_diff;
    data_mean_diff_sig(~cond) = nan;
end
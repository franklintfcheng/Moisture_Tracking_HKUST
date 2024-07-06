%%
function [contf_data, contf_levels, cbar_YTickLabel, cbar_YTick]= contfdata4uneven(varargin)

    %/ create a set of valid parameters and their default value
    pnames = {'contf_data', 'contf_levels', 'cbar_YTickLabel'};  
    dflts  = cell(1, length(pnames));

    [          contf_data,   contf_levels,   cbar_YTickLabel] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

%%
    format long g
    if length(unique(string(diff(contf_levels)))) > 1   %/ Trick: using string() to convert double to str
                                                        %/ helps remove trailing zeros that confuse unique()!
        fprintf('*** Detected uneven contf_levels! Replacing the values of contf_data with indices... ***\n')
        contf_data_bc = contf_data;             %/ replicate

        for k = 1:length(contf_levels)+1
            if k == 1
                cond                = (contf_data   < contf_levels(k));
                contf_data_bc(cond) = k-1;
            elseif k == length(contf_levels)+1
                cond                = (contf_levels(k-1) <= contf_data);
                contf_data_bc(cond) = k-1;
            else
                [n1, n2]            = size(contf_data);
                contf_data_2Dto1D   = reshape(contf_data, [], 1);
                cond                = (contf_levels(k-1) <= contf_data_2Dto1D & contf_data_2Dto1D < contf_levels(k));

                contf_data_bc_2Dto1D        = reshape(contf_data_bc, [], 1);
                contf_data_bc_2Dto1D(cond)  = rescale(contf_data_2Dto1D(cond), k-1, k);  %/ [IMPORTANT]: As shading does not work well with discrete data, we need to rescale the data.
                contf_data_bc               = reshape(contf_data_bc_2Dto1D, n1, n2);
            end
        end

        if isempty(cbar_YTickLabel)
            cbar_YTickLabel = contf_levels(2:1:end-1); %/ By default, show each uneven level
        end

        cbar_YTick = findismember_loop(contf_levels, cbar_YTickLabel);  %/ update
        contf_data   = contf_data_bc;                                   %/ update
        contf_levels = 1:length(contf_levels);                          %/ update
    else
        %/ Do nothing. Return original input.
        cbar_YTick = [];
    end
end
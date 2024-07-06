%%
function output = T_CI_2table(x, varnames, rownames, str_round_or_signif, digits)

    [nhs, nphase, nseason, nattr] = size(x);
    disp(size(x));

    %/ Round up L, NL and NLO, before summing.
    if isequal(str_round_or_signif, 'significant')
        x = round(x, digits, str_round_or_signif);
    elseif isequal(str_round_or_signif, 'round')
        x = round(x, digits);
    else
        error('Wrong input!')
    end
    
    nattr = nattr+1;                                                       %/ +1 for 'Tot' column.
    nrows_each_hotspot = nseason*nattr;
    
    output = nan(nhs*nrows_each_hotspot, nphase);
    for i = 1:nhs
        for s = 1:nseason
            idx = (nrows_each_hotspot*(i-1) + nattr*(s-1)+1):((i-1)*nrows_each_hotspot + nattr*s);
%             disp(idx)
            output(idx,:) = cat(1, squeeze(x(i,:,s,:))', nansum(squeeze(x(i,:,s,:))', 1) );
        end
    end
    output(output == 0) = nan;   % set 0 to nan since nansum([nan nan]) = 0.
    
    output = array2table(output, 'VariableNames', varnames, 'RowNames', rownames);
end
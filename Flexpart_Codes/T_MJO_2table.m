%%
function output = T_MJO_2table(x,varnames, rownames)

%     dim = size(x);
    [nhs, nphase, nattr] = size(x);
    disp(size(x));

    nattr  = nattr + 1;                                                    %/ +1 for 'Tot' column.
    output = nan(nhs*nattr, nphase);
    for i = 1:nhs
        output(nattr*(i-1)+1:nattr*i,:) = cat(1, squeeze(x(i,:,:))', nansum(squeeze(x(i,:,:))', 1) );
    end

    output = round(output, 2, 'significant');
    output = array2table(output, 'VariableNames', varnames, 'RowNames', rownames);
end
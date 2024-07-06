function vm = vertmean(varargin)

    %/ IMPORTANT: data has to be 4D, with the 3rd dim being the p-level dim.
    pnames = {'data', 'level'}; 
    dflts  = cell(length(pnames),1);
    [          data,    level] = ...
        internal.stats.parseArgs(pnames, dflts, varargin{:});

    %%

    %/==========================================================================================
    error('This function has been deprecated! Use the function vertinte(''take_vertmean'', 1) instead!');
    %/==========================================================================================

    % %/ NOTE:
    % %/      No unit conversion of pressure is needed. Will be canceled out.
    % 
    % dP = abs(diff(level));
    % s = 0;
    % for p = 1:length(dP)
    % %     s = s + squeeze(mean(data(:,:,p:p+1,:),3)) * dP(p);
    % 
    %     A = squeeze(mean(data(:,:,p:p+1,:),3)) * dP(p);
    %     A(isnan(A)) = 0;
    %     s = s + A;
    % end
    % vm = s/sum(dP);

end
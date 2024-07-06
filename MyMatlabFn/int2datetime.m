function dt = int2datetime(dates_int, format, timezone)
    
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 28 Jun 2024
    %/
    %/ NOTE: Make sure the input 'format' matches 'dates_int'!!
    %/=====================================================================

    switch nargin
        case 1
            format = 'yyyyMMddHHmm';
            timezone = 'UTC';
        case 2
            timezone = 'UTC';

        case 3
            if ~ischar(timezone) && ~iscell(timezone)
                error('timezone must be a char or a cell!');
            elseif isecell(timezone)
                timezone = timezone{:};
            end
        otherwise
            error('Invalid input of data! Only three input at most!')
    end

    if ~isinteger(dates_int)
        % type = class(dates_int);
        % warning('Detected that the input date is not an integer but %s. Convert it into integers...', type);
        dates_int = int64(dates_int);
    end
    % disp(dates_int(1:3))
    dt = datetime(string(dates_int), 'InputFormat', format, 'TimeZone', timezone);
end
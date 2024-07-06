function X_datetime_int = datetime2int(X_datetime, format)

    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 28 Jun 2024
    %/=====================================================================

    % %/ Why bother to do so? ==> It avoids using statistics toolbox!!!!
    % if nargin == 1    %/ if the number of input argument = 1, set default format
    %     format = 'yyyymmddHHMM';
    % end
    
    switch nargin
        case 1
            format = 'yyyyMMddHHmm';
        case 2
            %/ do nothing
        otherwise
            error('Invalid input of data! Only two input at most!')
    end

    if isdatetime(X_datetime) == 0
        class(X_datetime)
        error('The input array is not of a datetime type!');
    end
    
    if isequal(format, 'yyyymmddHHMM') || isequal(format, 'yyyyMMddHHmm') %/ since some Matlab functions regard 'MM' as months or minutes
        X_double       = fix(X_datetime.Year*1e8 + X_datetime.Month*1e6 + X_datetime.Day*1e4 + X_datetime.Hour*1e2  + X_datetime.Minute);
        X_datetime_int = int64(X_double);
        
    elseif isequal(format, 'yyyymmdd') || isequal(format, 'yyyyMMdd')     %/ since some Matlab functions regard 'MM' as months or minutes
        X_double       = fix(X_datetime.Year*1e4 + X_datetime.Month*1e2 + X_datetime.Day);
        X_datetime_int = int64(X_double);
        
    elseif isequal(format, 'yyyymm') || isequal(format, 'yyyyMM')         %/ since some Matlab functions regard 'MM' as months or minutes
        X_double       = fix(X_datetime.Year*1e2 + X_datetime.Month);
        X_datetime_int = int64(X_double);
        
    else
        error('Wrong input of format!')
    end
    
    %/ check if there is any NaT. Since int64(NaT) = 0. We correct it into missing value (-9999, cannot set nan)
    X_datetime_int(isnat(X_datetime)) = -9999;  
end

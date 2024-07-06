%%
function date_LST = UTCtoLST(date_UTC, lon_array)

    if ~isdatetime(date_UTC)
        error('Input date_UTC should be of datetime type!')
    end

    %/ NOTE: This function can also handle a 2D lon_array -> output a 2D date_LST
    
    lon_array(lon_array > 180) = lon_array(lon_array > 180) - 360;  %/ convert into -180 to 180.

    date_LST = date_UTC + hours(lon_array/15); %/ To the first order approximation, the solar noon move 15 deg every hour.
    
end

%--- alternatively ---

%     time = lon/15;
%     Hr   = floor(time);
%     Min  = floor((time - Hr)*60);                   %/ NOTE: if min = 59.9, it's still the 59th minute of an hour.
%     
%     date_LST = date_UTC + hours(Hr) + minutes(Min); %/ To the first order approximation, the solar noon move 15 deg every hour.

function dates = timesince2date(varargin)

    pnames = {'filename', 'date_format', 'time_shift', 'verbose'};
    dflts  =  cell(1, length(pnames));
    [          filename,   date_format,   time_shift,   verbose] ...
                            = internal.stats.parseArgs(pnames, dflts, varargin{:});

%% 
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 28 Jun 2024
    %/
    %/ Description: This simple function is designed to quickly read the time units of 
    %/              a NetCDF file in format of days since 1900-01-01 (for example)
    %/             
    %/      OUTPUT: human-readable dates in date_format
    %/
    %/=====================================================================
    if isempty(date_format)
        date_format = 'yyyymmddHHMM';
    end
    if isempty(time_shift)
        time_shift = 0;
    end
    time             = double(ncread(filename,'time'));
    time_units       = ncreadatt(filename,'time','units');
    time_units_parts = strsplit(time_units,' '); %/ split the string by ' '
    if contains(time_units_parts{1}, 'day')
        time_units_conv = 1;
    elseif contains(time_units_parts{1}, 'hour')
        time_units_conv = 1/24;
    elseif contains(time_units_parts{1}, 'minute')
        time_units_conv = 1/24/60;
    elseif contains(time_units_parts{1}, 'second')
        time_units_conv = 1/24/60/60;
    else
        error('code not set for %s!', time_units_parts);     
    end
    ref_date   = datetime(time_units_parts{3},'InputFormat','yyyy-MM-dd');
    % dates = int64(str2num(datestr(ref_date + (time+time_shift)*time_units_conv, date_format)));
    dates = int64(str2double(string(ref_date + (time+time_shift)*time_units_conv, date_format)));

    if verbose
        disp(dates([1:3,end-2:end]))  %/ Print the first and last three dates
    end
end
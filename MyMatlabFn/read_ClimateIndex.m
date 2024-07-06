function ClimateIndex = read_ClimateIndex(varargin)

    pnames = {  'type',   'datafolder', 'indexname', 'skipheadings', 'dataformat'};
    dflts  = {'Stndrd',             [],          [],              0,         '%f'};
    [type, datafolder, indexname, skipheadings, dataformat] = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%    
    
    %/ MAIN
    if iscell(datafolder)
        datafolder = datafolder{:};
    end

    filename = strcat(datafolder, indexname, '.txt');
    fprintf('reading %s...\n',filename)
    text_lines  = textread(filename, '%s','delimiter', '\n');
    
    data=[];
    for line = (1+skipheadings):length(text_lines) % select data during 1978 - 2017
    %     disp(sscanf(char(text_lines(line)), dataformat)')
        data = cat(1, data, sscanf(char(text_lines(line)), dataformat)'); %split the cell after converting it to char array
    end
    
    if ismember(type, {'MJO_format'})
        % datayear           = data(:,1);
        data_date_yyyymmdd = data(:,1)*1e4 + data(:,2)*1e2 + data(:,3);

        ClimateIndex.(indexname)(:, 1) = data_date_yyyymmdd;
        ClimateIndex.(indexname)(:, 2) = data(:,7);    %/ Amplitude
        ClimateIndex.(indexname)(:, 3) = data(:,6);    %/ Phase     phase == 999 means missing values!
        
    else
        datayear = unique(floor(data(:,1))); %can also handle the date array like 1948.0067  1948.013 etc.

        if numel(num2str(datayear(1))) == 6 %/ it may have a date format of yyyymm
            datayear = unique(floor(data(:,1)/100)); % correct to yyyy
            data_date_yyyymm = unique(floor(data(:,1))); 
        else
            stDate = datetime(datayear(1),1,1, 'Format','yyyy-MM-dd');
            edDate = datetime(datayear(end),12,31, 'Format','yyyy-MM-dd');
            data_date = stDate:edDate;
            data_date_yyyymm = unique((data_date.Year*100 + data_date.Month)');
        end
        %/ double check
        if numel(num2str(datayear(1))) ~= 4
            error('datayear is not in 4 digits, check the date format!')
        end
        
        if ismember(type, {'Stndrd'})
            data = data(:,2:end);
            ClimateIndex.(indexname) = nan(length(datayear)*12,1);
            for t = 1:length(datayear)
                ClimateIndex.(indexname)(:, 1) = data_date_yyyymm;
                ClimateIndex.(indexname)((1:12)+(t-1)*12, 2) = data(t,:);
            end

        elseif ismember(type, {'OneCol'})
            data = data(:,2);
            ClimateIndex.(indexname)(:, 1) = data_date_yyyymm;
            ClimateIndex.(indexname)(:, 2) = data;

        elseif ismember(type, {'TwoCol'})
            data = data(:,3);
            ClimateIndex.(indexname)(:, 1) = data_date_yyyymm;
            ClimateIndex.(indexname)(:, 2) = data;
        else
            error('Wrong input of type. Check the input!');
        end

        %/ check if all years have complete data record.
        if mod(length(ClimateIndex.(indexname)(:, 2)), 12) ~= 0
            error('Contains years with incomplete data!')
        end
        
    %     %/ check missing values (no need, deprecated)
    %     ind = find( data_date_yyyymm >= select_year(1)*100 & data_date_yyyymm <= (select_year(end)+1)*100);
    %     % disp(size(ind))
    %     if ~isempty(find(any(data(ind) == [-99.99, -9999])))
    %         error('The Subsetted data may contain missing values!')
    %     end
    end
    
    % stDate = datetime(select_year(1),1,1, 'Format','yyyy-MM-dd');
    % edDate = datetime(select_year(end),12,31, 'Format','yyyy-MM-dd');
    % slct_date = stDate:edDate;
    % slct_date_yyyymm = unique([slct_date.Year*100 + slct_date.Month]');

end




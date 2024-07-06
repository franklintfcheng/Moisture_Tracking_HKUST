%%
function IBTrACS = read_IBTrACS(varargin)

    pnames = {'data_path', 'agency', };
    dflts =  {        [],     'hko', };
    [          data_path,    agency, ] ...
                            = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %/==========================================
    %/ Author: Fandy Cheng (fandycheng@ust.hk)
    %/ Date of Last Update: 24 Jun 2023
    %/==========================================
    
    if isempty(data_path)  data_path = '/disk/r128/tfchengac/IBTrACS_data/IBTrACS.ALL.v04r00.nc';  end
    clc; ncdisp(data_path);

    IBTrACS_basic_varname = {'time', 'basin', 'subbasin', 'number', 'season', 'name', 'nature', 'wmo_wind'};
    if isequal(agency, 'hko')
        IBTrACS_varname = [IBTrACS_basic_varname, 'hko_lon', 'hko_lat', 'hko_cat', 'hko_wind', 'hko_pres'];
    else
        error('Invalid ''agency''. Check your input or function!');
    end
    
    %/ Loop over all the selected variables
    IBTrACS_raw = [];
    for i = 1:length(IBTrACS_varname)
        IBTrACS_raw.(IBTrACS_varname{i}) = ncread(data_path, IBTrACS_varname{i});
    end
    IBTrACS = IBTrACS_raw;
    
    %/ Restructure the data
    [~, n2, n3]      = size(IBTrACS_raw.basin);
    IBTrACS.basin    = cell(n2, n3); 
    IBTrACS.subbasin = cell(n2, n3); 
    IBTrACS.cat      = cell(n2, n3);
    IBTrACS.nature   = cell(n2, n3); 
    ind_cat          = find(contains(IBTrACS_varname, 'cat'));
    for i = 1:n2
        for j = 1:n3
            IBTrACS.basin(i,j)    = {IBTrACS_raw.basin(:,i,j)'};
            IBTrACS.subbasin(i,j) = {IBTrACS_raw.subbasin(:,i,j)'};
            IBTrACS.cat(i,j)      = {IBTrACS_raw.(IBTrACS_varname{ind_cat})(:,i,j)'};
            IBTrACS.nature(i,j)   = {IBTrACS_raw.nature(:,i,j)'};
        end
    end
    IBTrACS.cat_list    = unique(IBTrACS.cat);
    IBTrACS.nature_list = unique(IBTrACS.nature);
    IBTrACS.name        = string(squeeze(IBTrACS_raw.name(:, :))'); %/ as netcdf store char in row major order, but matlab in column major order, so need to transpose it.
    
    basedate                    = datenum(1858,11,17, 0, 0, 0);
    IBTrACS.date_yyyymmdd       = nan(size(IBTrACS_raw.time));
    IBTrACS.date_yyyymmddHHMM   = cell(size(IBTrACS_raw.time));

    for t = 1:size(IBTrACS_raw.time, 2)
    %     disp(t)
        b = IBTrACS_raw.time(:, t);
        b(isnan(b)) = []; %/ remove NaN

        IBTrACS.date_yyyymmdd(1:length(b), t)       = str2num(datestr(basedate + b,'yyyymmdd'));
        IBTrACS.date_yyyymmddHHMM(1:length(b), t)   = cellstr(datestr(basedate + b,'yyyymmdd HH:MM'));
    end


end
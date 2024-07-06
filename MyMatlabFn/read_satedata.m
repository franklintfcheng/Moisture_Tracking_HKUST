function satedata = read_satedata(varargin)

    pnames = {'sate_product', 'sate_path',  'slct_year',    'select_field', 'lon_range', 'lat_range',  'NumWorkers',  };
    dflts  = cell(length(pnames), 1);
    [          sate_product,    sate_path,    slct_year,      select_field,  lon_range,   lat_range,    NumWorkers,  ] ...
                            = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 26 Jan 2024
    %/
    %/   'select_field': The desired time scale field of the output data
    %/
    %/=====================================================================
    
    %/ Convert it into string from cell, otherwise isequal would not work.
    if isempty(select_field)
        select_field = 'daily';
    elseif iscell(select_field)
        select_field = select_field{:};
    end
    
    if isempty(NumWorkers)
        NumWorkers = 5;
    end
    %================================================================================
    %/ *** READ VARIOUS SATELLITE / GAUGE DATA ***
    %/ NOTE: GLASS data has a very different temporal reoslution (e.g.,
    %/       value in every few weeks). The current function thus has not
    %/       handled it yet.
    flag_WB       = 1;  %/ whether it is water budget
    longname      = [];
    ins_or_acc    = 'acc';  
    unit_multiply = 1;
    if isequal(sate_product, 'CMORPH_P')
        var    = 'cmorph';   
        unit_multiply = 0.5;  %/ rain rate (mm/hr) to rain accumulation (mm per 30 min)
    elseif isequal(sate_product, 'IMERG_P')
        var    = 'precipitationCal';   %/ Recommended by the IMERG manual
    elseif isequal(sate_product, 'GPCC_P')
        var    = 'precip';
    elseif isequal(sate_product, 'GPCP_P')
        var    = 'sat_gauge_precip';  %/ merged satellite-gauge precip.
        unit_multiply = 24; %/ mean mm/d to mm/month
    elseif isequal(sate_product, 'CRU_P')
        var    = 'pre';
    elseif isequal(sate_product, 'GLEAM_E')
        var    = 'E'; 
    elseif isequal(sate_product, 'GLEAM_SMsurf')
        var    = 'SMsurf'; 
    elseif isequal(sate_product, 'HARv2_P')
        var           = 'prcp'; 
        unit_multiply = 24; %/ mean mm/hr to mm/day
    elseif isequal(sate_product, 'HARv2_E')
        var           = 'et'; 
        unit_multiply = 24; %/ mean mm/hr to mm/day
    elseif isequal(sate_product, 'HARv2_T2m')
        var        = 't2'; 
        ins_or_acc = 'ins';
        flag_WB    = 0;
    elseif isequal(sate_product, 'TPR_P')
        % var        = 'Prep'; %/ for hourly
        var        = 'tp';  %/ for daily (For strange reason, the varname changes)
        longname   = 'Precip';
    elseif isequal(sate_product, 'TPR_E')
        var        = 'SFCEVP'; 
        longname   = 'SFCEVP';
    elseif isequal(sate_product, 'TPR_SRO')
        var        = 'sfcr'; 
        longname   = 'SurfaceRunoff';
    elseif isequal(sate_product, 'TPR_SSRO')
        var        = 'ssfcr'; 
        longname   = 'SubsurfaceRunoff';
    elseif isequal(sate_product, 'HadCRUT5_T2m')
        var        = 'tas_mean';   %/ CAVEAT: this data are already anomalies relative to 1961-1990, no raw data is found.
        ins_or_acc = 'ins';
        flag_WB    = 0;
    elseif isequal(sate_product, 'NOAA_CDR_OLR')
        var        = 'olr';   %/ units are already in 'W m-2'
        ins_or_acc = 'ins';   
        flag_WB    = 0;
    elseif isequal(sate_product, 'NOAA_Interp_OLR')
        var        = 'olr';   %/ units are already in 'W m-2'
        ins_or_acc = 'ins';   
        flag_WB    = 0;
    else
        error('invalid ''sate_product''!');
    end
    %================================================================================
    
    %/ Get basic information
    dt_slct_mo = []; dt_slct_hr = []; dt_slct_min = [];  %/ time interval in mo/hr/min
    satedata.(sate_product) = [];  
    if contains(sate_product, 'CMORPH') 
        if isempty(sate_path) 
            sate_path = '/disk/r128/tfchengac/CMORPH_data/';  %/ half-hourly, 
        end
        sate_filename_random = strcat(sate_path, '30min/8km/2017/09/28/CMORPH_V1.0_ADJ_8km-30min_2017092820.nc');
        dt_slct_min = 30;     
        ori_field = 'subdaily';

    elseif contains(sate_product, 'IMERG')
        if isempty(sate_path) 
            sate_path = '/disk/r128/tfchengac/IMERG_v6_data/';  %/ daily, 
        end
        sate_filename_random = strcat(sate_path, '3B-DAY-E.MS.MRG.3IMERG.20220623-S000000-E235959.V06.nc4');
        dt_slct_hr = 24;      
        ori_field = 'daily';

    elseif contains(sate_product, 'GPCC')
        if isempty(sate_path) 
            sate_path = '/disk/r128/tfchengac/GPCC_data/';  %/ daily
        end
        sate_filename_random = strcat(sate_path, 'opendata.dwd.de/climate_environment/GPCC/full_data_daily_v2022/full_data_daily_v2022_10_1982.nc');
        dt_slct_hr = 24; 
        ori_field = 'daily';
        
    elseif contains(sate_product, 'GPCP')
        if isempty(sate_path) 
            sate_path = '/disk/r128/tfchengac/GPCP_data/';  %/ monthly
        end
        sate_filename_random = strcat(sate_path, 'GPCPMON_L3_201912_V3.1.nc4');
        dt_slct_mo = 1;
        ori_field = 'monthly';

    elseif contains(sate_product, 'CRU') && ~contains(sate_product, 'HadCRUT5')  %/ avoid misunderstanding
        if isempty(sate_path) 
            sate_path = '/disk/r128/tfchengac/CRU_data/';  %/ monthly
        end
        sate_filename_random = strcat(sate_path, 'cru_ts4.07.1901.2022.pre.dat.nc');
        dt_slct_mo = 1;
        ori_field = 'monthly';
        
    elseif contains(sate_product, 'GLEAM')
        if isempty(sate_path) 
            sate_path = '/disk/r128/tfchengac/GLEAM_v3.6a_data/';  %/ daily
        end
        sate_product_suffix = strsplit(sate_product, '_');
        sate_filename_random = strcat(sate_path, sprintf('daily/1980/%s_1980_GLEAM_v3.6a.nc', sate_product_suffix{2}));
        dt_slct_hr = 24;        
        ori_field = 'daily';
        
    elseif contains(sate_product, 'HARv2')
        if isempty(sate_path) 
            sate_path = '/disk/r128/tfchengac/HARv2_data/';  %/ daily 
        end
        sate_filename_random = strcat(sate_path, sprintf('data.klima.tu-berlin.de/HAR/v2/d10km/d/2d/HARv2_d10km_d_2d_%s_1995.nc', var));
        dt_slct_hr = 24;        
        ori_field = 'daily';
    
    elseif contains(sate_product, 'TPR')
        if isempty(sate_path) 
            sate_path = '/disk/r128/tfchengac/TPR_data/';  %/ daily 
        end
        if ismember(sate_product, {'TPR_SSRO', 'TPR_SRO'})
            sate_filename_random = strcat(sate_path, sprintf('biggeo.gvc.gu.se/TPReanalysis/TP9km/Daily/%s/WRFOut_TP9km_Daily_Mean_%s_2016.nc', longname, longname));
            dt_slct_hr = 24;
            ori_field  = 'daily';
        elseif isequal(sate_product, 'TPR_P')
            % sate_filename_random = strcat(sate_path, sprintf('biggeo.gvc.gu.se/TPReanalysis/TP9km/Hourly/%s/WRFOut_TP9km_HourlyP_1999_03.nc', longname));
            % dt_slct_hr = 1;
            % ori_field  = 'subdaily';
            sate_filename_random = strcat(sate_path, sprintf('biggeo.gvc.gu.se/TPReanalysis/TP9km/Daily/%s/WRFOut_TP9km_Daily_Total_P_2011.nc', longname));
            dt_slct_hr = 24;
            ori_field  = 'daily';
        elseif isequal(sate_product, 'TPR_E')
            sate_filename_random = strcat(sate_path, sprintf('biggeo.gvc.gu.se/TPReanalysis/TP9km/Hourly/%s/WRFOut_TP9km_Hourly%s_1999_03.nc', longname, longname));
            dt_slct_hr = 1;
            ori_field  = 'subdaily';
        else
            error('Code not set!');
        end
        
    elseif contains(sate_product, 'HadCRUT5')
        if isempty(sate_path) 
            sate_path = '/disk/r128/tfchengac/HadCRUT5_Analysis_data/';  %/ monthly
        end
        sate_filename_random = strcat(sate_path, 'HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc');
        dt_slct_mo = 1;         
        ori_field = 'daily';
        
    elseif isequal(sate_product, 'NOAA_CDR_OLR')
        if isempty(sate_path)
            sate_path = '/disk/r128/tfchengac/NOAA_CDR_data/';  %/ daily
        end
        sate_filename_random = strcat(sate_path, 'www.ncei.noaa.gov/data/outgoing-longwave-radiation-daily/access/olr-daily_v01r02_20200101_20201231.nc');
        dt_slct_hr = 24;
        ori_field = 'daily';

    elseif isequal(sate_product, 'NOAA_Interp_OLR')
        if isempty(sate_path)
            sate_path = '/disk/r128/tfchengac/NOAA_data/';  %/ twice daily
        end
        sate_filename_random = strcat(sate_path, 'olr.2xdaily.1979-2022.nc');
        dt_slct_hr = 12;
        ori_field = 'subdaily';
    else
        error('code not set!')
    end
    ncdisp(sate_filename_random);
        
    if contains(sate_product, 'HadCRUT5')
        lon     = double(ncread(sate_filename_random, 'longitude')); 
        lat     = double(ncread(sate_filename_random, 'latitude'));
    else
        lon     = double(ncread(sate_filename_random, 'lon')); 
        lat     = double(ncread(sate_filename_random, 'lat'));
    end
    
    %/ For downscaling products (HARv2, TPR), lon, lat are 2D unevenly distributed matrices
    if contains(sate_product, {'HARv2', 'TPR'})
        %/ As it covers only Asian High Mountain, we do not subset lon lat.
        nlon    = size(lon,1);
        nlat    = size(lat,2);
        ind_lon = 1:nlon;
        ind_lat = 1:nlat;
    else
        if isempty(lon_range)    lon_range = [min(lon) max(lon)];  end
        if isempty(lat_range)    lat_range = [min(lat) max(lat)];  end
        ind_lon = find(lon >= lon_range(1) & lon <= lon_range(2));
        ind_lat = find(lat >= lat_range(1) & lat <= lat_range(2));
        nlon    = length(ind_lon);
        nlat    = length(ind_lat);
        lon     = lon(ind_lon);
        lat     = lat(ind_lat);
    end

    %/ Prepare slct_dates for double-checking
    % if isequal(ori_field, 'daily')
    %     slct_dates = date_array_gen('year_list', slct_year, 'dt_slct_hr', 24, 'output_date_format', 'yyyymmdd');          
    % elseif isequal(ori_field, 'monthly')
    %     slct_dates = date_array_gen('year_list', slct_year, 'dt_slct_mo', 1, 'output_date_format', 'yyyymmdd');      
    % elseif isequal(ori_field, 'subdaily')
    %     slct_dates = date_array_gen('year_list', slct_year, 'dt_slct_mo', dt_slct_mo, 'dt_slct_hr', dt_slct_hr, 'dt_slct_min', dt_slct_min, 'output_date_format', 'yyyymmddHHMM');  
    % else
    %     error('code not set!'); 
    % end

    output_data     = cell(length(slct_year), 1); 
    output_dates    = cell(length(slct_year), 1); 
    date_of_unreas  = [];
    TPR_E_thres     = 1e3;

    %/ Year Loop
    % if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
    %     parpool('Processes', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
    % end
    
    % parfor y = 1:length(slct_year)
    for y = 1:length(slct_year)   %/ testing
        %/ Broadcast variables
        ind_lon_bc   = ind_lon;
        ind_lat_bc   = ind_lat;
        year_list_bc = slct_year;
        year_bc      = year_list_bc(y);
        slct_dates_dt_eachyr  = date_array_gen('year_list', year_bc, 'dt_slct_mo', dt_slct_mo, 'dt_slct_hr', dt_slct_hr, 'dt_slct_min', dt_slct_min,...
                                              'output_date_format', 'datetime');
        slct_dates_UTC_eachyr = date_array_gen('year_list', year_bc, 'dt_slct_mo', dt_slct_mo, 'dt_slct_hr', dt_slct_hr, 'dt_slct_min', dt_slct_min,...
                                              'output_date_format', 'yyyymmddHHMM');                          

        if isequal(ori_field, 'daily')
            ind_EOD = [find(diff(datetime2int(slct_dates_dt_eachyr, 'yyyymmdd')) ~= 0); length(slct_dates_dt_eachyr)]; %/ index of the end of the day -> to calc daily value.
        elseif isequal(ori_field, 'monthly')
            ind_EOD = [find(diff(datetime2int(slct_dates_dt_eachyr, 'yyyymm')) ~= 0); length(slct_dates_dt_eachyr)]; %/ index of the end of the day -> to calc daily value.
        else
            ind_EOD = 1:length(slct_dates_UTC_eachyr); %/ just store every subdaily time step
        end
        ndate = length(ind_EOD);

        %/ Set t_list
%         t_list = 1:800;  %/ testing
        t_list = 1:length(slct_dates_UTC_eachyr);
        
        cnt = 0; aggr_data = 0; freq_data = 0; file_exist_flag = zeros(length(t_list), 1);
        for t = t_list
            casedate_yyyyMMddHH = floor(slct_dates_UTC_eachyr(t)/1e2);
            casedate_yyyyMMdd   = floor(slct_dates_UTC_eachyr(t)/1e4);
            casedate_year       = floor(slct_dates_UTC_eachyr(t)/1e8);
            casedate_month      = mod(floor(slct_dates_UTC_eachyr(t)/1e6), 1e2);
            casedate_day        = mod(floor(slct_dates_UTC_eachyr(t)/1e4), 1e2);
            casedate_min        = mod(slct_dates_UTC_eachyr(t), 1e2);  %/ minutes
            
            %/ Daily or subdaily products
            if ~isempty(dt_slct_hr) || ~isempty(dt_slct_min) 
                if contains(sate_product, 'CMORPH')
                    sate_folder   = strcat(sate_path, sprintf('30min/8km/%d/%02d/%02d/', casedate_year, casedate_month, casedate_day));
                    str_ver       = 'V1.0_ADJ';
                    sate_filename = sprintf('CMORPH_%s_8km-30min_%d.nc', str_ver, casedate_yyyyMMddHH);
                    sate_fullpath = [sate_folder, sate_filename];
                    if casedate_min == 0   ind_time = 1;     %/ odd
                    else                   ind_time = 2;            end
                    ntime         = 1;
                    all_dates_in_one  = 0;  %/ Whether the nc file contains all values for the target year
                
                elseif contains(sate_product, 'IMERG')
                    sate_folder   = sate_path;
                    sate_filename = sprintf('3B-DAY-E.MS.MRG.3IMERG.%d-S000000-E235959.V06.nc4', casedate_yyyyMMdd);
                    sate_fullpath = [sate_folder, sate_filename];
                    ind_time      = 1;
                    ntime         = 1;
                    all_dates_in_one  = 0;  %/ Whether the nc file contains all values for the target year
                
                elseif contains(sate_product, 'GPCC')
                    sate_folder   = strcat(sate_path, 'opendata.dwd.de/climate_environment/GPCC/full_data_daily_v2022/');
                    sate_filename = sprintf('full_data_daily_v2022_10_%d.nc', casedate_year);
                    sate_fullpath = [sate_folder, sate_filename];
                    ind_time      = 1;
                    ntime         = length(double(ncread(sate_fullpath, 'time'))); %/ Since the nc file contains all dates 
                    all_dates_in_one  = 1;  %/ Whether the nc file contains all values for the target year
                    
                elseif contains(sate_product, 'GPCP')
                    sate_folder   = sate_path;
                    sate_filename = sprintf('GPCPMON_L3_%d%02d_V3.1.nc4', casedate_year, casedate_month);
                    sate_fullpath = [sate_folder, sate_filename];
                    ind_time      = 1;
                    ntime         = length(double(ncread(sate_fullpath, 'time'))); %/ Since the nc file contains all dates 

                elseif contains(sate_product, 'GLEAM')
                    sate_folder   = strcat(sate_path, sprintf('daily/%d/', casedate_year));
                    sate_filename = sprintf('%s_%d_GLEAM_v3.6a.nc', var, casedate_year);
                    sate_fullpath = [sate_folder, sate_filename];
                    ind_time      = 1;
                    ntime         = length(double(ncread(sate_fullpath, 'time'))); %/ Since the nc file contains all dates 
                    all_dates_in_one  = 1;  %/ Whether the nc file contains all values for the target year
                    
                elseif contains(sate_product, 'HARv2')
                    sate_folder   = strcat(sate_path, 'data.klima.tu-berlin.de/HAR/v2/d10km/d/2d/');
                    sate_filename = sprintf('HARv2_d10km_d_2d_%s_%d.nc', var, casedate_year);
                    sate_fullpath = [sate_folder, sate_filename];
                    ind_time      = 1;
                    ntime         = length(double(ncread(sate_fullpath, 'time'))); %/ Since the nc file contains all dates 
                    all_dates_in_one  = 1;  %/ Whether the nc file contains all values for the target year
                    
                elseif contains(sate_product, 'TPR')
                    if ismember(sate_product, {'TPR_E'})  %/ monthly nc file
                        sate_folder   = strcat(sate_path, 'biggeo.gvc.gu.se/TPReanalysis/TP9km/Hourly/');
                        % if isequal(sate_product, 'TPR_P')
                            % sate_filename = sprintf('%s/WRFOut_TP9km_HourlyP_%d_%02d.nc', longname, casedate_year, casedate_month);
                        % elseif isequal(sate_product, 'TPR_E')
                        sate_filename = sprintf('%s/WRFOut_TP9km_HourlySFCEVP_%d_%02d.nc', longname, casedate_year, casedate_month);
                        % end
                        sate_fullpath = [sate_folder, sate_filename];
                        raw_dates = date_array_gen('year_list', year_bc, 'st_month', casedate_month, 'ed_month', casedate_month,...
                                                   'dt_slct_hr', dt_slct_hr, 'output_date_format', 'yyyymmddHHMM');  
                        ind_time      = findismember_loop(raw_dates, slct_dates_UTC_eachyr(t));
                        ntime         = 1;
                        all_dates_in_one  = 0;  %/ Whether the nc file contains all values for the target year
                        
                    elseif ismember(sate_product, {'TPR_P','TPR_SRO','TPR_SSRO'}) %/ yearly nc file
                        sate_folder   = strcat(sate_path, 'biggeo.gvc.gu.se/TPReanalysis/TP9km/Daily/');
                        if ismember(sate_product, {'TPR_P'})
                            sate_filename = sprintf('%s/WRFOut_TP9km_Daily_Total_P_%d.nc', longname, casedate_year);
                        else
                            sate_filename = sprintf('%s/WRFOut_TP9km_Daily_Mean_%s_%d.nc', longname, longname, casedate_year);
                        end
                        sate_fullpath = [sate_folder, sate_filename];
                        ind_time      = 1;
                        ntime         = length(double(ncread(sate_fullpath, 'time')));
                        all_dates_in_one  = 1;  %/ Whether the nc file contains all values for the target year
                    else
                        error('code not set!');
                    end
                    
                elseif contains(sate_product, 'NOAA_CDR')
                    if ismember(sate_product, {'NOAA_CDR_OLR'})
                        sate_folder   = strcat(sate_path, 'www.ncei.noaa.gov/data/outgoing-longwave-radiation-daily/access/');
                        sate_filename = sprintf('olr-daily_v01r02_%d0101_%d1231.nc', casedate_year, casedate_year);
                        sate_fullpath = [sate_folder, sate_filename];
                        ind_time      = 1;
                        ntime         = length(double(ncread(sate_fullpath, 'time'))); %/ Since the nc file contains all dates 
                        all_dates_in_one  = 1;  %/ Whether the nc file contains all values for the target year
                    else
                        error('code not set!');
                    end

                elseif contains(sate_product, 'NOAA_Interp')
                    if ismember(sate_product, {'NOAA_Interp_OLR'})
                        sate_folder   = sate_path;
                        sate_filename = 'olr.2xdaily.1979-2022.nc'; %/ one file only, specify ntime to read it in the year loop!
                        sate_fullpath = [sate_folder, sate_filename];

                        %/ First, convert time since XXXX-XX-XX into dates
                        raw_dates     = timesince2date('filename', sate_fullpath);
                        ind_date      = findismember_loop(floor(raw_dates/1e4), casedate_year);
                        ind_time      = ind_date(1);
                        ntime         = length(ind_date);
                        all_dates_in_one  = 1;  %/ Whether the nc file contains all values for the target year
                    else
                        error('code not set!');
                    end
                end

            elseif ~isempty(dt_slct_mo)   %/ Monthly products
                if contains(sate_product, 'GPCP')
                    time_shift    = -1;  %/ For some reason, the dates are shifted by one time unit. Hence, set time_shift = -1 to correct it.
                    sate_folder   = sate_path;
                    sate_filename = sprintf('GPCPMON_L3_%d%02d_V3.1.nc4', casedate_year, casedate_month);
                    sate_fullpath = [sate_folder, sate_filename];
                    all_dates_in_one  = 0;  %/ Whether the nc file contains all values for the target year
                elseif contains(sate_product, 'CRU') && ~contains(sate_product, 'HadCRUT5')
                    time_shift    = 0;
                    sate_folder   = sate_path;
                    sate_filename = 'cru_ts4.07.1901.2022.pre.dat.nc'; %/ one file only, specify ntime to read it in the year loop!
                    sate_fullpath = [sate_folder, sate_filename];
                    all_dates_in_one  = 1;  %/ Whether the nc file contains all values for the target year
                elseif contains(sate_product, 'HadCRUT5')
                    time_shift    = 0;
                    sate_folder   = sate_path;
                    sate_filename = 'HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc';
                    sate_fullpath = [sate_folder, sate_filename];
                    all_dates_in_one  = 1;  %/ Whether the nc file contains all values for the target year
                else
                    error('code not set!'); 
                end
                
                %/ First, convert time since XXXX-XX-XX into dates
                raw_dates = timesince2date('filename', sate_fullpath, 'time_shift', time_shift);
                ind_date  = findismember_loop(floor(raw_dates/1e8), casedate_year);
                ind_time  = ind_date(1);
                ntime     = length(ind_date);
            else
                error('Specify dt_slct_min, dt_slct_hr, or dt_slct_mo!'); 
            end
            
            if isfile(sate_fullpath)
                file_exist_flag(t) = 1;
                fprintf('*** [y = %d] [%02d/%02d] Reading ''%s'' from %s ... ***\n', year_bc, t, length(t_list), sate_product, sate_filename);
            else
                warning('!!! The file is not found: %s !!!\n Skipping', sate_fullpath);
                continue;
            end
            
            ind_lat_bc_parfor = ind_lat_bc;
            ind_lon_bc_parfor = ind_lon_bc;
            if contains(sate_product, 'IMERG')  %/ NOTE: only IMERG is in (lat, lon, time)
                raw_data = double(ncread(sate_fullpath, var,...
                                                [ind_lat_bc_parfor(1)  ind_lon_bc_parfor(1)  ind_time],... 
                                                [nlat                  nlon                  ntime]));
                raw_data = permute(raw_data, [2, 1, 3]);  %/ permute to (lon, lat, time) for consistency
                
            else %/ Otherwise, assume data in (lon, lat, time)
                raw_data = double(ncread(sate_fullpath, var,...
                                                [ind_lon_bc_parfor(1)  ind_lat_bc_parfor(1)  ind_time],... 
                                                [nlon                  nlat                  ntime]));
            end
            raw_data = squeeze(raw_data);
            cond_realzero = (raw_data == 0);
            
            %/ Fill missing values with zeros
            if flag_WB
                cond = (raw_data < 0 | isnan(raw_data)); %/ Check BOTH -ve values and NaNs in water budget data (e.g., P, E, RO)
                if ~isempty(find(cond, 1))
%                     fprintf('!!! [WARNING]: %.2f%% of the elements (%d) are negative or NaNs (i.e. missing values)! Correcting them to 0 !!!\n',...
%                             length(find(cond))/numel(raw_data)*100, length(find(cond)));
                    raw_data(cond) = 0;        %/ set -ve values to zero
                end
            else  
                cond = (isnan(raw_data));%/ Check NaNs in data     
                if ~isempty(find(cond, 1))
%                     fprintf('!!! [WARNING]: %.2f%% of the elements (%d) are NaNs (i.e. missing values)! Correcting them to 0 !!!\n',...
%                             length(find(cond))/numel(raw_data)*100, length(find(cond)));
                    raw_data(cond) = 0;        %/ set NaNs values to zero (for reduction to work properly)
                end
            end
            
            %/ Aggregation
            if isequal(sate_product, 'TPR_E')
                %/ [Known Issues]: Rectify the unreasonably huge values in TPR_E
                ind = find(mean(raw_data, [1,2], 'omitnan') > TPR_E_thres, 1); 
                
                if ~isempty(ind)
                    fprintf('!!! Unreasonably huge values in %s are detected on %s ***\n!!! Replace the field with that at the previous time step... ***\n',...
                            sate_product, string(casedate_yyyyMMddHH)) 
                    date_of_unreas = [date_of_unreas; casedate_yyyyMMddHH]; %/ record those dates
                    raw_data = raw_data_prev;
                end
                
                %/ Double check if the correction has been correctly done
                ind = find(mean(raw_data, [1,2], 'omitnan') > TPR_E_thres, 1); 
                if ~isempty(ind)
                    error('There are still unreasonably huge values in %s on %s! Check your code!', sate_product, string(casedate_yyyyMMddHH));
                end
                aggr_data     = aggr_data + raw_data;  
                raw_data_prev = raw_data;
            else
                aggr_data     = aggr_data + raw_data;  
            end
            
            %/ Compute frequency (ignore missing data)
            if isequal(ins_or_acc, 'ins')
                freq = logical(raw_data);  
                freq(cond_realzero) = 1;  %/ count those *real* zeros (e.g., Temp)
                freq_data = freq_data + freq; 
            end
            if t == 1                                                      %/ initialize
                output_eachyr = zeros(size(raw_data,1), size(raw_data,2), ndate);
                output_dates_eachyr = nan(ndate, 1);
            end
            
            %/ Whether the nc file contains all values for the target year?
            if all_dates_in_one
                output_eachyr = aggr_data;
                
                if ismember(ori_field, {'daily', 'monthly'})
                    output_dates_eachyr = floor(slct_dates_UTC_eachyr./1e4);
                else
                    output_dates_eachyr = slct_dates_UTC_eachyr;
                end
                break;  %/ Break the date loop of the year
            else
                if any(ismember(t, ind_EOD))                                   %/ store into daily data at the end of day
                    %/ Only the below WSV will have to take average.
                    cnt = cnt + 1;    %/ since we will skip the last time step (or should we?), cnt = 356 or 366
                    
                    if isequal(ins_or_acc, 'acc')
                        output_eachyr(:,:,cnt) = aggr_data;
                    elseif isequal(ins_or_acc, 'ins')
                        output_eachyr(:,:,cnt) = aggr_data./freq_data;
                    end
                    
                    if isequal(ori_field, 'daily')
                        output_dates_eachyr(cnt) = casedate_yyyyMMdd;
                    elseif isequal(ori_field, 'monthly')
                        output_dates_eachyr(cnt) = casedate_year*1e4+casedate_month*1e2+1; %/ First day of each month
                    else
                        output_dates_eachyr(cnt) = slct_dates_UTC_eachyr(cnt);
                    end
                    aggr_data = 0;    %/ reset
                end
            end
        end

        %/ If data has been loaded for this year, we remove the missing dates
        if any(file_exist_flag)
            %/ Remove NaN (missing dates) if any.
            ind_missing_dates = find(isnan(output_dates_eachyr));
            output_eachyr(:,:,ind_missing_dates)   = [];
            output_dates_eachyr(ind_missing_dates) = [];
        end

        %/ [IMPORTANT] Convert subdaily to daily (save memory!!)
        if isequal(ori_field, 'subdaily')
            [output_eachyr, output_dates_eachyr] = subdaily2daily('subdaily_data', output_eachyr, 'dates', output_dates_eachyr, 'ins_or_acc', ins_or_acc); 
        end

        % %/ If the data has not been aggregated before, see if we need to do it here
        % if flag_aggr == 0
        %     %/ If the product is daily/subdaily, while we quiry monthly
        %     %/ field, then we do daily2monthly
        %     if isequal(ori_field, 'monthly') && (~isempty(dt_slct_min) || ~isempty(dt_slct_hr))
        %         [output_eachyr, output_dates_eachyr] = daily2monthly('daily_data', output_eachyr, 'dates', output_dates_eachyr, 'ins_or_acc', ins_or_acc);
        %     end
        % end
        output_data{y}  = output_eachyr;
        output_dates{y} = output_dates_eachyr;
        output_eachyr   = [];
    end
    output_data         = cat(3, output_data{:});      %/ Avoid this as it will implicitly duplicate the matrix, causing memory issue
    output_dates = cat(1, output_dates{:});
    
    %====== Double Checking ======== 

    %/ Double-check if TPR_E has been corrected
    if isequal(sate_product, 'TPR_E')
        fprintf('!!! The %s values on %s have been corrected. !!!\n', sate_product, join(string(date_of_unreas), ', '))
        ind = find(mean(output_data, [1,2], 'omitnan') > TPR_E_thres); 
        if ~isempty(ind)
            error('!!! There are still unreasonably huge values in %s on %s! Check your code. !!!', sate_product, join(string(output_dates(ind)), ', '));
        end
    end

    %/ If monthly CRU_prcp a land grid cell will have zero prcp (or will it?)
    if isequal(sate_product, 'CRU_prcp')
        fprintf('!!! Correcting zeros to NaNs (since no oceanic values for %s data!). !!!\n', sate_product);
        output_data(output_data == 0) = nan;
    end
    
    %====== Convert ori_field into select_field if not daily =====%
    if isequal(select_field, 'monthly')
        [output_data, output_dates] = daily2monthly('daily_data', output_data, 'dates', output_dates, 'ins_or_acc', ins_or_acc);
    elseif isequal(select_field, 'yearly')
        [output_data, ~, output_dates, ~] = daily2any('data_daily', output_data, 'dates', output_dates, 'mth', 0, 'ins_or_acc', ins_or_acc);
    elseif ~isequal(select_field, 'daily')
        error('code not set for converting %s into %s!', ori_field, select_field);
    end

    % if ~isequal(ori_field, select_field)
        % fprintf('*** The queried field (%s) differs from the original field (%s) of the data. Post-processing will be performed. ***\n', select_field, ori_field)

        % if isequal(ori_field, 'subdaily')
        %     if isequal(select_field, 'daily')
        %         [output_data, output_dates] = subdaily2daily('subdaily_data', output_data, 'dates', output_dates, 'ins_or_acc', ins_or_acc); 
        % 
        %     elseif isequal(select_field, 'monthly')
        %         [output_data, output_dates] = daily2monthly('daily_data', output_data, 'dates', output_dates, 'ins_or_acc', ins_or_acc);
        % 
        %     elseif isequal(select_field, 'yearly')
        %         [output_data, ~, output_dates, ~] = daily2any('data_daily', output_data, 'dates', output_dates, 'mth', 0, 'ins_or_acc', ins_or_acc);
        % 
        %     else
        %         error('code not set for converting %s into %s!', ori_field, select_field);
        %     end

        % elseif isequal(ori_field, 'daily')
        % if isequal(select_field, 'monthly')
        %     [output_data, output_dates] = daily2monthly('daily_data', output_data, 'dates', output_dates, 'ins_or_acc', ins_or_acc);
        % elseif isequal(select_field, 'yearly')
        %     [output_data, ~, output_dates, ~] = daily2any('data_daily', output_data, 'dates', output_dates, 'mth', 0, 'ins_or_acc', ins_or_acc);
        % elseif ~isequal(select_field, 'daily')
        %     error('code not set for converting %s into %s!', ori_field, select_field);
        % end
        % else
            % error('code not set for converting %s into %s!', ori_field, select_field);
        % end
    % end

    %====== Store into satedata ======== 
    satedata.(sate_product).(select_field)  = output_data;
    satedata.(sate_product).lon             = lon;
    satedata.(sate_product).lat             = lat;
    clear output_data; %/ spare memory!
    
    %====== Units ========
    if any(contains(sate_product, {'HadCRUT5', 'HARv2_T2m'}))
        satedata.(sate_product).units = 'K';
    elseif any(contains(sate_product, {'OLR'}))
        satedata.(sate_product).units = 'W m^{-2}';
    else
        %/ Get the correct magnitude
        satedata.(sate_product).(select_field) = satedata.(sate_product).(select_field)*unit_multiply; 

        %/ Subdaily? Daily? Monthly?
        if isequal(select_field, 'daily')
            satedata.(sate_product).units = 'mm/day'; 
        elseif isequal(select_field, 'monthly')
            satedata.(sate_product).units = 'mm/month'; 
        else
            satedata.(sate_product).units  = 'mm';  
        end
    end

    %/ Convert dates to datetime
    if ismember(select_field, {'daily'})
        flds_date       = 'date_yyyymmdd_AllYr';
        output_dates_dt = int2datetime(output_dates, 'yyyyMMdd');

    elseif ismember(select_field, {'monthly'})
        output_dates    = floor(output_dates/100);  %/ yyyyMMdd -> yyyMM
        flds_date       = 'date_yyyymm';
        output_dates_dt = int2datetime(output_dates, 'yyyyMM');

    elseif ismember(select_field, {'subdaily'})
        flds_date       = 'date_yyyymmdd_AllYr';
        output_dates_dt = int2datetime(output_dates);
    else
        error('invalid ''select_field''!');
    end
    satedata.(sate_product).(flds_date)     = output_dates;       %/ int; yyyymmddHHMM (A redundant field for the convinience of coding)  
    % satedata.(sate_product).casedate_UTC    = output_dates;       %/ int; yyyymmddHHMM
    % satedata.(sate_product).casedate_UTC_dt = output_dates_dt;    %/ datetime in UTC;
    
%     %/ [Ac hoc treatment]: Rectify the unreasonably huge values in TPR_E
%     if isequal(sate_product, 'TPR_E')
%         if ismember(select_field, {'daily'})
%             ind = find(nanmean(satedata.(sate_product).(select_field), [1,2]) > 1e3); 
%         elseif ismember(select_field, {'monthly'})
%             ind = find(nanmean(satedata.(sate_product).(select_field), [1,2]) > 1e4); 
%         else
%             error('code not set!');
%         end
%         if ~isempty(ind)
%            fprintf('!!! Unreasonably huge values in %s are detected on %s ***\n!!! Rectifying them based on the value on the previous day... ***\n',...
%                     sate_product, join(strcat(string(output_dates(ind))), ', ')) 
%            for i = 1:length(ind)
%                satedata.(sate_product).(select_field)(:,:,ind(i)) = satedata.(sate_product).(select_field)(:,:,ind(i)-1); 
%            end
%         end
%         ind = find(nanmean(satedata.(sate_product).(select_field), [1,2]) > 1e3); 
%         if ~isempty(ind)
%             error('There are still unreasonably huge values in %s! Check your code!', sate_product);
%         end
%     end

    
end
%%
function [str_years, lon, lat, maxtraj_day, str_BLH_factor, str_optimal, str_traj_rm_jump, dt_slct,...
          data_folder, plotting_folder, WSV_dir, str_expmntinfo, basin_catalog, str_domain] = set_flexpart_fpath(varargin)
    
    % create a set of valid parameters and their default value
    pnames = {'expmnt', 'output_res', 'dt', 'year_list', 'ldirect', 'str_remark', 'str_RHc_dqc', 'optimal_tracking', 'optimal_rr', 'traj_rm_jump',...
              'BLH_factor', 'str_src_z', 'masterfolder','from_basin'};
    dflts  = cell(length(pnames), 1);
    [          expmnt,   output_res,   dt,   year_list,   ldirect,   str_remark,   str_RHc_dqc,   optimal_tracking,   optimal_rr,   traj_rm_jump,...
               BLH_factor,   str_src_z,   masterfolder,  from_basin] ...
                 = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    
    if contains(expmnt, '_EA_') 
        forcing     = 'EA';
        maxtraj_day = 20;                        %/ max traj days
        str_optimal = sprintf('_optimal%.2f',optimal_rr);

    elseif contains(expmnt, '_CERA_') 
        forcing = 'CERA';
        if optimal_tracking
            if ismember(ldirect, {'bwd'})
                maxtraj_day = 20;                        %/ max traj days
                if ~isempty(optimal_rr) && optimal_rr ~= 0.9
                    str_optimal = sprintf('_optimal%.2f',optimal_rr);
                else
                    str_optimal = '_optimal';
                end
            elseif ismember(ldirect, {'fwd'})
                maxtraj_day   = 20;                      %/ max traj days
                if ~isempty(optimal_rr) && optimal_rr ~= 0.1
                    str_optimal = sprintf('_optimal%.2f',optimal_rr);
                else
                    str_optimal = '_optimal';
                end
            end
        else
            maxtraj_day = 10;                        %/ max traj days
            str_optimal = [];
        end
    else
        error('Invalid input of ''expmnt'' = %s', expmnt);
    end

    %/ Load the parameters                 
    dt_slct     = 3;
    % leap        = dt_slct/dt;                       %/ Leap from FLEXPART time interval to user-defined interval
    lon         = 0:output_res:360-output_res;                    %/ this is the lon that I input to calcuate uptake map, although in FLexpart it is in [-179, 180].
    lat         = -90:output_res:90;
    
    %/ set the strings
    str_years = sprintf('_%d-%d', year_list(1), year_list(end));
    
    if traj_rm_jump == 0      str_traj_rm_jump = '_intact';                        else str_traj_rm_jump = [];   end
    if BLH_factor ~= 1        str_BLH_factor   = sprintf('_%.1fBLH', BLH_factor);  else  str_BLH_factor   = [];  end
    
    %/ indicate where to load the WSV data from
    WSV_dir     = cell(length(year_list),1);
    % flexpart_date_dt = cell(length(year_list),1);
    for y = 1:length(year_list)
        year = year_list(y);
        
        %========================= Set WSV folder =============================
        %/ Folder to save processed data (traj + WSV)
        if isequal(forcing, 'EA')
            WSV_folder = '/disk/r149/tfchengac/FLEXPART/';
                
        elseif isequal(forcing, 'CERA')
            if RHc_dqc_scheme == 1
                if ismember(ldirect, {'fwd'})
                    WSV_folder     = '/disk/r034/tfchengac/FLEXPART/';   
    
                elseif ismember(ldirect, {'bwd'})
                    if year >= 1991
                        WSV_folder = '/disk/r034/tfchengac/FLEXPART/'; 
                    else
                        WSV_folder = '/disk/r059/tfchengac/FLEXPART/'; 
                    end
                end
            else
                if ismember(ldirect, {'fwd'})
                    WSV_folder     = '/disk/r012/tfchengac/FLEXPART/';   
    %                 WSV_folder_hs3 = WSV_folder;  
    
                elseif ismember(ldirect, {'bwd'})
                    if year >= 2001
                        WSV_folder = '/disk/r149/tfchengac/FLEXPART/';  %/ 10 yrs -> 5 TB (~500 Gb a year)
                    elseif year >= 1981 && year < 2001
                        WSV_folder = '/disk/r012/tfchengac/FLEXPART/';  %/ 20 yrs -> 10 TB (~500 Gb a year)
                    elseif year >= 1971 && year < 1981
                        WSV_folder = '/disk/r148/tfchengac/FLEXPART/';  %/ 10 yrs -> 5 TB (~500 Gb a year)
                    end
                end
            end
        end
        WSV_dir{y}     = strcat(WSV_folder,     expmnt, '/prcssd_data_', num2str(year),'/');
        % %========================= Load the header ============================
        % if isequal(forcing, 'EA')
        %     flexpart_output_dir = strcat('/disk/r149/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
        % elseif isequal(forcing, 'CERA')
        %     if year >= 2001 && year <= 2010
        %         flexpart_output_dir = strcat('/disk/r037/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
        %     elseif year >= 1991 && year <= 2000
        %         flexpart_output_dir = strcat('/disk/r034/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
        %     elseif year >= 1981 && year <= 1990
        %         flexpart_output_dir = strcat('/disk/r012/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
        %     elseif year >= 1971 && year <= 1980
        %         flexpart_output_dir = strcat('/disk/r014/tfchengac/flexpart_output/', expmnt,'_',num2str(year),'/');
        %     end
        % end
    
    %     readp = 1; %/ read release points (0/1)
    %     [header, ~] = flex_header(flexpart_output_dir, 0, readp, 0);
    %     header.dates_str = string(importdata([flexpart_output_dir 'dates'])); %/ read dates file in string
    %     header.dates_dt = datetime(header.dates_str,'InputFormat','yyyyMMddHHmmss', 'Format', 'yyyy-MM-dd HH:mm:ss');
    % %     Area = calc_grid_area_header(header); %/ my function to calculating gird area.
    % 
    %     %/ Since we added a buffer month before each year's expmnt, 
    %     %/ we have to save the datetime in cells without concatenating them!!
    %     if leap == 1    
    %         flexpart_date_dt{y} = header.dates_dt(1:leap:end);
    %     elseif  leap == 2    
    %         flexpart_date_dt{y} = header.dates_dt(2:leap:end); %/ since it starts from year(t-1)-12-01 03:00
    %     end
    end
    
    %/ load basin_catalog based on 'from_basin'
    [basin_catalog, ~, str_from_basin] = load_from_basin('from_basin', from_basin);
    str_domain = str_from_basin;
    
    %/ set the experiment name and folder name.
    % if ismember(ldirect, {'fwd'})
    %     str_domain_traj     = '';  
    %     str_expmntinfo = strcat({' '}, str_RHc_dqc, {' '}, ldirect, {' '}, num2str(maxtraj_day),...
    %                             {'d'}, str_optimal, {' '}, num2str(dt_slct), 'h',  str_domain_traj, str_src_z, str_traj_rm_jump, str_BLH_factor, str_remark);
    % 
    %     foldername = strcat({' '}, str_RHc_dqc, {' '}, ldirect, {' '}, num2str(maxtraj_day),...
    %                             {'d'}, str_optimal, {' '}, num2str(dt_slct), 'h', str_src_z, str_traj_rm_jump, str_BLH_factor, str_remark);
    % 
    % elseif ismember(ldirect, {'bwd'})
    %     str_domain_traj = '_glbland';
    %     str_expmntinfo = strcat({' '}, str_RHc_dqc, {' '}, ldirect, {' '}, num2str(maxtraj_day),...
    %                             {'d'}, str_optimal, {' '}, num2str(dt_slct), 'h',  str_domain_traj, str_src_z, str_traj_rm_jump, str_BLH_factor, str_remark);
    % 
    %     foldername = str_expmntinfo;
    % end

    str_domain_traj = '_glbland';
    str_expmntinfo = strcat({' '}, str_RHc_dqc, {' '}, ldirect, {' '}, num2str(maxtraj_day),...
                            {'d'}, str_optimal, {' '}, num2str(dt_slct), 'h',  str_domain_traj, str_src_z, str_traj_rm_jump, str_BLH_factor, str_remark);

    foldername = str_expmntinfo;
    data_folder     = strcat(masterfolder, expmnt, '/prcssd_data_4plotting/');
    plotting_folder = strcat(masterfolder, expmnt, '/Plottings', strrep(foldername, ' ', '_'), '/');
    mkdir(char(data_folder));
    mkdir(char(plotting_folder));
    
    %/ how would dqc = 0.2 g/kg means in terms of mm/day?
    % A = mean(mean(Area));
    % mass = 5.095e18 / 4995000;
    % P_LA = dqc * 10e-3 * mass * 24/dt / A;
    % fprintf('*** dqc = %.3f g/kg corresponds to %.3f mm/day ***\n', dqc, P_LA);

end
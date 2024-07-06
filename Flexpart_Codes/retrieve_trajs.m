%%
function [trajs_catalog] = retrieve_trajs(varargin)

%/ create a set of valid parameters and their default value
pnames = {'WSV_dir', 'year_list', 'mth', 'st_month', 'st_day', 'ed_month', 'ed_day', 'str_RHc_dqc', 'ldirect', 'maxtraj_day',...
          'str_optimal', 'dt_slct', 'str_BLH_factor', 'str_remark', 'str_src_z',...
          'basin_name', 'basin_bndry', 'NumWorkers', 'test_mode'};  
dflts  = cell(1, length(pnames));

[         WSV_dir,    year_list,   mth,   st_month,   st_day,   ed_month,   ed_day,   str_RHc_dqc,   ldirect,   maxtraj_day,...
          str_optimal,    dt_slct,  str_BLH_factor,   str_remark,   str_src_z,...
          basin_name,    basin_bndry,    NumWorkers,  test_mode] ...
               = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
%%
%==========================================================================
%/      Author: Franklin Cheng
%/ Last Update: 30 Jan 2024
%/
%/ Description: This function is currently designd to retrieve backward
%/              trajectories (ldirect = 'bwd') from the 'basin_bndry'.
%/
%/              It will be updated to retrieve forward trajectories in the
%/              future.      
%==========================================================================
    
slct_dates_dt    = cell(length(year_list),1); 
WSV_name         = 'traj';   %/ by default

%/ Set dates
mth_bc            = mth;
st_month_bc       = st_month;
st_day_bc         = st_day;
ed_month_bc       = ed_month;
ed_day_bc         = ed_day;

%/ Retireve data based on 'mth' if the period is not given
if ~(isempty(st_month_bc) && isempty(st_day_bc) && isempty(ed_month_bc) && isempty(ed_day_bc))
    dt_target_dates = date_array_gen('year_list', year_list, 'st_month', st_month_bc, 'st_day', st_day_bc,...
                               'ed_month', ed_month_bc, 'ed_day', ed_day_bc,...   
                               'dt_slct_hr', dt_slct, 'output_date_format', 'datetime', 'skip_the_incomplete', 0);
else
    if mth == 0
        season = [];
    else
        season = mth_bc - 12;
    end
    dt_target_dates = date_array_gen('year_list', year_list, 'season', season,...   
                               'dt_slct_hr', dt_slct, 'output_date_format', 'datetime', 'skip_the_incomplete', 0);
end
disp(dt_target_dates([1:3,end-2:end]))


trajs_attr       = cell(length(year_list), 1);
trajs_date       = cell(length(year_list), 1);
trajs_date_ntraj = cell(length(year_list), 1);

for y = 1:length(year_list)
    %/ boardcast var
    season = []; st_month = []; ed_month = []; st_day = []; ed_day = [];
    year_bc = year_list(y);
    mth_bc = mth;
    
    if test_mode        %/ only consider 2010 Jan 1.
        if y ~= length(year_list)    continue;    end
        fprintf('==============================================================================\n')
        fprintf('============= Test Mode On (Only run Jan 1 2010, 8 time steps) ===============\n')
        fprintf('==============================================================================\n')
        st_month = 1; st_day = 1;   ed_month  = 1;  ed_day = 1;
    else
        if  mth_bc == 0    %/ mth_bc=0 --> whole period
            st_month = 1;     st_day = 1;   ed_month  = 12; 
            if year_bc == 2010 && ismember(ldirect, {'bwd'})  
                ed_day = 30;  %/ no 1231 in 2010 in our FLEXPART output.

            elseif year_bc == 2010 && ismember(ldirect, {'fwd'})  
                ed_day = 15;  %/ no data after 1215 in 2010 in our FLEXPART output.
            else
                ed_day = 31;
            end

        elseif mth_bc == 13   season = 1; %/ MAM
        elseif mth_bc == 14   season = 2; %/ JJA
        elseif mth_bc == 15   season = 3; %/ SON
        elseif mth_bc == 16   season = 4; %/ DJF
        else
            st_month  = mth_bc; st_day = 1;
            ed_month  = st_month; 

            if ed_month == 12 && year_bc == 2010 && ismember(ldirect, {'bwd'})  
                ed_day = 30;  %/ no 1231 in 2010 in our FLEXPART output.

            elseif ed_month == 12 && year_bc == 2010 && ismember(ldirect, {'fwd'})  
                ed_day = 15;  %/ no data after 1215 in 2010 in our FLEXPART output.
            else
                ed_day = eomday(year_bc,ed_month);
            end        
        end
    end
    slct_dates_dt{y} = date_array_gen('year_list', year_bc, 'st_month', st_month, 'st_day', st_day,...
                                      'ed_month', ed_month, 'ed_day', ed_day, 'season', season,...   %/ omit 20101231 if using season function.
                                      'dt_slct', dt_slct, 'output_date_format', []);

    slct_dates_dt_eachyr    = slct_dates_dt{y};
    ndate                   = length(cat(1, slct_dates_dt_eachyr));
    trajs_attr_eachyr       = cell(ndate, 1);
    trajs_date_eachyr       = cell(ndate, 1);
    trajs_date_ntraj_eachyr = nan(ndate, 1);
    
    if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
        parpool('Process', NumWorkers) 
    end
    parfor t = 1:length(slct_dates_dt_eachyr)   
%     for t = [1, 230:235] % 1:length(slct_dates_dt_eachyr)  
        slct_dates_dt_eachyr_bc = slct_dates_dt_eachyr;
        
        if isequal(datestr(slct_dates_dt_eachyr_bc(t,:), 'mmddHHMM'), '12312100')
            %/ NOTE: No 1231 21:00 data from flexpart expmnt, so we copy the 1231 18:00 instead.
            str_date = datestr(slct_dates_dt_eachyr_bc(t-1,:), 'yyyymmddHHMM'); 
            
        elseif year_bc == 2010 && isequal(datestr(slct_dates_dt_eachyr_bc(t,:), 'mmddHHMM'), '12302100') 
            %/ Since no 12302100 in year 2010, here we borrow 12301800 to compute daily value.
            str_date = datestr(slct_dates_dt_eachyr_bc(t-1,:), 'yyyymmddHHMM'); 
        else
            str_date = datestr(slct_dates_dt_eachyr_bc(t,:), 'yyyymmddHHMM');
        end

        if ~isempty(str_optimal)
            str_optimal_traj = '_optimal';
        else
            str_optimal_traj = [];
        end
        traj_suffix = strcat(str_RHc_dqc, '_', ldirect,'_',...
                                num2str(maxtraj_day), 'd', str_sharpcut, str_optimal_traj, '_', num2str(dt_slct), 'h_', str_domain_trajfile,'_', str_src_z, str_BLH_factor, str_remark,...
                                datestr(slct_dates_dt_eachyr_bc(t,:), 'yyyymmddHHMM'), '.mat');
           
    
        % %/ Since DJF data will use two-year data due to the 'date_array_gen' function. Need to correctly locate the dir.
        % b = cellfun(@(x) strsplit(x, {'_', '/'}), WSV_dir, 'UniformOutput', false);  %/ cellfun to strplit each cell.
        % b = cellfun(@(x)                x{end-1},       b, 'UniformOutput', false);  %/ cellfun to select year component from each cell.
        % ind = findismember_loop(b, {str_date(1:4)});
        
        %/ Loading file
        traj_fullpath = strcat(WSV_dir{ind}, WSV_name, traj_suffix);
        if isfile(traj_fullpath)
            fprintf('*** [year = %d, mth = %d, t = %d] Loading from it... %s *** \n', year_bc, mth_bc, t, strcat(WSV_name, traj_suffix))
            domfill = par_load(traj_fullpath, 'domfill');
        else
            error('!!! Data is not found: %s !!!',traj_fullpath);
        end
        
        %/ 1. Correct the lon range to strictly in (-180, 180]) (sometimes there are x = 182.xx)
        x = domfill.traj(:,1,:);
        x(x > 180) = x(x > 180) - 360;
        domfill.traj(:,1,:) = x;
        
        %/ 2. Subsetting traj with destination position within 'basin_bndry'
        traj_x_des = domfill.traj(:,1,1);
        traj_y_des = domfill.traj(:,2,1);
%         traj_x_des(traj_x_des > 180) = traj_x_des(traj_x_des > 180) - 360;
        
        try
            [in, ~] = inpoly2([traj_x_des, traj_y_des], basin_bndry);                          %/ inpoly2 is 600xx faster than inpolygon!! 
        catch 
%             warning(E.identifier, 'Error caught!: %s. Proceed with inpolygon().', E.message)
            [in, ~] = inpolygon(traj_x_des, traj_y_des, basin_bndry(:,1), basin_bndry(:,2));    %/ sometimes inpoly2 could run into a bug, then we use inpolygon.
        end
        ind = find(in == 1);

        if isempty(ind)
            %/ NOTE: do NOT skip the loop! Keep it running to save empty files! Or the program will always rerun..
            warning('!!! No trajs are found. !!!\n'); 
        end
        domfill.traj  = domfill.traj(ind,:,:);
        domfill.q     = domfill.q   (ind,:);
        domfill.T     = domfill.T   (ind,:);
        domfill.topo  = domfill.topo(ind,:);
        
        %/ Reshape (231,161) to (231,1,161), for example.
        domfill.q    = reshape(domfill.q,    size(domfill.q,1),    1, size(domfill.q,2));
        domfill.T    = reshape(domfill.T,    size(domfill.T,1),    1, size(domfill.T,2));
        domfill.topo = reshape(domfill.topo, size(domfill.topo,1), 1, size(domfill.topo,2));
        
        %/ Restore the true altitude of the particle.
        domfill.traj_real        = domfill.traj;
        domfill.traj_real(:,3,:) = domfill.traj(:,3,:) + domfill.topo;
        
        %/ Store into one matrix (ntraj, nvar, ntimestep)
        %/ [CAVEAT] Not every time step stores T along the full trajectory.
        %/          For now, store only x,y,z,q.
        trajs_attr_eachyr{t} = cat(2, domfill.traj_real, domfill.q, domfill.T);
        ntraj = size(domfill.traj_real, 1);
        
        trajs_date_eachyr{t} = str_date;                                   %/ record the date
        trajs_date_ntraj_eachyr(t) = ntraj;                                %/ record the ntraj in each date
    end

    %/ store the results for the year.
    trajs_attr{y}       = cat(1, trajs_attr_eachyr{:}); 
    trajs_date{y}       = cat(1, trajs_date_eachyr); 
    trajs_date_ntraj{y} = cat(1, trajs_date_ntraj_eachyr); 
    clear trajs_attr_eachyr     %/ spare memory
end

%/ concatenate
trajs_attr       = cat(1, trajs_attr{:});
trajs_date       = cat(1, trajs_date{:});
trajs_date_ntraj = cat(1, trajs_date_ntraj{:});

%/ combine cell arrays (trajs + date + ntraj of date)
trajs_catalog = [];
trajs_catalog.traj       = trajs_attr(:,1:3,:);
trajs_catalog.q          = trajs_attr(:,4,:);
trajs_catalog.date       = trajs_date;
trajs_catalog.date_ntraj = trajs_date_ntraj;

if length(size(trajs_attr, 2)) == 5
    trajs_catalog.T      = trajs_attr(:,5,:);
end

%/ Frequency map of trajs
% [ind_lon, ind_lat] = point2map('pt_x', reshape(trajs{1}(:,1,:), [], 1),...
%                                'pt_y', reshape(trajs{1}(:,2,:), [], 1),...
%                                'lon', lon, 'lat', lat, 'NumWorkers', NumWorkers);
% tic
% trajs_freq_map = 0;
% parfor i = 1:length(ind_lon)
%     xx = ind_lon(i);
%     yy = ind_lat(i);
%     
%     A = zeros(length(lon), length(lat));
%     A(xx, yy) = 1;
%     
%     trajs_freq_map = trajs_freq_map + A;  %/reduction variable
%     
% end
% fprintf('*** Finished computing traj freq map (%.1f seconds) *** \n', toc)




end
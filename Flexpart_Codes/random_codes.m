%% speed up performance techniquesA = rand(1000,1);

A = rand(1000,1);

tic
for idx=1:100000, b=mean(A); end
toc

tic
for idx=1:100000, b=sum(A,1)/size(A,1); end
toc
% ⇨ Elapsed time is 0.318459 seconds.

%%

data = randi(1000,1,1e6); % 1M elements between 1-1000
% Using find
tic
n1 = find(data==10);
toc

tic
n2 = strfind(data,10); %/ strind is much faster!!!
toc

isequal(n1,n2)

%%

    D = rand(1e6, 1);
    tic
    for i = 1:20
        a = 0;
        parfor j = 1:60
            a = a + sum(D);
        end
    end
    toc
    
    tic
    D = parallel.pool.Constant(D);
    for i = 1:20
        b = 0;
        parfor j = 1:60
            b = b + sum(D.Value);
        end
    end
    toc

    
    isequal(a, b)
    
%% DiffMap
% %/ NOTE: Run <flexpart_traj_freq_Dorling.m> first if want to visualize traj data here!!!
% 
% close all;
% slct_contf_str = ''; slct_cont_str = ''; slct_vector_str = ''; slct_hatch_str = '';
% contf_data = []; contf_lon = []; contf_lat = []; contf_levels = []; contf_unit = []; cbar_YTickLabel = []; cbar_YTick = []; colmap = []; cbar_interval = []; pcolor_mode = [];
% cont_data = []; cont_data_raw = []; cont_lon = []; cont_lat = []; cont_levels = []; cont_colmap = []; cont_linewi = []; cont_labelsize = [];
% Udata = []; Vdata = []; uv_lon = []; uv_lat = []; vec_step_lon = []; vec_step_lat = []; vector_color = []; 
% vecscale = []; vecscale2 = []; shaftwidth = []; headlength = []; headlength = []; vec_lbs = []; vec_mag_ref = []; bndry_data = [];
% traj_data = []; traj_levels = []; c_traj = []; traj_unit = []; plot_traj = 0; hatch_data = []; hatch_lon = []; hatch_lat = []; hatch_thres = []; hatch_mode = []; hatch_linewi =[]; hatch_intvl = [];
% str_mth = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'MAM', 'JJA', 'SON', 'DJF'};
% 
% savefig          = 0;
% draw_plot        = 0;
% draw_cbar_only   = 0;                                                     %/ also keep all_in_one == 0
% cbar_location    = 'southoutside';
% save_table       = 0;
% savemat          = 0;
% deseason_mode    = 0;  %/ I prefer not to deseason, esp. for MJO events. It does not reflect the pattern at all.
% detrend_mode     = 0;  %/ not to detrend, otherwise will affect L and NL computation
% SR_mode          = 1;
% LR_NLR_text_mode = 1;
% 
% if SR_mode           str_SR = '_SR';               else   str_SR = '';         end
% if LR_NLR_text_mode  str_text = '_LRNLR';          else   str_text = '';       end
% if detrend_mode      str_detrend  = '_detrend';    else   str_detrend = '';    end
% 
% sig_or_both   = 2;     %/ affect contour only
% alpha = 0.05;
% if sig_or_both    str_sig = char(sprintf('_sig%.2f', alpha));  else  str_sig = char(sprintf('_both%.2f', alpha)); end
% NumWorkers = 20;
% 
% global_map = 0;     %/ 1]: to check climate variability's global impact. 0]: to check regional influence on hotspots
% if global_map   
%     hs_list = 1;
%     slct_contf_str  = 'P';
%     slct_cont_str   = 'Z850';
%     slct_vector_str = 'IVT';
% else
%     hs_list = 1; %:length(AYR_hotspot_bc);
%     slct_contf_str  = 'BL_Pm';               %/ from each hotspots
% %     slct_contf_str  = 'contr_map';
%     slct_vector_str = '';
%     slct_hatch_str  = '';              
% %     slct_vector_str = 'IVT';
% %     slct_hatch_str  = 'P';                 %/ use hatched pattern/points to denote prcp!
% end
% 
% str_slct_event = 'MJO_event';
% % str_slct_event = 'complex_ENSO_event';
% % str_slct_event = 'ENSO_event';
% event_phase = ClimIndices.(str_slct_event);
% 
% rand_subset_portion = 0.01;    %/ 1%
% % plot_traj = 1; slct_traj_data = 'xyn';
% % plot_traj = 1; slct_traj_data = 'xyt';
% % plot_traj = 1; slct_traj_data = 'xyq';
% % plot_traj = 1; slct_traj_data = 'xyz';
% % slct_contf_str = 'traj_freq';
% % slct_contf_str = 'total_init_contr';
% % slct_contf_str = 'avg_init_contr';
% 
% if ismember(str_slct_event, {'ENSO_event'})
%     str_event_phase      = {'El', 'La', 'Neu'};
%     str_event_code       = [   1,   -1,     0];
%     skip_conv = 0;
%     
%     %========= Defaulted phase1 minus phase2 ==========%
%     slct_phase1 = { 'El',  'La'};
%     slct_phase2 = repmat({'Neu'}, 1, length(slct_phase1));                 %/ 2nd phase is always Neu
%     %===============================================%
%     
% elseif ismember(str_slct_event, {'complex_ENSO_event'})
%     str_event_phase = {'EP-El', 'Mixed-El', 'CP-El', 'EP-La', 'Mixed-La', 'CP-La', 'Neu'};
%     str_event_code  = [      1,          2,       3,      -1,         -2,      -3,    0];
%     skip_conv = 0;
%     
%     %========= Defaulted phase1 minus phase2 ==========%
%     slct_phase1 = {'EP-El', 'Mixed-El', 'CP-El', 'EP-La', 'Mixed-La', 'CP-La'};
%     slct_phase2 = repmat({'Neu'}, 1, length(slct_phase1));                 %/ 2nd phase is always Neu
%     %===============================================%
%     
% elseif ismember(str_slct_event, {'AMO_event'})
%     str_event_phase      = {'AMO+', 'AMO-', 'Neu'};
%     str_event_code       = [     1,     -1,     0];
%     skip_conv = 0;
%     
%     %========= Defaulted phase1 minus phase2 ==========%
%     slct_phase1 = {'AMO+', 'AMO-'};
%     slct_phase2 = repmat({'Neu'}, 1, length(slct_phase1));                 %/ 2nd phase is always Neu
%     %===============================================%
% 
% elseif ismember(str_slct_event, {'MJO_event'})
%     str_event_phase      = {'Neu', 'Phase 1', 'Phase 2', 'Phase 3','Phase 4','Phase 5','Phase 6','Phase 7','Phase 8'};
%     str_event_code       = 0:8;
%     skip_conv = 1;                                                         %/ do not convert data into monthly
%     
%     %========= Defaulted phase1 minus phase2 ==========%
%     slct_phase1 = {'Phase 1', 'Phase 2', 'Phase 3','Phase 4','Phase 5','Phase 6','Phase 7','Phase 8'};
%     slct_phase2 = repmat({'Neu'}, 1, length(slct_phase1));                 %/ 2nd phase is always Neu
%     %===============================================%
% end
% if plot_traj      str_traj = char(strcat({'traj '}, slct_traj_data));  else  str_traj = [];    end
% 
% % T       = cell(length(hs_list)*2, length(slct_phase1));  %/ two values for each hotspots
% % T_onfig = cell(length(hs_list)*2, length(slct_phase1));  %/ two values for each hotspots
% T       = {};  %/ two values for each hotspots
% T_onfig = {};  %/ two values for each hotspots
% T_double          = nan(length(hs_list), length(str_event_phase)-1, 2);  %/ hotspots * Indices (+ve/-ve) * Seasons * 2 (local vs non-local)
% T_double_fraction = nan(length(hs_list), length(str_event_phase)-1, 2);  %/ hotspots * Indices (+ve/-ve) * Seasons * 2 (local vs non-local)
% NLL_fldname_all   = cell(length(AYR_hotspot_bc), 1);
% 
% for top = hs_list
%     if isequal(AYR_hotspot_bc(top).name, 'Shebelle')   warning('!!! Skip Shebelle. !!!');   end
%     if any(ismember(slct_contf_str, {'contr_map', 'BL_Pm'}))
%         str_hs = sprintf('_AYR_hs#%d_%s', top, AYR_hotspot_bc(top).name);
% 
%         %/ 17 hotspots (regardless of hotspot_list_ver - wouldn't affect)
%         hs_ind{1} = [1, 4, 6, 7, 15];           %/ South America
%         hs_ind{2} = [2, 3, 5, 8, 10];           %/ Africa (skip Shebelle)
%         hs_ind{3} = [11, 12, 14, 16, 17, 13];   %/ Asia 
% 
%         flag_error = 1;
%         for ii = 1:3
%             if any(ismember(hs_ind{ii}, top))
%                 remaining_hs  = setdiff(hs_ind{ii}, top);
%                 color         = repmat([0 0 0]./255, length(hs_ind{ii}), 1);
%             end
%         end
%         bndry_data    = {AYR_hotspot_bc([remaining_hs, top]).bndry_simp}; 
%         
%         if ~isempty(remaining_hs)
%             color(end, :) = [255 0 51]./255;                                   %/ taret bndry color
%         else
%             color         = [255 0 51]./255;                                   %/ taret bndry color
%         end
%         marker = []; markersize = []; markerfacecolor = [0 0 0]; linewi = 1.5;
%     else
%         str_hs = '_global';
%         bndry_data = []; marker = []; markersize = []; markerfacecolor = []; linewi = []; color = [];
%     end
%     
%     contf_data_phase1 = {}; contf_data_phase2 = {};
%     cont_data_phase1  = {}; cont_data_phase2  = {};
%     hatch_data_phase1 = {}; hatch_data_phase2 = {};
%     Udata_phase1 = {};      Udata_phase2 = {};
%     Vdata_phase1 = {};      Vdata_phase2 = {};
%     
%     %====== data loading and pre-processing ======%
%     if ~isempty(slct_vector_str)
%         var_set = {slct_contf_str, slct_cont_str, slct_hatch_str, strcat('u',slct_vector_str), strcat('v',slct_vector_str)};
%         var_set(cellfun(@isempty, var_set)) = [];     %/ remove empty cell
% 
%         str_var_set = {slct_contf_str, slct_cont_str, slct_hatch_str, strcat('uv',slct_vector_str)}; %/ for titlename/figname
%         str_var_set(cellfun(@isempty, str_var_set)) = [];
%         str_var_set = strjoin(str_var_set, '_');
%     else
%         var_set = {slct_contf_str, slct_cont_str, slct_hatch_str};
%         var_set(cellfun(@isempty, var_set)) = [];
%         str_var_set = strjoin(var_set, '_');
%     end
% 
%     select_data = findismember_loop(dataname, var_set)';
%     if length(select_data) ~= length(var_set)       error('The input data is not all found in dataname!');  end
% 
%     for k = select_data
%         if any(ismember(dataname{k}, {'contr_map', 'BL_Pm'}))
%             if isfield(WSV, 'hs_name')
%                 if isequal(WSV.hs_name, AYR_hotspot_bc(top).name)    %/ make sure the contr_map is related to the requested hotspot!!!
%                     fprintf('!!! The request WSV %s data has been loaded. !!! \n', dataname{k});
%                     WSV_for_hs_exist = 1;
%                 else
%                     WSV_for_hs_exist = 0;
%                 end
%             else
%                 WSV_for_hs_exist = 0;
%             end
% 
%             if WSV_for_hs_exist == 0
%                 str_season = [];
%                 Pm_data_filename = string(strcat(masterfolder, expmnt, '/prcssd_data_4plotting/',...
%                                           'WSV_', dataname{k}, str_hs, strrep(str_expmntinfo, ' ', '_'), str_years, str_season, '.mat'));
%                 fprintf('*** Loading %s ... ***\n', Pm_data_filename)
%                 load(Pm_data_filename);
%             end
%             
%             %/ daily to monthly
%             X             = WSV.(strcat(slct_contf_str,'_daily'));
%             X_lon         = WSV.lon;
%             X_lat         = WSV.lat;
%             X_dates       = WSV.slct_dates;
%             
%         elseif ismember(dataname{k}, {'traj_freq', 'total_init_contr', 'avg_init_contr'})
%             X_dates = nan(length(year_list)*12, 1);
%             X       = nan(length(lon), length(lat), length(year_list)*12);
% 
%             for y = 1:length(year_list)
%                 years = year_list(y);    %/ since the output traj for each hotspot can still be very large -> Save data year by year.
% 
%                 for mth = 1:12
%                     if mth == 0      str_date = num2str(years);  
%                     else             str_date = num2str(years*1e2 + mth);    end
% 
%                     data_filename = string(strcat(masterfolder, expmnt, '/prcssd_data_4plotting/',...
%                                                 slct_contf_str, str_hs, strrep(str_expmntinfo, ' ', '_'), '_', str_date, '.mat'));
% 
%                     tic; fprintf('*** Loading: %s *** \n', data_filename); toc;
%                     L = par_load(data_filename, strcat(slct_contf_str, '_map'));           %/ output is trajs_freq_map
% 
%                     X_dates(mth+(y-1)*12) = years*1e2 + mth;
%                     X(:,:,mth+(y-1)*12) = L;
%                 end
%             end
%             X_lon = lon;
%             X_lat = lat;
%         else
%             select_fields = 'daily';
% 
%             if isfield(dataset, dataname{k})     
%                 fprintf('!!! The requested variable %s exists! Skip loading. !!!\n', dataname{k});
%             else
%                 if ismember(dataname{k}, {'E', 'P', 'EmPSST', 'E_minus_P', 'SST', 'T2m', 'uIVT', 'vIVT', 'IVTdiv', 'Z850'})
%                     filepath = sprintf('/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/%s_%s_global_1970-2010.mat', dataname{k}, select_fields);
%                 else
%                     filepath = sprintf('/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/%s_%s_global%s.mat', dataname{k}, select_fields, str_years);
%                 end
%                 fprintf('*** Loading %s ... ***\n', filepath)
%                 load(filepath)
%                 flds = {select_fields, 'lon', 'lat', 'date_yyyymmdd_AllYr'};
%                 for i = 1:length(flds)
%                     dataset.(dataname{k}).(flds{i}) = S.(flds{i});
%                 end
%                 S = [];
%             end
%             X       = dataset.(dataname{k}).(select_fields);
%             X_lon   = dataset.(dataname{k}).lon;
%             X_lat   = dataset.(dataname{k}).lat;
%             X_dates = dataset.(dataname{k}).date_yyyymmdd_AllYr;
%         end
%         
%         %/ daily to monthly & detrend
%         [X, X_dates] = conv2monthly('data', X, 'dates', X_dates, 'select_year', year_list, 'skip_conv', skip_conv, 'detrend_mode', detrend_mode);
%         
%         %/ deseason
%         if deseason_mode    X = deseasonalize(X, X_dates);            end
%         
%         %/ Subset the events for the given years
%         ind = findismember_loop(event_phase(:,1), X_dates);
%         if isempty(ind)     error('ind is empty! check the dates!');  end
%         event_phase_bc = event_phase(ind,:);
% 
%         L_allphase     = zeros(length(slct_phase1), 1);
%         NLLO_allphase  = zeros(length(slct_phase1), 1);
%         str_plot_phase = {};
%         
%         for i = 1:length(slct_phase1)
%             code_phase1 = str_event_code(findismember_loop(str_event_phase, slct_phase1(i)));
%             code_phase2 = str_event_code(findismember_loop(str_event_phase, slct_phase2(i)));
%             
%             %/ time indices for event_phase
%             ind1 = find(event_phase_bc(:,2) == code_phase1);
%             ind2 = find(event_phase_bc(:,2) == code_phase2);
% 
%             %/ Correct time indices for X (often screw this up!)
%             ind1_for_X = findismember_loop(X_dates, event_phase_bc(ind1,1));
%             ind2_for_X = findismember_loop(X_dates, event_phase_bc(ind2,1));
% 
%             if ismember(dataname{k}, slct_contf_str)
%                 contf_data_phase1{i} = X(:,:,ind1_for_X);
%                 contf_data_phase2    = X(:,:,ind2_for_X);
%                 contf_lon  = X_lon;
%                 contf_lat  = X_lat;
%                 
%             elseif ismember(dataname{k}, slct_cont_str)
%                 cont_data_phase1{i} = X(:,:,ind1_for_X);
%                 cont_data_phase2    = X(:,:,ind2_for_X);
%                 cont_lon  = X_lon;
%                 cont_lat  = X_lat;
%                 
%             elseif ismember(dataname{k}, slct_hatch_str)
%                 hatch_data_phase1{i} = X(:,:,ind1_for_X);
%                 hatch_data_phase2    = X(:,:,ind2_for_X);
%                 hatch_lon  = X_lon;
%                 hatch_lat  = X_lat;
%                 
%             elseif ismember(dataname{k}, {strcat('u',slct_vector_str)})
%                 Udata_phase1{i} = X(:,:,ind1_for_X);
%                 Udata_phase2    = X(:,:,ind2_for_X);
%                 uv_lon  = X_lon;
%                 uv_lat  = X_lat;
%                 
%             elseif ismember(dataname{k}, {strcat('v',slct_vector_str)})
%                 Vdata_phase1{i} = X(:,:,ind1_for_X);
%                 Vdata_phase2    = X(:,:,ind2_for_X);
%             end
%             str_plot_phase{i} = sprintf('%s(%d) - %s(%d)', slct_phase1{i}, length(ind1_for_X), slct_phase2{i}, length(ind2_for_X));
%             
%             %/ obtain the table for basin-level difference (fwd contr_map)
%             if ismember(slct_contf_str, {'contr_map'})
%                 if SR_mode
%                     %/ 17 hotspots (regardless of hotspot_list_ver)
%                     if hotspot_list_ver == 2
%                         hs_ind{1} = [1, 4, 6, 7, 15];           %/ South America
%                         hs_ind{2} = [2, 3, 5, 8, 10];           %/ Africa (skipped Shebelle)
%                         hs_ind{3} = [11, 12, 14, 16, 17];       %/ Asia   (skipped New Guinea)
%                     else
%                         hs_ind{1} = [1, 4, 6, 7, 15];           %/ South America
%                         hs_ind{2} = [2, 3, 5, 8, 9, 10];        %/ Africa
%                         hs_ind{3} = [11, 12, 14, 16, 17];       %/ Asia (skipped New Guinea)
%                     end
%                     flag_to_skip = 1;
%                     for ii = 1:3
%                         if any(ismember(hs_ind{ii}, top))
%                             remaining_hs  = setdiff(hs_ind{ii}, top);
%                             flag_to_skip = 0;
%                         end
%                     end
%                     if flag_to_skip
%                        remaining_hs  = [];
%                        warning('The SR network of %s is skipped!', AYR_hotspot_bc(top).name)
%                     end
% 
%                     if ~isempty(remaining_hs)
%                         NLL_fldname = {AYR_hotspot_bc(remaining_hs).name};
%                         NLL_fldname = strrep(strrep(strrep(NLL_fldname, ' ', '_'), '.', ''), '-', '_');
%                     else
%                         NLL_fldname  = {'NLL'};                        %/ NLL for New Guinea
%                     end
%                 else
%                     NLL_fldname = {'NLL'};
%                 end
%                 NLL_fldname_all{top} = NLL_fldname;                    %/ store for later use
% 
%                 %/ loop over regression of contr. to NLL
%                 L_diff_sig_char = {}; NLL_diff_sig_char = {}; L_diff_sig_char_onfig = {}; NLL_diff_sig_char_onfig = {};
%                 L_diff_sig_double    = nan; 
%                 NLLO_diff_sig_double = nan(1, length(NLL_fldname)); 
%                 for ff = 1:length(NLL_fldname)
% 
%                     %/ forward contr. from the target hotspot
%                     L   = AYR_hotspot_bc(top).L_amount_mm_wgt;
%                     NLL = AYR_hotspot_bc(top).(strcat(NLL_fldname{ff}, '_amount_mm_wgt'));  %/ non-local land
%                     
%                     %/ NOT to convert into monthly & detrend
%                     date = AYR_hotspot_bc(top).date;
%                     [L,  ~]  = conv2monthly('data', L,   'dates', date, 'select_year', year_list, 'skip_conv', skip_conv, 'detrend_mode', detrend_mode);
%                     [NLL, ~] = conv2monthly('data', NLL, 'dates', date, 'select_year', year_list, 'skip_conv', skip_conv, 'detrend_mode', detrend_mode);
% 
%                     %/ forward contr. from the nearby hotspots to the target hotspot
%                     hs_name_strrep = strrep(strrep(strrep({AYR_hotspot_bc.name}, ' ', '_'), '.', ''), '-', '_');
%                     top_other = findismember_loop(hs_name_strrep, NLL_fldname(ff));
%                     if isempty(top_other)
%                         warning('top_other is empty for %s! Setting NLLO to zeros!!', AYR_hotspot_bc(top).name); 
%                         NLLO = zeros(length(date),1);
%                     else
%                         NLLO = AYR_hotspot_bc(top_other).(strcat(hs_name_strrep{top}, '_amount_mm_wgt'));  %/ NLLO: from other hotspots
%                         
%                         %/ NOT to convert into monthly & detrend
%                         [NLLO, ~] = conv2monthly('data', NLLO, 'dates',  date, 'select_year', year_list, 'skip_conv', skip_conv, 'detrend_mode', detrend_mode);  %/ 1971 - 2010 (extract one additional year to do lag-corr).
%                     end                        
% 
%                     %/ Correct time indices for L and NL (often screw this up!)
%                     ind1_for_LNL = findismember_loop(date, event_phase_bc(ind1,1));
%                     ind2_for_LNL = findismember_loop(date, event_phase_bc(ind2,1));
% 
%                     L_1 = L(ind1_for_LNL);
%                     L_2 = L(ind2_for_LNL);
% 
%                     NLL_1 = NLL(ind1_for_LNL);
%                     NLL_2 = NLL(ind2_for_LNL);
% 
%                     NLLO_1 = NLLO(ind1_for_LNL);
%                     NLLO_2 = NLLO(ind2_for_LNL);
%                     
%                     L_allphase(i)    = mean(L_1,'omitnan');                         %/ daily mean of Local recycling in phase x 
%                     NLLO_allphase(i) = NLLO_allphase(i) + mean(NLLO_1,'omitnan');   %/ daily mean of nonLocal recycling from nearby hotspots in phase x 
% 
%                     alpha_table = [0.05, 0.01]; %/ 95% 99% of confidence
%                     str_alpha = {'*', '**'};
% 
%                     L_diff_sig_char = '';       NLL_diff_sig_char{ff} = ''; 
%                     L_diff_sig_char_onfig = ''; NLL_diff_sig_char_onfig{ff} = '';
%                     for ii = 1:length(alpha_table)
%                         [L_diff_sig,    L_diff]    = ttest2_sig_fn(L_1,    L_2,    alpha_table(ii), 1);
%                         [NLL_diff_sig,  NLL_diff]  = ttest2_sig_fn(NLL_1,  NLL_2,  alpha_table(ii), 1);
%                         [NLLO_diff_sig, NLLO_diff] = ttest2_sig_fn(NLLO_1, NLLO_2, alpha_table(ii), 1);
% 
%                         if ~isnan(L_diff_sig)
%                             L_diff_sig_char       = sprintf('%.2g%% (%.2g)%s', L_diff_sig/mean(L_2)*100, L_diff_sig, str_alpha{ii});
%                             L_diff_sig_char_onfig = sprintf('%.2g%% (%.2g)%s', L_diff_sig/mean(L_2)*100, L_diff_sig, str_alpha{ii});    %/ add the unit (%) to the text on fig
%                             L_diff_sig_double     = L_diff_sig;
%                         end
%                         if ~isnan(NLL_diff_sig)
%                             NLL_diff_sig_char{ff}       = sprintf('%.2g%% (%.2g)%s', NLL_diff_sig/mean(NLL_2)*100, NLL_diff_sig, str_alpha{ii});
%                             NLL_diff_sig_char_onfig{ff} = '';              %/ not show NLL on fig.
%                         end
%                         if ~isnan(NLLO_diff_sig)
%                             NLLO_diff_sig_double(ff) = NLLO_diff_sig;
%                         end
%                     end
%                 end
% 
%                 if top == 1     st_ind = 1;
%                 else            st_ind = length([NLL_fldname_all{1:top-1}])+top-1+1;        end
%                 ind = st_ind:st_ind+length(NLL_fldname_all{top});
% 
%                 T       (ind, i)    = [{L_diff_sig_char}, NLL_diff_sig_char]';
%                 T_onfig (top, i)    = {L_diff_sig_char_onfig};
%                 T_double(top, i, 1) = L_diff_sig_double;
%                 T_double(top, i, 2) = nansum(NLLO_diff_sig_double);
%             end
%             
%             %/ obtain the table for basin-level difference (bwd BL_Pm)
%             if ismember(slct_contf_str, {'BL_Pm'})
% 
%                 %/ loop over regression of contr. to NLL
%                 L_diff_sig_char = {}; NLL_diff_sig_char = {}; L_diff_sig_char_onfig = {}; NLL_diff_sig_char_onfig = {};
%                 L_diff_sig_double    = nan; 
%                 NLLO_diff_sig_double = nan(1, length(NLL_fldname)); 
%                 for ff = 1:length(NLL_fldname)
% 
%                     %/ backward contribution from the target hotspot
%                     L   = AYR_hotspot_bc(top).L_amount;    %/ local
%                     NLL = AYR_hotspot_bc(top).NLL_amount;  %/ non-local land
%                     NLO = AYR_hotspot_bc(top).NLO_amount;  %/ non-local ocean
%                     
%                     %/ NOT to convert into monthly & detrend
%                     date = AYR_hotspot_bc(top).date;
%                     [L,  ~]  = conv2monthly('data', L,   'dates', date, 'select_year', year_list, 'skip_conv', skip_conv, 'detrend_mode', detrend_mode);
%                     [NLL, ~] = conv2monthly('data', NLL, 'dates', date, 'select_year', year_list, 'skip_conv', skip_conv, 'detrend_mode', detrend_mode);
% 
%                     %/ forward contr. from the nearby hotspots to the target hotspot
%                     hs_name_strrep = strrep(strrep(strrep({AYR_hotspot_bc.name}, ' ', '_'), '.', ''), '-', '_');
%                     top_other = findismember_loop(hs_name_strrep, NLL_fldname(ff));
%                     if isempty(top_other)
%                         warning('top_other is empty for %s! Setting NLLO to zeros!!', AYR_hotspot_bc(top).name); 
%                         NLLO = zeros(length(date),1);
%                     else
%                         NLLO = AYR_hotspot_bc(top_other).(strcat(hs_name_strrep{top}, '_amount_mm_wgt'));  %/ NLLO: from other hotspots
%                         
%                         %/ NOT to convert into monthly & detrend
%                         [NLLO, ~] = conv2monthly('data', NLLO, 'dates',  date, 'select_year', year_list, 'skip_conv', skip_conv, 'detrend_mode', detrend_mode);  %/ 1971 - 2010 (extract one additional year to do lag-corr).
%                     end                        
% 
%                     %/ Correct time indices for L and NL (often screw this up!)
%                     ind1_for_LNL = findismember_loop(date, event_phase_bc(ind1,1));
%                     ind2_for_LNL = findismember_loop(date, event_phase_bc(ind2,1));
% 
%                     L_1 = L(ind1_for_LNL);
%                     L_2 = L(ind2_for_LNL);
% 
%                     NLL_1 = NLL(ind1_for_LNL);
%                     NLL_2 = NLL(ind2_for_LNL);
% 
%                     NLLO_1 = NLLO(ind1_for_LNL);
%                     NLLO_2 = NLLO(ind2_for_LNL);
%                     
%                     L_allphase(i)    = mean(L_1,'omitnan');                         %/ daily mean of Local recycling in phase x 
%                     NLLO_allphase(i) = NLLO_allphase(i) + mean(NLLO_1,'omitnan');   %/ daily mean of nonLocal recycling from nearby hotspots in phase x 
% 
%                     alpha_table = [0.05, 0.01]; %/ 95% 99% of confidence
%                     str_alpha = {'*', '**'};
% 
%                     L_diff_sig_char = '';       NLL_diff_sig_char{ff} = ''; 
%                     L_diff_sig_char_onfig = ''; NLL_diff_sig_char_onfig{ff} = '';
%                     for ii = 1:length(alpha_table)
%                         [L_diff_sig,    L_diff]    = ttest2_sig_fn(L_1,    L_2,    alpha_table(ii), 1);
%                         [NLL_diff_sig,  NLL_diff]  = ttest2_sig_fn(NLL_1,  NLL_2,  alpha_table(ii), 1);
%                         [NLLO_diff_sig, NLLO_diff] = ttest2_sig_fn(NLLO_1, NLLO_2, alpha_table(ii), 1);
% 
%                         if ~isnan(L_diff_sig)
%                             L_diff_sig_char       = sprintf('%.2g%% (%.2g)%s', L_diff_sig/mean(L_2)*100, L_diff_sig, str_alpha{ii});
%                             L_diff_sig_char_onfig = sprintf('%.2g%% (%.2g)%s', L_diff_sig/mean(L_2)*100, L_diff_sig, str_alpha{ii});    %/ add the unit (%) to the text on fig
%                             L_diff_sig_double     = L_diff_sig;
%                         end
%                         if ~isnan(NLL_diff_sig)
%                             NLL_diff_sig_char{ff}       = sprintf('%.2g%% (%.2g)%s', NLL_diff_sig/mean(NLL_2)*100, NLL_diff_sig, str_alpha{ii});
%                             NLL_diff_sig_char_onfig{ff} = '';              %/ not show NLL on fig.
%                         end
%                         if ~isnan(NLLO_diff_sig)
%                             NLLO_diff_sig_double(ff) = NLLO_diff_sig;
%                         end
%                     end
%                 end
% 
%                 if top == 1     st_ind = 1;
%                 else            st_ind = length([NLL_fldname_all{1:top-1}])+top-1+1;        end
%                 ind = st_ind:st_ind+length(NLL_fldname_all{top});
% 
%                 T       (ind, i)    = [{L_diff_sig_char}, NLL_diff_sig_char]';
%                 T_onfig (top, i)    = {L_diff_sig_char_onfig};
%                 T_double(top, i, 1) = L_diff_sig_double;
%                 T_double(top, i, 2) = nansum(NLLO_diff_sig_double);
%             end
%         end
%         
%         if ismember(slct_contf_str, {'contr_map'})
%             %/ Convert into fractional change by the sum of contr. from local & non-local hotspots
%             total_L_NLLO = L_allphase + NLLO_allphase;        %/ 8 x 1
%             total_L_NLLO = reshape(total_L_NLLO, 1, [], 1);   %/ reshape to  1 x 8 x 1
% 
%             T_double_fraction(top,:,:) = T_double(top,:,:)./(total_L_NLLO)*100; %/ (1 x 8 x 2 ) ./ (1 x 8 x 1)
%         end
%         X = [];  %/ release memory
%     end
% 
%     if plot_traj
%         rand_traj_filename = string(strcat(masterfolder, expmnt, '/prcssd_data_4plotting/',...
%                                         'randtraj', num2str(rand_subset_portion), str_hs, strrep(str_expmntinfo, ' ', '_'),...
%                                         str_years, '.mat'));
%                                     
%         if isfile(rand_traj_filename)
%             fprintf('*** Required traj data already exists (rand: %.4f, %s), loading from it... ***\n', rand_subset_portion, str_years);
%             rand_traj_catalog = par_load(rand_traj_filename, 'rand_traj_catalog');
%             
%         else
%             trajs               = cell(length(year_list)*12, 1);
%             trajs_dates_alltraj = cell(length(year_list)*12, 1);
% 
%             for y = 1:length(year_list)
% 
%                 years = year_list(y);    %/ since the output traj for each hotspot can still be very large -> Save data year by year.
%                 trajs_y               = cell(12, 1);
%                 trajs_dates_alltraj_y = cell(12, 1);
% 
%                 if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
%                     parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
%                 end
% 
%                 tic;
%                 parfor mth = 1:12
%                     if mth == 0      str_date = num2str(years);  
%                     else             str_date = num2str(years*1e2 + mth);    end
% 
%                     data_filename = string(strcat(masterfolder, expmnt, '/prcssd_data_4plotting/',...
%                                                 'traj', str_hs, strrep(str_expmntinfo, ' ', '_'), '_', str_date, '.mat'));
% 
%                     fprintf('*** Loading: %s *** \n', data_filename)
%                     trajs_catalog = par_load(data_filename, 'trajs_catalog');           %/ output is trajs_catalog
% 
%                     %/ randomly select only 1% of trajs for processing. (lessen memory burden)
%                     ntraj = size(trajs_catalog{1}, 1);
% 
%                     trajs_dates = fix(str2double(trajs_catalog{2}));
%                     trajs_dates_ntraj = cell2mat(trajs_catalog{3});
% 
%                     trajs_dates_alltraj_eachmth = nan(size(trajs_catalog{1},1),1);
%                     for i = 1:length(trajs_dates)
%                         if i == 1               ind = 1:trajs_dates_ntraj(i);
%                         else                    ind = (sum(trajs_dates_ntraj(1:i-1))+1):sum(trajs_dates_ntraj(1:i));
%                         end
%                         trajs_dates_alltraj_eachmth(ind) = trajs_dates(i);
%                     end
% 
%                     if ~isempty(rand_subset_portion)
% 
%                         rng(seed);                                                              %/ rng() has to run right before the randi().
%                         ind_rand   = randperm(ntraj, fix(ntraj*rand_subset_portion))';          %/ use randperm() instead of randi() to obtain unique random number!       
%                         trajs_y{mth}                = trajs_catalog{1}(ind_rand, :, :);
%                         trajs_dates_alltraj_y{mth}  = trajs_dates_alltraj_eachmth(ind_rand);
% 
%                         fprintf('*** Randomly subset %.2f%% (%d) of the traj of %d ... ***\n', rand_subset_portion*100, length(ind_rand), years*1e2+mth)
%                     end
%                 end
%                 toc;
% 
%                 trajs((1:12)+12*(y-1))               = trajs_y;
%                 trajs_dates_alltraj((1:12)+12*(y-1)) = trajs_dates_alltraj_y;
% 
%             end
%             trajs               = cat(1, trajs{:});
%             trajs               = permute(trajs, [1 3 2]);
%             trajs_dates_alltraj = cat(1, trajs_dates_alltraj{:});
%             
%             rand_traj_catalog = {trajs, trajs_dates_alltraj};
%             if savemat
%                 fprintf('** Saving clustering results into %s ***\n', rand_traj_filename);
%                 tic; save(rand_traj_filename, 'rand_traj_catalog', '-v7.3'); toc;
%             end
%         end
%         
%         trajs               = rand_traj_catalog{1};
%         trajs_dates_alltraj = rand_traj_catalog{2};
%         clear rand_traj_catalog
%         
%         %/ check nan
%         [ind1, ind2, ind3] = find(isnan(trajs));
%         if ~isempty(ind1)
%             nan_date = unique(trajs_dates_alltraj(ind1));
%             for i = 1:length(nan_date)
%                 if mod(nan_date(i), 1e8) == 12162100 || nan_date(i) == 201012152100
%                     warning('NaN detected in trajs on %s, likely due to the NaNs on 12312100 or 12302100 in 2010 by FLEXPART experiment', string(nan_date(i)))
%                 else
%                     error('NaN detected in trajs on %s. Why??', string(nan_date(i)))
%                 end
%             end
%         end
% 
%         %/ subset traj conditioned on climate events
%         event_dates = event_phase(:,1);
%         
%         dates_pve = event_dates(event_phase(:,2) == 1);
%         dates_nve = event_dates(event_phase(:,2) == -1);
%         dates_neu = event_dates(event_phase(:,2) == 0);
%                 
%         ind_traj_pve = findismember_loop(floor(trajs_dates_alltraj/1e6), dates_pve); 
%         ind_traj_nve = findismember_loop(floor(trajs_dates_alltraj/1e6), dates_nve); 
%         ind_traj_neu = findismember_loop(floor(trajs_dates_alltraj/1e6), dates_neu); 
% 
%         dates_pve = unique(floor(trajs_dates_alltraj(ind_traj_pve)/1e6));
%         dates_nve = unique(floor(trajs_dates_alltraj(ind_traj_nve)/1e6));
%         dates_neu = unique(floor(trajs_dates_alltraj(ind_traj_neu)/1e6));
%         
%         %/ Dorling's CC
%         ncluster = 50;
% %         ncluster = lip(50:50:500);
%         seed     = 129;                                                    %/ do not change the seed.
%         max_iter = 10000;
% %         loss_threshold = 0.01;                                            %/ 1%
%         loss_threshold = 0.005;                                            %/ 0.5%
%         
%         str_event_phase      = {'pve', 'nve', 'neu'};
%         ind_traj_phase = {ind_traj_pve, ind_traj_nve, ind_traj_neu};
%         dates_phase    = {dates_pve,    dates_nve,    dates_neu};
%         
%         if length(ncluster) == 1     str_ncluster = num2str(ncluster);
%         else                         str_ncluster = strjoin(string([ncluster(1), abs(diff(ncluster(1:2))), ncluster(end)]), '_'); end
%         
%         Dorling_filename = string(strcat(masterfolder, expmnt, '/prcssd_data_4plotting/',...
%                                         'Dorling_CC_',  str_event_phase, '_n', str_ncluster, '_loss', num2str(loss_threshold),...
%                                         '_randtraj', num2str(rand_subset_portion),...
%                                         str_hs, strrep(str_expmntinfo, ' ', '_'), str_years, '.mat'));
% 
%         for i = 1:length(str_event_phase)
%             if isfile(Dorling_filename{i})
%                 fprintf('** Data already exists. Skipping ***\n');         %/ will read it in later lines.
%             else
%                 fprintf('*** Performing Dorling CC ... ***\n')
%                 [seed_trajs, trajs_clu, RMSD] = DorlingCC('ncluster', ncluster, 'trajs', trajs(ind_traj_phase{i},:,:), 'seed', seed, 'max_iter', max_iter,...
%                                                           'loss_threshold', loss_threshold, 'NumWorkers', NumWorkers);
%                 
%                 Dorling_catalog = {seed_trajs, trajs_clu, RMSD};
%                 if savemat
%                     fprintf('** Saving clustering results into %s ***\n', Dorling_filename{i});
%                     tic; save(Dorling_filename{i}, 'Dorling_catalog', '-v7.3'); toc;
%                 end
%             end
%         end
%     end
% 
%     % -------------- Plot different phases --------------
%     if draw_plot
%         th = []; point_data = []; traj_data = []; cnt = 0;
%         for i = 1:length(slct_phase1)
%             cnt = cnt + 1;
%             titlename     = [];
%             FigName = sprintf('%s %s%s%s%s %d-%d%s', str_slct_event, str_var_set, str_traj, str_hs, str_sig, floor(year_list(1)), floor(year_list(end)), str_text);
%             FigName_underscore = strrep(FigName, ' ', '_');
%             if draw_cbar_only
%                 create_fig = 1;
%                 savepath   = char(strcat(plotting_folder, FigName_underscore)); 
%                 if cnt ~= 1   break;   end
%             else
%                 create_fig = 0;
%                 savepath   = []; 
%             end
%             
%             grid_mode = -1;
%             fontsize = 8;
%             panel_row = 8; panel_col = 1;
%             p1_size       = 0.02;
%             p22_size      = 1e-10;
%             map_proj      = 'Miller Cylindrical'; 
%             if global_map
%                 map_lon_lower = -179;
%                 map_lon_upper = 180;
%                 map_lat_lower = -40;
%                 map_lat_upper = 40;
% 
%                 marginleft    = 15;
%                 marginright   = 15;
%                 margintop     = 3;
%                 marginbottom  = 3;
% 
%             else
%                 cen = AYR_hotspot(top).centroid;
%                 map_lon_lower = cen(1) - 50;
%                 map_lon_upper = cen(1) + 50;
%                 map_lat_lower = cen(2) - 25;
%                 map_lat_upper = cen(2) + 15;
% 
%                 marginleft    = 15;
%                 marginright   = 15;
%                 margintop     = 3;
%                 marginbottom  = 3;
%             end
% 
%             if cnt == 1                                                          %/ initiate figure
%                 close all;
%                 fig = figure;
%                 p = panel(fig);
%                 p.pack(2, 2);
%                 p(1).repack(p1_size);                                          %/ placeholder to avoid over-cutting at the top.
%                 p(2,1).pack(panel_row, panel_col);
%                 p(2,2).repack(p22_size);                                          %/ placeholder to avoid over-cutting on the right.
%                 set(gcf, 'Position', get(0, 'Screensize'));                    %/ set the figure position & size before drawing! Otherwise vector plot is incorrect!
%                 set(gcf, 'color', 'w');                                        %/ set figure bg color to be white
%                 m = 0; n = 0;
%             end
%             if mod(m, panel_row) == 0                  m = 1; n = n + 1;
%             else                                       m = m + 1;              end
%             p(2, 1, m, n).select();
% 
%             %======= contf =======%
%             if ~isempty(slct_contf_str)
%                 %/ ttest2
%                 X1 = contf_data_phase1{i};
%                 X2 = contf_data_phase2;
%                 [data_diff_sig, data_diff] = ttest2_sig_fn(X1, X2, alpha, 3);
% 
%                 contf_data = data_diff_sig;
%                 point_data = [];
% 
%     %             if sig_or_both == 2
%     %                 contf_data = data_diff;
%     %                 
%     %                 [ind_lon_sig, ind_lat_sig] = find(~isnan(data_diff_sig));
%     %                 point_data = nan(length(ind_lon_sig), 2);
%     %                 for k = 1:length(ind_lon_sig)
%     %                     point_data(k,1) = X_lon(ind_lon_sig(k));
%     %                     point_data(k,2) = X_lat(ind_lat_sig(k));
%     %                 end
%     %             end
% 
%                 colmap_name = '*RdYlBu';
%                 if ismember(slct_contf_str, {'contr_map'})
%                     contf_levels = -1.1:0.1:1.1;
%                     contf_unit   = '';
% %                     contf_unit   = 'mm/day';
%                     cbar_YTick      = -1:0.5:1;
%                     cbar_YTickLabel = cbar_YTick;
% 
%                 elseif ismember(slct_contf_str, {'P', 'E'})
%                     contf_levels = -3.3:0.3:3.3;
%                     contf_unit   = 'mm/day';
%                     cbar_YTick      = -3:1:3;
%                     cbar_YTickLabel = cbar_YTick;
%                     colmap_name = 'BrBG';
% 
%                 elseif ismember(slct_contf_str, {'traj_freq'})
%                     contf_levels = [-10:1:10]*5e2;
%                     contf_unit   = 'Freq.';
% 
%                 elseif ismember(slct_contf_str, {'total_init_contr'})
%                     contf_levels = [-10:1:10]*5e1;
%                     contf_unit   = 'g/kg';
% 
%                 elseif ismember(slct_contf_str, {'avg_init_contr'})
%                     contf_levels = [-1:0.1:1];
%                     contf_unit   = 'g/kg';
%                 end
%                 NoOfColor = length(contf_levels)-1;
%                 colmap = brewermap(NoOfColor+2, colmap_name);
%                 colmap(NoOfColor/2+1:NoOfColor/2+2,:) = [];  
% %                 colmap = brewermap(NoOfColor, colmap_name);
% %                 colmap(NoOfColor/2:NoOfColor/2+1,:) = [1 1 1; 1 1 1];  
%             end
% 
%             %======= cont =======%
%             if ~isempty(slct_cont_str)
%                 %/ttest2
%                 X1 = cont_data_phase1{i};
%                 X2 = cont_data_phase2;
%                 [data_diff_sig, data_diff] = ttest2_sig_fn(X1, X2, alpha, 3);
% 
%                 if sig_or_both == 1
%                     cont_data = data_diff_sig;
%                 else
%                     cont_data     = data_diff_sig;
%                     cont_data_raw = data_diff;
%                 end
% 
%                 if ismember(slct_cont_str, {'Z850'})
%                     cont_levels   = [-24:4:24];          %/ in m
%     %                 cont_levels   = [-20:2:20];          %/ in m
%                     cont_levels(cont_levels == 0) = [];
% 
%                 elseif ismember(slct_cont_str, {'P', 'E_minus_P', 'E'})
%                     cont_levels   = -5:1:5;            %/ in mm/day
%                     cont_levels(cont_levels == 0) = [];
%                 else
%                     error('Not yet set the code!')
%                 end
%                 cont_colmap    = [repmat([153 0 255], length(cont_levels)/2, 1); repmat([0 204 102], length(cont_levels)/2, 1)]./255;
%                 cont_labelsize = fontsize;
%                 cont_linewi  = 1;
% %                 cont_linewi  = 1.5;
%             end
% 
%             %======= hatch (w/ +ve and -ve indication) ======%
%             if ~isempty(slct_hatch_str)
%                 %/ ttest2
%                 X1 = hatch_data_phase1{i};
%                 X2 = hatch_data_phase2;
%                 [data_diff_sig, data_diff] = ttest2_sig_fn(X1, X2, alpha, 3);
% 
%                 hatch_data = data_diff_sig;
%                 
%                 hatch_mode = 1;
%                 hatch_linewi = 1; %1.5;
%                 hatch_intvl = 10; 
%                 if ismember(slct_hatch_str, {'P'})
%                     hatch_thres = 0;  %/ mm/day
% %                     hatch_thres = 0.5;  %/ mm/day
%                     hatch_unit = 'mm/day';
% 
%                 elseif ismember(slct_hatch_str, {'E'})
%                     hatch_thres = 0;  %/ mm/day       %/ for testing
% %                     hatch_thres = 0;  %/ mm/day       %/ for testing
%                     hatch_unit = 'mm/day';
%                 else
%                     error('Not yet set the code!')
%                 end
%                 
%                 fprintf('*** Mean magnitude of hatch data: %.2f %s. hatch_thres = %.2f %s ***\n',...
%                             mean(reshape(abs(hatch_data),[],1), 'omitnan'), hatch_unit, hatch_thres, hatch_unit)
%             end
%             
%             %======= vector =======%
%             if ~isempty(slct_vector_str)
%                 %/ ttest2
%                 X1 = Udata_phase1{i};
%                 X2 = Udata_phase2;
%                 [Udata_sig, Udata] = ttest2_sig_fn(X1, X2, alpha, 3);
% 
%                 X1 = Vdata_phase1{i};
%                 X2 = Vdata_phase2;
%                 [Vdata_sig, Vdata] = ttest2_sig_fn(X1, X2, alpha, 3);
% 
%                 %/ Restore vector data if one of the component is sig (IMPORTANT!!)
%                 [nlon, nlat] = size(Udata);
%                 Udata_sig_restore = zeros(nlon, nlat);
%                 Vdata_sig_restore = zeros(nlon, nlat);
%                 for xx = 1:nlon
%                 for yy = 1:nlat
%                     if ~isnan(Udata_sig(xx, yy)) || ~isnan(Vdata_sig(xx, yy)) % either U or V is sig, then retore its values in UV
%                         Udata_sig_restore(xx, yy) = Udata(xx, yy);
%                         Vdata_sig_restore(xx, yy) = Vdata(xx, yy);
%                     end
%                 end
%                 end
% 
%                 %/ final UV data 
%                 Udata = Udata_sig_restore;
%                 Vdata = Vdata_sig_restore;
% 
%                 vector_color = [0 102 255]./255; %blue 
%                 vector_edgecolor = 'w';
% 
%                 if global_map
%                     vec_step_lon = 20;
%                     vec_step_lat = 10;
%                     vecscale     = 200;            % the smaller value the bigger vector. for winds (keep it below 300 to have nice-looking vectors)
%                     vecscale2    = 0.5;            % control shaft length. (adjust to have nice-looking vectors)
%                     shaftwidth   = 1.5;            % control width of the shaft, the larger value the thicker
%                     headlength   = shaftwidth*3;   % control length of the head, the larger value the longer
%                 else
%                     vec_step_lon = 12;
%                     vec_step_lat = 6;
%                     vecscale     = 300;            % the smaller value the bigger vector. for winds
%                     vecscale2    = 1;              % control shaft length.
%                     shaftwidth   = 1.2;            % control width of the shaft, the larger value the thicker
%                     headlength   = shaftwidth*4;   % control length of the head, the larger value the longer
%                 end
% 
%                 if cnt ~= length(slct_phase1)
%                     vec_mag_ref = [];
%                 else
%                     %/ vector reference
%                     if ismember(slct_vector_str, {'IVT'})
%                         vec_mag_ref = 100;
%                         vec_lbs = strcat(num2str(vec_mag_ref), {' kg/m/s'});
%                     else
%                         vec_mag_ref = 5; 
%                         vec_lbs = strcat(num2str(vec_mag_ref), {' m/s'});
%                     end
%                 end
%                 slct_vector_str2 = strcat('uv', slct_vector_str);
%             end
%             X1 = [];  X2 = [];
%             cbar_mode = 0;
%             cbar_interval = 3;
%             pcolor_mode = 1;
%             glb_plateau_mode = 1;
%             plateau_hgt = 1500;
%             plateau_col = [255 51 204]./255;
%             coast_col  = [.4 .4 .4];
%             coast_wi   = 1.5;
% 
%             if ismember(AYR_hotspot(top).name, {'New Guinea'})
%                 glb_data_mode = 0;  %/ in order to show data from 0 to 360
%             else
%                 glb_data_mode = 1;  %/ in order to show data from -179 to 180
%             end
%             
%             plot_contfmap('contf_data', contf_data, 'contf_lon', contf_lon, 'contf_lat', contf_lat, 'contf_levels', contf_levels,...
%                       'contf_unit', contf_unit, 'colmap', colmap, 'cbar_interval', cbar_interval, 'pcolor_mode', pcolor_mode,...
%                       'point_data', point_data, 'marker', marker, 'markersize', markersize, 'markerfacecolor', markerfacecolor, 'linewi', linewi, 'color', color,...
%                       'hatch_data', hatch_data, 'hatch_lon', hatch_lon, 'hatch_lat', hatch_lat, 'hatch_thres', hatch_thres, 'hatch_mode', hatch_mode, 'hatch_linewi', hatch_linewi, 'hatch_intvl', hatch_intvl,...
%                       'cont_data', cont_data, 'cont_data_raw', cont_data_raw, 'cont_lon', cont_lon, 'cont_lat', cont_lat, 'cont_levels', cont_levels, 'cont_colmap', cont_colmap, 'cont_linewi', cont_linewi, 'cont_labelsize', cont_labelsize,...
%                       'Udata', Udata, 'Vdata', Vdata, 'uv_lon', uv_lon, 'uv_lat', uv_lat, 'vec_step_lon', vec_step_lon, 'vec_step_lat', vec_step_lat,...
%                       'vector_color', vector_color, 'vecscale', vecscale, 'vecscale2', vecscale2, 'shaftwidth', shaftwidth, 'headlength', headlength, 'vec_lbs', vec_lbs, 'vec_mag_ref', vec_mag_ref,...
%                       'bndry_data', bndry_data, 'titlename', titlename, 'savepath', savepath,...
%                       'glb_data_mode', glb_data_mode, 'glb_plateau_mode', glb_plateau_mode, 'plateau_hgt', plateau_hgt, 'plateau_col', plateau_col,...
%                       'map_proj', map_proj, 'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper,...
%                       'fontsize', fontsize, 'create_fig', create_fig, 'grid_mode', grid_mode, 'cbar_mode', cbar_mode, 'cbar_YTick', cbar_YTick, 'cbar_YTickLabel', cbar_YTickLabel, ...
%                       'draw_cbar_only', draw_cbar_only, 'cbar_location', cbar_location)
%             fprintf('done \n')
%             
%             %/ Add text on map!
%             if LR_NLR_text_mode && ismember(slct_contf_str, {'contr_map'}) && draw_cbar_only == 0
%                 text_fontsize = fontsize + 4;
%                 str_change = '';
%                 if ~isempty(T_onfig{top, i})
%                     str_change   = sprintf(' L: %s', T_onfig{top, i});
%                 end
%                     
%                 if ~isempty(str_change)
%                     h = m_text(map_lon_lower+0.5, map_lat_lower+0.5, str_change,'fontsize',text_fontsize,'color','k', 'verticalalignment','bottom',...
%                             'horizontalalignment','left', 'backgroundcolor', 'w', 'edgecolor', 'k', 'margin', 1);
%                     drawnow;
%                 end
%             end
%             
%             th(cnt) = title(titlename, 'fontsize', fontsize+5, 'fontweight', 'normal');
%             titlePos = get(th(cnt), 'position');
%             titlePos(2) = titlePos(2);   %/ change y position of title
%             set(th(cnt), 'position', titlePos)
%             drawnow; pause(0.05);
%         end
% 
%         p(2).de.marginleft   = marginleft;    
%         p(2).de.marginright  = marginright;
%         p(2).de.marginbottom = marginbottom;
%         p(2).de.margintop    = margintop;
%         drawnow(); pause(0.05);
%         set(th(1:end), 'fontweight', 'normal', 'fontsize', fontsize+3);
%         drawnow(); pause(1);
% 
%         if savefig 
% %             export_fig(char(strcat(plotting_folder, FigName_underscore,'.png')),'-r300','-png','-opengl', '-nocrop');
%             export_fig(char(strcat(plotting_folder, FigName_underscore,'.pdf')),'-pdf','-painters', '-c[inf, inf, inf, inf]', '-nocrop'); %'-transparent');
%         end
%     end
% end
% 
% if save_table
%     
%     %/ Set 0 to nan in case nansum([NaN NaN]) = 0.
%     T_double(T_double == 0)                         = nan;
%     T_double_fraction(T_double_fraction == 0)       = nan;
%     
%     %/ Save T_double
%     matfilename = char(strcat(plotting_folder, 'hs_MJO_diff_double', str_detrend, str_SR, '.mat'));
%     save(matfilename, 'T_double'); %/ do NOT save in '-v7.3'! R cannot read it.
% 
%     %/ Save T_double_fraction
%     matfilename = char(strcat(plotting_folder, 'hs_MJO_diff_double_fraction', str_detrend, str_SR, '.mat'));
%     save(matfilename, 'T_double_fraction'); %/ do NOT save in '-v7.3'! R cannot read it.
%     
% %-------------------------------------------------------------------------%
%     hs_old_ordering = {AYR_hotspot_bc(:).name};
%     hs_new_ordering = {'Amazon', 'Orinoco', 'La Plata', 'NE. South America', 'Tocantins',...
%                        'Congo', 'Nile', 'SE. Africa', 'Zambezi', 'Shebelle', 'Lake Turkana-Swash',...
%                        'Ganges', 'Yangtze', 'Mekong', 'Pearl', 'Irrawaddy-Salween', 'New Guinea'};
% 
%     ind = findismember_loop(hs_old_ordering, hs_new_ordering);   %/ reorder sources by continent.
% 
%     rownames_1 = repmat(hs_new_ordering, 2, 1);
%     rownames_1 = reshape(rownames_1, [], 1);
% 
%     attr_labels = {'Local'; 'Non-local'};
%     rownames_2  = repmat(attr_labels, length(AYR_hotspot_bc), 1);
%     rownames    = strcat(rownames_1, '_', rownames_2);
%     
%     varnames    = strrep(slct_phase1, ' ', '_');
%     
%     %/ Save T_double as csv
%     matfilename = char(strcat(plotting_folder, 'hs_MJO_diff_double', str_detrend, str_SR, '.mat'));
%     load(matfilename);              %/ output is T_double
%     T_double_bc = T_double(ind,:,:);   %/ reorder sources by continent.
% 
%     T_bc = T_MJO_2table(T_double_bc, varnames, rownames);
%     csvfilename = char(strcat(plotting_folder, 'hs_MJO_diff_double', str_detrend, str_SR, '.csv'));
%     fprintf('*** Writing into %s ***\n', csvfilename);
%     writetable(T_bc, csvfilename, 'WriteRowNames',true);
%     
%     %/ Save T_double_fraction as csv
%     matfilename = char(strcat(plotting_folder, 'hs_MJO_diff_double_fraction', str_detrend, str_SR, '.mat'));
%     load(matfilename);                                   %/ output is T_double_fraction
%     T_double_fraction_bc = T_double_fraction(ind,:,:);   %/ reorder sources by continent.
% 
%     T_bc = T_MJO_2table(T_double_fraction_bc, varnames, rownames);
%     csvfilename = char(strcat(plotting_folder, 'hs_MJO_diff_double_fraction', str_detrend, str_SR, '.csv'));
%     fprintf('*** Writing into %s ***\n', csvfilename);
%     writetable(T_bc, csvfilename, 'WriteRowNames',true);
%     
%   
% end
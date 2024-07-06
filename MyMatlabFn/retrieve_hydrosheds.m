function [bndry_data, id_list, basin_name_list, simp_bndry_data] = retrieve_hydrosheds(varargin)

    pnames = {'slct_conti', 'simplify_bndry', 'verbose'};
    dflts  = cell(length(pnames), 1);
    [           slct_conti,  simplify_bndry,   verbose] = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%

    hydrosheds_folder = '/disk/r059/tfchengac/hydrosheds/';
    contin_id = {'Africa', 'Europe', 'Siberia', 'Asia', 'Australia', 'South America', 'North America', 'Arctic', 'Greenland'};
    file_code = {    'af',     'eu',      'si',   'as',       'au' ,            'sa',            'na',     'ar',        'gr'}'; %/ arbitrary ordering, 'na_mex' is made by me.
%     file_code = {'na_mex'); %/ ignore it for now.

%     if isequal(slct_conti, 'global')
%         ind = 1:length(contin_id);
%     else
    ind = findismember_loop(contin_id, slct_conti);
%     end
    
    bndry_data = {};        simp_bndry_data = {}; point_data = []; cnt = 0;
    id_list    = int64([]); basin_name_list = {}; text_data = {};
    for k = ind
        if isequal(file_code{k}, 'na_mex') 
            level = 4;
            basin_file = sprintf('hybas_%s_lev%02d_v1c.shp', 'na', level);
        else
            level = 3;
            basin_file = sprintf('hybas_%s_lev%02d_v1c.shp', file_code{k}, level);
        end
        if verbose  fprintf('*** Reading watershed shapefile: %s *** \n', basin_file);  end
        
        %/ CAVEAT: shapread() uses Map_Toolbox licenses! Avoid calling it repeatedly if
        %/         multiple programs are to be launched!
        [S,A]      = shaperead(strcat(hydrosheds_folder, basin_file),'UseGeoCoords',true);
        basin_name = hydrosheds_name('file_code', file_code{k}, 'level', level);
%         length(S)
        
        basin_name_list = [basin_name_list; basin_name(:,1)]; 
        for i = 1:length(S) 
            cnt = cnt + 1;
            box_reg           = [S(i).Lon; S(i).Lat]';     
            bndry_data(cnt,1) = {box_reg};                                 %/ use the river basin's boundary
            id_list(cnt,1)    = int64(A(i).HYBAS_ID);
            contin{cnt,1}     = contin_id{floor(A(i).HYBAS_ID/1e9)};
            
            %/ Simplify the boundary data (only for those with a length > 3000)!!
            if simplify_bndry
                a = bndry_data{cnt,1};
                
                if length(a) < 3000
                    simp_bndry_data(cnt,1) = {a};
                else
                    ind_nan = find(isnan(a(:,1))); %/ since nans are mean to cut the boundary, we need to keep those nans when simplifying the boundary!
                    ind_nan = [0; ind_nan; length(a)];
                    b = [];
                    for j = 1:length(ind_nan)-1
                        a_subset = a(ind_nan(j)+1:ind_nan(j+1),:);
                        if isempty(a_subset)  continue;   end

                        intvl = fix(length(a_subset)*0.01);  %/ simplied to a size of just ~100
                        if intvl ~= 0
                            b = [b; a_subset(1,:); a_subset(2:intvl:end-1,:); a_subset(end,:)];
                        else
                            b = [b; a_subset];  %/ otherwise just copy it.
                        end
                    end
                    simp_bndry_data(cnt,1) = {b};
                end
            end
            
%             %/ skip the unwanted basins
%             if (i == 1 || i == 23)     && isequal(file_code{k}, 'na')      continue; end   %/ skip the level3 mexico.
%             if ~ismember(i, [1, 2, 5]) && isequal(file_code{k}, 'na_mex')  continue; end   %/ skip those not from mexico.
%             if  ismember(i, [11:13])   && isequal(file_code{k}, 'sa')      continue; end   %/ skip those from South America.

            %/ Get box_reg (i.e., boundaries)
%             if i == 10 && isequal(file_code{k}, 'sa')              
%                 box_reg = box_region({'NESA'});           flag_outline = 1;  %/ Replace those river basins with the customized one
%                 basin_name{i} = 'NE. South America';                                      %/ update basin name
% 
%             elseif i == 1 && isequal(file_code{k}, 'na_mex')
%                 box_reg = box_region({'SMexico'}); 	  flag_outline = 1;
%                 basin_name{i} = 'S. Mexico';                         %/ update basin name
% 
%             elseif i == 2 && isequal(file_code{k}, 'na_mex')
%                 box_reg = box_region({'CMexico'});      flag_outline = 1;
%                 basin_name{i} = 'C. Mexico';                         %/ update basin name
% 
%             elseif i == 5 && isequal(file_code{k}, 'na_mex')
%                 box_reg = box_region({'WNMexico'}); 	  flag_outline = 1;
%                 basin_name{i} = 'WN. Mexico';                        %/ update basin name
% 
%             elseif i == 17 && isequal(file_code{k}, 'au')
%                 %/ original one is too heavy)
%                 box_reg = box_region({'NewGuinea'}); 	  flag_outline = 0;  %/ do not auto-outline.
% 
%             elseif i == 10 && isequal(file_code{k}, 'as')
%                 %/ original one is too heavy)
%                 box_reg = box_region({'Pearl'}); 	      flag_outline = 0;  %/ do not auto-outline.
%             else
%             box_reg = [S(i).Lon; S(i).Lat]';          flag_outline = 0;
%             end
            
            %/ obtain the points from bndry in land-only hotspot var.
%             mask_2D = inpolygon(lon_2D, lat_2D, box_reg(:,1), box_reg(:,2));
%             reg_2D(~mask_2D) = nan;                                  %/ set ocean grids to be nan if outside the box region.

    %         [reg_lat_ind, reg_lon_ind] = find(~isnan(hotspot_var_land));
%             point_data = [point_data; [lon_2D(1,reg_lon_ind)', lat_2D(reg_lat_ind, 1)]];

            %/ create the final logical matrix for outlining boudaries
%             logical_2D = reg_2D;
%             logical_2D(~isnan(logical_2D)) = 1;                                %/ first set nonnan to be 1 (as values can be zeros even in the target region)
%             logical_2D(isnan(logical_2D))  = 0;                                %/ then set nan to be 0

            %/ outline boundaries
%             cnt = cnt + 1;
% %             if flag_outline == 1                                               %/ use the outlined boundary 
% %                 bndry_data(cnt,1) = get_bndry_from_logi('logical_2D', logical_2D', 'bndry_lon', lon_m179_180, 'bndry_lat', WSV.(slct_var).lat, 'outputcell', 1);
% %             else
%             bndry_data(cnt,1) = {box_reg};                                 %/ use the river basin's boundary
% %             end
%             id_list(cnt,1)    = int64(A(i).HYBAS_ID);
%             contin{cnt,1}     = contin_id{floor(A(i).HYBAS_ID/1e9)};
            

        end
    end


end
function [basin_catalog, str_from_basin_traj, str_from_basin] = load_from_basin(varargin)

    % create a set of valid parameters and their default value
    pnames = {'from_basin'};
    dflts  = cell(length(pnames), 1);
    [          from_basin,   ] = internal.stats.parseArgs(pnames, dflts, varargin{:});

%%
    %======================================================================
    %/      Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: 28 Jun 2024
    %======================================================================

    %/ Obtain 'basin_catalog'
    if from_basin == -1             %/ global (-1)
        basin_catalog       = [];
        str_from_basin      = 'glb';
        str_from_basin_traj = 'glb';

    elseif from_basin == 0          %/ global land (0) 
        basin_catalog       = [];
        str_from_basin      = 'glbland';
        str_from_basin_traj = 'glbland';

    elseif from_basin == 1   %/ load the 17 original hotspots 
        catalog_filename = string(strcat('/disk/r059/tfchengac/FLEXPART/', expmnt, '/prcssd_data_4plotting/',...
                                      'AYR_hs_p95_area_sum_BL_Pm_RHc85_dqc0.1_bwd_20d_optimal_3h_glbland_nozbound_intact_RH2_1971-2010.mat'));
        fprintf('*** Loading: %s *** \n', catalog_filename)
        load(catalog_filename, 'AYR_hotspot');   
        basin_catalog = AYR_hotspot;
        str_from_basin = sprintf('AYR_hs%d', length(basin_catalog));
        str_from_basin_traj = 'glbland';
        
    elseif from_basin == 2  %/ load the 3 new hotspots only
        catalog_filename = string(strcat('/disk/r059/tfchengac/FLEXPART/', expmnt, '/prcssd_data_4plotting/',...
                                      'AYR_hs_p95_area_sum_Pm_noOcns_RHc85_dqc0.1_bwd_20d_optimal_3h_glbland_nozbound_intact_RH2_1971-2010_v3.mat'));
        fprintf('*** Loading: %s *** \n', catalog_filename)
        load(catalog_filename, 'AYR_hotspot');   

        ind = find(~ismember({AYR_hotspot.name}, {'Guianas', 'WC. Africa', 'Kalimantan'})); %/ remove those that are not the 3 new ones.
        AYR_hotspot(ind) = [];
        basin_catalog = AYR_hotspot;
        str_from_basin = sprintf('AYR_hs%d', length(basin_catalog));
        str_from_basin_traj = 'glbland';
        
    elseif from_basin == 3  %/ load the 17 most updated hotspots (v3, Cheng & Lu 2022)
        catalog_filename = string(strcat('/disk/r059/tfchengac/FLEXPART/', expmnt, '/prcssd_data_4plotting/',...
                                      'AYR_hs_p95_area_sum_Pm_noOcns_RHc85_dqc0.1_bwd_20d_optimal_3h_glbland_nozbound_intact_RH2_1971-2010_v3.mat'));
        fprintf('*** Loading: %s *** \n', catalog_filename)
        load(catalog_filename, 'AYR_hotspot');   
        basin_catalog = AYR_hotspot;
        str_from_basin = sprintf('AYR_hs%d_v3', length(basin_catalog));
        str_from_basin_traj = 'glbland';
        
    elseif from_basin == 4
        catalog_name = 'TP_basins';
        catalog_filename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/catalog_', catalog_name, '.mat'];
        fprintf('*** Loading: %s *** \n', catalog_filename)
        load(catalog_filename, 'basin_catalog');
        str_from_basin      = sprintf('TP_basins%d', length(basin_catalog));
        str_from_basin_traj = str_from_basin;

    elseif from_basin == 5
        slct_reg = strcat({'TP_'}, {'Grass', 'Semidesert', 'Tundra'});
        catalog_name = strjoin(slct_reg,'_');
        catalog_filename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/catalog_', catalog_name, '.mat'];
        fprintf('*** Loading: %s *** \n', catalog_filename)
        load(catalog_filename, 'basin_catalog');
        str_from_basin      = 'TP_main_LoVeg';
        str_from_basin_traj = 'TP_basins15';   %/ adopt the traj file saved in from_basin == 4 to save time and space!

    elseif from_basin == 6
    %     catalog_name = 'TP_grids_0.5x0.5';
        catalog_name = 'TP_grids_1x1';
        catalog_filename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/catalog_', catalog_name, '.mat'];
        fprintf('*** Loading: %s *** \n', catalog_filename)
        load(catalog_filename, 'basin_catalog');
        str_from_basin      = catalog_name;
        str_from_basin_traj = 'TP_basins15';   %/ Since a grid box may cover areas outside TP boundary, 
                                               %/ we use the traj file saved in from_basin == 4 for consistency and convenience.

    %=========================== Pakistan Flood Case Study =====================%
    elseif from_basin == 7   
        catalog_name = 'Pakistan_box';
        catalog_filename = ['/disk/r059/tfchengac/FLEXPART/domfill_EA_0.5deg_MPI_202207-08_PakCase/prcssd_data_4plotting/catalog_', catalog_name, '.mat'];
        fprintf('*** Loading: %s *** \n', catalog_filename)
        load(catalog_filename, 'basin_catalog');
        str_from_basin      = catalog_name;
        str_from_basin_traj = catalog_name;   

    elseif from_basin == 8
        catalog_name = 'IPCC-PAK';
        catalog_filename = ['/disk/r059/tfchengac/FLEXPART/domfill_EA_0.5deg_MPI_202207-08_PakCase/prcssd_data_4plotting/catalog_', catalog_name, '.mat'];
        fprintf('*** Loading: %s *** \n', catalog_filename)
        load(catalog_filename, 'basin_catalog');
        str_from_basin      = catalog_name;
        str_from_basin_traj = catalog_name;   

    %=========================== Australia Flood Case Study =====================%
    elseif from_basin == 9
        catalog_name = 'Australia_box';
        catalog_filename = ['/disk/r059/tfchengac/FLEXPART/domfill_EA_0.5deg_MPI_202202_AusCase/prcssd_data_4plotting/catalog_', catalog_name, '.mat'];
        fprintf('*** Loading: %s *** \n', catalog_filename)
        load(catalog_filename, 'basin_catalog');
        str_from_basin      = catalog_name;
        str_from_basin_traj = catalog_name;   

    %=========================== Scotland AR Case Study =====================%
    elseif from_basin == 10
        catalog_name = 'Scotland_box';
        catalog_filename = ['/disk/r059/tfchengac/FLEXPART/domfill_EA_0.5deg_MPI_202309-10_ScotCase/prcssd_data_4plotting/catalog_', catalog_name, '.mat'];
        fprintf('*** Loading: %s *** \n', catalog_filename)
        load(catalog_filename, 'basin_catalog');
        str_from_basin      = catalog_name;
        str_from_basin_traj = catalog_name;   
    
    else
        error('Invalid input of ''from_basin == %d''!', from_basin)
    end

end
% if from_basin == -1                                       
%     basin_catalog   = [];
%     str_domain      = 'glb';
% 
% elseif any(ismember(from_basin, [0:3]))
%     AYR_hs_data_filename = strcat(masterfolder, expmnt, '/prcssd_data_4plotting/',...
%                                  'AYR_hs_p95_area_sum_', slct_var, str_noOceans,'_RHc85_dqc0.1_bwd_20d_optimal_3h_glbland_nozbound_intact_RH2_1971-2010', str_figver,'.mat');
%     fprintf('*** The queried data are found. *** \n*** Loading: %s *** \n', AYR_hs_data_filename)
%     load(AYR_hs_data_filename, 'AYR_hotspot');                                            %/ output is AYR_hotspot (based on glbland)
%     basin_catalog = AYR_hotspot;
% 
%     if from_basin == 0
%         str_domain = 'glbland';
%     elseif from_basin == 1 || from_basin == 3
%         str_domain = 'AYR_hs17';
%     elseif from_basin == 2 
%         str_domain = 'AYR_hs3';
%     end
% 
% elseif from_basin == 4
%     slct_reg = {'TP_basins'};
%     catalog_filename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/catalog_', strjoin(slct_reg,'_'), '.mat'];
%     fprintf('*** The queried data are found. *** \n*** Loading: %s *** \n', catalog_filename)
%     load(catalog_filename, 'basin_catalog');
%     str_domain = sprintf('TP_basins%d', length(basin_catalog));
% 
% elseif from_basin == 5
%     slct_reg = strcat({'TP_'}, {'Grass', 'Semidesert', 'Tundra'});
%     catalog_filename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/catalog_', strjoin(slct_reg,'_'), '.mat'];
%     fprintf('*** The queried data are found. *** \n*** Loading: %s *** \n', catalog_filename)
%     load(catalog_filename, 'basin_catalog');
%     str_domain = 'TP_main_LoVeg';
% 
% elseif from_basin == 6
% %     catalog_name = 'TP_grids_0.5x0.5';
%     catalog_name = 'TP_grids_1x1';
%     catalog_filename = ['/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/catalog_', catalog_name, '.mat'];
%     fprintf('*** Loading: %s *** \n', catalog_filename)
%     load(catalog_filename, 'basin_catalog');
%     str_domain      = catalog_name;
% 
% else
%     error('Invalid input of ''from_basin''!')
% end




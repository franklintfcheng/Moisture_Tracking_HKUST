function slct_ens = read_CMIP_ens(varargin)

    pnames = {'project_name', 'model_list', 'var', 'exp', 'ori_field', 'remove_nest_cell'};

    dflts =  cell(1, length(pnames));
    [          project_name,   model_list,    var,  exp,     ~,         remove_nest_cell] ...
                            = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 29 Feb 2024
    %/
    %/ NOTE: If the model has no ensemble for a certain variable, 
    %/       it is not a good practice to input an empty cell for slct_ens{m}.
    %/       Instead, you just shouldn't call this model in the first place.
    %/
    %/=====================================================================
    
    if ischar(model_list)               
        model_list = {model_list};   
    end
    % if ~isequal(ori_field, 'monthly')   error('[read_CMIP_ens]: code not ready for %s CMIP6 data!', ori_field);   end
    
    slct_ens = cell(length(model_list), 1); 
    if isequal(project_name, 'TP')
        %========================= Relevant to DAMIP ======================
        %/ [IMORPOTANT] Omit CESM2 for missing uva for r1i1p1f1 to r3i1p1f1
        %/              Omit E3SM models for missing ua, va
        %/              Omit HadGEM models for their 360-day calendar
        %==================================================================
        for m = 1:length(model_list)
            %/ By default
            if contains(model_list{m}, 'MME')
                slct_ens{m} = {'3EM'};  %/ For TP project, use three ensembles (can be different across models) for MME
            else
                slct_ens{m} = {'r1i1p1f1','r2i1p1f1','r3i1p1f1'}; 
            end

            %/ Handle the exception
            if ismember(model_list{m}, {'BCC-CSM2-MR'})
                if ismember(exp, {'ssp245','ssp585'})
                    slct_ens{m} = {'r1i1p1f1'};
                end

            elseif ismember(model_list{m}, {'CNRM-CM6-1'}) 
                if ismember(var, {'od550aer'})
                    slct_ens{m} = {'r1i1p1f2','r2i1p1f2','r3i1p1f2'};  %/ only these three ensembles available for od550aer
                else
                    slct_ens{m} = {'r1i1p1f2','r3i1p1f2','r4i1p1f2'};  
                end

            elseif ismember(model_list{m}, {'GISS-E2-1-G'}) 
                if ismember(var, {'od550aer'})
                    if ismember(exp, {'historical', 'hist-aer','ssp245','ssp585'})
                        slct_ens{m} = {'r1i1p3f1','r2i1p3f1','r3i1p3f1'};
                    end
                else
                    slct_ens{m} = {'r1i1p1f2','r2i1p1f2','r3i1p1f2'};
                end

            elseif ismember(model_list{m}, {'GFDL-ESM4'})
                if ismember(exp, {'hist-aer','hist-GHG','ssp585'}) || ...
                    (ismember(exp, {'historical'}) && ismember(var, {'evspsbl', 'wap'})) || ...
                    (ismember(exp, {'ssp245'}) && ismember(var, {'wap'})) || ...
                    ismember(exp, {'ssp370'}) || ...
                    (ismember(exp, {'historical', 'hist-aer','ssp245','ssp585'}) && ismember(var, {'od550aer'}))

                    slct_ens{m} = {'r1i1p1f1'};
                end

            elseif ismember(model_list{m}, {'NorESM2-LM'})
                if ismember(exp, {'ssp585'})
                    slct_ens{m} = {'r1i1p1f1'};
                end
            end
        end

    elseif ismember(project_name, {'CMIP6_MJO', 'CMIP6_pr_diagnosis'})
        %==================================================================
        %/ [IMORPOTANT] 
        %/              Omit AWI-CM-1-1-MR for missing rlut for historical
        %/              Omit AWI-ESM-1-1-LR for missing rlut for ssp245, ssp585
        %/              Omit BCC-CSM2-MR for missing rlut for historical
        %/              Omit BCC-ESM1 for missing rlut for ssp245, ssp585
        %/              Omit CESM2-FV for missing uva for r1i1p1f1 to r3i1p1f1
        %/              Omit CMCC-CM2-HR4 for missing rlut for ssp245, ssp585
        %/              Omit all E3SM models for missing ua, va
        %/              Omit EC-Earth3-AerCHEM for missing ssp585
        %/              Omit FGOAL-f3-L for missing ssp245, ssp585
        %/              Omit GISS-E2-1-G for missing wap
        %/              Omit GISS-E2-2-G for missing ssp245, ssp585
        %/              Omit all HadGEM models for their 360-day calendar
        %/              Omit ICON-ESM-LR for missing ssp245, ssp585
        %/              Omit IPSL-CM5A2-INCA, IPSL-CM6A-LR-INCA for missing sspw245, ssp585
        %/              Omit KIOST-ESM for missing ua, va
        %/              Omit MPI-ESM-1-2-HAM for missing ssp245, ssp585
        %/              Omit NorCPM1 for missing ssp245, ssp585
        %/              Omit all NorESM2-LM models for missing rlut
        %/              Omit SAM0-UNICON for missing ssp245, ssp585 
        %==================================================================

        for m = 1:length(model_list)
            %/ By default
            if contains(model_list{m}, 'MME')
                slct_ens{m} = {'1EM'};   %/ For CMIP6_MJO project, use one ensemble for MME
            else
                slct_ens{m} = {'r1i1p1f1'}; 
            end

            %/ Handle the exception
            if ismember(model_list{m}, {'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'CNRM-ESM2-1',...
                                        'MIROC-ES2L', 'UKESM1-0-LL', 'GISS-E2-1-G'}) 
                slct_ens{m} = {'r1i1p1f2'}; 
            elseif ismember(model_list{m}, {'CAMS-CSM1-0'}) 
                slct_ens{m} = {'r2i1p1f1'};  
            elseif ismember(model_list{m}, {'CESM2'}) 
                % slct_ens{m} = {'r10i1p1f1'};  %/ since r1i1p1f1 does not contain ua, va
                slct_ens{m} = {'r4i1p1f1'};  %/ since r10i1p1f1 does not contain historical pr
            elseif ismember(model_list{m}, {'MIROC-ES2H'}) 
                slct_ens{m} = {'r1i1p4f2'};   
            end
        end
    else
        error('code not set!');
    end

    %/ Whether to remove nested cell (will also remove the empty ones)
    if remove_nest_cell   
        slct_ens = slct_ens{:};
        slct_ens(cellfun(@isempty,slct_ens)) = []; 
    end
end
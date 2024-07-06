%%
function [CC_trajs_prop, mem, RMSD] = DorlingCC(varargin)

%/ create a set of valid parameters and their default value
pnames = {'ncluster', 'trajs_prop',   'dim_of_ntraj', 'dim_of_nvar', 'dim_of_ntrajtime',  'str_trajs_prop', 'seed',  'max_iter',  'loss_thres_perc',  'NumWorkers'};  
dflts  = {        [],           [],               [],            [],                 [],               [],    129,       10000,                  1,            []};
[           ncluster,   trajs_prop,     dim_of_ntraj,   dim_of_nvar,   dim_of_ntrajtime,   str_trajs_prop,   seed,    max_iter,    loss_thres_perc,   NumWorkers] = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

%%
%/ Description
%/          'ncluster':  A value or an array (sensitive test).
if isempty(ncluster)           error('ncluster is empty! Check your input!');     end
if any(diff(ncluster) > 0)     error('Input the ncluster in descending order!');  end %/ must start from the largest # of clusters.
if isempty(dim_of_ntraj)       error('Specify ''dim_of_ntraj''!');                end
if isempty(dim_of_nvar)        error('Specify ''dim_of_nvar''!');                 end
if isempty(dim_of_ntrajtime)   error('Specify ''dim_of_ntrajtime''!');            end
if isempty(str_trajs_prop)     error('Specify ''str_trajs_prop''!');              end

%/ (a1) Permute the input 'trajs_prop' into (ntraj, trajtime, nvar)
dim4permute = [dim_of_ntraj, dim_of_ntrajtime, dim_of_nvar];   %/ e.g.,  [1, 3, 2]
trajs_prop  = permute(trajs_prop, dim4permute);
[ntraj, ntrajtime, nvar] = size(trajs_prop);
fprintf('*** Input ''trajs_prop'': ntraj = %d, ntrajtime = %d, nvar = %d ***\n', ntraj, ntrajtime, nvar);
%%
%/ (a2) Convert lon into geodistance/arclen (rhumbline) to a ref point at equator.
% ind_x  = findismember_loop(str_trajs_prop, {'x'}); %/ As distance() takes lat lon.
% ref_pt = [0, 0];  %/ a reference point at the equator
% arclen = nan(ntraj, ntrajtime);
% for t = 1:ntrajtime
%     arclen(:,t) = distance('rh',[zeros(ntraj, 1), squeeze(trajs_prop(:,t,ind_x))],ref_pt); %/ arc length in degree (by default)
% end

ind_yx  = findismember_loop(str_trajs_prop, {'y','x'}); %/ As distance() takes lat lon.
ref_pt = [0, 0];  %/ a reference point at the equator
D = nan(ntraj, ntrajtime);
for t = 1:ntrajtime
    %/ [IMPORTANT]: 'gc' helps to correctly measure the distance of points
    %/                   around the north pole. 
    %/              'rh' will fail.
    D(:,t) = distance('gc',squeeze(trajs_prop(:,t,ind_yx)),ref_pt); %/ arc length in degree (by default)
end

r = 6371;       %/ Earth's raidus (km)
D = D/180*pi*r; %/ convert from degree to distance (km)

%/ replace lon with geodistance [But keeping the lat as one of the vars]
% ind_x  = findismember_loop(str_trajs_prop, {'x'});
% trajs_prop(:,:,ind_x) = D;  

%/ D not replace lon with D, but to append it.
%/ without lon it is hard to distinguish zonally mirrored trajectories.
trajs_prop(:,:,end+1) = D; 
nvar = nvar + 1;

%/ double checking
% delta_D = abs(diff(D', 1));
% Maxdelta_D = max(delta_D, [], 'all')
% 
% [ind_row, ind_col] = find(delta_D == Maxdelta_D)
% delta_D_test = delta_D(:,ind_col(1));
% traj_prop_test = squeeze(trajs_prop(ind_col(1), :, :));
% 
% %/ id 199961 -> a traj passing across the north pole
% ind_Case1 = 199961;
% delta_D_Case1 = delta_D(:, ind_Case1);
% traj_prop_Case1 = squeeze(trajs_prop(ind_Case1, :, :));


%/ (b) Standardization 
mean_array = mean(reshape(trajs_prop, [], nvar), 1, 'omitnan');             %/ (1, nvar)
std_array  = std(reshape(trajs_prop, [], nvar), 0, 1, 'omitnan');           %/ (1, nvar)

% There are boundary problems [unsolved]
% (180 - mean_array(2))/std_array(2)   %/ 202
% (-179 - mean_array(2))/std_array(2) %/ -201.57

for i = 1:nvar
    trajs_prop(:,:,i) = (trajs_prop(:,:,i) - mean_array(i))/std_array(i);
end

RMSD          = nan(length(ncluster),  2);
mem           = cell(length(ncluster), 1);
CC_trajs_prop = cell(length(ncluster), 1);

for c = 1:length(ncluster)
    nclu   = ncluster(c);   
    
    if nclu == ncluster(1)
        %/ (c) Generate Seed Trajs
        rng(seed);                                      %/ rng() has to run right before the randi().
        ind_rand   = randperm(ntraj, nclu)';            %/ use randperm() instead of randi() to obtain unique random number!       
        seed_trajs = trajs_prop(ind_rand, :, :);  
        
    else
        %/ (g) Merge the Pair of Closest Seed Trajs 
        tic
        ntimeOfmerge = ncluster(c-1) - ncluster(c);     %/ since ncluster is in descending order
        for m = 1:ntimeOfmerge                          %/ merge recursively until reaching the desired # of clusters (Agglomerative) 

            D_seed_2D = nan(ncluster(c), ncluster(c));
            for k1 = 1:ncluster(c)
            for k2 = 1:ncluster(c)

                if k2 <= k1   continue;   end           %/ keep only the upper triangle of the matrix (speedup by 2x!)
                D_seed_2D(k1, k2) = EulerDis(seed_trajs(k1,:,:), seed_trajs(k2,:,:));

            end
            end
            
            [ind_k1, ind_k2] = find(D_seed_2D == min(min(D_seed_2D)));

            merged_seed = mean(seed_trajs([ind_k1, ind_k2], :, :), 1);

            seed_trajs(ind_k1,:,:) = merged_seed;       %/ replace the 1st old seed
            seed_trajs(ind_k2,:,:) = [];                %/ remove  the 2nd old seed
            
            fprintf('*** m = %d, nclu = %d ***\n', m, size(seed_trajs,1))
        end
        fprintf('!!! Finished merging (%.1f seconds) !!!\n', toc);
    end
    
    %/ (d) Assign trajs to each seed traj based on the given properties
    trajs_clu_old = nan(ntraj, 1);
    trajs_clu_new = nan(ntraj, 1);
    for iter = 1:max_iter
        tic
        if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
            parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
        end
        parfor n = 1:ntraj
%         for n = 1:ntraj     %/ testing
            %/ boardcast var (avoid overhead, will speedup by ~30%)
            seed_trajs_bc = seed_trajs;  
            trajs_j_bc    = trajs_prop(n,:,:);    %/ 1 x ntrajtime x nvar
            
            %/ compute the distance of a traj to all seed trajs_prop
            D = nan(nclu, 1);
            for k = 1:nclu
                D(k) = EulerDis(trajs_j_bc, seed_trajs_bc(k,:,:)); 
            end

            [~, I] = min(D);
            trajs_clu_new(n) = I;
        end
        loss = length(find(trajs_clu_new - trajs_clu_old ~= 0))/ntraj*100;  %/ in %
        fprintf('*** nclu = %d, iter = %d/%d, loss (mismatch rate) = %.4f%% (%.1f seconds) ***\n', nclu, iter, max_iter, loss, toc);
        
        if iter ~= 1 && loss < loss_thres_perc   %/ it must go with >1 iteration to see if clustering (largely) converges.
            fprintf('!!! Clustering stops at loss = %.4f%% (< %.4f%%) !!!\n', loss, loss_thres_perc);
            break;
        elseif iter == max_iter
            warning('!!! Maximum iteration has reached (%d/%d). The membership has not yet converged !!!\n', iter, max_iter);
        else
            trajs_clu_old = trajs_clu_new;
        end

        %/ (e) Recalculate the seed trajs_prop by averaging trajs_prop with the new clustering.
        for k = 1:nclu                                      %/ use for-loop since trajs_prop is too large for parfor to run efficiently.
            ind = find(trajs_clu_new == k);
            seed_trajs(k,:,:) = mean(trajs_prop(ind, :, :), 1);
        end
    end
    
    if iter == max_iter && loss > loss_thres_perc*100
       error('!!! Iteration reaches the max. Results still do not converge. Terminating the program... !!!'); 
    end
    
    %/ (f) Compute total RMSD for the given # of clusters
    MS = nan(nclu, 1);
    for k = 1:nclu
        ind = find(trajs_clu_new == k);

        square_deviation = squeeze((trajs_prop(ind,:,:) - seed_trajs(k,:,:)).^2);
        MS(k) = sum(sum(sum(square_deviation, 'omitnan'), 'omitnan'), 'omitnan');
    end

    RMSD(c,1) = nclu;
    RMSD(c,2) = sqrt(sum(MS)/ntraj/nvar);
    fprintf('*** For nclu = %d, RMSD = %.4f ***\n', RMSD(c,1), RMSD(c,2));

    mem{c} = trajs_clu_new;   %/ store the converged membership for the given cluster #.
    
    %/ Restore the true values (inverse standardization), create a new
    %/ variable, since we need the standardized seed_trajs for the next loop of nclu.
    seed_trajs_restore = nan(size(seed_trajs));
    for i = 1:nvar
        seed_trajs_restore(:,:,i) = seed_trajs(:,:,i)*std_array(i) + mean_array(i);
    end

    CC_trajs_prop{c} = seed_trajs_restore; %/ store the converged seed_trajs for the given cluster #.
end

%/ Remove cell format if ncluster is a value.
if length(ncluster) == 1
    CC_trajs_prop = CC_trajs_prop{:};
    mem = mem{:};
end

end
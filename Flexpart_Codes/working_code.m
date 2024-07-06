



%/ IF perform interpolation of data to the level of sp
            if isequal(interp_mode, 'linear')
                ind_k = k-1:k;                                  %/ For ECMWF data, p-level is always from small (TOA) to large (sfc).
                x   = level_bc_old_reshape_ij(ind_k);           %/ The old level
                xq  = level_bc_new_reshape_ij(ind_k);           %/ The new level
                v   = data_bc_reshape_ij(ind_k,t);
                data_bc_reshape_ij(ind_k,t) = interp1(x,v,xq);  %/ Using 1D interp and update 
            end
            %%
            
            
            level_bc_old_reshape_ij = dataset.q.level;
            level_bc_new_reshape_ij = dataset.q.level;
            k = 25; level_bc_new_reshape_ij(k) = 930;       %/ case 1
%             k = 27; level_bc_new_reshape_ij(k) = 1015;    %/ case 2
            ind_k = k-1:k;                                  %/ For ECMWF data, p-level is always from small (TOA) to large (sfc).
            
            x   = level_bc_old_reshape_ij(ind_k);           %/ The old level
            xq  = level_bc_new_reshape_ij(ind_k);           %/ The new level
            v = squeeze(dataset.q.daily(1,1,ind_k,1))*1000;
            vq = interp1(x,v,xq, 'linear', 'extrap');       %/ Using 1D interp (with extrapolation!) and update 
            
            
            
            
            
            %%
function f_reconst = my_fourier(varargin)
    
    %/ create a set of valid parameters and their default value
    pnames = {'f',     'n_list'};
    dflts  = cell(length(pnames), 1);
               [f,      n_list] = ...
                    internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
    %%
    %/============================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Feb 28, 2024
    %/
    %/ OUTPUT: A reconstructed matrix with the same size as f, based on the
    %/           Fourier Harmonics indicated in n_list.
    %/
    %/============================================================================
    
    fprintf('*** Running my_fourier... ***\n')

    if ~isnumeric(n_list)
        error('n_list should only be a numeric array!')
    end

    %/ Get the time length (Always assume the last dim is time dim!)
    sz = size(f);
    time_dim = length(sz);
    T  = sz(end);
    if T == 1   
        error('T is 1. Check the code!');    
    end
    
    t              = 0:1:T-1;
    singleton_dims = [ones(1,length(sz)-1),T];
    t              = reshape(t, singleton_dims); %/ If input is 2D, size(t) = [1, T]; If input is 3D, size(t) = [1, 1, T]; etc.

    %/ NOTE: n == how many cycles per time period T, if n > # of time grid,
    %/       it becomes impossible to resolve n cycles in the period
    %/       the results will become strange --> need to avoid it.
    if max(n_list) > T-1
        error('Input n_list should not exceed the total # of time grids!!');
    end
    
    %/=====================================================================
    %/ NOTE: Do NOT use parfor here (whether Threads or Processes), because
    %/       they will create multiple copies of data that potentially blow
    %/       up the loop. More importantly, a simple for-loop is the fatest
    %/       in this case since matrix operation is inherently multithreaded.
    %/=====================================================================
    dt = 1; f_reconst = 0;
    for n = n_list
        %/=====================================================================
        %/ Calculate the Fourier Series with the given list of N harmonics
        %/ Show the reconstructed signal.                                 
        %/ NOTE: include const term (n=0) to restore the magnitude.      
        %/=====================================================================
        if n == 0
            a0 = 2/T*sum(f.*cos(n*2*pi.*t/T)*dt, time_dim); %/ integrate numerically. 
            f_reconst = f_reconst + repmat(a0/2, singleton_dims); %/ make it an const array (for visualization)
        else
            an = 2/T*sum(f.*cos(n*2*pi.*t/T)*dt, time_dim); %/ integrate numerically. E.g. (360 x 91).*(1 x 1 x 73) = (360 x 91 x 73) by elementwise multiplication.
            bn = 2/T*sum(f.*sin(n*2*pi.*t/T)*dt, time_dim); %/ integrate numerically

            f_reconst = f_reconst + an.*cos(n*2*pi.*t/T) + bn.*sin(n*2*pi.*t/T);
        end
    end
end
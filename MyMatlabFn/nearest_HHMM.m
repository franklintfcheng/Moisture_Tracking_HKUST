%%
function [N_HHMM, ind_N_HHMM] = nearest_HHMM(date_HHMM, diurnal_time)

    %/ this will take the last 4 digits as HHMM.
    date_HHMM_db = mod(date_HHMM, 1e4);  
    
    %/ double() is needed. Otherwise it will not output decimals!!
    date_HHMM_db = double(date_HHMM_db);
    date_HHMM_db(date_HHMM < 0) = nan;   %/ set any -ve dates (ie., missing values) into nan.
    
    %/ append 2400 to diurnal_time_db to correctly assign 23:59 to 0:00 (for example).
    diurnal_time_db = double(diurnal_time);
    diurnal_time_db(end+1) = 2400;
    
    date_HHMM_db_1D = reshape(date_HHMM_db, [], 1);  %/ reshape to 1D.
    
    K           = length(date_HHMM_db_1D);
    N_HHMM      = nan(size(date_HHMM_db_1D));
    ind_N_HHMM  = nan(size(date_HHMM_db_1D));
    
    for k = 1:K
        date_HHMM_k = date_HHMM_db_1D(k);

        if isnan(date_HHMM_k)    continue;    end  %/ skip if the date is nan.
        
        date_HHMM_decimal    = floor(date_HHMM_k/1e2)     + mod(date_HHMM_k, 1e2)./60;     
        diurnal_time_decimal = floor(diurnal_time_db/1e2) + mod(diurnal_time_db, 1e2)./60;

        %/ WARNING: min() will automatically output the value at the smallest
        %           index position if multiple minima exists.
        D = abs(date_HHMM_decimal - diurnal_time_decimal);
        I = find(D == min(D)); 
        
        %/ if at 12:45 between 12:30 and 13:00, then choose 13:00 --> [T-15m, T+15m)
        if length(I) == 2
            I = I(2);           
        elseif length(I) > 2
            error('something wrong with I.')
        end
        
        N_HHMM(k)     = diurnal_time_db(I);
        ind_N_HHMM(k) = I;  
    end
    N_HHMM(N_HHMM == 2400) = 0;                             %/ restore 24:00 back to 0:00.
    ind_N_HHMM(ind_N_HHMM == length(diurnal_time_db)) = 1;  %/ restore the last index position to the 1st.
    
%     N_HHMM = int64(N_HHMM); %/ make use it's an integer type. %/ do NOT set it to int, since it will auto convert nan into 0.
    
    %/ reshape back to the original dimension.
    N_HHMM     = reshape(N_HHMM,     size(date_HHMM_db, 1), size(date_HHMM_db, 2));
    ind_N_HHMM = reshape(ind_N_HHMM, size(date_HHMM_db, 1), size(date_HHMM_db, 2));
    
end




function idx = findismember_loop(A, B, row_mode) 
%%
    %/=====================================================================
    %/      Author: Franklin Cheng
    %/ Last update: 10 Aug 2022
    %/
    %/ Description: Find the indices of B in A (without auto-sorting and
    %/              with duplicates)
    %/
    %/=====================================================================
    
    if nargin == 2
        if ischar(A)    
            A = {A};  
        end
        if ischar(B)    
            B = {B};  
        end
    
        if length(unique(B)) == length(B) %/ No duplicated elements in B -> use the *much faster* routine here.
            [tf,loc] = ismember(A, B);
            [~,p] = sort(loc(tf));
            idx = find(tf);
            idx = idx(p);
        else                              %/ Otherwise, use the dump for-loop.
            idx = [];
            for i = 1:length(B)
                ind = find(ismember(A, B(i)));
                
                if ~isempty(ind)
                    if length(ind) > 1
                        warning('duplicated items is found in array A! Returning the first index by default.');
                    end
                    idx = cat(1, idx, ind(1));
                end
            end
        end
        
    elseif isequal(row_mode, 'rows')   %/ Assume it is to match rows. But this is SLOW!!!
        idx = nan(size(B, 1), 1);
        for i = 1:size(B, 1)
            idx(i) = find(ismember(A, B(i,:), 'rows'));
        end
        
%         idx = [];
%         for i = 1:size(B, 1)
%             ind = find(ismember(A, B(i,:), 'rows'));
% 
%             if ~isempty(ind)
%                 idx = cat(1, idx, ind);
%             end
%         end
        
    else
        error('Wrong input argument!!')
    end
    
%     idx = nan(1, length(B));
%     for i = 1:length(B)
%         ind = find(ismember(A, B(i)));
%         
%         if ~isempty(ind)
%             if length(ind) ~= 1
%                 if length(B) == 1 %/ then we output the indices of the redundant value
%                     idx = ind;
%                 else
%                     error('This function could not handle redundant values in both A and B!')
%                 end
%             else
%                 idx(i) = ind;
%             end
%         end
%     end
%     idx(isnan(idx)) = [];

end

function [data_trans, lon_2D, lat_2D] = conv_to_circular_data(data, lon, lat, map_lon_upper)

    if ~isvector(lon) || ~isvector(lat)
        error('lon, lat must be vectors!');
    end
    
    if find(lon == min(abs(lon))) ~= 1 && ~isempty(find(lon < 0, 1))    
        %/ If the smallest lon value (e.g., 0) is not the first element of
        %/ lon array, and the lon array contains -ve values,
        %/ e.g., -179.75:0.5:179.75,
        %/ then no need to translate the indices in lon and data.
        [lon_2D, lat_2D] = meshgrid(lon, lat);
        data_trans = data;
    else
        if nargin ~= 4
            map_lon_upper = 180;  %/ default
        end
        % map_lon_upper = 180
        ind = find(lon > map_lon_upper);
        lon(ind) = lon(ind) - 360;                %/ change to [0:180,-179:0] for a correct contf plot (if map_lon_upper == 180)

        [lon_2D, lat_2D] = meshgrid(lon, lat);  

        ind_translate = [ind', 1:ind(1)-1];       %/ put the left lon dim to the right. --> [-179:0, 0:180]
        lon_2D  = lon_2D(:, ind_translate);     
        lat_2D  = lat_2D(:, ind_translate);      
        
        if size(data, 1) == length(lon)   
            data = data'; 
        end
        data_trans = data(:, ind_translate); 
    end
%     disp(lon)
end

% function [data_trans, lon_2D, lat_2D] = conv_to_circular_data(data, lon, lat)
% %%
%     if find(lon == min(abs(lon))) ~= 1 && ~isempty(find(lon < 0))    %/ e.g., -179.75:0.5:179.75 (without zero)      
%         [lon_2D, lat_2D] = meshgrid(lon, lat);
%         
%     else
%         ind = find(lon > 180);
%         lon(ind) = lon(ind) - 360;                 %/ change to [0:180,-179:0] for a correct contf plot.
% 
%         [lon_2D, lat_2D] = meshgrid(lon, lat);  
% 
%         ind_translate = [ind', 1:ind(1)-1];       %/ put the left lon dim to the right. --> [-179:0, 0:180]
%         lon_2D  = lon_2D(:, ind_translate);     
%         lat_2D  = lat_2D(:, ind_translate);      
%         
%         if size(data, 1) == length(lon)   data = data'; end
%         data_trans = data(:, ind_translate); 
%     end           
% %     disp(lon)
% end
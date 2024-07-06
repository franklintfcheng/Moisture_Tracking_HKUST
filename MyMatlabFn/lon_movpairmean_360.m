%%
function mean_x = lon_movpairmean_360(x, dx_thres)

    %/ x must be a column vector. If not, then transpose it.

    if size(x,2) ~= 1 
        x = x';
    end
    
    x_pairwise = cat(2, x(1:end-1), x(2:end));

    for xi = 1:size(x_pairwise,1) %/ dont use length(), can be problematic when there is only one pair: size = (1,2)
        
        if abs(diff(x_pairwise(xi, :))) > dx_thres          %/ dx_thres can set to 300, which means an unusual value which indicates the traj is crossing the lon boundary.
            mean_x(xi,1) = (sum(x_pairwise(xi,:)) + 360)/2; %/ [x(t) + x(t+1) + 360] / 2
        else
            mean_x(xi,1) = mean(x_pairwise(xi,:));
        end
        
    end
    mean_x(mean_x >= 360) = mean_x(mean_x >= 360) - 360; %/ reset to the correct lon range [0, 360)

end
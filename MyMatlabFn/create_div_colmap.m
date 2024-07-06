%%
function div_colmap = create_div_colmap(nve_color, pve_color, cont_levels)
    
    %/ A quick creation of colmap for diverging colmap

    if any(ismember(cont_levels, 0))  %/ if cont_levels contains zero-level
    
        div_colmap = [repmat(nve_color, (length(cont_levels)-1)/2, 1);
                       [0 0 0];
                       repmat(pve_color, (length(cont_levels)-1)/2, 1);];
    else
        div_colmap = [repmat(nve_color, length(cont_levels)/2, 1);
                      repmat(pve_color, length(cont_levels)/2, 1);];
    end

end
%%
function bndry_data = add_pts_to_bndry(varargin)

    pnames = {'bndry_data', 'suppl_Dx', 'suppl_Dy'};
    dflts  = {          [],        2,       2};
    [           bndry_data,   suppl_Dx,   suppl_Dy] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %%
    %/      Author: Franklin Cheng
    %/ Last Update: Jul 4, 2023

    %/ Description: This function is designed to supplement points in a
    %/              long straight line on the map. This will allow for curvature when
    %/              plotting in spherical ('ortho') projection.

    %/ Search for any straight line along a latitude
    ind              = find(abs(diff(bndry_data(:,2))) < 1e-10);
    added_pts_cell = cell(length(ind), 1);
    for ii = 1:length(ind)
        if bndry_data(ind(ii),1) > bndry_data(ind(ii)+1,1)
            Dx = -1 * suppl_Dx;
        else
            Dx = suppl_Dx;
        end
        added_pts_lon      = [bndry_data(ind(ii),1):Dx:bndry_data(ind(ii)+1,1)]';
        added_pts_lonlat   = [added_pts_lon, repmat(bndry_data(ind(ii),2), length(added_pts_lon), 1)];
        added_pts_cell{ii} = added_pts_lonlat(2:end,:);  %/ the first index is needless (avoid redundant points)
        if ii == length(ind)
            shift = 0;
            for jj = 1:length(ind)
                bndry_data = [bndry_data(1:ind(jj)+shift,:); added_pts_cell{jj}; bndry_data(ind(jj)+shift+1:end,:)];
                shift = shift + size(added_pts_cell{jj}, 1);
            end
        end
    end

    %/ Likewise, search for any straight line along a longitude
    ind              = find(abs(diff(bndry_data(:,1))) < 1e-10);  %/ to handle round-off randomness
    added_pts_cell   = cell(length(ind), 1);
    for ii = 1:length(ind)
        if bndry_data(ind(ii),2) > bndry_data(ind(ii)+1,2)
            Dy = -1 * suppl_Dy;
        else
            Dy = suppl_Dy;
        end
        added_pts_lat      = [bndry_data(ind(ii),2):Dy:bndry_data(ind(ii)+1,2)]';
        added_pts_lonlat   = [repmat(bndry_data(ind(ii),1), length(added_pts_lat), 1), added_pts_lat];
        added_pts_cell{ii} = added_pts_lonlat(2:end,:);  %/ the first index is needless (avoid redundant points)
        if ii == length(ind)
            shift = 0;
            for jj = 1:length(ind)
                bndry_data = [bndry_data(1:ind(jj)+shift,:); added_pts_cell{jj}; bndry_data(ind(jj)+shift+1:end,:)];
                shift = shift + size(added_pts_cell{jj}, 1);
            end
        end
    end
end




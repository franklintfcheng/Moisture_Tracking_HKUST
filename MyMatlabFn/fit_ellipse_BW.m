%%
function [Centr, Orient, MajLen, MinLen, Eccen,  bndry_data] = fit_ellipse_BW(varargin)

    %/ Author: Fandy
    %/   Date: 7 Dec 2021
    
    %/ *** Some original codes from Akira Agata ***
    %/ https://www.mathworks.com/matlabcentral/answers/495720-how-to-fit-an-ellipse-to-an-image-in-matlab

    %/ NOTE: it sets binary and fit it with an ellipse, and output E with
    %/       important attributes of the object.

    %/ if data_2D = an RGB image, then we have to executes these codes to binarize it
    %     Igray = rgb2gray(data_2D);
    %     BW = imbinarize(Igray);

    pnames = {'field_data', 'x', 'y', 'lonlat_or_not', 'draw_ellipse',  'largest_n', 'create_fig'};
    dflts  = {          [],  [],  [],               0,              0,            1,            0};
    [          field_data,   x,   y,   lonlat_or_not,    draw_ellipse,    largest_n,   create_fig ] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    %-------
    
    %%
    if length(size(field_data)) ~= 2    error('The function ''fit_ellipse_BW'' only handles 2D data!');   end
    [nx, ny] = size(field_data);  %/ assume lon lat time, or just lon lat.
    
    if ~isempty(x) && ~isempty(y)
        if nx ~= length(x) || ny ~= length(y)
            error('Make sure field_data''s dimensions fit x and y!');
        end
    end
    
    if ~islogical(field_data)
        BW = ~isnan(field_data);  %/ then convert into logical based on NaN = 0, non-NaN = 1.
    else
        BW = field_data;
    end
    
    %/ BE CAREFUL: If input field is geofield,  need to transpose it from (lon,lat) to (lat,lon) for a correct fitting.
    flag_to_transpose = 0;
    if lonlat_or_not == 1
        flag_to_transpose = 1;
        BW = BW';
    end

    Centr  = nan(largest_n, 2);
    Orient = nan(largest_n, 1);
    MajLen = nan(largest_n, 1);
    MinLen = nan(largest_n, 1);
    Eccen  = nan(largest_n, 1);

    for n = 1:largest_n
        BW_nth = bwareafilt(BW, n);      %/ IMPORTANT: n = 1 (keep the largest object), n = 2 (keep the two largest objects)
        
        if n > 1
            BW_n_minus_1 = bwareafilt(BW, n-1);
            BW_nth       = BW_nth & ~BW_n_minus_1;     % Logical Operation to get the n-th largest object
%             figure
%             imshow(BW_nth);
        end
        
        % Calculate centroid, orientation and major/minor axis length of the ellipse
        E = regionprops(BW_nth, {'Centroid','Orientation','MajorAxisLength','MinorAxisLength', 'Eccentricity', 'BoundingBox'});
        
        if isempty(E)
            warning('No ellipse is fitted with the %d-th largest object from the given binary_array. Output nans.', n); 
            continue;
        else
            Centr(n,:)    = E.Centroid;  %/ x, y
            Orient(n)     = E.Orientation;
            MajLen(n)     = E.MajorAxisLength;
            MinLen(n)     = E.MinorAxisLength;
            Eccen(n)      = E.Eccentricity;
        end
    end
    
    %/ remove nan elements
    ind_nan = find(isnan(Eccen));
    Centr(ind_nan,:) = [];
    Orient(ind_nan)  = [];
    MajLen(ind_nan)  = [];
    MinLen(ind_nan)  = [];
    Eccen(ind_nan)   = [];
    
    if flag_to_transpose
        %/ In this case, BW_2D has been transposed. When drawing the ellipse on geomap,
        %/ the orientation angle is in fact the angle betn minor axis and
        %/ vertical axis. Multiply by -1 to correct it.
        Orient = Orient * -1;    
    end
    phi = Orient;

    bndry_data = cell(length(Eccen), 1);
    for n = 1:length(Eccen)
%         if any(isnan(Centr(n,:)))   continue;    end   %/ skip if nan.
        
        %/The properties are based on lon lat if given, otherwise based on the pixel size (NOTE: it starts from 1 !!!).
        if ~isempty(x) && ~isempty(y)
            dlon = abs(x(2) - x(1));
            dlat = abs(y(2) - y(1));

            %/ NOTE: regionprops define centroids by pixels. 

            %/ Check https://www.mathworks.com/help/images/ref/regionprops.html
            Centr(n,1) = x(floor(Centr(n,1))) + dlon*mod(Centr(n,1), 1);   
            Centr(n,2) = y(floor(Centr(n,2))) + dlat*mod(Centr(n,2), 1);   

            MajLen(n) = sqrt((MajLen(n).*cosd(Orient(n)) * dlon).^2 + (MajLen(n).*sind(Orient(n)) * dlat).^2 );
            MinLen(n) = sqrt((MinLen(n).*sind(Orient(n)) * dlon).^2 + (MinLen(n).*cosd(Orient(n)) * dlat).^2 );
        end

        theta         = linspace(0,2*pi);
        col           = (MajLen(n)/2)*cos(theta);
        row           = (MinLen(n)/2)*sin(theta);
        M             = makehgtform('translate',[Centr(n,:), 0],'zrotate',deg2rad(phi(n)));
        D             = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];    
        bndry_data{n} = D(1:2,:)';
        
        if draw_ellipse
            fontsize = 35;
            %/ Overlaying the fitted ellipse
            if create_fig
                figure
                set(gcf, 'color','w');
                set(gcf, 'Position', get(0, 'Screensize'));
                res_x = unique(round(diff(x), 3));  %/ round up to 3 decimal points to avoid numerical errors
                res_y = unique(round(diff(y), 3));

                if flag_to_transpose
                    h = pcolor(x-res_x/2, y-res_y/2, field_data');     %/ shift by half grid to correct the bug in pcolor.
                else
                    h = pcolor(x-res_x/2, y-res_y/2, field_data);     %/ shift by half grid to correct the bug in pcolor.
                end
                set(h, 'edgecolor','none');
                shading flat

                min_val      = min(field_data, [], 'all', 'omitnan');
                max_val      = max(field_data, [], 'all', 'omitnan');
                contf_levels = linspace(min_val, max_val, 21);
                colmap       = brewermap(length(contf_levels)-1, '*RdYlBu');
                caxis([min(contf_levels) max(contf_levels)]);
                colormap(gca, colmap)
                box on
                set(gca, 'linewidth', 2, 'fontsize', fontsize)
                axis equal
            end

            linewidth = 2;
            hold on;
            plot(bndry_data{n}(:,1), bndry_data{n}(:,2),'r','LineWidth',linewidth)

            xMinor = Centr(n,1) + [ -1,  1 ] * MinLen(n)/2*sind(phi(n));  %/ X-axis of two vertices of Minor axis
            yMinor = Centr(n,2) + [  1, -1 ] * MinLen(n)/2*cosd(phi(n));  %/ Y-axis of two vertices of Minor axis
            line(xMinor,yMinor, 'linewidth', linewidth, 'Color', 'g');

            xMajor = Centr(n,1) + [ -1, 1 ] * MajLen(n)/2*cosd(phi(n));  %/ X-axis of two vertices of Major axis
            yMajor = Centr(n,2) + [ -1, 1 ] * MajLen(n)/2*sind(phi(n));  %/ Y-axis of two vertices of Major axis
            line(xMajor,yMajor, 'linewidth', linewidth, 'Color', 'b');

            plot(Centr(n,1), Centr(n,2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'm');

            %/ horizontal/vertical lines through the centroid
            xline(Centr(n,1), 'linewidth', linewidth*0.5, 'Color', 'k', 'linestyle', '--');
            yline(Centr(n,2), 'linewidth', linewidth*0.5, 'Color', 'k', 'linestyle', '--');
            hold on;

            txt_ellipse = {sprintf('e=%.2f', Eccen(n)), strcat('\phi=', sprintf('%.1f', Orient(n)), '\circ')};
            text(3, -3, txt_ellipse, 'fontsize', fontsize*0.8, 'BackgroundColor', 'w', 'Color', 'k', 'EdgeColor', 'k');
            drawnow; pause(0.05);

            set(gca, 'XAxisLocation', 'bottom')
            set(gca, 'YAxisLocation', 'left')
            set(gca, 'YDir', 'normal')       %/ since imshow will auto reverse YDir!!
            set(gca, 'XDir', 'normal') 
        end
    end
    

    
    
    
    
    
end
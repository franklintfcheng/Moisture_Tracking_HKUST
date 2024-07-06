%%
function [Centr, Orient, MajLen, MinLen, Eccen, bndry_data, Ellip_Area, Grid_Area, MajLen_geo, output_BW] = fit_ellipse_BW_geo(varargin)

    pnames = {'field_data', 'lon', 'lat', 'seek_N', 'largest_n', 'draw_ellipse',  'create_fig'};
    dflts  = {          [],    [],    [],       20,           1,              0,            0};
    [           field_data,   lon,   lat,   seek_N,   largest_n,   draw_ellipse,   create_fig ] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %==================================================================================================
    %/       Author: Franklin Cheng (fandycheng@ust.hk)
    %/  Last Update: October 10, 2023
    %/
    %/     'seek_N': to find N black-and-white patches
    %/  'largest_n': to select the n largest patches from N
    %/
    %/      *** Some original codes from Akira Agata ***
    %/      https://www.mathworks.com/matlabcentral/answers/495720-how-to-fit-an-ellipse-to-an-image-in-matlab
    %/
    %/      NOTE: it sets binary and fit it with an ellipse, and output E with
    %/              important attributes of the object.
    %/
    %/            if data_2D = an RGB image, then we have to executes these codes to binarize it
    %/              Igray = rgb2gray(data_2D);
    %/              BW = imbinarize(Igray);
    %==================================================================================================
    
%     if seek_N < largest_n
%         warning('The input largest_n is larger than seek_N. Now setting seek_N = largest_n.')
%     end
    
    if length(size(field_data)) ~= 2    error('The function ''fit_ellipse_BW'' only handles 2D data!');   end
    [nx, ny] = size(field_data);  %/ assume lon lat time, or just lon lat.
    
    if ~isempty(lon) && ~isempty(lat)
        if nx ~= length(lon) || ny ~= length(lat)
            error('Make sure field_data''s dimensions fit lon and lat!');
        end
    end
    
    if ~islogical(field_data)
        BW = ~isnan(field_data);  %/ then convert into logical based on NaN = 0, non-NaN = 1.
    else
        BW = field_data;
    end
    
    %/ BE CAREFUL: If input field is geofield,  need to transpose it from (lon,lat) to (lat,lon) for a correct fitting.
    BW = BW';
    
    Centr      = nan(seek_N, 2);
    Orient     = nan(seek_N, 1);
    MajLen     = nan(seek_N, 1);
    MinLen     = nan(seek_N, 1);
    Eccen      = nan(seek_N, 1);
    Ellip_Area = nan(seek_N, 1);
    Grid_Area  = nan(seek_N, 1);
    bndry_data = cell(seek_N, 1);
    MajLen_geo = nan(seek_N, 1);
    output_BW  = nan(nx,ny,seek_N);
    xMinor     = nan(seek_N, 2);
    yMinor     = nan(seek_N, 2);
    xMajor     = nan(seek_N, 2);
    yMajor     = nan(seek_N, 2);
    
    BW_nth_tot_area = zeros(seek_N,1);
    grid_area_2D    = calc_grid_area('lon', lon, 'lat', lat)*1e-6; %/ m^2 -> km^2
    BW_nth          = false(size(BW,1),size(BW,2), seek_N);
    
    for n = 1:seek_N
        BW_nth(:,:,n) = bwareafilt(BW, n);      %/ n = 1 (keep the largest), n = 2 (keep the two largest)

        if n > 1
            BW_n_minus_1  = bwareafilt(BW, n-1);
            BW_nth(:,:,n) = BW_nth(:,:,n) & ~BW_n_minus_1;     % Logical Operation to get the n-th largest object
        end
        Ngrid = length(find(BW_nth(:,:,n) == 1));
%         fprintf('*** n = %d, no. of grid cells = %d ***\n', n, Ngrid)
%         figure
%         imshow(BW_nth(:,:,n));
        if Ngrid == 0       break;          end
        BW_nth_tot_area(n) = sum(grid_area_2D(BW_nth(:,:,n)), 'all');

        %/ Calculate centroid, orientation and major/minor axis length of the ellipse
        X = BW_nth(:,:,n);
        E = regionprops(X, {'Centroid','Orientation','MajorAxisLength','MinorAxisLength', 'Eccentricity', 'BoundingBox'});

        %/ Derive the boundary data (and other attributes) from regionprops.
        E = get_geo_attrs_regionprops(lon, lat, E);

        Ellip_Area(n) = E.Area;
        Grid_Area(n)  = BW_nth_tot_area(n);
        Centr(n,:)    = E.Centroid;         %/ lon, lat
        Orient(n)     = E.Orientation;
        MajLen(n)     = E.MajorAxisLength;
        MinLen(n)     = E.MinorAxisLength;
        Eccen(n)      = E.Eccentricity;
        bndry_data{n} = E.bndry_data;
        MajLen_geo(n) = E.MajLen_geo;
        xMinor(n,:)   = E.xMinor;
        yMinor(n,:)   = E.yMinor;
        xMajor(n,:)   = E.xMajor;
        yMajor(n,:)   = E.yMajor;                               
        output_BW(:,:,n) = BW_nth(:,:,n)';   %/ Store the BW logical matrix for the fitted ellpses. 
                                             %/ Tranpose from (lat, lon) to (lon, lat).                            
    end

    output_BW = flip(output_BW, 2);          %/ NOTE: Restore the descending lat dim.
    
    %/ Sort by elliptical area
    a = Ellip_Area;
    a(isnan(a)) = 0;  %/set nan to zeros
    [~, I] = sort(a, 1, 'descend');
    if length(a) > largest_n 
        I = I(1:largest_n);
    end
    
    %/ Sorting & remove nan;
    Ellip_Area = Ellip_Area(I);
    Grid_Area  = Grid_Area(I);
    Centr      = Centr(I,:);
    Orient     = Orient(I);
    MajLen     = MajLen(I);
    MinLen     = MinLen(I);
    Eccen      = Eccen(I);
    bndry_data = bndry_data(I);
    MajLen_geo = MajLen_geo(I);
    xMinor     = xMinor(I);
    yMinor     = yMinor(I);
    xMajor     = xMajor(I);
    yMajor     = yMajor(I);
    output_BW  = output_BW(:,:,I);

    %/ remove nan elements
    ind_nan                = find(isnan(Eccen));
    Ellip_Area(ind_nan)    = [];
    Grid_Area(ind_nan)     = [];
    Centr(ind_nan,:)       = [];
    Orient(ind_nan)        = [];
    MajLen(ind_nan)        = [];
    MinLen(ind_nan)        = [];
    Eccen(ind_nan)         = [];
    bndry_data(ind_nan)    = [];
    MajLen_geo(ind_nan)    = [];
    xMinor(ind_nan)        = [];
    yMinor(ind_nan)        = [];
    xMajor(ind_nan)        = [];
    yMajor(ind_nan)        = [];
    output_BW(:,:,ind_nan) = [];
    
    %/ Draw the ellipses one by one.
    for n = 1:length(Eccen)
        if draw_ellipse
            fontsize = 35;
            %/ Overlaying the fitted ellipse
            if create_fig
                figure
                set(gcf, 'color','w');
                set(gcf, 'Position', get(0, 'Screensize'));
                res_x = unique(round(diff(lon), 3));  %/ round up to 3 decimal points to avoid numerical errors
                res_y = unique(round(diff(lat), 3));

%                 if flag_to_transpose
                h = pcolor(lon-res_x/2, lat-res_y/2, field_data');     %/ shift by half grid to correct the bug in pcolor.
%                 else
%                 h = pcolor(lon-res_x/2, lat-res_y/2, field_data);     %/ shift by half grid to correct the bug in pcolor.
%                 end
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

            line(xMinor(n),yMinor(n), 'linewidth', linewidth, 'Color', 'g');
            line(xMajor(n),yMajor(n), 'linewidth', linewidth, 'Color', 'b');

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

%============================ Unused Codes ===============================%
        
%     if inclusive  %/ WARNING: code not ready.
%         for n = 1:seek_N
%             BW_nth(:,:,n) = bwareafilt(BW, n);      %/ n = 1 (keep the largest), n = 2 (keep the two largest)
% 
%             if n > 1
%                 BW_n_minus_1  = bwareafilt(BW, n-1);
%                 BW_nth(:,:,n) = BW_nth(:,:,n) & ~BW_n_minus_1;     % Logical Operation to get the n-th largest object
%             end
%             Ngrid = length(find(BW_nth(:,:,n) == 1));
%     %         fprintf('*** n = %d, no. of grid cells = %d ***\n', n, Ngrid)
%     %         figure
%     %         imshow(BW_nth(:,:,n));
%             if Ngrid == 0       break;          end
%             
%             if n <= largest_n %/ Check if it already exceeds largest_n.
%                 X = BW_nth(:,:,n);
%         %         figure
%         %         imshow(X);
% 
%                 % Calculate centroid, orientation and major/minor axis length of the ellipse
%                 E = regionprops(X, {'Centroid','Orientation','MajorAxisLength','MinorAxisLength', 'Eccentricity', 'BoundingBox'});
% 
%                 if isempty(E)
%                     warning('No ellipse is fitted with the %d-th largest object from the given binary_array. Output nans.', n); 
%                     continue;
%                 else
% %                     Area(n)       = BW_nth_tot_area(I(n));
%                     Centr(n,:)    = E.Centroid;  %/ lon, lat
%                     Orient(n)     = E.Orientation;
%                     MajLen(n)     = E.MajorAxisLength;
%                     MinLen(n)     = E.MinorAxisLength;
%                     Eccen(n)      = E.Eccentricity;
%                 end
%             end
% %             BW_nth_tot_area(n) = sum(grid_area_2D(BW_nth(:,:,n)), 'all');
%         end
%     else
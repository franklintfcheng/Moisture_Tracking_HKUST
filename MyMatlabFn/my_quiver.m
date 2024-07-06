%%

function my_quiver(varargin)
    
    pnames       = {'X', 'Y', 'U', 'V', 'S', 'quiver_step_lat', 'quiver_step_lon', 'vector_colmap', 'headwidth', 'headlength', 'quiver_linewidth', 'linelength'};
    dflts        = { [],  [],  [],  [],  [],                [],                [],              [],           5,            8,                  1,         0.08};
     
    [                X,   Y,   U,   V,    S,  quiver_step_lat,    quiver_step_lon,   vector_colmap,   headwidth,   headlength,   quiver_linewidth,  linelength] = ...
            internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
    
    fprintf('*** Running my_quiver... ***\n')
    %/ X and Y should be 2D lon lat matrices from meshgrid().
    
    %/ Check and do not plot outside xlim and ylim (otherwise a bug will show up)
    if numel(X) ~= 1 %/ only when X Y U V are matrices. 
        xL = xlim;
        yL = ylim;

        cond = (X >= xL(1)) & (X <= yL(2)) & (Y >= yL(1)) & (Y <= yL(2));

        X(~cond) = nan;
        Y(~cond) = nan;
        U(~cond) = nan;
        V(~cond) = nan;

        if ~isempty(quiver_step_lat) && ~isempty(quiver_step_lon)
            X = X(2:quiver_step_lat:end-1, 2:quiver_step_lon:end-1);
            Y = Y(2:quiver_step_lat:end-1, 2:quiver_step_lon:end-1);
            U = U(2:quiver_step_lat:end-1, 2:quiver_step_lon:end-1);
            V = V(2:quiver_step_lat:end-1, 2:quiver_step_lon:end-1);
        end
    end
    
    %/ reshape to 1D vectors before arrow3()
    x1_reshp = reshape(X, [], 1);
    x2_reshp = reshape(linelength*U, [], 1);
    
    y1_reshp = reshape(Y, [], 1);
    y2_reshp = reshape(linelength*V, [], 1);
    
    arrow3([x1_reshp, y1_reshp], [x1_reshp + x2_reshp, y1_reshp + y2_reshp], S, headwidth, headlength);
    
%     for ii = 1:size(X, 1)
%         for jj = 1:size(X, 2)
%             x1 = X(ii,jj);
%             y1 = Y(ii,jj);
%             x2 = linelength*U(ii,jj);
%             y2 = linelength*V(ii,jj);
%             
%             
%             arrow3([x1, y1], [x1+x2, y1+y2], S, headwidth, headlength);
% 
%         end
%     end
%                 
%     h = arrow3(p1,p2,[],p(:,4),[],'cone');
%     % Change truecolor CData to colormap values.
%     c = p(:,5);
%     for n = 1:length(c)
%     set(h(n),'CData',c(n)*ones(2,21))
%     end
            
end
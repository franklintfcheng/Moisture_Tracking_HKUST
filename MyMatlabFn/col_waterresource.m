function colmap = col_waterresource(NoOfColors)
    %/ Colors *manually* captured from https://tn.water.usgs.gov/cgi-bin/grapher/graph_colormap_setup.pl
    %/ website setting:
    %/                  custom three-color gradient: CC9900 CCFFCC 0099FF
   
    if NoOfColors == 9
        colmap = [203 153 0;
                  204 179 51;
                  204 203 102;
                  204 230 153;
                  204 255 204;
                  153 230 217;
                  102 204 230;
                  51  179 242;
                  26  153 255]./255;
              
    elseif NoOfColors == 11
        colmap = [203 153 0;
                  204 173 41;
                  203 194 82;
                  204 214 122;
                  204 235 163;
                  204 255 204;
                  163 235 214;
                  122 214 224;
                  82  194 235;
                  41  173 245;
                  26  153 255]./255;

    else
        error('No result for the requested NoOfColors!')
    end



end

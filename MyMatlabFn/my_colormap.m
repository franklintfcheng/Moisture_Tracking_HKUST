%%
function colmap = my_colormap(varargin)

    %/ Author:      Fandy
    %/ Last Update: Feb 7, 2022
    
    %=========================================================================
    %/ Visit https://pratiman-91.github.io/colormaps/docs/collections/ncar_ncl
    %/    to collect more colmap!
    %=========================================================================
    
    switch nargin
        case 2   %/ given 2 inputs
            NoOfColors  = varargin{1};
            name        = varargin{2};
            flip_or_not = [];
            
        case 3   %/ given 3 inputs
            NoOfColors  = varargin{1};
            name        = varargin{2};
            flip_or_not = varargin{3};

        otherwise
            error('Unexpected inputs')
    end
    
    if nargin == 1      name = [];          end
    
    %/ database of colormap that I collected before.
    colmap = [];
    
    if isequal(name, 'precip')
        colmap = jet(NoOfColors-1); 
        colmap(end+1, :) = [244 78 224]./255; %/ append purple color
    
    elseif isequal(name, 'precip2')
        colmap = brewermap(NoOfColors-1, '*spectral'); 
        colmap(end+1, :) = [137 117 241]./255; %/ append purple color 
        
    elseif NoOfColors == 6 && isequal(name, 'precip3_6lev')
        
        colmap = [255 255 255;
                  182 202 255;
                  98  113 246;
                  58  189  58;
                  254 249  13;
                  230   0   0;]./255;
              
    elseif NoOfColors == 12 && isequal(name, 'radar_12lev')
        
        colmap = [171 255 255;
                   84 161 255;
                   23   0 255;
                  127 230  91;
                   77 205  65;
                   42 179  54;
                   24 154  58;
                  255 255 177;
                  255 205 102;
                  255 137  75;
                  205  58  58;
                  138 103 206;]./255;
              
    elseif NoOfColors == 10 && isequal(name, 'perc2_10lev')
        colmap = [218 227 237;
                    186 203 251
                    152 179 249;
                    133 152 248;
                    180 207 113;
                    235 245 169;
                    254 250  85;
                    248 212  77;
                    242 169  60;
                    235  89   41;]./255;

    elseif NoOfColors == 11 && isequal(name, 'precip3_11lev') %/ a scary colorbar for precip :)
        colmap = [255 255 255;
                  182 202 255;
                  128 151 255;
                  98  113 246;
                  47  159  24;
                  186 249 110;
                  254 249  13;
                  253 199   1;
                  251 163   1;
                  247   0   0;
                  101   0 101;]./255;
       
    elseif NoOfColors == 11 && isequal(name, 'precip3_11lev_nowhite') %/ a scary colorbar for precip :)
        colmap = [235 245 255;
                  182 202 255;
                  128 151 255;
                  98  113 246;
                  47  159  24;
                  186 249 110;
                  254 249  13;
                  253 199   1;
                  251 163   1;
                  247   0   0;
                  101   0 101;]./255;
              
    elseif NoOfColors == 11 && isequal(name, 'radar_11lev')  
        colmap = [192 138 247;
                   74   8 140;
                  128 151 255;
                  98  113 246;
                  47  159  24;
                  186 249 110;
                  254 249  13;
                  253 199   1;
                  251 163   1;
                  247   0   0;
                  153   0  51;]./255;
              
    elseif NoOfColors == 6 && isequal(name, 'rainbow1_6lev')
        colmap = [245  45  37;
                  250 195  31;
                  243 235  64;
                  163 221 110;
                   65 230 239;
                  203 158 194;]./255;
        
    elseif NoOfColors == 6 && isequal(name, 'set1_6lev')
        
        colmap = [249 128 114;
                  128 117 211;
                  251 180  98;
                  179 222 105;
                  188 128 189;
                  251 232 109]./255;

    elseif NoOfColors == 9 && isequal(name, 'GnRdPu_9lev')
        
        colmap = [158  202   58;
                   53  179   84;
                   38  135   56;
                  255  239    4;
                  228  190   26;
                  242  140   45;
                  233   43   44;
                  193   28   30;
                  173   80  154;
                  ]./255;      
              
    elseif NoOfColors == 10 && isequal(name, 'GnRdPu_10lev')
        
        colmap = [158  202   58;
                   53  179   84;
                   38  135   56;
                  255  239    4;
                  228  190   26;
                  242  140   45;
                  233   43   44;
                  217   30   35;
                  193   28   30;
                  173   80  154;
                  ]./255;
    
    elseif NoOfColors == 27 && isequal(name, 'topo_land_27lev')
%         42  109 181  %/ dark blue for oceans
%         30 160 244;  %/ light blue for oceans
        colmap = [30  160 244;
                  79  123  49;
                  108 153  44;
                  146 185  51;
                  164 199  72;
                  174 205  92;
                  185 210 108;
                  205 220 139;
                  226 234 170;
                  244 236 196;
                  241 230 183;
                  238 225 170;
                  232 217 151;
                  222 203 125;
                  212 191 103;
                  206 177  75;
                  203 161  59;
                  198 152  45;
                  185 134  37;
                  164 112  34;
                  140  90  38;
                  120  77  35;
                  103  67  31;
                   86  55  25;
                  135 110  98;
                  184 173 168;
                  246 246 250;
                  ]./255;
        
    elseif NoOfColors == 12 && isequal(name, 'sunshine_12lev')   %/ try this!
            colmap = [  73   8  112;
                       149  14  223;
                       183  75  243;
                       203 126  246;
                       225 180  250;
                       236 208  252;
                       255 245  204;
                       255 230  112;
                       253 204   51;
                       251 175   52;
                       249 110    1;
                       230  40   31;]./255;

   elseif NoOfColors == 12 && isequal(name, 'sunshine_12lev_white')   %/ try this!
            colmap = [  73   8  112;
                       149  14  223;
                       183  75  243;
                       203 126  246;
                       225 180  250;
                       255 255  255;
                       255 255  255;
                       255 230  112;
                       253 204   51;
                       251 175   52;
                       249 110    1;
                       230  40   31;]./255;
                   
    elseif NoOfColors == 12 && isequal(name, 'BuYlRd_12lev')   %/ try this!
        colmap = [9    50  121;
                  30  110  200;
                  60  160  240;
                  130 210  255;
                  160 239  254;
                  230 255  255;
                  255 250  220;
                  255 232  120;
                  251 160    0;
                  248  96    1;
                  224  20    0;
                  165   0    0]./255;

    elseif NoOfColors == 16
        if isequal(name, 'amwg_blueyellowred_16lev')
            colmap = [ 130 32  240;
                       0    0  182;
                       43   57 255;
                       65  106 225;
                       30  144 255;
                       46  191 255;
                       160 210 255;
                       210 245 255;
                       255 255 200;
                       254 225  50;
                       251 170   0;
                       249 109   0;
                       247   0   0;
                       160  35  35;
                       153   0 106;
                       248 105 180; ]./255;
                   
        elseif isequal(name, 'amwg_blueyellowred_16lev_white')
            colmap = [ 130 32  240;
                       0    0  182;
                       43   57 255;
                       65  106 225;
                       30  144 255;
                       46  191 255;
                       160 210 255;
                       255 255 255;
                       255 255 255;
                       254 225  50;
                       251 170   0;
                       249 109   0;
                       247   0   0;
                       160  35  35;
                       153   0 106;
                       248 105 180; ]./255;
        end
    end
    
    if NoOfColors == 4
        if isequal(name, 'BuPi') || isempty(name) %/ by default
            colmap = [236  184  215;
                 237  120  185;
                 122  151  185;
                  91  127  183]./255;
        end
    end
    
    if isempty(colmap)
        error('No colormap fits your request!!');
    end
    
    if isequal(flip_or_not, 'flip')
        colmap = flip(colmap, 1);
    end
    
    
end
function output = box_region(region_name)
    
    %======================================================================
    %/      Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: 28 Jun 2024
    %======================================================================

    if ischar(region_name)  
        region_name = {region_name};
    end
    shift = 0.01;  %/ to avoid overlapping.
    
    %----------------------- Ocean (with land mask) --------------------------%
    %/ Global (lon: -180, 180, lat: -90, 90)
    vertices.global = [-180+shift  -90
                       -180+shift   90
                        180         90
                        180        -90
                       -180+shift  -90];
    
    %/ Niño 4 (5N-5S, 160E-150W)
    vertices.Nino4  = [ 160  -5
                        160   5
                       -150   5
                       -150  -5
                       -160  -5];
    
    %/ Niño 3.4 (5N-5S, 170W-120W): 
    vertices.Nino34 = [-170  -5
                       -170   5
                       -120   5
                       -120  -5
                       -170  -5];
    
    %/ Niño 3 (5N-5S, 150W-90W)
    vertices.Nino3  = [-150  -5
                       -150   5
                        -90   5
                        -90  -5
                       -150  -5];
    
    %/ Gulf of Guniea (outline by hand)
    vertices.GoG    = [-8    4
                       -8    3
                       12  -17
                       16  -17
                       16    8
                       -8    8
                       -8    4];
                   
    %/ Northern Tropical Atlantic Ocean (not overlapped with GoG)
    vertices.NTAO    = [ -68    12
                         -8.1   12
                         -8.1    3
                         -68     3
                         -68    12];
    
    %/ Southern Tropical Atlantic Ocean
    vertices.STAO    =   [-60         3-shift
                          -8-shift    3-shift
                           2-shift   -7
                          -36        -7
                          -60         3-shift];
                     
    %/ Eastern South Atlantic Ocean
    vertices.ESAO   = [-18        -7-shift
                        2-shift   -7-shift
                        14        -19
                        20        -35
                        20        -60+shift
                       -18        -60+shift
                       -18        -7-shift];
                   
    %/ Western South Atlantic Ocean
    vertices.WSAO   = [-36        -7-shift
                       -18-shift  -7-shift
                       -18-shift  -60+shift
                       -70        -60+shift  %< 
                       -70        -42
                       -36        -7-shift];
                   
    %/ Gulf of Mexico
    vertices.GoM    = [-81  26
                       -81  22.7
                       
                       -88  21
                       -92  17
                       -99  17
                       -99  31
                       -82  31
                       -81  26];
    
    %/ Caribean Sea
    vertices.CaribeanSea  = [-88        21
                             -90        16
                             -84        10
                             -82        8.5
                             -80        9
                             -78        9
                             -77        7
                             -76        7
                             -68-shift  7
                             -68-shift  12+shift
                             -61        12+shift
                             -61        16
                             -66        18
                             -74        20
                             -81        22.7-shift
                             -88        21];
    
                         
    %/ Hudson Bay
    vertices.HudsonBay    = [-88     66
                             -86     66
                             -84     66
                             -76     62
                             -76     51
                             -80     51
                             -95     58
                             -95     61
                             -88     66];
                         
    %/ Western North Atlantic Ocean
    vertices.WNAO         = [-44        12+shift
                             -61+shift  12+shift
                             -61+shift  16+shift
                             -66+shift  18+shift
                             -74+shift  20+shift
                             -81+shift  22.7
                             -81        26
                             -82        31
                             -76        40
                             -76        51
                             -76        62
                             -84        66
                             -68        66
                             -44        66
                             -44        12+shift];
    
    %/ Eastern North Atlantic Ocean
    vertices.ENAO         = [-44+shift  12+shift
                             -6         12+shift
                             -6         43
                             -2         43
                              20        54
                              30        60
                              26        66
                             -44+shift  66
                             -44+shift  12+shift];
    
                         
    %/ Arctic Ocean
    vertices.AO    = [-180+shift  66+shift
                      -180+shift  90
                       180        90
                       180        66+shift
                        45        66+shift
                        38        62.5
                        26        66+shift
                      -180+shift  66+shift];
    
    %/ MeditSea
    vertices.MeditSea     = [-6   37
                              4   44
                              14  46
                              32  40
                              38  36
                              34  31
                              32  30
                              18  30
                             -6   35
                             -6   37];
    
    %/ BlackSea      
    vertices.BlackSea     = [30   41
                             26   42
                             32   48
                             40   47
                             42   41
                             30   41];
                              
    %/ CaspianSea
    vertices.CaspianSea   = [48   46
                             51   47
                             53   47
                             54   46
                             55   42
                             54   36
                             48   37
                             46   44
                             48   46];
    
    %/ Arabian Sea
    vertices.AS    = [59.5     26
                      62       26
                      78       26
                      78       8
                      42.75    8
                      42.75    12
                      59.5     22
                      59.5     26];  
                  
    
    %/ PersianOmanGulf
    vertices.PersianOmanGulf  = [48          29
                                 52          23
                                 59.5-shift  23
                                 59.5-shift  26
                                 57          27
                                 49          31
                                 48          29];
                             
    %/ Bay of Bengal
    vertices.BoB   = [78-shift 8
                      78-shift 21
                      91.5     22.75
                      99       17
                      99       8
                      78-shift 8];
                  
    %/ Red Sea
    vertices.RedSea  = [32     30
                        36     30
                        44     14
                        44     13
                        42.75  12
                        38     17
                        32     30];
    
                  
    %/ Western Indian Ocean
    vertices.WIO  = [78        8-shift
                     78       -60+shift
                     20+shift -60+shift
                     20+shift -35
                     27       -33
                     30       -30
                     34       -20
                     38       -17
                     38       -5
                     50        8-shift
                     78        8-shift];
                  
    %/ Eastern Indian Ocean
    vertices.EIO  = [78+shift   8-shift
                     99         8-shift
                     103-shift  4
                     103-shift -1
                     104       -5
                     107       -7
                     116       -9
                     126       -9
                     133       -12
                     130       -15
                     124       -17
                     121       -20
                     116       -22
                     116       -31
                     132       -31
                     138       -33
                     147       -38
                     147       -60+shift
                     78+shift  -60+shift
                     78+shift   8-shift];
    
    %/ South China Sea
    vertices.SCS   = [118      25
                      121      25
                      121      22
                      122      21
                      122      18
                      120      12
                      117       8
                      117       6
                      111       1
                      111      -3
                      106      -3
                      103      -1
                      103       4
                      102       6
                      105       9
                      108      11
                      108      16
                      105      19
                      108      22
                      112      22
                      118      25];
    
    %/ Gulf of Thailand
    vertices.GoT   = [99         14
                      102        13
                      108-shift  11
                      105-shift  9
                      102-shift  6
                      99         8
                      99         14];
    
    %/ Western Tropical Pacific
    vertices.WTPO = [150        -10
                     128+shift  -10
                     126        -9+shift
                     116        -9+shift
                     107        -7+shift
                     104        -5+shift
                     105        -3-shift
                     111+shift  -3-shift
                     111+shift   1
                     117+shift   6
                     117+shift   8
                     119        10
                     180        10
                     180       -10
                     150       -10];
    
    %/ Yellow Sea & East China Sea
    vertices.YSECS = [127       25
                      121+shift 25
                      121+shift 25+shift
                      118       25+shift
                      120       27
                      120       30
                      121       32
                      117       39
                      121       41
                      126.5     41
                      126.5     35
                      126.5     33
                      131       33
                      130       29
                      127       25];     
    
    %/ Sea of Japan
    vertices.SoJ   = [131          33+shift
                      126.5+shift  33+shift
                      126.5+shift  41
                      138          47
                      140          52
                      143          52
                      142          45
                      140          39
                      138          36
                      131          34
                      131          33+shift];
    
    %/ Western North Pacific
    vertices.WNPO = [180        10+shift
                     119        10+shift
                     120+shift  12
                     122+shift  18
                     122+shift  21
                     121+shift  22
                     121+shift  25-shift
                     127+shift  25-shift
                     130+shift  29
                     131+shift  33
                     131+shift  34
                     138        36
                     140+shift  39
                     142+shift  45
                     143        52+shift
                     135        52+shift
                     135        55
                     143        60
                     154        60
                     156        62
                     164        63
                     180        66
                     180        10+shift];
                  
                  
    %/ Eastern North Pacific
    vertices.ENPO = [-180+shift  66
                     -160        66
                     -160        62
                     -150        62
                     -134        59
                     -122        49
                     -122        38
                     -112        31
                     -99         17
                     -90         16
                     -85         10+shift
                     -180+shift  10+shift
                     -180+shift  66];
    
    %/ Eastern Tropical Pacific 
    vertices.ETPO   = [-85         10
                       -82         8.5
                       -80         9
                       -78         9
                       -77         7
                       -77        -10
                       -180+shift -10
                       -180+shift  10
                       -85         10];
    
    %/ Eastern South Pacific 
    vertices.ESPO = [-77        -10-shift
                     -65-shift  -18
                     -70-shift  -42+shift
                     -70-shift  -60+shift   %< 
                     -180+shift -60+shift
                     -180+shift -10-shift
                     -77        -10-shift];
    
    %/ Western South Pacific
    vertices.WSPO = [129       -10-shift
                     133       -12
                     136       -16
                     140       -18
                     145       -18
                     152       -27
                     147+shift -38
                     147+shift -60+shift
                     180       -60+shift
                     180       -10-shift
                     129       -10-shift];
                 
                 
    %/ Southern Ocean
    vertices.SO     = [-180+shift     -90
                       -180+shift     -60
                        180           -60
                        180           -90
                       -180+shift     -90];
       
                  
    % vertices.WestIO = [38     8
    %                    77.75  8
    %                    77.75  -20
    %                    30     -20
    %                    38     -8
    %                    38     8];
    %             
    % vertices.EastIO = [77.75  8
    %                    100    8
    %                    103    0
    %                    105   -3
    %                    116   -3
    %                    120   -5.5
    %                    133   -10
    %                    133   -20
    %                    77.75 -20
    %                    77.75 8]; 
    
    % vertices.NSCS = [120 24
    %                 121 22
    %                 121 12
    %                 119 11
    %                 105 11
    %                 105 24
    %                 120 24];  
    %             
    % vertices.SSCS = [119.25 11
    %                 117.75  9
    %                 117     8
    %                 117     6
    %                 113     2
    %                 111     -3
    %                 99      -3
    %                 99      14
    %                 119.25  11];    
    %                    
    % vertices.PS = [131 33
    %                131 34
    %                140 36.75
    %                140 7
    %                117 7
    %                117 24
    %                131 33];
    % 
    % vertices.NP = [ 136   7
    %                 136   65
    %                 190   65
    %                 190   7
    %                 116   7];
    % 
    % vertices.TropPO = [116 7
    %                    190 7
    %                    190 -10
    %                    116 -10
    %                    116 7];
    %             
    % vertices.SouthPO = [133 -10
    %                     190 -10
    %                     190 -20
    %                     133 -20
    %                     133 -10];
    
                
    % vertices.SoJ   = [130    33
    %                   127.5  35.5      
    %                   127.5  40
    %                   141.5  52.25
    %                   142    51.5
    %                   142    43
    %                   140.25 43
    %                   140.25 38
    %                   136    35.5
    %                   130    33];
    % 
    %                    
    % vertices.YSECS = [117    41
    %                   118    33
    %                   119    33
    %                   120    31.75
    %                   120    30.25
    %                   118    24
    %                   122    24
    %                   123    24
    %                   125    25
    %                   126    25
    %                   130    29
    %                   130    30
    %                   130.75 30.5
    %                   130.75 33
    %                   130    33
    %                   128    35
    %                   126    41
    %                   117    41];     
                  
    %--------------------------------- Land ----------------------------------%
    
    %/ All Land (for double-checking)
    vertices.AllLand = [-180+shift  -90
                        -180+shift   90
                         180         90
                         180        -90
                        -180+shift  -90];
    
    %/ Pearl
    vertices.('Pearl')     = [106    20
                              105.5  20.5
                              105    20.7
                              104.5  20.7
                              103.5  21.7
                              102.6  21.3
                              102.4  22
                              102    22.4
                              101.5  22.5
                              101    23.5
                              100.9  24
                              100    25.5
                              100.5  25.5
                              101.5  24.8
                              102    25.5
                              102.3  25.5
                              102.3  24.5
                              102.8  24.5
                              103    25
                              103.3  25.3
                              103.7  25.3
                              103.5  25.5
                              103.8  26
                              103.8  26.2
                              104.1  26.3
                              103.8  26.4
                              104    26.7
                              104.5  26.7
                              105    26.5
                              105.5  26.3
                              109    26.3
                              109    26
                              110    26
                              110    26.3
                              110.5  26
                              110.5  25.3
                              111.5  25
                              111.5  24.6
                              112.1  24.7
                              112.1  25.2
                              113    25.6
                              114.7  25.4
                              114.2  24.7
                              114.7  24.6
                              115.5  25.2
                              115.9  25
                              116    25.7
                              116.5  26
                              116.5  25.7
                              117    25.5
                              117.5  25.8
                              118    25.5
                              119.3  25.5
                              118.8  25
                              118.6  24.55
                              118.2  24.55
                              118.2  24.25
                              117.5  23.7
                              116    22.8
                              115    22.7
                              114.5  22.7
                              114.2  22.3
                              113.6  22.7
                              113.5  22.25
                              113.1  22.5
                              113    22
                              111.5  21.5
                              110.5  21.3
                              110.3  21
                              110.5  20.5
                              110    20.4
                              109.7  21
                              109.8  21.4
                              109.4  21.4
                              109    21.6
                              108.5  21.6
                              108    21.5
                              107    21
                              106    20];
                          
    %/ New Guinea
    vertices.('NewGuinea') = [131   -1
                              132   -0.5
                              134   -0.8
                              135   -3
                              138   -1.5
                              145   -4
                              148   -6
                              148.3 -8
                              150.7 -10
                              150   -10.5
                              147.5 -10
                              146   -8
                              144   -8
                              143   -9.3
                              141   -9
                              140   -8
                              137.8 -8.2
                              138.2 -7.6
                              138.6 -7.5
                              138   -5.5
                              134   -4
                              132.7 -4
                              132   -2.8
                              133.7 -2.5
                              133.9 -2
                              133   -2.3
                              132   -2.3
                              131, -1];
    
    %/ Eastern north South America
    vertices.('NESA') = [-48  -0.5
                     -49   -2
                     -49   -4
                     -47   -5
                     -46   -6
                     -46  -10
                     -44  -10
                     -42   -9
                     -40   -7
                     -37   -8
                     -36  -10
                     -34  -10
                     -34 -0.5
                     -48  -0.5];
    
    %/ Western North Mexico
    vertices.('WNMexico') = [-114   32
                                -112   29
                                -109   26
                                -105   22
                                -104   24
                                -106   25
                                -109   31
                                -110   32
                                -110   31
                                -112   32
                                -114   32];
    
    %/ Western Central Mexico
    vertices.('CMexico') = [-104  23
                            -103  23
                            -102  22
                            -101  25
                            -100  24.9
                            -99   25.5
                            -97   25.5
                            -98   22
                            -97   20
                            -95   18
                            -96   17
                            -98   18
                            -102  18
                            -104  19
                            -105  20
                            -105  21
                            -104  23];
    
    %/ Western South Mexico
    vertices.('SMexico') = [-100  17
                            -97   17
                            -96   16.9
                            -94   18
                            -90   21
                            -87   21
                            -89   16
                            -84   15
                            -84   11
                            -83    9
                            -100  17];
                  
    %/ Pan-Tibetan Plateau 
    vertices.('Pan_TP_rough') = [  62   28
                                   62   35
                                   70   43
                                   80   46
                                   94   44.1
                                  109   40
                                  105   22  
                                   62   28];
    
    %/ Pan-East Asian Monsoon Region (to TP)
    vertices.('Pan_EAM') = [ 100+shift   45;
                             140         45;
                             140         33;  %/ add mid-points for the boundary to curve under the 'robin' projection
                             140         20;
                             100+shift   20;  
                             100+shift   33;  %/ add mid-points for the boundary to curve under the 'robin' projection
                             100+shift   45;];  
                         
    %/ Pan-Indian Monsoon Region (5S the southmost)
    vertices.('Pan_IM') = [  80    33;
                             78   32.5;
                             78    32;
                             76    32;
                             70    30;
                             57    23;
                             44    13;
                             35    -5;
                             100   -5;
                             100    5;
                             100   33;
                             80    33];
    
    %/ Easterly Region in north India 
    vertices.('Easterly_North_India') =  [65  30;
                                          95  30;
                                          95  20;
                                          65   20;
                                          65  30]; 
    
    vertices.('Pan_NH_MidLat_Westerly') = [ 180        15+shift;
                                            -180+shift 15+shift;  
                                            -180+shift 60;
                                            180        60;
                                            180        15+shift;]; 
    
    vertices.('Arunachal_Pradesh') = [91  28;
                                      91  30;
                                      98  30;
                                      98  28;
                                      91  28;];     
    
    vertices.('TP_southernflank') = [ 82   27;
                                      82   29;
                                      85   29;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      90   29;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      95   29;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      102  29;
                                      102  27;
                                      95   27;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      90   27;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      85   27;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      82   27;];
    
    vertices.('southern_TP') =      [ 82   27;
                                      82   30;
                                      85   30;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      90   30;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      95   30;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      102  30;
                                      102  27;
                                      95   27;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      90   27;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      85   27;  %/ add mid-points for the boundary to curve under the 'robin' projection
                                      82   27;];
    
    %/ Indonesian-Australian monsoon index (CMJO study)
    vertices.('IAM_index') = get_box_bndry([120, 150], [-15, -5]);
                                  
    %/ SCS summer monsoon index (CMJO study)
    vertices.('SCSSM_index') = get_box_bndry([110, 120], [10, 20]);
       
    %/ Monsoon onset over Kerala (MOK) index (CMJO study)
    vertices.('MOK_index') = get_box_bndry([60, 80], [5, 15]);
    
    %/ Yurong's AR research region 
    vertices.('ALA') = get_box_bndry([-170, -140], [55, 65]);
    
    
    %/ Tibetan Plateau -> retrieve from contour() using C2xyz()!
    if any(ismember(region_name, {'TP', 'Pan_TP'}))
        res = 0.25;   %/ use a higher resolution to outline the boundary.
        map_lon_lower = 55;
        map_lon_upper = 121;
        map_lat_lower = 17;
        map_lat_upper = 47;
        lon_grids = map_lon_lower:res:map_lon_upper;
        lat_grids = map_lat_lower:res:map_lat_upper;
    %     [lon_grids_2D, lat_grids_2D] = meshgrid(lon_grids,lat_grids);
    
        if ismember(region_name, {'Pan_TP'})
            plateau_hgt = 1500;
        elseif ismember(region_name, {'TP'})
            plateau_hgt = 3000;
        else
            error('code not set!');
        end
        avg_topo_map = interp_topo('lon_grids', lon_grids, 'lat_grids', lat_grids);
    %     figure(1)
    %     C = contour(lon_grids_2D, lat_grids_2D, avg_topo_map',...
    %                       [plateau_hgt, plateau_hgt]);
    %     close(figure(1))
        C = contourc(lon_grids, lat_grids, avg_topo_map', [plateau_hgt, plateau_hgt]);
    
        [Lon, Lat, ~] = C2xyz(C);  %/ magic function!
        
        %/ NOTE: By default, we select the polygon with the most # of vertices.
        [~, I] = max(cellfun(@(x) length(x), Lon));
        vertices.(region_name{:}) = [Lon{I}', Lat{I}'];
        
    end
                               
    %/ Antarctica
    vertices.('Antarctica') = [180       -60;
                              -180+shift -60;
                              -180+shift -90;
                               180       -90;
                               180       -60;];
                               
    output = [];
    for i = 1:length(region_name)
        output = [output; [vertices.(region_name{i}); [nan nan]]];         %/ use [nan nan] to separate boundary data.
    end
end

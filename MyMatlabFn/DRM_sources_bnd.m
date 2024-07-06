function output = DRM_sources_bnd(source_name)

    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: May 24, 2024
    %/
    %/ DESCRIPTION:
    %/       This function is designed for outputing the boundary of the
    %/       30 source regions defined in Cheng and Lu 2020 and Cheng et
    %/       al. (2021)
    %/
    %/ REFERENCES:
    %/  Cheng, T.F., & Lu, M. (2020): Moisture Source–Receptor Network of the East Asian Summer Monsoon Land Regions and the Associated Atmospheric Steerings. J. Climate, 33, 9213–9231, https://doi.org/10.1175/JCLI-D-19-0868.1.
    %/  Cheng, T.F., Lu, M. & Dai, L. (2021): Moisture channels and pre-existing weather systems for East Asian rain belts. npj Clim Atmos Sci 4, 32. https://doi.org/10.1038/s41612-021-00187-6
    %/
    %/=====================================================================
    
    vertices.EASM_110140E_2040N  = [110    20
                                    110    40
                                    140    40
                                    140    20
                                    110    20]; % rows are vertices, columns are lat/lon
    
    vertices.E1_SC  = [ 110    21 
                        110    27.5
                        120.75  27.5
                        120    26.75
                        119.75    26.5
                        119.5    25.5
                        119.25    25.25
                        118.75    24.75
                        118.25    24.5
                        118    24
                        117.75    24
                        117.25    23.5
                        116.75    23.25
                        116.5    23
                        116      22.75
                        114.75    22.75
                        114.25    22.25
                        113.75    22.75
                        113.5     22.25
                        113       22
                        112.5     21.75
                        112.25    21.75
                        111.75    21.75
                        111       21.5
                        110    21]; % rows are vertices, columns are lat/lon
    
    vertices.E1_TW     = [120.25  23.75
                          120.5   24.25
                          121     25
                          121.5   25.25
                          122     25
                          121.75  24.75
                          121.75  24
                          121.5   23.75
                          121.25   23
                          121     22.5
                          120.75   22
                          120.5   22.5
                          120.25   22.75
                          120.25  23.75]; % rows are vertices, columns are lat/lon
    
    vertices.E2   = [   110    27.5 
                        110    33
                        120.75    33
                        121.75    32
                        121.75    31.75
                        121.25   31.75
                        121.5    31.5
                        122      31
                        120.5    30.25
                        121.25   30.25
                        122      30
                        121.75    29
                        121.5    28.25
                        121      28
                        120.5    27.5
                        110    27.5]; % rows are vertices, columns are lat/lon
    
    
    %/ Outline the E3 region
    vertices.E3    = [ 129.25, 33.5
                        130.25, 33.75
                        130.75, 34
                        131, 34.5
                        131.25, 34.5
                        131.75, 34.75
                        131.75, 34.75
                        132.75, 35.5
                        135.25, 35.75
                        135.75, 35.75
                        136, 36.25
                        136.5, 36.75
                        136.75, 37.25
                        137.25, 37.5
                        137.25, 37
                        137.5, 37
                        138.5, 37.5
                        139, 38
                        139.25, 38.25
                        139.75, 39
                        140,    39
                        140,    35.75
                        139.5,    35.25
                        139.25,    35.25
                        138.75,    34.5
                        138.5,    35
                        138,    34.5
                        137.25,    34.5
                        136.75,    35
                        136.5,     34.75
                        136.75,    34.5
                        136.75,    34.25
                        136.25,    34
                        135.75,    33.5
                        135.75,    33.5
                        135.25,    33.75
                        135,       34.25
                        134.75,    34
                        134.75,    33.75
                        134.25,    33.5
                        134,    33.25
                        133.75,    33.5
                        133.5,    33.5
                        132.75,    32.75
                        132.25,    33
                        131.75,    32.75
                        131.5,    32
                        131.5,    31.5
                        130.75,    31
                        130.25,    31.25
                        130,    31.5
                        130,    32
                        129.75,    32.25
                        129.75,    32.75
                        129.25, 33.5];
    
    vertices.E4    = [ 110     33
                        110     39.25
                        117.75  39.25
                        117.5  38.75
                        117.5  38.75
                        118     38.25
                        118.75  38
                        119     37.75
                        118.75  37.5
                        119     37
                        119.75  37
                        120.5   37.75
                        121     37.75
                        122.5     37.25
                        122.5     37
                        121.5     36.75
                        120.75    36.25
                        120       35.75
                        119.25    35
                        119.5    34.75
                        120.5    34.25
                        120.5    34.25
                        120.75    33
                        110    33]; % rows are vertices, columns are lat/lon
                    
    %/ Outline the E5 region
    vertices.E5    = [
                        124.25,  40;
                        124.5,   40.25;
                        126,     41;
                        126.5,   41.5;
                        127,     41.75;
                        127.25,  41.5;
                        128,     41.5;
                        128.25,  41.5;
                        128,     42;
                        129,     42.25;
                        129.5,    42.5;
                        130,      43;
                        130.25,   42.75;
                        130.75,   42.5;
                        130	    42
                        129.75	41.75
                        129.75	41
                        129.25	40.5
                        128.25	40
                        128     40
                        127.5   39.75
                        127.5   39.75
                        127.5   39.25
                        128.5   38.75
                        128.5   38.5
                        128.75  38.25
                        128.75  38
                        129.5   37.25
                        129.5   36.25
                        129.75   36
                        129.75   36
                        129.5   35.5
                        128.75   35
                        126.5   34.25
                        126.25   35
                        126.5   35.75
                        126.5   36.25
                        126.25   36.75
                        126.5     37
                        126.5     37.25
                        126     37.75
                        125     37.75
                        124.75     38
                        125.25     39.25 
                        124.5     39.75 
                        124.25	40]; % rows are vertices, columns are lat/lon
                    
    %/ Outline the E6 region
    vertices.E6    = [ 110,     39.25;
                        110,     42.5;
                        110.75,  43.25;
                        111.75,  43.75;
                        111.75,  44;
                        111.25,  44.5;
                        111.75,  45;
                        131.25    45;
                        131.25,   43;
                        130.25,   42.75;
                        130,      43;
                        129.5,    42.5;
                        129,     42.25;
                        128,     42;
                        128.25,  41.5;
                        128,     41.5;
                        127.25,  41.5;
                        127,     41.75;
                        126.5,   41.5;
                        126,     41;
                        124.5,   40.25;
                        124.25,  40;
                        123.5,   39.75;
                        123.25,  39.75;
                        122.5,   39.5;
                        121.5,   38.75;
                        121.25,   39;
                        121.75,   39.25;
                        121.5,   39.5;
                        121.5,   39.75;
                        122,     40;
                        122.5,   40.5;
                        122,     40.75;
                        121.25,  41;
                        121,     40.75;
                        120.25,     40.25;
                        119.75,     40;
                        119.25,     39.5;
                        119,        39.25;
                        110,     39.25]; % rows are vertices, columns are lat/lon
                    
    vertices.EastAfrica = [ 30    31.5 
                            32    31.5
                            32.5  31
                            34.25 31.25
                            34.75 29.75
                            34.75 29.5
                            34.25 28
                            42    15
                            44    11.75
                            55    13
                            55    -40
                            30    -40];
                      
    vertices.Indochina = [92    21
                          93.25 22.25
                          93.25 23
                          93.5  23
                          93.5  23.75
                          94    23.75
                          94.5  24.5
                          94.75 25.25
                          95.25 25.75
                          95.25 26.5
                          96.25 27.25
                          97    27.5
                          97.5  28.25
                          98    28.25
                          98.5  27.5
                          98.5 26
                          97.5 25
                          97.5  24
                          98.75 24
                          98.75 23.25
                          99.25 23.25
                          99.25 22.25
                          100   22
                          100   21.5
                          101    21.5
                          101.75 21.25
                          101.75 22.25
                          102.5  22.5
                          104   22.5
                          104.25 22.75
                          105    23
                          105.25 23.25
                          106.25 22.75
                          106.75 22.75
                          106.5  22.5
                          106.75 22
                          107    21.75
                          107.75 21.5
                          107.75 17
                          110    17
                          110    8
                          105    8
                          105    1
                          104    1
                          102    2
                          94     10
                          94     16
                          92     21];
      
    vertices.MC =[92 -11
                  92 12.5
                  105 0
                  108 0
                  108 0
                  109 2
                  116 6.75
                  
                  120 19
                  165 19
                  165 -11
                  92 -11]; 
                  
    vertices.Australia = [120 -11
                          150 -11
                          150 -40
                          120 -40
                          120 -11];
    
    vertices.BoB = [77.75 8
                    77.75 21
                    91.5  22.75
                    99    17
                    99    8
                    77.75 8];
    
    vertices.ASnOthers = [32    30
                          78    30
                          78    8
                          42.75 8
                          42.75 13
                          38.25 16
                          32    30];    
                
    vertices.WestIO = [38     8
                       77.75  8
                       77.75  -20
                       30     -20
                       38     -8
                       38     8];
                
    vertices.EastIO = [77.75  8
                       100    8
                       103    0
                       105   -3
                       116   -3
                       120   -5.5
                       133   -10
                       133   -20
                       77.75 -20
                       77.75 8]; 
    
    
    vertices.NSCS = [120 24
                    121 22
                    121 12
                    119 11
                    105 11
                    105 24
                    120 24];  
                
    vertices.SSCS = [119.25 11
                    117.75  9
                    117     8
                    117     6
                    113     2
                    111     -3
                    99      -3
                    99      14
                    119.25  11];    
                
    vertices.SeaOfJapan = [130    33
                           127.5  35.5      
                           127.5  40
                           141.5  52.25
                           142    51.5
                           142    43
                           140.25    43
                           140.25    38
                           136    35.5
                           130    33];
    
                       
    vertices.YSECS = [117    41
                      118    33
                      119    33
                      120    31.75
                      120    30.25
                      118    24
                      122    24
                      123    24
                      125    25
                      126    25
                      130    29
                      130    30
                      130.75 30.5
                      130.75 33
                      130    33
                      128    35
                      126    41
                      117    41];                   
                       
    %                    
    vertices.PS = [131 33
                   131 34
                   140 36.75
                   140 7
                   117 7
                   117 24
                   131 33];
    
    vertices.NP = [ 136   7
                    136   65
                    190   65
                    190   7
                    116   7];
                
    
    vertices.TropPO = [116 7
                       190 7
                       190 -10
                       116 -10
                       116 7];
                
    vertices.SouthPO = [133 -10
                        190 -10
                        190 -20
                        133 -20
                        133 -10];
    
    vertices.SWChina = [97.5  31
                        111 31
                        111 31
                        111 18
                        97.5  18
                        97.5  31];
    
    vertices.CChina = [100      42.5
                       101.75   42.5
                       102.25   42
                       103      42
                       103.5    41.75
                       104.5    41.75
                       105.25   41.5
                       105.5    41.75
                       106      42
                       106.75   42.25
                       109.75   42.5
                       110      42.5
                       110      31
                       100      31
                       100      42];
                    
    vertices.IndianSC = [71.75  36.5
                         71.75  35
                         71.25  34.5
                         71     34
                         70.25  33.75
                         70.5   33.25
                         69.75  33
                         69.5   32.75
                         69.5   32
                         69     31.5
                         68.5   31.75
                         68.25  31.75
                         67.75  31.5
                         66.75  31.25
                         66.5   31
                         66.25  29.75
                         66     29.75
                         64.25  29.5
                         63     29.25
                         61     29.75
                         62     28.5
                         63     28.25
                         63     27.25
                         63.5   27
                         62.75  26.5
                         62.5   26.5
                         62     26.25
                         61.75  25.75
                         61.75  25
                         78     5
                         92     6
                         92     20
                         99     30
                         94     30
                         74     35
                         71.75  36.5];
    
                     
                     
    vertices.MiddleEast = [30    11
                           30    41
                           31    41
                           32.5  41.75
                           35    41.75
                           37.5  41
                           40.5  41
                           41.5  41.5
                           45    41.25
                           49    37.5
                           52    36.5
                           54    37
                           55.5  38
                           57.25 38
                           61    36.5
                           61.25 35.25
                           62.75 35.25
                           65.75 37.5
                           68.5  37
                           72.5  39
                           72.5  11
                           30    11];    
               
    output = vertices.(source_name);
end

%%
function basin_name = hydrosheds_name(varargin)
    
    pnames = {'file_code', 'level'};
    dflts  = cell(1, length(pnames));
    [          file_code,  level] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %%
    if isempty(file_code)
        file_code = {'sa', 'na', 'af', 'as', 'si', 'au', 'eu', 'ar', 'gr'};
    else
        if ischar(file_code)
            file_code = {file_code};
        end
    end
    
%     str_noname = {'noname'};
    str_noname = {''};
    basin_name = cell(length(file_code), 2);  
    contin     = cell(length(file_code), 1);
    
    for k = 1:length(file_code)
        if level == 3
            if ismember(file_code(k), {'sa'})
                contin{k} = 'South America';
                
                %/ ten items in each row
                a = {'Atraro', 'Magdalena', 'N. SA', 'Orinoco', 'Guianas', '', 'Amazon', '', 'Tocantins', '',...
                           '',          '', '', 'Sao Francisco', 'E. SA', 'Uruguay', 'La Plata', '', 'Salado', 'Colorado SA',...
                        'S. SA',          '', 'Chile', 'Peru', '', '', '', '', '', '',...
                   'Titicaca'}';

               
            elseif ismember(file_code(k), {'na'})
                contin{k} = 'North America';
                
                %/ ten items in each row
                a = {'', 'Colorado', 'California', '', 'Columbia', 'Fraser',  'Churchill', 'Nelson', '', '',...
                     '', '',   '',   'St.Lawrence', '', '', 'SE. US', '',  'Mississippi', '',...
                     '',  'Bravo', '', '', '', '', '', '', ''}';

            elseif ismember(file_code(k), {'af'})
                contin{k} = 'Africa';
                
                a = {'',      'E. Africa', 'Shebelle', 'SE. Africa',  '', 'Zambezi', 'Sabi-Buzi', 'Limpopo', 'S. Africa', '',...
                    'Orange',    '',     'SW. Africa', 'Congo', 'WC. Africa', '', 'Niger',  '',  'Volta', '',...
                    'Senegal',  '',   '',      '',       '', '', 'Nile',   '',  '', '',...
                    '',  'Madagascar', '', '', '', '', '', '', '', 'S. Lake Chad',...
                    'Okavango', 'N. Lake Chad', '', '', 'Lake Turkana-Awash', '', '', '', ''}';

            elseif ismember(file_code(k), {'as'})
                contin{k} = 'Asia';
                
                % http://www.wepa-db.net/policies/state/china/river.htm
                a = {'', 'Amur', 'SE. Russia Coast', 'Liao River-Korea', 'Haihe', 'Yellow', 'Huai', 'Yangtze',...
                     'SE. China Coast', 'Pearl', 'Vietnam', 'Mekong', 'Chao Phraya', 'S. Indochina',...
                     'Irrawaddy', '', 'Ganges', 'S. India', '', 'NW. India',...
                     'Indus', '', '', 'Hokkaido', 'Mainland Japan', '', '', 'Taiwan', '', 'Hainan',...
                     '', 'Sri Lanka', 'Tarim', 'Amu Darya', 'Balkhash', 'Syr Darya',...
                     'Qaidam', 'Gobi3', 'Gobi1','Gobi5', 'Gobi2', 'Gobi4', 'Gobi6', 'S. Inner TP', 'N. Inner TP'}';

            elseif ismember(file_code(k), {'si'})
                contin{k} = 'Siberia';
                
                a = {'', 'Ob', '', '', 'Yenisei', '', 'N. Siberia', '', 'Lena', '',...
                     'E. Siberia', '', '', '', '', '', '', ''}';

            elseif ismember(file_code(k), {'au'})
                contin{k} = 'Australia';
                
                noOfbasin = 47;
                a = repmat(str_noname, noOfbasin, 1); 
                a(8)  = {'Kalimantan'};
                a(17) = {'New Guinea'};
                a(18) = {'N. Australia'};
                
            elseif ismember(file_code(k), {'eu'})
                contin{k} = 'Europe';
                noOfbasin = 50;
%                 a = repmat(str_noname, noOfbasin, 1); 
                a = cellstr(strcat('eu', string(1:noOfbasin)))';
                a(38) = {'CapsianSea East Coast'};
                a(42) = {'Persian Gulf Coast'};
                a(46) = {'Rub al Khali'};
                a(47) = {'Helmand'};
                a(48) = {'Central Iran'};
                
            elseif ismember(file_code(k), {'ar'})
                contin{k} = 'Arctic';
                noOfbasin = 19;
                a = repmat(str_noname, noOfbasin, 1); 
                
            elseif ismember(file_code(k), {'gr'})
                contin{k} = 'Greenland';
                noOfbasin = 4;
                a = repmat(str_noname, noOfbasin, 1); 
            else
                error('Wrong input of file_code %s for a given level %d!', file_code{k}, level);
            end

        elseif level == 4
            if ismember(file_code(k), {'na_mex'})
                %/ ten items in each row
                a = {'WSM', 'WCM', 'WNM'};
            end
        end

        %/ Name the empty cell with a temporary code name 
        ind = find(cellfun(@isempty, a));
        for i = 1:length(ind)
            a{ind(i)} = char(strcat(file_code(k), num2str(ind(i))));  %/ char() -> write char (instead of a string) into cell.
                                                                      %/ it will avoid bug in ismember() in reg_extractor.m  
        end
        
%         a(ismember(a, {''})) = {'na'};   %/ convert '' to 'na' for the convenience of using functions like cellstr.
%         a(ismember(a, {''})) = str_noname;   %/ convert '' to 'na' for the convenience of using functions like cellstr.
%         basin_name(1:length(a)) = a;
        basin_name{k,1} = a;
        basin_name{k,2} = repmat(contin(k), length(a), 1);
    end
    
    %/ Col 1: basin name; Col 2: the continent of the basin 
    basin_name = [cat(1, basin_name{:,1}), cat(1, basin_name{:,2})]; %/ concatenate all nested cells.
    
end
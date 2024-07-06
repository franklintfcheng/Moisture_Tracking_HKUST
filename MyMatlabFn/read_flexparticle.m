function [partoutput, numpart_esti] = read_flexparticle(varargin)
    
    pnames = {'outputfolder', 'date_flag', 'numpart_esti'}; 
    dflts  = {[], [], 'unknown'};
    %/ parse function arguments
    [outputfolder, date_flag, numpart_esti] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:});

    if ismember(numpart_esti, 'unknown')             
        %/ If numpart is unknown, then input any string (e.g., 'unknown'), it will estimate it.
        numpart_esti = readpart10(outputfolder, date_flag, numpart_esti);
    end

    partoutput = readpart10(outputfolder, date_flag, numpart_esti);
end


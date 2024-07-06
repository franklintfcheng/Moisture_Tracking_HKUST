function fid = Cell2CsvTxt(data2D,filename1)
    fid = fopen(filename1,'wt');
    if fid>0
        [n1, n2]= size(data2D);

        for k=1:n1
            if isnumeric(data2D{k,1}) %Assume each row has only one type of data.
                ty = '%f,';
            else
                if ismember(filename1(end-3:end),{'.txt'})
                    ty = '%s ';
                else
                    ty = '%s,';
                end
            end
            fmt = repmat(ty,1,n2);
            fmt = fmt(1:end-1); %remove the last char that may cause an additional unwanted element.
            
            fprintf(fid, [fmt,'\n'], data2D{k,:});
        end
        fclose(fid);
    end
end     

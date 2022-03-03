function [] = replaceLineInFile(filename,num, str)
    fid = fopen(filename,'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);
    
    % Change cell A
    for i=1:length(num)
        A{num(i)} = convertStringsToChars(str(i));
    end
    
    % Write cell A into txt
    fid = fopen(filename, 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
end


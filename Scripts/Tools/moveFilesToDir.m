function moveFilesToDir(dirName, dateStamp)
    % This function create dirName in D:\Results,
    % and then move all the directories in D:\Results with the same
    % time stamp (date) to the new directory.
    
    dirPath = 'D:\Results\';
    newDirPath = [dirPath, dirName];
    mkdir([dirPath, dirName]);
    files = dir(dirPath);

    for id = 1:length(files)
        % Get the file name (minus the extension)
        [~, name] = fileparts(files(id).name);
        if length(name) < 11
            continue;
        end
        timeStamp = name(1:11);

        if files(id).isdir && strcmp(timeStamp, dateStamp)
            newName = name(22:end);
            movefile([dirPath, name], [newDirPath,'\', newName]);
        end
    end

end
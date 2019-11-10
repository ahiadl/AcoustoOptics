function moveFilesToDir(dirName, dateStamp)

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
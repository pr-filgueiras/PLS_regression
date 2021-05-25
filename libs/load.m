function [] = load(path)
    files = dir(path);
    dirFlags = [files.isdir];
    subFolders = files(dirFlags);
    
    for k = 1 : length(subFolders)
        dirName = subFolders(k).name;
        if strcmp(dirName, '.') || strcmp(dirName, '..')
            continue;
        end
      addpath(strcat(path, dirName));
      disp(dirName);
    end
end


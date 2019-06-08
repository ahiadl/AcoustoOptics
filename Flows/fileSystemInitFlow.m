function [measDir, rawDataDir, figsDir] = fileSystemInitFlow(general)

dateStr = strrep(['Results-',datestr(datetime('now'))], ':', '-');
measDir = [general.path, '\', dateStr,'-',general.measName];
rawDataDir = [measDir, '\rawData'];
figsDir = [measDir, '\Figures'];
mkdir(measDir);
mkdir(rawDataDir);
mkdir(figsDir);

end


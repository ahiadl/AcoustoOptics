
dllStr = 'C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\Other\VS\MathLibrary\x64\Debug\MathLibrary.dll';
libStr = 'C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\Other\VS\MathLibrary\x64\Debug\MathLibrary.lib';
hdrStr = 'C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\Other\VS\MathLibrary\MathLibrary\MathLibrary.h';

loadlibrary(dllStr, hdrStr)
% loadlibrary(dllStr, @testHeader)
calllib('MathLibrary', 'fibonacci_init', 1, 1);

arg  = calllib('MathLibrary', 'fibonacci_current')
arg2 = calllib('MathLibrary', 'fibonacci_next')

fibonacci_current
fibonacci_index
%%
addpath(fullfile(matlabroot,'extern','examples','shrlib'))

if not(libisloaded('shrlibsample'))
    loadlibrary('shrlibsample')
end
libfunctions('shrlibsample')
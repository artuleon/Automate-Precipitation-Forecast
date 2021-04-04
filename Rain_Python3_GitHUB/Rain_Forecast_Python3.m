home_dir = pwd; %curent folder
tic;
%string_temp4 = ['python.exe '  home_dir '\Precipitation_Map.py']
string_temp4 = ['python.exe '  home_dir '\Forecast_GFS.py']


[status,cmdout] = system(string_temp4)
toc;

tic;
path_asc2dssGrid_executable =  [ home_dir '\Forecast_GFS\'];
cd(path_asc2dssGrid_executable);
string_temp = [path_asc2dssGrid_executable 'ASCIIToDSS.bat'];
[status,cmdout] = system(string_temp);
disp('time for converting to DSS')
toc;












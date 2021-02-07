close all;
if ~isfolder('pastJobs')
  mkdir('pastJobs');
end
addpath('./pastJobs');
addpath('./library');

courantAlert;
loadApp;

currentFolder = pwd;
startSimulationDatetime = datetime;
startSimulationDatetime.Format = 'yyyy-MM-dd''T''HHmmss';
datetimePath = char(startSimulationDatetime);
startSimulationDatetime.Format = 'uuuu/MM/dd HH:mm:ss';
startSimulationDatetime = char(startSimulationDatetime);
cd('./pastJobs');
mkdir(datetimePath);
newDirAbsolutePath = [currentFolder '/pastJobs/' datetimePath];
addpath(newDirAbsolutePath);
cd(newDirAbsolutePath);
mkdir('figures');
figPath = [newDirAbsolutePath '/figures'];
addpath(figPath);
cd(currentFolder);

inputParamMatFilename = 'inputParam.mat';
save(inputParamMatFilename, 'inputParam', '-v7.3');
movefile(inputParamMatFilename, newDirAbsolutePath);

disp('Simulation started'); 
addpath('./colormap');
pltColor.map = colormapTurbo;
pltColor.mapEJ = colormapRdYlBu;
EBstring = {'Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz'};
EBtex = {'$E_x$', '$E_y$', '$E_z$', '$B_x$', '$B_y$', '$B_z$'};
EBtexCat = {'E_{x}', 'E_{y}', 'E_{z}', 'B_{x}', 'B_{y}', 'B_{z}'};
% GUIから実行する場合は次の行をコメントアウト
% addpath('./params/'); [filename,path] = uigetfile('./params/*.mat'); load(filename);
figNum = 1;
figObjNum = floor(ntime/ndskip/8);
f_xyObj = gobjects(1, ndskip*8);
% ndskip = 8;
nkmax = 50;
endTime = startTime + ntime;

if(jobnumber == 1) 
  renormalization;
  initialization;
  position;
  if inputParam.poisson
    charge_with_c;
    potential;
  end
  energy_with_c;
end

openVideos;

num_v = 40;
dv = cv/num_v;
div = 2*[1:num_v+1]-1;
veloDistAxis.para = [-1:dv/cv:1];
veloDistAxis.perp = [0:dv/cv:1];
veloDistAxis.paraLabel = '$v_{\parallel}c^{-1}$';
veloDistAxis.perpLabel = '$v_{\perp}c^{-1}$';
global tmp_editedHist

parameterFileForContiniusJob = ['_jobnum' num2str(jobnumber) '_' datetimePath '.mat'];
% paramEJ.spName = {'All Species', 'Species 1', 'Species 2', 'Species 3', 'Species 4'};
% paramEJ.direction = {'Parallel', 'Perpendicular', 'All directions'};

divide_k = 2;

if inputParam.dispersionEx
  kxkytEx = zeros(nx/divide_k+1, ny*2/divide_k, ntime);
end
if inputParam.dispersionEy
  kxkytEy = zeros(nx/divide_k+1, ny*2/divide_k, ntime);
end
if inputParam.dispersionEz
  kxkytEz = zeros(nx/divide_k+1, ny*2/divide_k, ntime);
end
if inputParam.dispersionBx
  kxkytBx = zeros(nx/divide_k+1, ny*2/divide_k, ntime);
end
if inputParam.dispersionBy
  kxkytBy = zeros(nx/divide_k+1, ny*2/divide_k, ntime);
end
if inputParam.dispersionBz
  kxkytBz = zeros(nx/divide_k+1, ny*2/divide_k, ntime);
end

enableCAccel = true;
tic;
for itime = 1:ntime
  jtime = jtime +1;
  bfield;
  rvelocity_with_c;
  position;
  current_with_c;
  if inputParam.check_EJ && mod(jtime, ndskip) == 0
    % calcEJ_with_c(inputParam, nxp1, nxp2, nyp1, nyp2, np, ex, ey, ez, bx0, by0, bz0,  x, y, q, vx, vy, vz, jtime, X2, Y2, pltColor, v_EJ, figPath);
    calcEJ(inputParam, nxp1, nxp2, nyp1, nyp2, np, ex, ey, ez, bx0, by0, bz0,  x, y, q, vx, vy, vz, jtime, X2, Y2, pltColor, v_EJ, figPath);
  end
%   if inputParam.check_BJ
%     calcBJ;
%   end
  position;
  bfield;
  efield;
  if inputParam.poisson
    charge_with_c;
    potential;
  end
  diagnostics;
  timeDisp = sprintf('Time = %10.3f / %10.3f', jtime*dt, ntime*dt);
  disp(timeDisp);
end 
toc;
jobnumber = jobnumber + 1;
startTime = endTime;

closeVideos;
moveVideos;


clear v_velocitydist  ignoreAxis;
save(parameterFileForContiniusJob, '-v7.3');
movefile(parameterFileForContiniusJob, newDirAbsolutePath);

% if check.wkxky
%   disp('Culculating For Dispersion Plot');
%   dispersionPlot;
% end
dispersionFigNum = 11;
if inputParam.dispersionEx
dispersionFigNum = dispersionPlot(kxkytEx, cv, wc, newDirAbsolutePath, inputParam, pltColor, EBtex, EBtexCat, startSimulationDatetime, 1, dispersionFigNum);
end
if inputParam.dispersionEy
dispersionFigNum = dispersionPlot(kxkytEy, cv, wc, newDirAbsolutePath, inputParam, pltColor, EBtex, EBtexCat, startSimulationDatetime, 2, dispersionFigNum);
end
if inputParam.dispersionEz
dispersionFigNum = dispersionPlot(kxkytEz, cv, wc, newDirAbsolutePath, inputParam, pltColor, EBtex, EBtexCat, startSimulationDatetime, 3, dispersionFigNum);
end
if inputParam.dispersionBx;
dispersionFigNum = dispersionPlot(kxkytBx, cv, wc, newDirAbsolutePath, inputParam, pltColor, EBtex, EBtexCat, startSimulationDatetime, 4, dispersionFigNum);
end
if inputParam.dispersionBy
dispersionFigNum = dispersionPlot(kxkytBy, cv, wc, newDirAbsolutePath, inputParam, pltColor, EBtex, EBtexCat, startSimulationDatetime, 5, dispersionFigNum);
end
if inputParam.dispersionBz
dispersionFigNum = dispersionPlot(kxkytBz, cv, wc, newDirAbsolutePath, inputParam, pltColor, EBtex, EBtexCat, startSimulationDatetime, 6, dispersionFigNum);
end

disp('Successfully Finished. Files are stored to:');
disp(newDirAbsolutePath);
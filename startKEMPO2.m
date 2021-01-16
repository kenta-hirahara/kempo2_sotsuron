close all; clc;
addpath('./pastJobs');

courantAlert;
loadApp;

currentFolder = pwd;
startSimulationDatetime = char(datetime('now','Format','yyyy-MM-dd''T''HHmmss'));
cd('./pastJobs');
mkdir(startSimulationDatetime);
newDirAbsolutePath = [currentFolder '/pastJobs/' startSimulationDatetime];
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
% GUIから実行する場合はこれ以降5行をコメントアウト
% addpath('./params/');
% [filename,path] = uigetfile('./params/*.mat');
% load(filename);
% fileWithoutDotM = strtok(filename, '.');
% fileOnlyAlphabet = strtok(filename, '_');

if(jobnumber == 1) 
  renormalization;
  initialization;
  position; 
  charge;
  potential;
  energy;
end  

nkmax = 50; 
endTime = startTime + ntime;

openVideos;

num_v = 40;
dv = cv/num_v;
div = 2*[1:num_v+1]-1;
veloDistAxis.para = [-1:dv/cv:1];
veloDistAxis.perp = [0:dv/cv:1];
veloDistAxis.paraLabel = texlabel('v_{para}*c^(-1)');
veloDistAxis.perpLabel = texlabel('v_{perp}*c^(-1)');
global tmp_editedHist

parameterFileForContiniusJob = ['_jobnum' num2str(jobnumber) '_' startSimulationDatetime '.mat'];
paramEJ.spName = {'All Species', 'Species 1', 'Species 2', 'Species 3', 'Species 4'};
paramEJ.direction = {'Parallel', 'Perpendicular', 'All directions'};

divide_k = 2;
kxkyt = zeros(nx/divide_k+1, ny*2/divide_k, ntime);

tic;  
for itime = 1:ntime
  jtime = jtime +1;
  bfield;
  rvelocity;
  position;
  current;
  if check.EJ
    calcEJ;
  end
  position;
  bfield;
  efield;
  charge;
  potential;
  diagnostics;
  timeDisp = sprintf('Time = %10.3f / %10.3f', jtime*dt, ntime*dt);
  disp(timeDisp);
end 
kxkytMatFilename = ['kxkyt_', cell2mat(EBstring(EB.number)), '.mat';
save(kxkytMatFilename, 'kxkyt', '-v7.3');
movefile(kxkytMatFilename, newDirAbsolutePath);

jobnumber = jobnumber + 1;
startTime = endTime;

closeVideos;
moveVideos;
toc;

clear app event im fig ax v_velocitydist f_velocitydist ignoreAxis;
save(parameterFileForContiniusJob, '-v7.3');
movefile(parameterFileForContiniusJob, newDirAbsolutePath);

if check.wkxky
  disp('Culculating For Dispersion Plot');
  dispersionPlot;
end

disp('Successfully Finished. Files are stored to:');
disp(newDirAbsolutePath);
disp('Continuous job started');

close all;
currentFolder = pwd;

EJplot = 0;
nkmax = 150; 
surfORimagesc = 2;

endTime = startTime + ntime;

openVideos;

num_v = 40;
dv = cv/num_v;
div = 2*[1:num_v+1]-1;
para = -1*cv:dv:cv;
perp = 0:dv:cv;
global zmax
zmax = zeros(1,ns);
% skipNumber = 4;
asterisk = '*';
cd(currentFolder);

tic;  
for itime = 1:ntime
  jtime = jtime +1;
  bfield;
  rvelocity;
  position;
  current;
  position;
  bfield;
  efield;
  charge;
  potential;
  diagnostics;
  if mod(itime, floor(ntime/10))==0
    disp(asterisk)
    asterisk = [asterisk, '*'];
    if size(asterisk, 2)==6
      asterisk = '*';
    end
  end
end

jobnumber = jobnumber + 1;
startTime = endTime;

closeVideos;
toc;
cd(currentFolder);
clear app event im fig ax kxkyt;
save(['jobnum' num2str(jobnumber) 'tori.mat']);
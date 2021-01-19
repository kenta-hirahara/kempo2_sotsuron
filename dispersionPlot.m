% 以降で時間方向にfft
kxkyw = abs(fft(kxkyt, ntime, 3)/ntime*2);
clear kxkyt;
kxkyw = log10(kxkyw);
kxkyw = fftshift(kxkyw, 2);
kxkyw = fftshift(kxkyw, 3);
kxkywMatFilename = 'kxkyw.mat';
save(kxkywMatFilename, 'kxkyw', '-v7.3');
movefile(kxkywMatFilename, newDirAbsolutePath);

wmax = 0.5*pi;  
kmax = pi;
dkx = 2*kmax/nx;
dky = 2*kmax/ny;
dw = 2*wmax/ntime;
w_axis = [-size(kxkyw, 3)/2+1:size(kxkyw, 3)/2] * dw / abs(wc);
kx_axis = [0:size(kxkyw, 1)-1] * dkx * cv / abs(wc);
ky_axis = [-size(kxkyw, 2)/2+1:size(kxkyw, 2)/2] * dky * cv / abs(wc);
dispersionMatFilename = 'dispersion_params.mat';
save(dispersionMatFilename,'wmax','kmax','dkx','dky',...
'dw','w_axis','kx_axis','ky_axis');
movefile(dispersionMatFilename, newDirAbsolutePath);
fig = figure(100);
fig.Name = 'Dispersion Relation';
fig.Position = [0, 100, 1000, 350];

ax(1) = subplot(1,2,1);
im = imagesc(kx_axis, w_axis, squeeze(kxkyw(:,size(kxkyw, 2)/2,:))');
colormap(pltColor.map); c(1) = colorbar; shading flat;
caxis([-4, -1]);
c(1).Label.Interpreter = 'latex';
c(1).Label.String = ['$\log{10}|' cell2mat(EBstring(EB.number)) '|$'];
ax(1).YDir='normal';

ax(1).Title.String = EBstring(EB.number);
ax(1).XLabel.Interpreter = 'latex';
ax(1).XLabel.FontSize = inputParam.Fontsize;
ax(1).XLabel.String = '$k_{x}$';
ax(1).YLabel.Interpreter = 'latex';
ax(1).YLabel.FontSize = inputParam.Fontsize;
ax(1).YLabel.String = '$\omega$';
ax(1).XLim = [0,20];
ax(1).YLim = [0,20];

ax(2) = subplot(1,2,2); im = imagesc(ky_axis, w_axis, squeeze(kxkyw(1,:,:))');
c(2) = colorbar; shading flat;
colormap(pltColor.map);
caxis([-4, -1]);
c(2).Label.Interpreter = 'latex';
c(2).Label.String = ['$\log{10}|' cell2mat(EBstring(EB.number)) '|$'];
ax(2).YDir='normal';
ax(2).Title.String = EBstring(EB.number);
ax(2).XLabel.Interpreter = 'latex';
ax(2).XLabel.String = '$k_{y}$';
ax(2).XLabel.FontSize = inputParam.Fontsize;
ax(2).YLabel.Interpreter = 'latex';
ax(2).YLabel.String = '$\omega$';
ax(2).YLabel.FontSize = inputParam.Fontsize;
ax(2).XLim = [0,20];
ax(2).YLim = [0,20];

dispersionParallelFigName = strcat('dispersionParallel_', fileWithoutDotM, '.fig');
savefig(fig, dispersionParallelFigName);
movefile(dispersionParallelFigName, newDirAbsolutePath);
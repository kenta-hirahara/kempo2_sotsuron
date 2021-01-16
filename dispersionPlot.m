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

fig  = figure(100);
fig.Name = 'Dispersion Relation';
fig.Position = [0, 100, 1000, 350];
ax(1) = subplot(1,2,1);
im = imagesc(kx_axis, w_axis, squeeze(kxkyw(:,ny/2,:))');
colormap(pltColor.map); colorbar; shading flat;
caxis([-5, -2]);

ax(1).YDir='normal';
ax(1).Title.String = EBstring(EB.number);
ax(1).XLabel.String = 'kx';
ax(1).YLabel.String = '\omega';
ax(1).XLim = [0,30];
ax(1).YLim = [0,30];

ax(2) = subplot(1,2,2); im = imagesc(ky_axis, w_axis, squeeze(kxkyw(1,:,:))');
colorbar; shading flat;
colormap(pltColor.map);
caxis([-5, -2]);

ax(2).YDir='normal';
ax(2).Title.String = EBstring(EB.number);
ax(2).XLabel.String = 'ky';
ax(2).YLabel.String = '\omega';
ax(2).XLim = [0,30];
ax(2).YLim = [0,30];

dispersionParallelFigName = strcat('dispersionParallel_', fileWithoutDotM, '.fig');
savefig(fig, dispersionParallelFigName);
movefile(dispersionParallelFigName, newDirAbsolutePath);
% % 斜め伝搬ここから
% krw = zeros(21, size(kxkyw, 3));
% [kyDispersion, kxDispersion] = meshgrid(ky_axis, kx_axis);
thetaDispersion = 10;
phi_kSpace = atan(dkx*tan(deg2rad(thetaDispersion))/dky);

nkx = 1;
nky = 1;
kx_div_ky = nkx/nky; % int/int
ky_div_kx = nky/nkx; % int/int
krw = zeros(size(kxkyw, 1), size(kxkyw, 3)/2+1); %wの正の部分のみ
for w=1:size(kxkyw, 3)/2+1
  for i=1:size(kxkyw, 1)
    krw(i, w) = kxkyw(i, nx/divide_k-1+i, w+ntime/2-1);
  end
end

fig  = figure(10);
fig.Name = 'Oblique Dispersion Relation';
fig.Position = [0, 100, 800, 600];
ax = axes();
im_krw = imagesc(krw');
colormap(pltColor.map); colorbar; shading flat;
ax.XLabel.String = 'k(mode)'; ax.YLabel.String = '\omega';
ax.YDir = 'normal';
ax.Title.String = sprintf('%s \n %3.3f degree oblique mode', cell2mat(EBstring(EB.number)), rad2deg(atan(kx_div_ky)));
% xDispersion = [0:dkx:20]*cos(deg2rad(thetaDispersion));
% yDispersion = [0:dkx:20]*sin(deg2rad(thetaDispersion));
% for i=1:size(kxkyw, 3)
%   krw(:, i) = interp2(kxDispersion, kyDispersion, kxkyw(:, :, i), xDispersion, yDispersion);
% end
% fig  = figure(101);
% fig.Name = 'Dispersion Relation oblique';
% fig.Position = [0, 100, 1000, 350];

% im = imagesc(krw);
% colormap(pltColor.map);
% colorbar; shading flat;
% caxis([-5, -2]);
% % 斜め伝搬ここまで
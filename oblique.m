function  figNum = oblique(nkx, nky, kxkyw, kx_axis, w_axis, inputParam, pltColor, EB, figNum)
  kx_div_ky = nkx/abs(nky); % int/int
  ky_div_kx = nky/nkx; % int/int
  EBstring = {'Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz'};
  % krw = zeros(size(kxkyw, 1), size(kxkyw, 3)/2+1); %wの正の部分のみ
  numOfPoints = min(ceil(size(kxkyw,1)/nkx), ceil(size(kxkyw,1)/abs(nky)));
  krw = zeros(numOfPoints, size(kxkyw, 3)); %wの正の部分のみ
  for w=1:size(kxkyw, 3)
    for i=1:numOfPoints
      krw(i, w) = kxkyw(1+(i-1)*nkx, size(kxkyw,1)+(i-1)*nky-1, w);
    end
  end
  
  k_norm = sqrt(nkx^2+nky^2)/nkx;
  fig  = figure(figNum);
  fig.Name = 'Oblique Dispersion Relation';
  fig.Position = [0, 100, 520, 400];
  ax = axes();
  im_krw = imagesc(kx_axis(1:nkx:numOfPoints*nkx)*k_norm, w_axis, krw');

  hold on
  w = [0:0.001:20];
  ksq = (w.^2+(w*inputParam.wp(1))./(1-w))/(inputParam.cv^2);
  p = plot(sqrt(ksq)*inputParam.cv,w);
  hold off

  colormap(pltColor.map); c = colorbar; shading flat;
  c.Label.Interpreter = 'latex';
  c.Label.FontSize = inputParam.Fontsize;
  c.Label.String = ['$\log_{10}|' cell2mat(EBstring(EB.number)) '|$'];
  c.LineWidth = 2;

  ax.LineWidth = 2;
  ax.XLabel.Interpreter = 'latex';
  ax.XLabel.String = '$kc\Omega_e^{-1}$';
  ax.XLabel.FontSize = inputParam.Fontsize;
  ax.YLabel.Interpreter = 'latex';
  ax.YLabel.String = '$\omega\Omega_e^{-1}$';
  ax.YLabel.FontSize = inputParam.Fontsize;
  ax.YDir = 'normal';
  ax.PlotBoxAspectRatio = [100, 100, 1];
  ax.Title.FontSize = inputParam.Fontsize;
  ax.Title.String = sprintf('%3.3f degrees', rad2deg(atan(ky_div_kx))-inputParam.phi);
  ax.XLim = [0,10];
  ax.YLim = [0,10];

  p.LineWidth = 1.4;
  figNum = figNum + 1;
end
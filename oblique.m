function  oblique(nkx, nky, kxkyw, kx_axis, w_axis, inputParam, pltColor, EBstring, EB)
  nkx = 1;
  nky = 1;
  kx_div_ky = nkx/nky; % int/int
  ky_div_kx = nky/nkx; % int/int
  % krw = zeros(size(kxkyw, 1), size(kxkyw, 3)/2+1); %wの正の部分のみ
  numOfPoints = min(ceil(size(kxkyw,1)/nkx), ceil(size(kxkyw,1)/nky));
  krw = zeros(numOfPoints, size(kxkyw, 3)); %wの正の部分のみ
  for w=1:size(kxkyw, 3)
    for i=1:numOfPoints
      krw(i, w) = kxkyw(1+(i-1)*nkx, size(kxkyw,1)+(i-1)*nky-1, w);
    end
  end
  
  fig  = figure(10);
  fig.Name = 'Oblique Dispersion Relation';
  fig.Position = [0, 100, 800, 600];
  ax = axes();
  im_krw = imagesc(kx_axis(1:max(nkx,nky):end)*sqrt(nkx^2+nky^2), w_axis, krw');
  colormap(pltColor.map); colorbar; shading flat;
  ax.XLabel.Interpreter = 'latex';
  ax.XLabel.String = '$k$';
  ax.XLabel.FontSize = inputParam.Fontsize;
  ax.YLabel.Interpreter = 'latex';
  ax.YLabel.String = '$\omega$';
  ax.YDir = 'normal';
  ax.Title.String = sprintf('%s \n %3.3f degree oblique mode', cell2mat(EBstring(EB.number)), rad2deg(atan(kx_div_ky)));
end
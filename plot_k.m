function  plot_k(xyEB, k, prm, nkmax, pltColor, EBstring, jtime)
  Z = zeros(prm.nx/2, prm.ny/2);
  Y = abs(fft2(cell2mat(xyEB(k)), prm.nx,prm.ny)) / (prm.nx*prm.ny) * 4; % Yはnx行, ny列の行列
  i = 2:prm.nx/2; j = 2:prm.ny/2;
  Z(i,j)= Y(i,j) + Y(prm.nx-i+2, prm.ny-j+2);
  Z(1,j) = Y(1,j) + Y(1,prm.ny-j+2);
  Z(i,1) = Y(i,1) + Y(prm.nx-i+2,1);
  Z(1,1) = Y(1,1); 
  KK = 1:nkmax+1;
  ax(k) = subplot(2, 3, k); sc = imagesc(Z(KK,KK)');
  colormap(pltColor.map); colorbar; shading flat;
  switch k
  case {1,2}
    caxis([0, 0.4]);
  case {3}
    caxis([0, 0.3]);
  case {4,5,6}
    caxis([0, 8e-3]);
  end
  ax(k).DataAspectRatio = [100, 100, 1];
  sc(1).XData = 0:nkmax; sc(1).YData = 0:nkmax;
  ax(k).XLim = [0, 20]; ax(k).YLim = [0, 20]; 
  ax(k).YDir = 'normal';
  ax(k).XLabel.String = 'kx(mode)'; ax(k).YLabel.String = 'ky(mode)';
  % ax(k).ZLabel.String = cell2mat(EBstring(k));
  ax(k).Title.String = sprintf('%s \n Time = %10.3f / %10.3f', ...
  cell2mat(EBstring(k)), jtime*prm.dt, prm.ntime*prm.dt);
end
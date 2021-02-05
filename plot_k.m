function  plot_k(xyEB, k, inputParam, nkmax, pltColor, jtime)
  Z = zeros(inputParam.nx/2, inputParam.ny/2);
  Y = abs(fft2(cell2mat(xyEB(k)), inputParam.nx,inputParam.ny)) / (inputParam.nx*inputParam.ny) * 4; % Yはnx行, ny列の行列
  i = 2:inputParam.nx/2; j = 2:inputParam.ny/2;
  Z(i,j)= Y(i,j) + Y(inputParam.nx-i+2, inputParam.ny-j+2);
  Z(1,j) = Y(1,j) + Y(1,inputParam.ny-j+2);
  Z(i,1) = Y(i,1) + Y(inputParam.nx-i+2,1);
  Z(1,1) = Y(1,1); 
  KK = 1:nkmax+1;
  ax(k) = subplot(2, 3, k); sc = imagesc(Z(KK,KK)');
  colormap(pltColor.map); 
  c(k) = colorbar; shading flat;
  switch k
  case {1,2}
    caxis([0, 0.4]);
  case {3}
    caxis([0, 0.3]);
  case {4,5,6}
    caxis([0, 8e-3]);
  end
  forTitle = {'E_x', 'E_y', 'E_z', 'B_x', 'B_y', 'B_z'};

  ax(k).Title.Interpreter = 'latex';
  ax(k).FontSize = inputParam.Fontsize;
  ax(k).DataAspectRatio = [100, 100, 1];
  sc(1).XData = 0:nkmax; sc(1).YData = 0:nkmax;
  ax(k).XLim = [0, 20]; ax(k).YLim = [0, 20]; 
  ax(k).YDir = 'normal';
  ax(k).XLabel.Interpreter = 'latex';
  ax(k).XLabel.FontSize = inputParam.Fontsize;
  ax(k).XLabel.String = '$k_{x}$';
  ax(k).YLabel.Interpreter = 'latex';
  ax(k).YLabel.FontSize = inputParam.Fontsize;
  ax(k).YLabel.String = '$k_{y}$';
  ax(k).Title.String = sprintf('$%s$ \n Time = %10.3f / %10.3f', ...
  cell2mat(forTitle(k)), jtime*inputParam.dt, inputParam.ntime*inputParam.dt);
  ax(k).LineWidth = 2;
  c(k).LineWidth = 2;
end
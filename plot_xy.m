function plot_xy(xyEB, pltColor, k, EBstring, jtime, prm)
  ax(k) = subplot(2, 3, k);imagesc(cell2mat(xyEB(k))');
  colormap(pltColor.map); 
  c(k) = colorbar; shading flat;
  switch k 
  case {1,2}
    caxis([-0.6, 0.6]);
  case {3}
    caxis([-0.3, 0.3]);
  case {4}
    caxis([-2e-2, 2e-2]);
  case {5,6}
    caxis([-1.2e-2, 1.2e-2]);
  end

  % c(k).Label.Interpreter = 'latex';
  % c(k).Label.FontSize = prm.Fontsize;
  % c(k).Label.String = ['$' cell2mat(EBstring(k)) '$'];
  % c(k).Location = 'northoutside';
  ax(k).DataAspectRatio = [100, 100, 1];
  ax(k).YDir = 'normal';
  ax(k).XLabel.Interpreter = 'latex';
  ax(k).XLabel.FontSize = prm.Fontsize;
  ax(k).XLabel.String = '$x$';
  ax(k).YLabel.Interpreter = 'latex';
  ax(k).YLabel.FontSize = prm.Fontsize;
  ax(k).YLabel.String = '$y$';
  ax(k).Title.String = sprintf('%s \n Time = %10.3f / %10.3f',...
   cell2mat(EBstring(k)), jtime*prm.dt, prm.ntime*prm.dt);
end
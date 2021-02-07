function plot_xy(xyEB, pltColor, k, jtime, inputParam)
  ax(k) = subplot(2, 3, k);imagesc(cell2mat(xyEB(k))');
  colormap(pltColor.map); 
  c(k) = colorbar; shading flat;
  switch k 
  case {1,2,3}
    caxis([-0.6, 0.6]);
  case {4,5,6}
    caxis([-2e-2, 2e-2]);
  end
  
  c(k).LineWidth = 2;
  % c(k).Label.Interpreter = 'latex';
  % c(k).Label.FontSize = inputParam.Fontsize;
  % c(k).Label.String = ['$' cell2mat(EBstring(k)) '$'];
  % c(k).Location = 'northoutside';
  ax(k).FontSize = inputParam.Fontsize;
  ax(k).LineWidth = 2;
  ax(k).DataAspectRatio = [100, 100, 1];
  ax(k).Title.Interpreter = 'latex';
  ax(k).YDir = 'normal';
  ax(k).XLabel.Interpreter = 'latex';
  ax(k).XLabel.FontSize = inputParam.Fontsize;
  ax(k).XLabel.String = '$x$';
  ax(k).YLabel.Interpreter = 'latex';
  ax(k).YLabel.FontSize = inputParam.Fontsize;
  ax(k).YLabel.String = '$y$';
  forTitle = {'E_x', 'E_y', 'E_z', 'B_x', 'B_y', 'B_z'};
  if inputParam.ajamp
    if jtime < inputParam.ctime
      ax(k).Title.String = sprintf('$%s$ \n Time = %10.3f / %10.3f \n $J_{ext}$ direction: %s',...
          cell2mat(forTitle(k)), jtime*inputParam.dt, inputParam.ntime*inputParam.dt, inputParam.directionJ);
    else
      ax(k).Title.String = sprintf('$%s$ \n Time = %10.3f / %10.3f \n $J_{ext}$ %s stopped',...
          cell2mat(forTitle(k)), jtime*inputParam.dt, inputParam.ntime*inputParam.dt, inputParam.directionJ);
    end
  else
  ax(k).Title.String = sprintf('$%s$ \n Time = %10.3f / %10.3f',...
    cell2mat(forTitle(k)), jtime*inputParam.dt, inputParam.ntime*inputParam.dt);
  end
end
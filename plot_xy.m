function plot_xy(xyEB, pltColor, k, EBstring, jtime, prm)
  ax(k) = subplot(2, 3, k);imagesc(cell2mat(xyEB(k))');
  colormap(pltColor.map);  colorbar; shading flat;
  switch k 
  case {1,2}
    caxis([-0.6, 0.6]);
  case {3}
    caxis([-0.3, 0.3]);
  case {4}
    caxis([-1e-1, 1e-1]);
  case {5,6}
    caxis([-8e-3, 8e-3]);
  end

  ax(k).DataAspectRatio = [100, 100, 1];
  ax(k).XLabel.String = 'X'; ax(k).YLabel.String = 'Y'; ax(k).ZLabel.String = cell2mat(EBstring(k));
  ax(k).YDir = 'normal';
  ax(k).Title.String = sprintf('%s \n Time = %10.3f / %10.3f',...
   cell2mat(EBstring(k)), jtime*prm.dt, prm.ntime*prm.dt);
end

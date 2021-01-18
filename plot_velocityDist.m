function plot_velocityDist(editableHistogram, prm, veloDistAxis, pltColor, spName, jtime)
  % contour plot of velocity distribution
  global tmp_editedHist
  for k = 1:prm.ns
    ax(k) = subplot(prm.ns, 2, 2*k-1);
    editedHist = cell2mat(editableHistogram(1, k))';
    imagesc(veloDistAxis.para, veloDistAxis.perp, editedHist);
    colormap(ax(k), pltColor.map); colorbar;
    caxis([0, 1e6]/prm.vPara(k)/prm.vPerp(k));
    ax(k).YDir = 'normal';
    ax(k).XLabel.Interpreter = 'latex';
    ax(k).XLabel.FontSize = prm.Fontsize;
    ax(k).XLabel.String = veloDistAxis.paraLabel;
    ax(k).YLabel.Interpreter = 'latex';
    ax(k).YLabel.FontSize = prm.Fontsize;
    ax(k).YLabel.String = veloDistAxis.perpLabel;
    ax(k).DataAspectRatio = [100, 100, 1];
    ax(k).Title.String = sprintf('%s \n Time = %10.3f / %10.3f', ...
    cell2mat(spName(1+k)), jtime*prm.dt, prm.ntime*prm.dt);
    
    diff_ax(k) = subplot(prm.ns, 2, 2*k);
    imagesc(veloDistAxis.para, veloDistAxis.perp, editedHist-cell2mat(tmp_editedHist(1, k))');
    colormap(diff_ax(k), pltColor.mapEJ); colorbar;
    caxis([-1e4, 1e4]/prm.vPara(k)/prm.vPerp(k));
    diff_ax(k).YDir = 'normal';
    diff_ax(k).XLabel.Interpreter = 'latex';
    diff_ax(k).XLabel.FontSize = prm.Fontsize;
    diff_ax(k).XLabel.String = veloDistAxis.paraLabel;
    diff_ax(k).YLabel.Interpreter = 'latex';
    diff_ax(k).YLabel.FontSize = prm.Fontsize;
    diff_ax(k).YLabel.String = veloDistAxis.perpLabel;
    diff_ax(k).DataAspectRatio = [100, 100, 1];
    diff_ax(k).Title.String = sprintf('diff %s \n Time = %10.3f / %10.3f', ...
    cell2mat(spName(1+k)), jtime*prm.dt, prm.ntime*prm.dt);
  end
end
function plot_velocityDist(editableHistogram, inputParam, veloDistAxis, pltColor, jtime)
  % contour plot of velocity distribution
  global tmp_editedHist
  spName = {'All Species', 'Species 1', 'Species 2', 'Species 3', 'Species 4'};

  thetaForCirc = [0:0.001:pi];
  circ_x = cos(thetaForCirc)*0.2;
  circ_y = sin(thetaForCirc)*0.2;

  for k = 1:inputParam.ns
    ax(k) = subplot(inputParam.ns, 2, 2*k-1);
    editedHist = cell2mat(editableHistogram(1, k))';
    imagesc(veloDistAxis.para, veloDistAxis.perp, editedHist);
    hold on;
    for circIter = 1:4
      p(circIter) = plot(circ_x*circIter, circ_y*circIter);
      p(circIter).LineWidth = 1;
      p(circIter).LineStyle = '-.';     
      p(circIter).Color = [1, 1, 1];
      hold on;
    end
    hold off;
    colormap(ax(k), pltColor.map); c(k) = colorbar;
    caxis([0, 2e5]/inputParam.vPara(k)/(inputParam.vPerp(k)^2));
    c(k).LineWidth = 2;
    ax(k).Title.Interpreter = 'latex';
    ax(k).FontSize = inputParam.Fontsize;
    ax(k).LineWidth = 2;
    % ax(k).Title.Interpreter = 'latex';
    ax(k).YDir = 'normal';
    ax(k).XLabel.Interpreter = 'latex';
    ax(k).XLabel.FontSize = inputParam.Fontsize;
    ax(k).XLabel.String = veloDistAxis.paraLabel;
    ax(k).YLabel.Interpreter = 'latex';
    ax(k).YLabel.FontSize = inputParam.Fontsize;
    ax(k).YLabel.String = veloDistAxis.perpLabel;
    ax(k).DataAspectRatio = [100, 100, 1];
    tmpTitle = sprintf('%s \n Time = %10.3f / %10.3f', ...
    cell2mat(spName(1+k)), jtime*inputParam.dt, inputParam.ntime*inputParam.dt);
    ax(k).Title.String = replace(tmpTitle, 'Time', '$t\Omega_e$');
    diff_ax(k) = subplot(inputParam.ns, 2, 2*k);
    imagesc(veloDistAxis.para, veloDistAxis.perp, editedHist-cell2mat(tmp_editedHist(1, k))');
    hold on;
    for diff_circIter = 1:4
      diff_p(diff_circIter) = plot(circ_x*diff_circIter, circ_y*diff_circIter);
      diff_p(diff_circIter).LineWidth = 1;
      diff_p(diff_circIter).LineStyle = '-.';     
      diff_p(diff_circIter).Color = [0, 0, 0];
      hold on;
    end
    hold off;
    colormap(diff_ax(k), pltColor.mapEJ); diff_c(k) = colorbar;
    caxis([-5e3, 5e3]/inputParam.vPara(k)/inputParam.vPerp(k));
    diff_c(k).LineWidth = 2;
    diff_ax(k).XLabel.Interpreter = 'latex';
    diff_ax(k).LineWidth = 2;
    diff_ax(k).Title.Interpreter = 'latex';
    diff_ax(k).FontSize = inputParam.Fontsize;
    diff_ax(k).YDir = 'normal';
    diff_ax(k).XLabel.Interpreter = 'latex';
    diff_ax(k).XLabel.FontSize = inputParam.Fontsize;
    diff_ax(k).XLabel.String = veloDistAxis.paraLabel;
    diff_ax(k).YLabel.Interpreter = 'latex';
    diff_ax(k).YLabel.FontSize = inputParam.Fontsize;
    diff_ax(k).YLabel.String = veloDistAxis.perpLabel;  
    diff_ax(k).DataAspectRatio = [100, 100, 1];
    diff_tmpTitle = sprintf('diff %s \n Time = %10.3f / %10.3f', ...
    cell2mat(spName(1+k)), jtime*inputParam.dt, inputParam.ntime*inputParam.dt);
    diff_ax(k).Title.String = replace(diff_tmpTitle , 'Time', '$t\Omega_e$');
  end
end
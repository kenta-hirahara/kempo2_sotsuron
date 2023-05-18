function plotOnlyTempAnis(it, inputParam, ndskip, At, ax_TempAnis)
    IT=(1:it);    
    pt = IT*inputParam.dt*ndskip;
    speciesNames = {'$K_{\mathrm{p}{1}}$', '$K_{\mathrm{p}{2}}$', '$K_{\mathrm{p}{3}}$', '$K_{\mathrm{p}{4}}$', '$K_{\mathrm{p}{5}}$', '$K_{\mathrm{p}{6}}$', '$K_{\mathrm{p}{7}}$'};
    

    for j = 1:inputParam.ns
      pltAnis(j) = plot(pt,At(:,j)');
      pltAnis(j).LineWidth = 1.4;
      hold on;
    end
    hold off;
    xlabel('Time')
    lgd_TempAnis = legend(speciesNames(1:inputParam.ns));
    lgd_TempAnis.Box = 'on';
    lgd_TempAnis.Interpreter = 'latex';
    lgd_TempAnis.FontSize = inputParam.Fontsize*0.8;
    lgd_TempAnis.LineWidth = 2;

    ax_TempAnis.FontSize = inputParam.Fontsize*0.8;
    ax_TempAnis.FontWeight = 'bold';
    ax_TempAnis.LineWidth = 2;
    ax_TempAnis.Box = 'on';
    ax_TempAnis.Title.String = 'Temperature Anisotropy';
    ax_TempAnis.Title.FontSize = inputParam.Fontsize;
    ax_TempAnis.XLabel.Interpreter = 'latex';
    ax_TempAnis.XLabel.FontSize = inputParam.Fontsize*1.2;
    ax_TempAnis.XLabel.String = '$t\Omega_e$';
    ax_TempAnis.YLabel.Interpreter = 'latex';
    ax_TempAnis.YLabel.FontSize = inputParam.Fontsize*1.2;
    ax_TempAnis.YLabel.String = '$A$'; 
  end
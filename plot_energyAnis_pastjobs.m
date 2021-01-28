function fig_energyTempAnis = plot_energyAnis_pastjobs(it, inputParam, engt, eepara, eeperp, ebpara, ebperp, ke, At)
  
  energyCatenated = cat(1, engt, eepara, eeperp, ebpara, ebperp, ke');

  fig_energyTempAnis = figure(8);
  fig_energyTempAnis.Name = 'Energy History and Temperature Anisotropy';
  fig_energyTempAnis.Position = [0,0,1200,500];
  frame = getframe(gcf);
  
  ax_energy = subplot(1,2,1);
  IT=(1:it);    
  pt = IT*inputParam.dt*inputParam.ndskip;
  for i = 1:5+inputParam.ns
    s(i) = semilogy(energyCatenated(i,:));
    s(i).LineWidth = 1;
    hold on;
  end
  hold off;
  legendNames = {'$K_\mathrm{T}$', '$K_{\mathrm{E}\parallel}$', '$K_{\mathrm{E}\perp}$', '$K_{\mathrm{B}\parallel}$', '$K_{\mathrm{B}\perp}$'};
  speciesNames = {'$K_{\mathrm{p}{1}}$', '$K_{\mathrm{p}{2}}$', '$K_{\mathrm{p}{3}}$', '$K_{\mathrm{p}{4}}$', '$K_{\mathrm{p}{5}}$', '$K_{\mathrm{p}{6}}$', '$K_{\mathrm{p}{7}}$'};
  legendNames = cat(2, legendNames, speciesNames(1:inputParam.ns));
  lgd_energy = legend(legendNames);
  lgd_energy.Box = 'on';
  lgd_energy.Interpreter = 'latex';
  lgd_energy.Orientation = 'horizontal';
  lgd_energy.Location = 'southeast';
  lgd_energy.FontSize = 14 ;
  lgd_energy.NumColumns = 3;
  lgd_energy.LineWidth = 1;

  ax_TempAnis = subplot(1,2,2);
  for j = 1:inputParam.ns
    pltAnis(j) = plot(pt,At(:,j)');
    pltAnis(j).LineWidth = 1;
    hold on;
  end
  hold off;
  xlabel('Time')
  lgd_TempAnis = legend(speciesNames(1:inputParam.ns));
  lgd_TempAnis.Box = 'on';
  lgd_TempAnis.Interpreter = 'latex';
  lgd_TempAnis.FontSize = 14;
  lgd_TempAnis.LineWidth = 1;

  ax_energy.LineWidth = 1;
  ax_energy.Box = 'on';
  ax_energy.Title.String = 'Energy History';
  ax_energy.Title.FontSize = inputParam.Fontsize*0.8;
  ax_energy.XLabel.Interpreter = 'latex';
  ax_energy.XLabel.FontSize = inputParam.Fontsize;
  ax_energy.XLabel.String = '$t\Omega_e$'; 
  ax_energy.YLabel.Interpreter = 'latex';
  ax_energy.YLabel.FontSize = inputParam.Fontsize;
  ax_energy.YLabel.String = '$K$'; 

  ax_TempAnis.LineWidth = 1;
  ax_TempAnis.Box = 'on';
  ax_TempAnis.Title.String = 'Temperature Anisotropy';
  ax_TempAnis.Title.FontSize = inputParam.Fontsize*0.8;
  ax_TempAnis.XLabel.Interpreter = 'latex';
  ax_TempAnis.XLabel.FontSize = inputParam.Fontsize;
  ax_TempAnis.XLabel.String = '$t\Omega_e$';
  ax_TempAnis.YLabel.Interpreter = 'latex';
  ax_TempAnis.YLabel.FontSize = inputParam.Fontsize;
  ax_TempAnis.YLabel.String = '$A$'; 
end
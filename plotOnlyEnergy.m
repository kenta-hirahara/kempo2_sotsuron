function plotOnlyEnergy(it, inputParam, ndskip, energyCatenated, ax_energy)
  
    % energyCatenated = cat(1, engt, eepara, eeperp, ebpara, ebperp, ke');
  
    IT=(1:it);    
    pt = IT*inputParam.dt*ndskip;
    for i = 1:5+inputParam.ns
      s(i) = semilogy(pt, energyCatenated(i,:));
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
    lgd_energy.FontSize = inputParam.Fontsize*0.8;
    lgd_energy.NumColumns = 3;
    lgd_energy.LineWidth = 2;
  
    ax_energy.FontSize = inputParam.Fontsize*0.8;
    ax_energy.FontWeight = 'bold';
    % ax_energy.TickLabelInterpreter = 'latex'
    ax_energy.LineWidth = 2;
    ax_energy.Box = 'on';
    ax_energy.Title.String = 'Energy';
    ax_energy.Title.FontSize = inputParam.Fontsize;
    ax_energy.XLabel.Interpreter = 'latex';
    ax_energy.XLabel.FontSize = inputParam.Fontsize*1.2;
    ax_energy.XLabel.String = '$t\Omega_e$'; 
    ax_energy.YLabel.Interpreter = 'latex';
    ax_energy.YLabel.FontSize = inputParam.Fontsize*1.2;
    ax_energy.PlotBoxAspectRatio = [100, 100, 1];
    ax_energy.YLabel.String = '$K$'; 
  end
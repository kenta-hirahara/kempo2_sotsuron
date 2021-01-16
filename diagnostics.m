%diagnostics at an interval of dt*ndskip
xyEB = {ex(X2,Y2)*rne, ey(X2,Y2)*rne, ez(X2,Y2)*rne, (bx(X2,Y2)-bx0)*rnb, (by(X2,Y2)-by0)*rnb, (bz(X2,Y2)-bz0)*rnb}; %これは全プロットで使うのでここに固定
if mod(jtime, ndskip) == 0 %ndskipつまり8の倍数の時だけ描画
  energy;

  % figure of xy plot
  if inputParam.xyPlot
    f_xy = figure(1);
    f_xy.Name = 'xy plot';
    f_xy.Position = [0, 0, 1400, 750];
    frames_xy = getframe(f_xy);
    for k =1:6
      plot_xy(xyEB, pltColor, k, EBstring, jtime, inputParam);
    end
    writeVideo(v_xyEorB, frames_xy);
    if mod(jtime, floor(ntime/ndskip/20)) == 0
      figFilename = ['xy', num2str(jtime), '.fig'];
      savefig(f_xy, figFilename);
      movefile(figFilename, figPath);
    end
  end

  % figure of kxky plot
  if inputParam.kxkyPlot
    f_kxky = figure(2);
    f_kxky.Name = 'kxky plot';
    f_kxky.Position = [0, 0, 1400, 750];
    frames_kxky = getframe(f_kxky);
    for k =1:6 %E, Bの全成分プロット
      plot_k(xyEB, k, inputParam, nkmax, pltColor, EBstring, jtime);
    end
    writeVideo(v_kxkyEorB, frames_kxky);
    if mod(jtime, floor(ntime/ndskip/20)) == 0
      figFilename = ['kxky', num2str(jtime), '.fig'];
      savefig(f_kxky, figFilename);
      movefile(figFilename, figPath);
    end
  end

  % figure of velocity distribution plot
  if inputParam.veloDistPlot
    f_velocitydist = figure(3);
    f_velocitydist.Position = [0, 0, 0, 0];
    
    n2 = 0;
    for  k=1:ns
      n1 = n2 + 1;
      n2 = n2 + np(k);

      ignoreAxis(k) = subplot(2,2,k);     
      i = n1:n2;
      h(k) = histogram2(vx(i), sqrt(vy(i).*vy(i) + vz(i).*vz(i)));
      h(k).XBinLimits = [-1*cv, cv];
      h(k).YBinLimits = [0, cv];
      h(k).NumBins = [2*num_v+1, num_v+1];
      editableHistogram(1, k) = {h(k).Values};
      editableHistogram(1, k) = ...
      {cell2mat(editableHistogram(1, k)) ./ (pi*((cv/num_v)^3) * div) * abs(q(k))}; 
    end
    close(f_velocitydist);

    editedVelocityDist = figure(4);
    editedVelocityDist.Name = 'velocity distributions';
    editedVelocityDist.Position = [0, 0, 1400, 900];
    frames_velocitydist = getframe(editedVelocityDist);
    
    if jtime == ndskip
      tmp_editedHist = editableHistogram;
    end
    plot_velocityDist(editableHistogram, inputParam, veloDistAxis, pltColor, paramEJ.spName, jtime);
    
    writeVideo(v_velocitydist, frames_velocitydist);
    if mod(jtime, floor(ntime/ndskip/20)) == 0
      figFilename = ['veloDist', num2str(jtime), '.fig'];
      savefig(editedVelocityDist, figFilename);
      movefile(figFilename, figPath);
    end
  end
end

% figure of w-kx-ky diagram
% if check.wkxky
  % Z = zeros(nx/2, ny/2);
  kxky = fft2(cell2mat(xyEB(app.EBnumber)), nx, ny) / (nx*ny) * 4; 
  %ここでk空間の行列の要素数を一気に減らす、これをnplot(ntimeの約数, とりあえずntime)回して格納した. 
  kxkyt(:, :, itime) = cat(2, kxky(1:nx/divide_k+1, 1:nx/divide_k+1), kxky(1:nx/divide_k+1, end-nx/divide_k+2:end));
% end 

% Diagnostics at the end of the job 
% plotting time history of energies
if itime == ntime
  % fign_xy = fign_xy+1;
  fig = figure(8);
  fig.Name = 'Energy History and Temperature Anisotropy';
  fig.Position = [0,0,1200,500];
  frame = getframe(gcf);
  ax(1) = subplot(1,2,1);
  IT=(1:it);
  pt = IT*dt*ndskip;
  if ns==2
    semilogy(pt, engt, pt, eepara,pt, eeperp,pt,ebpara,pt,ebperp,pt,ke(IT,1),pt,ke(IT,2));
    legend('total','e-para','e-perp','b-para','b-perp','sp1','sp2');
  elseif ns==3
    semilogy(pt, engt, pt, eepara,pt, eeperp,pt,ebpara,pt,ebperp,pt,ke(IT,1),pt,ke(IT,2),pt,ke(IT,3));
    legend('total','e-para','e-perp','b-para','b-perp','sp1','sp2','sp3');
  end
  title('Energy History');
  %plotting time history of temperature anisotropy
  % fign_xy = fign_xy+1;
  ax(2) = subplot(1,2,2);
  plot(pt,At), xlabel('Time'),ylabel('Temperature Anisotropy');
  if ns==2
      legend ('sp1','sp2');
  elseif ns ==3
      legend('sp1','sp2', 'sp3');
  end
  energyAndAnisFigName = strcat('energy_anisotropy_', fileWithoutDotM, '.fig');
  savefig(fig, energyAndAnisFigName);
  movefile(energyAndAnisFigName, newDirAbsolutePath);
  % fign_xy = fign_xy-1;
end; 
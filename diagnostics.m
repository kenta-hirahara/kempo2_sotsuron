%diagnostics at an interval of dt*ndskip
xyEB = {ex(X2,Y2)*rne, ey(X2,Y2)*rne, ez(X2,Y2)*rne, (bx(X2,Y2)-bx0)*rnb, (by(X2,Y2)-by0)*rnb, (bz(X2,Y2)-bz0)*rnb}; %これは全プロットで使うのでここに固定
if mod(jtime, ndskip) == 0 %ndskipつまり8の倍数の時だけ描画
  energy_with_c;

  % figure of xy plot
  if inputParam.check_xyPlot
    f_xy = figure(1);
    f_xy.Name = 'xy plot';
    f_xy.Position = [0, 0, 1400, 750];
    frames_xy = getframe(f_xy);
    for k =1:6
      plot_xy(xyEB, pltColor, k, jtime, inputParam);
    end
    writeVideo(v_xyEorB, frames_xy);
    if mod(jtime, figObjNum) == 0
      % f_xyObj(figObjNum) = f_xy;
      figFilename = ['xy', num2str(jtime), '.fig'];
      savefig(f_xy, figFilename);
      movefile(figFilename, figPath, 'f');
    end
  end

  % figure of kxky plot
  if inputParam.check_kxkyPlot
    f_kxky = figure(2);
    f_kxky.Name = 'kxky plot';
    f_kxky.Position = [0, 0, 1400, 900];
    frames_kxky = getframe(f_kxky);
    for k =1:6 %E, Bの全成分プロット
      plot_k(xyEB, k, inputParam, nkmax, pltColor, jtime);
    end 
    writeVideo(v_kxkyEorB, frames_kxky);
    if mod(jtime, figObjNum) == 0
      figFilename = ['kxky', num2str(jtime), '.fig'];
      savefig(f_kxky, figFilename);
      movefile(figFilename, figPath, 'f');
    end
  end

  % figure of velocity distribution plot
  if inputParam.check_veloDistPlot
    f_velocitydist = figure(3);
    f_velocitydist.Position = [0, 0, 0, 0];
    
    n2 = 0;
    for  k=1:ns
      n1 = n2 + 1;
      n2 = n2 + np(k);

      ignoreAxis(k) = subplot(2,2,k);
      i = n1:n2;
      % 静磁場の方向をparaとする, paraは正負あり
      % v_xyz = [vx(i); vy(i); vz(i)];
      % v_para_vector = b0_xyz'*(b0_xyz*v_xyz)/(b0^b0);
      % v_para_norm = sqrt(sum(v_para_vector.^2, 1));
      % logicalArray = b0_xyz*v_xyz > 0;
      % v_para_direction = v_para_norm .* logicalArray + v_para_norm .* (logicalArray-1);
      % v_perp_norm = sqrt(sum((v_xyz - v_para_vector).^2, 1));

      % h(k) = histogram2(v_para_direction, v_perp_norm);
      h(k) = histogram2(vx(i), sqrt(vy(i).^2+vz(i).^2));

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
    plot_velocityDist(editableHistogram, inputParam, veloDistAxis, pltColor, jtime);
    
    writeVideo(v_velocitydist, frames_velocitydist);
    if mod(jtime, figObjNum) == 0
      figFilename = ['veloDist', num2str(jtime), '.fig'];
      savefig(editedVelocityDist, figFilename);
      movefile(figFilename, figPath, 'f');
    end
  end
end

if inputParam.dispersionEx 
  kxkyEx = fft2(cell2mat(xyEB(1)), nx, ny) / (nx*ny) * 4; 
  kxkytEx(:, :, itime) = cat(2, kxkyEx(1:nx/divide_k+1, 1:nx/divide_k+1), kxkyEx(1:nx/divide_k+1, end-nx/divide_k+2:end));
end
if inputParam.dispersionEy 
  kxkyEy = fft2(cell2mat(xyEB(2)), nx, ny) / (nx*ny) * 4; 
  kxkytEy(:, :, itime) = cat(2, kxkyEy(1:nx/divide_k+1, 1:nx/divide_k+1), kxkyEy(1:nx/divide_k+1, end-nx/divide_k+2:end));
end
if inputParam.dispersionEz 
  kxkyEz = fft2(cell2mat(xyEB(3)), nx, ny) / (nx*ny) * 4; 
  kxkytEz(:, :, itime) = cat(2, kxkyEz(1:nx/divide_k+1, 1:nx/divide_k+1), kxkyEz(1:nx/divide_k+1, end-nx/divide_k+2:end));
end
if inputParam.dispersionBx 
  kxkyBx = fft2(cell2mat(xyEB(4)), nx, ny) / (nx*ny) * 4; 
  kxkytBx(:, :, itime) = cat(2, kxkyBx(1:nx/divide_k+1, 1:nx/divide_k+1), kxkyBx(1:nx/divide_k+1, end-nx/divide_k+2:end));
end
if inputParam.dispersionBy 
  kxkyBy = fft2(cell2mat(xyEB(5)), nx, ny) / (nx*ny) * 4; 
  kxkytBy(:, :, itime) = cat(2, kxkyBy(1:nx/divide_k+1, 1:nx/divide_k+1), kxkyBy(1:nx/divide_k+1, end-nx/divide_k+2:end));
end
if inputParam.dispersionBz 
  kxkyBz = fft2(cell2mat(xyEB(6)), nx, ny) / (nx*ny) * 4; 
  kxkytBz(:, :, itime) = cat(2, kxkyBz(1:nx/divide_k+1, 1:nx/divide_k+1), kxkyBz(1:nx/divide_k+1, end-nx/divide_k+2:end));
end

% plotting time history of energies
if itime == ntime
  energyCatenated = cat(1, engt, eepara, eeperp, ebpara, ebperp, ke');
  % fig_energyTempAnis = plot_energyAnis_pastjobs(it, inputParam, energyCatenated, At);
  fig_energyTempAnis = figure(8);
  fig_energyTempAnis.Name = 'Energy History and Temperature Anisotropy';
  fig_energyTempAnis.Position = [0,0,1200,500];  
  ax_energy = subplot(1,2,1); plotOnlyEnergy(it, inputParam,ndskip, energyCatenated, ax_energy);
  ax_TempAnis = subplot(1,2,2); plotOnlyTempAnis(it, inputParam,ndskip, At, ax_TempAnis);

  energyAnisFilename = 'fig_energyTempAnis.fig';
  savefig(fig_energyTempAnis, energyAnisFilename);
  movefile(energyAnisFilename, newDirAbsolutePath);
  clear fig_energyTempAnis;
end; 
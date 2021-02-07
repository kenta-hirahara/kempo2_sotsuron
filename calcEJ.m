function calcEJ(inputParam, nxp1, nxp2, nyp1, nyp2, np, ex, ey, ez, bx0, by0, bz0,  x, y, q, vx, vy, vz, jtime, X2, Y2, pltColor, v_EJ, figPath)
  b0_xyz = [bx0, by0, bz0];
  % 0. E*Jを描画する行列を初期化(1+ns, 3, nx+2行ny+2列の行列)
  EJplot = zeros(1+inputParam.ns, 3, nxp2, nyp2);
  % 1. E*Jを格納する行列を初期化(3, nx+2行ny+2列の行列)
  EJ = zeros(3, nxp2, nyp2);
  % 2. grid のEを空間方向に平均化
  ExFullGrid = zeros(nxp2, nyp2);
  EyFullGrid = zeros(nxp2, nyp2);
  EzFullGrid = zeros(nxp2, nyp2);
  
  ExFullGrid(2:nxp2, 2:nyp2) = (ex(1:nxp1, 2:nyp2) + ex(2:nxp2, 2:nyp2)) * 0.5;
  EyFullGrid(2:nxp2, 2:nyp2) = (ey(1:nxp1, 2:nyp2) + ey(2:nxp2, 2:nyp2)) * 0.5;
  EzFullGrid(2:nxp2, 2:nyp2) = ...
  (ez(1:nxp1, 1:nyp1) + ez(2:nxp2, 1:nyp1) + ez(1:nxp1, 2:nxp2) + ez(2:nxp2, 2:nxp2)) * 0.25;
  
  % 3. ある超粒子に対するEをsf1, sf2, sf3, sf4を用いて導出
  n2 = 0;
  for  k=1:inputParam.ns
  	n1  = n2+1;
    n2  = n2 + np(k);
    tmp = EJ;
  	for m = n1:n2
      xmf = x(m) + 2.0; %配列xの要素は0からnx(グリッド数)までの数
      ymf = y(m) + 2.0;
      i = floor(xmf); %iは2からnx+2までの数
      j = floor(ymf);
      i1 = i+1;
      j1 = j+1;
      x1 = xmf - i;
      y1 = ymf - j;
      sf3 = x1*y1;
      sf2 = x1-sf3;
      sf4 = y1-sf3;
      sf1 = 1.0-x1-sf4;
      E_particle = zeros(1, 3); %ある超粒子mの位置におけるEを初期化
  
      %ある超粒子mの位置におけるEをxyz3成分算出
      E_particle(1) = ...
      ExFullGrid(i, j)*sf1 + ExFullGrid(i1, j)*sf2 + ExFullGrid(i1, j1)*sf3 + ExFullGrid(i, j1)*sf4;
      E_particle(2) = ...
      EyFullGrid(i, j)*sf1 + EyFullGrid(i1, j)*sf2 + EyFullGrid(i1, j1)*sf3 + EyFullGrid(i, j1)*sf4;
      E_particle(3) = ...
      EzFullGrid(i, j)*sf1 + EzFullGrid(i1, j)*sf2 + EzFullGrid(i1, j1)*sf3 + EzFullGrid(i, j1)*sf4;
  
      % 4. q*v*Eを計算
      EJp = q(k) .* [vx(m), vy(m), vz(m)] .* E_particle;
      % ここ修正の必要あり. vのpara求めて、EのparaとかけてE·Jpara求める. これをE·Jから引いたらE·Jperpもとまる
      % EJp_para = (dot(EJp, b0_xyz) / norm(b0_xyz)) * (b0_xyz/norm(b0_xyz));
      % EJp_perp = EJp - EJp_para;
      % 5. sf1,sf2,sf3,sf4を用いてE*Jを格納する行列へ配分
      
      EJ(:, i ,j ) = EJ(:, i ,j ) + EJp(:)*sf1;
      EJ(:, i1,j ) = EJ(:, i1,j ) + EJp(:)*sf2;
      EJ(:, i1,j1) = EJ(:, i1,j1) + EJp(:)*sf3;
      EJ(:, i,j1 ) = EJ(:, i ,j1) + EJp(:)*sf4;
    end
    % Species kについてもとまったので、これをEJplotに代入
    diffEJ = EJ - tmp;
    EJplot(1+k, 1, :, :) = diffEJ(1, :, :);
    EJplot(1+k, 2, :, :) = diffEJ(2, :, :) + diffEJ(3, :, :);
    EJplot(1+k, 3, :, :) = diffEJ(1, :, :) + diffEJ(2, :, :) + diffEJ(3, :, :);
  end
  % 6. 境界条件の処理
  EJ(:, 2,Y2) = EJ(:, 2 ,Y2) + EJ(:, nxp2,Y2  ); %Y2 = 2:ny+1;
  EJ(:, X2,2) = EJ(:, X2,2 ) + EJ(:, X2  ,nyp2); %X2 = 2:nx+1;
  EJ(:, 2,2)  = EJ(:, 2 ,2 ) + EJ(:, nxp2,nyp2);

  % perp方向を求めてEJ(2, :, :)に代入してしまう
  EJplot(1, 1, :, :) = EJ(1, :, :); % parallel
  EJplot(1, 2, :, :) = EJ(2, :, :) + EJ(3, :, :); % perpendicular
  EJplot(1, 3, :, :) = EJ(1, :, :) + EJ(2, :, :) + EJ(3, :, :); % all directions
  
  % figure of EJ plot
  paramEJ.spName = {'All Species', 'Species 1', 'Species 2', 'Species 3', 'Species 4'};
  paramEJ.direction = {'$E_{\parallel}\cdot J_{\parallel}$', '$E_{\perp}\cdot J_{\perp}$', '$E\cdot J$'};
  % if mod(jtime, inputParam.ndskip) == 0 %inputParam.ndskipの倍数の時だけ描画
    f_EJ = figure(88);
    f_EJ.Position = [0, 0, 1400, 900];
    f_EJ.Name = 'EJ plot';
    frames_EJ = getframe(f_EJ);
    for k=1:(inputParam.ns+1)*3
      switch mod(k,3)
        case 1
          EJcol = 1;
        case 2
          EJcol = 2;
        case 0
          EJcol = 3;
      end
      EJrow = (k-EJcol)/3+1;
      ax(k) = subplot(inputParam.ns+1, 3, k); imagesc(squeeze(EJplot(EJrow, EJcol, X2, Y2))');
      colormap(pltColor.mapEJ); 
      c(k) = colorbar; shading flat;
      ax(k).Title.FontSize = inputParam.Fontsize*0.8;
      ax(k).Title.Interpreter = 'latex';
      ax(k).DataAspectRatio = [100, 100, 1];
      tmpTitle = sprintf('%s, %s \n Time = %10.3f / %10.3f', ...
      cell2mat(paramEJ.spName(EJrow)), ...
      cell2mat(paramEJ.direction(EJcol)), jtime*inputParam.dt, inputParam.ntime*inputParam.dt);
      ax(k).Title.String = replace(tmpTitle, 'Time', '$t\Omega_e$');
      ax(k).LineWidth = 2;
      ax(k).XLabel.Interpreter = 'latex';
      ax(k).XLabel.FontSize = inputParam.Fontsize*0.8;
      ax(k).XLabel.String = '$x$';
      ax(k).YLabel.Interpreter = 'latex';
      ax(k).YLabel.FontSize = inputParam.Fontsize*0.8;
      ax(k).YLabel.String = '$y$';
      ax(k).YDir='normal';
      c(k).LineWidth = 2;
      caxis([-3e-11, 3e-11]);
    end
    writeVideo(v_EJ, frames_EJ);
    if mod(jtime, floor(inputParam.ntime/inputParam.ndskip/20)) == 0
      figFilename = ['EJ', num2str(jtime), '.fig'];
      savefig(f_EJ, figFilename);
      movefile(figFilename, figPath);
    end
  % end
end
function calcEJ(inputParam, nxp1, nxp2, nyp1, nyp2, np, ex, ey, ez, bx0, by0, bz0,  x, y, q, vx, vy, vz, jtime, X2, Y2, pltColor, v_EJ, figPath)
  [ejx,ejy,ejz] = c_calcEJ(x,y,vx,vy,vz,inputParam.ns,np,q,ex,ey,ez);

  ejx(2,Y2,:) = ejx(2 ,Y2,:) + ejx(nxp2,Y2  ,:); %Y2 = 2:ny+1;
  ejx(X2,2,:) = ejx(X2,2 ,:) + ejx(X2  ,nyp2,:); %X2 = 2:nx+1;
  ejx(2,2, :) = ejx(2 ,2 ,:) + ejx(nxp2,nyp2,:);

  ejy(2,Y2,:) = ejy(2 ,Y2,:) + ejy(nxp2,Y2  ,:); %Y2 = 2:ny+1;
  ejy(X2,2,:) = ejy(X2,2 ,:) + ejy(X2  ,nyp2,:); %X2 = 2:nx+1;
  ejy(2,2, :) = ejy(2 ,2 ,:) + ejy(nxp2,nyp2,:);

  ejz(2,Y2,:) = ejz(2 ,Y2,:) + ejz(nxp2,Y2  ,:); %Y2 = 2:ny+1;
  ejz(X2,2,:) = ejz(X2,2 ,:) + ejz(X2  ,nyp2,:); %X2 = 2:nx+1;
  ejz(2,2, :) = ejz(2 ,2 ,:) + ejz(nxp2,nyp2,:);

  EJplot = zeros(nxp2, nyp2, 3, 1+inputParam.ns);
  EJplot(:,:,1,1) = sum(ejx,3);              % parallel
  EJplot(:,:,2,1) = sum(ejy,3) + sum(ejz,3); % perpendicular
  EJplot(:,:,3,1) = EJplot(:,:,1,1) + EJplot(:,:,2,1);

  for is = 2:inputParam.ns+1
      EJplot(:,:,1,is) = ejx(:,:,is-1);                       % parallel
      EJplot(:,:,2,is) = ejy(:,:,is-1) + ejz(:,:,is-1);       %perpendicular
      EJplot(:,:,3,is) = EJplot(:,:,1,is) + EJplot(:,:,2,is); %All
  end
  % figure of EJ plot
  paramEJ.spName = {'All Species', 'Species 1', 'Species 2', 'Species 3', 'Species 4'};
  paramEJ.direction = {'$E_{\parallel}\cdot J_{\parallel}$', '$E_{\perp}\cdot J_{\perp}$', '$E\cdot J$'};
  %inputParam.ndskipつまり8の倍数の時だけ描画
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
      ax(k) = subplot(inputParam.ns+1, 3, k); imagesc(squeeze(EJplot(X2, Y2, EJcol, EJrow))');
      colormap(pltColor.mapEJ); 
      c(k) = colorbar; shading flat;
      ax(k).Title.FontSize = inputParam.Fontsize*0.8;
      ax(k).Title.Interpreter = 'latex';
      ax(k).DataAspectRatio = [100, 100, 1];
      ax(k).Title.String = sprintf('%s, %s \n Time = %10.3f / %10.3f', ...
      cell2mat(paramEJ.spName(EJrow)), ...
      cell2mat(paramEJ.direction(EJcol)), jtime*inputParam.dt, inputParam.ntime*inputParam.dt);
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
  end
end
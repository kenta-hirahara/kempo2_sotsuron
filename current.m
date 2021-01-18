
% ajxはnx+2行ny+2列の行列
ajx = zeros(nxp2, nyp2);
ajy = zeros(nxp2, nyp2);
ajz = zeros(nxp2, nyp2);

n2 = 0;
for  k=1:ns
	n1  = n2+1;
	n2  = n2 + np(k);
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
     
     ajxp = q(k)*vx(m); %ある超粒子のx方向電流密度
     ajyp = q(k)*vy(m);
     ajzp = q(k)*vz(m);
     
     
     ajx(i ,j ) = ajx(i ,j ) + ajxp*sf1;
     ajx(i1,j ) = ajx(i1,j ) + ajxp*sf2;
     ajx(i1,j1) = ajx(i1,j1) + ajxp*sf3;
     ajx(i,j1 ) = ajx(i ,j1) + ajxp*sf4;
     
     ajy(i ,j ) = ajy(i ,j ) + ajyp*sf1;
     ajy(i1,j ) = ajy(i1,j ) + ajyp*sf2;
     ajy(i1,j1) = ajy(i1,j1) + ajyp*sf3;
     ajy(i,j1 ) = ajy(i ,j1) + ajyp*sf4;
     
     ajz(i ,j ) = ajz(i ,j ) + ajzp*sf1;
     ajz(i1,j ) = ajz(i1,j ) + ajzp*sf2;
     ajz(i1,j1) = ajz(i1,j1) + ajzp*sf3;
     ajz(i,j1 ) = ajz(i ,j1) + ajzp*sf4;
  end
end
  % periodic boundary conditions
  ajx(2,Y2) = ajx(2 ,Y2) + ajx(nxp2,Y2  ); %Y2 = 2:ny+1;
  ajx(X2,2) = ajx(X2,2 ) + ajx(X2  ,nyp2); %X2 = 2:nx+1;
  ajx(2,2)  = ajx(2 ,2 ) + ajx(nxp2,nyp2);

  ajy(2,Y2) = ajy(2 ,Y2) + ajy(nxp2,Y2  );
  ajy(X2,2) = ajy(X2,2 ) + ajy(X2  ,nyp2);
  ajy(2,2)  = ajy(2 ,2 ) + ajy(nxp2,nyp2);
  
  ajz(2,Y2) = ajz(2 ,Y2) + ajz(nxp2,Y2  );
  ajz(X2,2) = ajz(X2,2 ) + ajz(X2  ,nyp2);
  ajz(2,2)  = ajz(2 ,2 ) + ajz(nxp2,nyp2);

% relocation of current densities from (Full grid, Full grid)
  ajx(nxp2,  Y2) = ajx(2,Y2);    
  ajx(X2  ,  Y2) = (ajx(X2,Y2)+ajx(X3,Y2))*0.5;
  
  ajy(X2  ,nyp2) = ajy(X2,2);
  ajy(X2  ,  Y2) = (ajy(X2,Y2)+ajy(X2,Y3))*0.5;
  
  ajz(nxp2,  Y2) = ajz(2,Y2);
  ajz(V2  ,nyp2) = ajz(V2,2);
  ajz(X2  ,  Y2) =(ajz(X2,Y2)+ajz(X3,Y2)+ajz(X2,Y3)+ajz(X3,Y3))*0.25;
% cancellation of uniform currents
if ajamp == 0
 ajxu = sum(ajx(X2,Y2), 'all')/nxny;
 ajx(X2,Y2) = ajx(X2,Y2) - ajxu;
 ajyu = sum(ajy(X2,Y2), 'all')/nxny;
 ajy(X2,Y2) = ajy(X2,Y2) - ajyu;
 ajzu = sum(ajz(X2,Y2), 'all')/nxny;
 ajz(X2,Y2) = ajz(X2,Y2) - ajzu;
end
% external current   
% if jtime < ctime
% switch check.externalCurrentLocation
%   case 0
    ajz(nxc , nyc) = ajamp*sin(omega*2*itime); 
  % case 1
%     ajz(nxc , :) = ajamp*sin(omega*2*itime); 
% end
% % ajx(nxc , nyc) = ajamp*sin(omega*2*itime); 
% % end

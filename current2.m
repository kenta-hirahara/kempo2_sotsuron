
ajx = zeros(nxp3, nyp3);
ajy = zeros(nxp3, nyp3);
ajz = zeros(nxp2, nyp2);

n2 = 0;
for  k=1:ns
	n1  = n2+1;
	n2  = n2 + np(k);
	for m = n1:n2
      x1 = x(m) - vx(m);
      x2 = x(m) + vx(m);
      y1 = y(m) - vy(m);
      y2 = y(m) + vy(m);
      i1 = floor(x1);
      i2 = floor(x2);
      j1 = floor(y1);
      j2 = floor(y2);
      xr = min(min(i1,i2)+1,max(max(i1,i2),0.5*(x1+x2)));
      yr = min(min(j1,j2)+1,max(max(j1,j2),0.5*(y1+y2)));
      Fx1 = q(k)*(xr - x1)*0.5;
      Fx2 = q(k)*vx(m) - Fx1;
      Fy1 = q(k)*(yr - y1)*0.5;
      Fy2 = q(k)*vy(m) - Fy1;
      Wx1 = (x1 + xr)*0.5 - i1;
      Wx2 = (xr + x2)*0.5 - i2;
      Wy1 = (y1 + yr)*0.5 - j1;
      Wy2 = (yr + y2)*0.5 - j2;
      i1p2 = i1 + 2;
      i1p3 = i1 + 3;
      j1p2 = j1 + 2;
      j1p3 = j1 + 3;
      i2p2 = i2 + 2;
      i2p3 = i2 + 3;
      j2p2 = j2 + 2;
      j2p3 = j2 + 3;
      ajx(i1p2,j1p2) = ajx(i1p2,j1p2) + Fx1*(1 - Wy1);
      ajx(i1p2,j1p3) = ajx(i1p2,j1p3) + Fx1*Wy1;
      ajy(i1p2,j1p2) = ajy(i1p2,j1p2) + Fy1*(1 - Wx1);
      ajy(i1p3,j1p2) = ajy(i1p3,j1p2) + Fy1*Wx1;
   
      ajx(i2p2,j2p2) = ajx(i2p2,j2p2) + Fx2*(1 - Wy2);
      ajx(i2p2,j2p3) = ajx(i2p2,j2p3) + Fx2*Wy2;
      ajy(i2p2,j2p2) = ajy(i2p2,j2p2) + Fy2*(1 - Wx2);
      ajy(i2p3,j2p2) = ajy(i2p3,j2p2) + Fy2*Wx2;
      
     xmf = x(m) + 2.0;
     ymf = y(m) + 2.0;
     i = floor(xmf);
     j = floor(ymf);
     ip1 = i+1;
     jp1 = j+1;
     x1 = xmf - i;
     y1 = ymf - j;
     sf3 = x1*y1;
     sf2 = x1-sf3;
     sf4 = y1-sf3;
     sf1 = 1.0-x1-sf4;
     
     ajzp = q(k)*vz(m);
      
     ajz(i ,  j ) = ajz(i  , j )  + ajzp*sf1;
     ajz(ip1, j ) = ajz(ip1, j )  + ajzp*sf2;
     ajz(ip1,jp1) = ajz(ip1,jp1)  + ajzp*sf3;
     ajz(i,  jp1) = ajz(i  ,jp1)  + ajzp*sf4;
  end
end
  % periodic boundary conditions
  ajx(2,WA)    = ajx(2 ,WA)   + ajx(nxp2,WA);
  ajx(nxp1,WA) = ajx(nxp1,WA) + ajx(1   ,WA);
  ajx(X2,2)    = ajx(X2,2 )   + ajx(X2,nyp2);
  ajx(X2,3)    = ajx(X2,3 )   + ajx(X2,nyp3);
  ajx(X2,nyp1) = ajx(X2,nyp1) + ajx(X2,   1);
  
  ajy(VA,2)    = ajy(VA,2)    + ajy(VA,nyp2);
  ajy(VA,nyp1) = ajy(VA,nyp1) + ajy(VA,   1);
  ajy(2,Y2)    = ajy(2,Y2)    + ajy(nxp2,Y2);
  ajy(3,Y2)    = ajy(3,Y2)    + ajy(nxp3,Y2);
  ajy(nxp1,Y2) = ajy(nxp1,Y2) + ajy(1   ,Y2);
  
  ajz(2,Y2) = ajz(2 ,Y2) + ajz(nxp2,Y2  );
  ajz(X2,2) = ajz(X2,2 ) + ajz(X2  ,nyp2);
  ajz(2,2)  = ajz(2 ,2 ) + ajz(nxp2,nyp2);

% relocation of current densities from (F,F)to (H, H)  
  ajz(nxp2,Y2) = ajz(2,Y2);
  ajz(V2,nyp2) = ajz(V2,2);
  ajz(X2,Y2)=(ajz(X2,Y2)+ajz(X3,Y2)+ajz(X2,Y3)+ajz(X3,Y3))*0.25;
% cancellation of uniform currents
if ajamp == 0
 ajxu = sum(sum(ajx(X2,Y2)))/nxny;
 ajx(X2,Y2) = ajx(X2,Y2) - ajxu;
 ajyu = sum(sum(ajy(X2,Y2)))/nxny;
 ajy(X2,Y2) = ajy(X2,Y2) - ajyu;
 ajzu = sum(sum(ajz(X2,Y2)))/nxny;
 ajz(X2,Y2) = ajz(X2,Y2) - ajzu;
end
% external current   
if jtime < ctime
ajz(nxc , nyc) = ajamp*sin(omega*2*itime); 
end

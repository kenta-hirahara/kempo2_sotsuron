% computing charge density defined at (F,F) grids

rho = rho0; %rho0がイオンの電荷に相当

n2 = 0;
for  k=1:ns
	n1  = n2+1;
	n2  = n2 + np(k);
	for m = n1:n2
     xmf = x(m) + 2.0;
     ymf = y(m) + 2.0;
     i = floor(xmf); %少数切り捨て
     j = floor(ymf);
     i1 = i+1;
     j1 = j+1;
     x1 = xmf - i;
     y1 = ymf - j;
     sf3 = x1*y1;
     sf2 = x1-sf3;
     sf4 = y1-sf3;
     sf1 = 1.0-x1-sf4;
     
     rho(i ,j ) = rho(i ,j ) + q(k)*sf1;
     rho(i1,j ) = rho(i1,j ) + q(k)*sf2;
     rho(i1,j1) = rho(i1,j1) + q(k)*sf3;
     rho(i,j1 ) = rho(i ,j1) + q(k)*sf4;
  end
end

  % periodic boundary conditions
  rho(2,Y2) = rho(2 ,Y2) + rho(nxp2,Y2  ) - rho0(nxp2,Y2); %初期値を引いてる
  rho(X2,2) = rho(X2,2 ) + rho(X2  ,nyp2) - rho0(X2,nyp2);
  rho(2,2)  = rho(2 ,2 ) + rho(nxp2,nyp2) - rho0(nxp2,nyp2);
  
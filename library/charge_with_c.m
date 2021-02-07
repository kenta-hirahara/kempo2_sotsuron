% computing charge density defined at (F,F) grids
rho = c_charge(x,y,ns,np,q,nx,ny);
rho = rho + rho0;
% periodic boundary conditions
rho(2,Y2) = rho(2 ,Y2) + rho(nxp2,Y2  ) - rho0(nxp2,Y2);
rho(X2,2) = rho(X2,2 ) + rho(X2  ,nyp2) - rho0(X2,nyp2);
rho(2,2)  = rho(2 ,2 ) + rho(nxp2,nyp2) - rho0(nxp2,nyp2);
  
  [ajx, ajy, ajz] = c_current(x,y,vx,vy,vz,ns,np,q,nx,ny);
  
  % periodic boundary conditions
  ajx(2,Y2) = ajx(2 ,Y2) + ajx(nxp2,Y2  );
  ajx(X2,2) = ajx(X2,2 ) + ajx(X2  ,nyp2);
  ajx(2,2)  = ajx(2 ,2 ) + ajx(nxp2,nyp2);
  ajy(2,Y2) = ajy(2 ,Y2) + ajy(nxp2,Y2  );
  ajy(X2,2) = ajy(X2,2 ) + ajy(X2  ,nyp2);
  ajy(2,2)  = ajy(2 ,2 ) + ajy(nxp2,nyp2);
  ajz(2,Y2) = ajz(2 ,Y2) + ajz(nxp2,Y2  );
  ajz(X2,2) = ajz(X2,2 ) + ajz(X2  ,nyp2);
  ajz(2,2)  = ajz(2 ,2 ) + ajz(nxp2,nyp2);

% relocation of current densities from (F,F)  
  ajx(nxp2,  Y2) = ajx(2,Y2);    
  ajx(X2  ,  Y2) = (ajx(X2,Y2)+ajx(X3,Y2))*0.5;
  ajy(X2  ,nyp2) = ajy(X2,2);
  ajy(X2  ,  Y2) = (ajy(X2,Y2)+ajy(X2,Y3))*0.5;
  ajz(nxp2,  Y2) = ajz(2,Y2);
  ajz(V2  ,nyp2) = ajz(V2,2);
  ajz(X2  ,  Y2) =(ajz(X2,Y2)+ajz(X3,Y2)+ajz(X2,Y3)+ajz(X3,Y3))*0.25;
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
% if jtime < ctime
%   ajy(nxc , 2:ny+1) = ajy(nxc , 2:ny+1) + ajamp*cos(omega*2*itime-(2*pi)*(1:ny)/ny); 
%   ajz(nxc , 2:ny+1) = ajz(nxc , 2:ny+1) + ajamp*sin(omega*2*itime-(2*pi)*(1:ny)/ny); 
% end



if jtime %< ctime
  switch inputParam.sigmoidJ
    case 0
      externalCurrent = ajamp*sin(omega*2*itime);
    case 1
      externalCurrent = 0.5*(tanh((itime-600)*0.005)+1)*ajamp*sin(omega*itime);
  end
  switch inputParam.location
    case 'central_point'
      switch inputParam.directionJ
        case 'x'
          ajx(nxc , nyc) = externalCurrent;  
        case 'y'
          ajy(nxc , nyc) = externalCurrent;  
        case 'z'
          ajz(nxc , nyc) = externalCurrent;  
      end
    case 'center_x'
      switch inputParam.directionJ
        case 'x'
          ajx(nxc , :) = externalCurrent;  
        case 'y'
          ajy(nxc , :) = externalCurrent;  
        case 'z'
          ajz(nxc , :) = externalCurrent;  
      end
    case 'center_y'
      switch inputParam.directionJ
        case 'x'
          ajx(: , nyc) = externalCurrent;  
        case 'y'
          ajy(: , nyc) = externalCurrent;  
        case 'z'
          ajz(: , nyc) = externalCurrent; 
      end
  end
end
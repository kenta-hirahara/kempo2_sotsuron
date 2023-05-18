% solving Maxwell's equation
ex(X2,Y2)=ex(X2,Y2)+(bz(X2,Y2) - bz(X2,Y1))*tcs - ajx(X2,Y2)*2; %X2 = 2:nxp1, Y2 = 2:nyp1;
ey(X2,Y2)=ey(X2,Y2)-(bz(X2,Y2) - bz(X1,Y2))*tcs - ajy(X2,Y2)*2;
ez(X2,Y2)=ez(X2,Y2)+(bx(X2,Y2)-bx(X2,Y3)+by(X3,Y2)-by(X2,Y2))*tcs - ajz(X2,Y2)*2;

% periodic boundary conditions 
ex(nxp2,Y2)=ex(2,   Y2);
ex(V2,nyp2)=ex(V2,   2);
ex(1 ,  W2)=ex(nxp1,W2);

ey(X2,nyp2)=ey(X2,   2);
ey(nxp2,W2)=ey(2,   W2);
ey(V2  , 1)=ey(V2,nyp1);

ez(nxp2,Y2)=ez(2,Y2);
ez(V2,nyp2)=ez(V2,2);
ez(1,W2)=ez(nxp1,W2);
ez(V2,1)=ez(V2,nyp1);

%  reflecting boundary conditions
%  ez(nxp2,Y2)=zeros(1,ny);
%  ez(2,Y2) = zeros(1,ny);
%  ez(V2,nyp2)= zeros(nxp1,1);
%  ez(V2,2) = zeros(nxp1,1);
%  ez(1,W2)= zeros(1,nyp1);
%  ez(nxp1,W2) = zeros(1,nyp1);
%  ez(V2,1)= zeros(nxp1,1);
%  ez(V2,nyp1) = zeros(nxp1,1);
% 

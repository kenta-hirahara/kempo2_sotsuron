% solving Maxwell's equation
% Lecture4のp3解いてる
bx(X2,Y2)=bx(X2,Y2)-ez(X2,Y2)+ez(X2,Y1); %bxはnx+2行 ny+2列の行列, X2 = 2:nxp1; Y2 = 2:nyp1; nx*nyの行列を更新
by(X2,Y2)=by(X2,Y2)+ez(X2,Y2)-ez(X1,Y2);
bz(X2,Y2)=bz(X2,Y2)-ex(X2,Y2)+ex(X2,Y3)-ey(X3,Y2)+ey(X2,Y2);

% periodic boundary conditions
bx(nxp2,Y2)=bx(2,Y2);
bx(V2,nyp2)=bx(V2,2);
bx(1,W2)=bx(nxp1,W2);

by(X2,nyp2)=by(X2,2);
by(nxp2,W2)=by(2,W2);
by(V2,1)=by(V2,nyp1);

bz(nxp2,Y2)=bz(2,Y2);
bz(V2,nyp2)=bz(V2,2);
bz(1,W2)=bz(nxp1,W2);
bz(V2,1)=bz(V2,nyp1);

% Computing Electrostatic Potential
%                Y. Omura

% Charge density for correction of electric field
% ex = zeros(nxp2,nyp2);
% ey = zeros(nxp2,nyp2);

psi(X2,Y2) =rho(X2,Y2)-ex(X2,Y2)+ex(X1,Y2)-ey(X2,Y2)+ey(X2,Y1);

% Fourier Transform of charge density
psi_k = fft2(psi(X2,Y2));

% Converting charge denisty spectra to potential spectra
psi_k = psi_k.*k_fft;

% Inverse Frourier Transform of potential spectra
psi(X2,Y2) = ifft2(psi_k);

psi(nxp2,Y2) = psi(2,Y2);
psi(X2,nyp2) = psi(X2,2);

% Correcting the electric field
ex(X2,Y2) = ex(X2,Y2) + psi(X2,Y2) - psi(X3,Y2);
ey(X2,Y2) = ey(X2,Y2) + psi(X2,Y2) - psi(X2,Y3);

% Boundary Condition
ex(nxp2,Y2)=ex(2,Y2);
ex(V2,nyp2)=ex(V2,2);
ex(1,W2)=ex(nxp1,W2);
ey(X2,nyp2)=ey(X2,2);
ey(nxp2,W2)=ey(2,W2);
ey(V2,1)=ey(V2,nyp1);

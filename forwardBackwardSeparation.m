function [Fx0, Fy0, Ffx, Ffy, Ffz, Fbx, Fby, Fbz] = forwardBackwardSeparation(Fx, Fy, Fz)
%forwardBackwardSeparation separate forward backward waves by means of
%helicity of waves
  dim = size(Fx);
  NX = dim(1)-2;
  NY = dim(2)-2;

  % 1. FFT2
  Fxk = fft2(Fx(2:NX+1,2:NY+1));
  Fyk = fft2(Fy(2:NX+1,2:NY+1));
  Fzk = fft2(Fz(2:NX+1,2:NY+1));

  % 2. Rotate B(Bxk, Byk, Bzk) to B'(0, Byk1, Bzk1) along wave vector 
  dkx = 2.0/double(NX); 
  dky = 2.0/double(NY);
  kx_axis = horzcat(0:NX/2, -NX/2+1:-1).*dkx;
  ky_axis = horzcat(0:NY/2, -NY/2+1:-1).*dky;
  [KY, KX] = meshgrid(ky_axis, kx_axis);

  cos_phi = abs(KX)./sqrt(KX.^2+KY.^2);
  sin_phi = abs(KY)./sqrt(KX.^2+KY.^2);
  %(KX, -KY), (-KX, KY)ÇÃóÃàÊÇÕâÒì]ï˚å¸Ç™ãtÇ»ÇÃÇ≈SIN(-É∆) -> -SIN(É∆)ÇæÇØãtï˚å¸Ç…âÒì]Ç∑ÇÈ
  sin_phi(1:end/2+1  ,end/2+2:end) = -abs(KY(1:end/2+1  ,end/2+2:end))./sqrt(KX(1:end/2+1  ,end/2+2:end).^2+KY(1:end/2+1  ,end/2+2:end).^2);
  sin_phi(end/2+2:end,1:end/2+1  ) = -abs(KY(end/2+2:end,1:end/2+1))./sqrt(KX(end/2+2:end,1:end/2+1).^2+KY(end/2+2:end,1:end/2+1).^2);

  cos_phi(1,1) = 0;
  sin_phi(1,1) = 0;

  Fx_kn =  Fxk.*cos_phi + Fyk.*sin_phi;
  Fy_kn = -Fxk.*sin_phi + Fyk.*cos_phi;
  Fz_kn = Fzk;

  % 3. ëOêiå„ëﬁîgÇÃèàóù 
  Ffy_kn = zeros(NX, NY);
  Ffz_kn = zeros(NX, NY);
  Fby_kn = zeros(NX, NY);
  Fbz_kn = zeros(NX, NY);
    
  Ffz_kn(1:end/2+1,:)   = 0.5*( 1i  * Fy_kn(1:end/2+1,  :) + Fz_kn(1:end/2+1,  :));
  Fbz_kn(1:end/2+1,:)   = 0.5*( -1i * Fy_kn(1:end/2+1,  :) + Fz_kn(1:end/2+1,  :));
  Ffz_kn(end/2+2:end,:) = 0.5*( -1i * Fy_kn(end/2+2:end,:) + Fz_kn(end/2+2:end,:));
  Fbz_kn(end/2+2:end,:) = 0.5*( 1i  * Fy_kn(end/2+2:end,:) + Fz_kn(end/2+2:end,:));
    
  Ffy_kn(1:end/2+1,  :) = Ffz_kn(1:end/2+1  ,:) * -1i;
  Ffy_kn(end/2+2:end,:) = Ffz_kn(end/2+2:end,:) * 1i;
  Fby_kn(1:end/2+1,  :) = Fbz_kn(1:end/2+1  ,:) * 1i;
  Fby_kn(end/2+2:end,:) = Fbz_kn(end/2+2:end,:) * -1i;
    
  % 4. âÒì]Çå≥Ç…ñﬂÇ∑ÅB
  Fx0_k = Fx_kn.*cos_phi;
  Fy0_k = Fx_kn.*sin_phi;
  
  Ffx_k = Fx_kn.*cos_phi - Ffy_kn.*sin_phi;
  Ffy_k = Fx_kn.*sin_phi + Ffy_kn.*cos_phi;
  Ffz_k = Ffz_kn;

  Fbx_k = Fx_kn.*cos_phi - Fby_kn.*sin_phi;
  Fby_k = Fx_kn.*sin_phi + Fby_kn.*cos_phi;  
  Fbz_k = Fbz_kn;  

  % 5. ãtFFT
  Fx0 = zeros(NX+2, NY+2);
  Fy0 = zeros(NX+2, NY+2);
  Ffx = zeros(NX+2, NY+2);
  Ffy = zeros(NX+2, NY+2);
  Ffz = zeros(NX+2, NY+2);
  Fbx = zeros(NX+2, NY+2);
  Fby = zeros(NX+2, NY+2);
  Fbz = zeros(NX+2, NY+2);

  Fx0(2:NX+1,2:NY+1) = real(ifft2(Fx0_k));
  Fy0(2:NX+1,2:NY+1) = real(ifft2(Fy0_k));
  Ffx(2:NX+1,2:NY+1) = real(ifft2(Ffx_k));
  Ffy(2:NX+1,2:NY+1) = real(ifft2(Ffy_k));
  Ffz(2:NX+1,2:NY+1) = real(ifft2(Ffz_k));
  Fbx(2:NX+1,2:NY+1) = real(ifft2(Fbx_k));
  Fby(2:NX+1,2:NY+1) = real(ifft2(Fby_k));
  Fbz(2:NX+1,2:NY+1) = real(ifft2(Fbz_k));
  
  % for periodic boundary
  Fx0(NX+2  , 2:NY+1) = Fx0(2     , 2:NY+1);
  Fx0(2:NX+1, NY+2  ) = Fx0(2:NX+1, 2     );
  Fx0(1     , 2:NY+2) = Fx0(NX+1  , 2:NY+2);
  Fy0(2:NX+1, NX+2  ) = Fy0(2:NX+1, 2     );
  Fy0(NX+2  , 2:NY+2) = Fy0(2     , 2:NY+2);
  Fy0(2:NX+2, 1     ) = Fy0(2:NX+2, NY+1  );
  
  Ffx(NX+2  , 2:NY+1) = Ffx(2     , 2:NY+1);
  Ffx(2:NX+1, NY+2  ) = Ffx(2:NX+1, 2     );
  Ffx(1     , 2:NY+2) = Ffx(NX+1  , 2:NY+2);
  Fbx(NX+2  , 2:NY+1) = Fbx(2     , 2:NY+1);
  Fbx(2:NX+1, NY+2  ) = Fbx(2:NX+1, 2     );
  Fbx(1     , 2:NY+2) = Fbx(NX+1  , 2:NY+2);

  Ffy(2:NX+1, NX+2  ) = Ffy(2:NX+1, 2     );
  Ffy(NX+2  , 2:NY+2) = Ffy(2     , 2:NY+2);
  Ffy(2:NX+2, 1     ) = Ffy(2:NX+2, NY+1  );
  Fby(2:NX+1, NX+2  ) = Fby(2:NX+1, 2     );
  Fby(NX+2  , 2:NY+2) = Fby(2     , 2:NY+2);
  Fby(2:NX+2, 1     ) = Fby(2:NX+2, NY+1  );

  Ffz(NX+2  , 2:NY+1) = Ffz(2     , 2:NY+1);
  Ffz(2:NX+2, NX+2  ) = Ffz(2:NX+2, 2     );
  Ffz(1     , 2:NY+2) = Ffz(NX+1  , 2:NY+2);
  Ffz(2:NX+2, 1     ) = Ffz(2:NX+2, NY+1  );
  Fbz(NX+2  , 2:NY+1) = Fbz(2     , 2:NY+1);
  Fbz(2:NX+2, NX+2  ) = Fbz(2:NX+2, 2     );
  Fbz(1     , 2:NY+2) = Fbz(NX+1  , 2:NY+2);
  Fbz(2:NX+2, 1     ) = Fbz(2:NX+2, NY+1  );
end


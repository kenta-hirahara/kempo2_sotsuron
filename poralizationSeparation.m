function [FRy, FRz, FLy, FLz] = poralizationSeparation(Fy, Fz)
%poralizationSeparation separate R and L waves from either forward or
%backward waves.
  
dim = size(Fy);
NX = dim(1)-2;
NY = dim(2)-2;

% 1. FFT2
Fyk = fft2(Fy(2:NX+1,2:NY+1));
Fzk = fft2(Fz(2:NX+1,2:NY+1));

% 2. Determine L and R mode;
FRy_k = zeros(NX, NY);
FRz_k = zeros(NX, NY);
FLy_k = zeros(NX, NY);
FLz_k = zeros(NX, NY);

FRz_k(1:end/2+1,:)   = 0.5*( 1i  * Fyk(1:end/2+1,:)   + Fzk(1:end/2+1,:));
FLz_k(1:end/2+1,:)   = 0.5*( -1i * Fyk(1:end/2+1,:)   + Fzk(1:end/2+1,:));
FRz_k(end/2+2:end,:) = 0.5*( -1i * Fyk(end/2+2:end,:) + Fzk(end/2+2:end,:));
FLz_k(end/2+2:end,:) = 0.5*( 1i  * Fyk(end/2+2:end,:) + Fzk(end/2+2:end,:));

FRy_k(1:end/2+1,  :) = FRz_k(1:end/2+1  ,:) * -1i;
FRy_k(end/2+2:end,:) = FRz_k(end/2+2:end,:) * 1i;
FLy_k(1:end/2+1,  :) = FLz_k(1:end/2+1  ,:) * 1i;
FLy_k(end/2+2:end,:) = FLz_k(end/2+2:end,:) * -1i;

FRy = zeros(NX+2, NY+2);
FRz = zeros(NX+2, NY+2);
FLy = zeros(NX+2, NY+2);
FLz = zeros(NX+2, NY+2);

FRy(2:NX+1,2:NY+1) = real(ifft2(FRy_k));
FRz(2:NX+1,2:NY+1) = real(ifft2(FRz_k));
FLy(2:NX+1,2:NY+1) = real(ifft2(FLy_k));
FLz(2:NX+1,2:NY+1) = real(ifft2(FLz_k));


% for periodic boundary
FRy(2:NX+1, NX+2  ) = FRy(2:NX+1, 2     );
FRy(NX+2  , 2:NY+2) = FRy(2     , 2:NY+2);
FRy(2:NX+2, 1     ) = FRy(2:NX+2, NY+1  );
FLy(2:NX+1, NX+2  ) = FLy(2:NX+1, 2     );
FLy(NX+2  , 2:NY+2) = FLy(2     , 2:NY+2);
FLy(2:NX+2, 1     ) = FLy(2:NX+2, NY+1  );

FRz(NX+2  , 2:NY+1) = FRz(2     , 2:NY+1);
FRz(2:NX+2, NX+2  ) = FRz(2:NX+2, 2     );
FRz(1     , 2:NY+2) = FRz(NX+1  , 2:NY+2);
FRz(2:NX+2, 1     ) = FRz(2:NX+2, NY+1  );
FLz(NX+2  , 2:NY+1) = FLz(2     , 2:NY+1);
FLz(2:NX+2, NX+2  ) = FLz(2:NX+2, 2     );
FLz(1     , 2:NY+2) = FLz(NX+1  , 2:NY+2);
FLz(2:NX+2, 1     ) = FLz(2:NX+2, NY+1  );
end


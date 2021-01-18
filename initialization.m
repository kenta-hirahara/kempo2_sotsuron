%  Initialization Variables
%
xlen = nx;
ylen = ny;
nxny = nx*ny;
g = 1.0;
npt =sum(np(1:ns));
b0 = wc/qm(1);

costh = cos(deg2rad(theta));
sinth = sin(deg2rad(theta));
cosph = cos(deg2rad(phi));
sinph = sin(deg2rad(phi));

bx0 = b0*sinth*cosph;
by0 = b0*sinth*sinph;
bz0 = b0*costh;
csq = cv*cv;
tcs = 2.0*csq;
it = 0;
  nxp1 = nx+1;
  nyp1 = ny+1;
  nxp2 = nx+2;
  nyp2 = ny+2;
  nxp3 = nx+3;
  nyp3 = ny+3;
  nxc = nx/2+1;
  nyc = ny/2+1;
  % nxny = nx*ny;
  bx = ones(nx+2, ny+2)*bx0;
  by = ones(nx+2, ny+2)*by0;
  bz = ones(nx+2, ny+2)*bz0;
  ex = zeros(nx+2, ny+2);
  ey = zeros(nx+2, ny+2);
  ez = zeros(nx+2, ny+2);
  rho = zeros(nx+2, ny+2);
  psi = zeros(nx+2, ny+2);
  
% definition of inductive electric field
%  exind = zeros(nxp2, nyp2);
%  exind(nx/4+1:nx/4*3, ny/8*3+1:ny/8*5) = ones(nx/4*2, ny/8*2)*exinductive;

  
% Indiex Arry for Matrix Operation
  X1 = 1:nx;
  Y1 = 1:ny;
  X2 = 2:nxp1;
  Y2 = 2:nyp1;
  X3 = 3:nxp2;
  Y3 = 3:nyp2;
  V1 = 1:nxp1;
  V2 = 2:nxp2;
  VA = 1:nxp3;
  W1 = 1:nyp1;
  W2 = 2:nyp2;
  WA = 1:nyp3;
  A2 = 2:nx;
  B2 = 2:ny;
  A3 = 3:nxp1; 
  B3 = 3:nyp1;
% Paritlce Initialization
 x = zeros(1,npt);
 y = zeros(1,npt);
 vx = zeros(1,npt);
 vy = zeros(1,npt);
 vz = zeros(1,npt);

n2 = 0;
for  k=1:ns
  n1 = n2 + 1;
  n2 = n2 + np(k);
  pchr = deg2rad(pch(k));
  vdpa = vd(k)*cos(pchr);
  vdpe = vd(k)*sin(pchr);

  x(n1:n2) = real(xlen*rand(1,np(k))); % xlen = nx
  % y(n1:n2) = ylen*rand(1,np(k));
  y(n1:n2) = real(ylen/np(k)*(1:np(k)));

  %x(n1:n2) = xlen*0.5+0.1*xlen*randn(1,np(k));
  %y(n1:n2) = ylen*0.5+0.1*ylen*randn(1,np(k));
  for i = n1:n2
    phase = 2*pi*rand;
    uxi = vPerp(k)*randn + vdpe*cos(phase);
    uyi = vPerp(k)*randn + vdpe*sin(phase);
    uzi = vPara(k)*randn + vdpa; 
    
    % rotation to the direction of the magnetic field
    % 座標変換しているので、uxが磁場に平行なx成分に必ずなるはず
    
    ux = uxi*costh*cosph - uyi*sinph + uzi*sinth*cosph;
    uy = uxi*costh*sinph + uyi*cosph + uzi*sinth*sinph;
    uz = -uxi*sinth                  + uzi*costh;
   
    g = cv /sqrt(csq +ux*ux +uy*uy +uz*uz);
   
    vx(i) = ux*g;
    vy(i) = uy*g;
    vz(i) = uz*g;

  end
end
xcen = xlen*0.5;
ycen = ylen*0.5;
%for m = 1:np(1)
%  r = sqrt((x(m)-xcen)^2+(y(m)-ycen)^2);
%   if (r < 4.0) & (r >0.1) 
%      vx(m) = (x(m) - xcen)*vd(1)/r;
%      vy(m) = (y(m) - ycen)*vd(1)/r;
%   end
%end   

% Coefficiendts for Poisson Solver with FFTW
k_fft = zeros(nx,ny);
kx_min = pi/nx;
ky_min = pi/ny;

% Higher than first order
% x = (1:kxmax), y = (1:kymax)
for iy = 2:ny/2+1
    rky  = sin(ky_min*(iy-1))*2.0;
    for ix = 2:nx/2+1
        rkx  = sin(kx_min*(ix-1))*2.0;
        k_fft(ix,     iy     ) = 1.0/(rkx.^2+rky.^2);
        k_fft(ix,     ny+2-iy) = 1.0/(rkx.^2+rky.^2);
        k_fft(nx+2-ix,iy     ) = 1.0/(rkx.^2+rky.^2);
        k_fft(nx+2-ix,ny+2-iy) = 1.0/(rkx.^2+rky.^2);
    end
end

% x = 0, y =(1:kymax)
for iy=2:ny/2+1
  rky  = sin(ky_min*(iy-1))*2.0;
  k_fft(1,iy     ) = 1.0/rky.^2;
  k_fft(1,ny+2-iy) = 1.0/rky.^2;
end

% x = (1:kxmax), y = 0
for ix= 2:nx/2+1
  rkx  = sin(kx_min*(ix-1))*2.0;
  k_fft(ix,     1) = 1.0/rkx.^2;
  k_fft(nx+2-ix,1) = 1.0/rkx.^2;
end
k_fft(1,1) = 0.0;
%
q = (xlen*ylen)./np(1:ns).*(wp(1:ns).^2)./qm(1:ns); %超粒子、1個あたりの粒子が持ってるweightがここで決まる, 
% このqを使えばいい. ここはプラズマ周波数から計算してるので重み付けになっている. 
% 例えば超粒子の数を2倍にしたとしても、その結果出てくる図の縦軸が変わって欲しくない, 分布関数として値を変えたくない.
% ここで取る方法は、粒子数を増やした分だけqが小さくなる
% density n = np/((nx*dx)*(ny*dy))  /// in 2D code
% q = rho/n
mass = q./qm(1:ns);
rho0 = - sum(q(1:ns).*np(1:ns))/nxny*ones(nxp2,nyp2);

color = get(gca,'ColorOrder');
fign_xy = 1;
fign_kxky = 2;
itime = 0;
jtime = 0;
rmax = dx*nx;
rmin = 0.0;
%ke=zeros(ntime,ns);
isFirstTime = true;
ke = zeros(floor(ntime/ndskip),ns);
At = zeros(floor(ntime/ndskip),ns);
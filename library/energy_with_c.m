it = it + 1;
kf = csq/nxny;
ef = 0.5/nxny;
bf = 0.5*csq/nxny;
ef2 = 0.5/(b0^2*nxny);
bf2 = 0.5*csq/(b0^2*nxny);
% kinetic energy
[ketmp, Attmp] = c_energy(vx, vy, vz, ns, np, mass, cv, bx0, by0, bz0);
ke(it,:) = ketmp'*kf;
engp(it) = sum(ke(it,:));
At(it,:) = Attmp';
% electric field energy
enge(it) = sum(sum(ex(X2,Y2).^2+ey(X2,Y2).^2+ez(X2,Y2).^2))*ef;
% magnetic field energy
engb(it) = sum(sum((bx(X2,Y2)-bx0).^2 + (by(X2,Y2)-by0).^2 + (bz(X2,Y2)-bz0).^2))*bf;
% totalenergy
engt(it) = engp(it) +enge(it)+engb(it);
% parallel and perpendicular energies
eepara(it) = sum(sum((ex(X2,Y2)*bx0 + ey(X2,Y2)*by0 + ez(X2,Y2)*bz0).^2))*ef2;
eeperp(it) = enge(it) - eepara(it); 
ebpara(it) = sum(sum(((bx(X2,Y2)-bx0)*bx0 + (by(X2,Y2)-by0)*by0 + (bz(X2,Y2)-bz0)*bz0).^2))*bf2;
ebperp(it) = engb(it) - ebpara(it); 
  
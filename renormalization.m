%  Renomalization of the Input Parameters
%
% distance
rnx = dx;
% time
rnt = dt/2;
% velocity
rnv = rnx/rnt;
% electric field
rne = rnx/(rnt*rnt);
% electric potential
rnp = rnx*rnx/(rnt*rnt);
% magnetic field
rnb = 1/rnt;
% current density
rna = rnx/(rnt*rnt*rnt);
% charge density
rnr = 1/(rnt*rnt);
% charge
rnq = rnx*rnx/(rnt*rnt);
% mass
rnm = rnx*rnx/(rnt*rnt);
%
wp = wp * rnt;
wc = wc * rnt;
omega = omega * rnt;
cv = cv/rnv;
vPara = vPara/rnv;
vPerp = vPerp/rnv;
vd  = vd /rnv;
ajamp = ajamp/rna;  
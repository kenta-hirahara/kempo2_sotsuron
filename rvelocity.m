pex(V2,W2)=(ex(V1,W2)+ex(V2,W2)).*0.5;
pey(V2,W2)=(ey(V2,W1)+ey(V2,W2)).*0.5;
pez(V2,W2)=(ez(V1,W1)+ez(V2,W1)+ez(V1,W2)+ez(V2,W2)).*0.25;

pbx(V2,W2)=(bx(V1,W2)+bx(V2,W2)).*0.5;
pby(V2,W2)=(by(V2,W1)+by(V2,W2)).*0.5;
pbz(V2,W2)=(bz(V1,W1)+bz(V2,W1)+bz(V1,W2)+bz(V2,W2)).*0.25;

n2 = 0;
for k=1:ns
  n1  = n2+1;
  n2  = n2 + np(k);
  for m = n1:n2
    xmf = x(m) + 2.0;
    ymf = y(m) + 2.0;
    i = floor(xmf);
    j = floor(ymf);
    i1 = i+1;
    j1 = j+1;
    x1 = xmf - i;
    y1 = ymf - j;
    sf3 = x1*y1;
    sf2 = x1-sf3;
    sf4 = y1-sf3;
    sf1 = 1.0-x1-sf4;
    
    sf1 = sf1*qm(k);
    sf2 = sf2*qm(k);
    sf3 = sf3*qm(k);
    sf4 = sf4*qm(k);
    
    ex1=sf1*pex(i,j)+sf2*pex(i1,j)+sf3*pex(i1,j1)+sf4*pex(i,j1);
    ey1=sf1*pey(i,j)+sf2*pey(i1,j)+sf3*pey(i1,j1)+sf4*pey(i,j1);
    ez1=sf1*pez(i,j)+sf2*pez(i1,j)+sf3*pez(i1,j1)+sf4*pez(i,j1);
    bx1=sf1*pbx(i,j)+sf2*pbx(i1,j)+sf3*pbx(i1,j1)+sf4*pbx(i,j1);
    by1=sf1*pby(i,j)+sf2*pby(i1,j)+sf3*pby(i1,j1)+sf4*pby(i,j1);
    bz1=sf1*pbz(i,j)+sf2*pbz(i1,j)+sf3*pbz(i1,j1)+sf4*pbz(i,j1);
    g = cv/sqrt(csq-vx(m)*vx(m)-vy(m)*vy(m)-vz(m)*vz(m));
    ux = vx(m)*g + ex1;
    uy = vy(m)*g + ey1;
    uz = vz(m)*g + ez1;
    g = cv/sqrt(csq + ux*ux + uy*uy + uz*uz);
    bx1 = bx1*g;
    by1 = by1*g;
    bz1 = bz1*g;
    boris = 2./(1. + bx1*bx1 + by1*by1 + bz1*bz1);
    uxt = ux + uy*bz1 - uz*by1;
    uyt = uy + uz*bx1 - ux*bz1;
    uzt = uz + ux*by1 - uy*bx1;
    ux = ux + boris*(uyt*bz1 - uzt*by1) + ex1;
    uy = uy + boris*(uzt*bx1 - uxt*bz1) + ey1;
    uz = uz + boris*(uxt*by1 - uyt*bx1) + ez1;
    g = cv /sqrt(csq +ux*ux +uy*uy +uz*uz);
    
    vx(m) = ux*g;  
    vy(m) = uy*g; 
    vz(m) = uz*g;
  end
end


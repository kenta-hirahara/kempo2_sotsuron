 % advancing position by half time step
  x = x + vx;
  y = y + vy;
 % periodic boundary condition 
  x = mod(real(x),xlen);
  y = mod(real(y),ylen);
  

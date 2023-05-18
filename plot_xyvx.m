%plotting a 3D phase diagram of x-y-vx
n2 = 0;    
for  k=1:ns
  n1 = n2 + 1;
  n2 = n2 + np(k);
   h= plot3(x(n1:n2)*rnx, y(n1:n2)*rnx,vx(n1:n2)*rnv,'color',color(mod(k-1,7)+1,:)); 
   h.Marker = 'o';  
   h.MarkerSize = 2
   hold on;
end
    xlabel('X'), ylabel('Y'), zlabel('Vx'),
    grid on, axis([rmin,rmax, rmin,rmax,vmin,vmax]),
    title(['time = ',num2str(jtime*dt)]),axis('square'),
    hold off;
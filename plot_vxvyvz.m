%plotting a 3D phase diagram of vx-vy-vz
n2 = 0;    
for  k=1:ns
  n1 = n2 + 1;
  n2 = n2 + np(k);
  h = plot3(vx(n1:n2)*rnv, vy(n1:n2)*rnv,vz(n1:n2)*rnv,'color',color(mod(k-1,7)+1,:));
  h.Marker = 'o';
  h.MarkerSize = 2;
  hold on;
end
    %  ylabel('Vy'), zlabel('Vz'),
    grid on
    % axis([vmin,vmax, vmin,vmax,vmin,vmax])
    % axis('square')
    hold off;

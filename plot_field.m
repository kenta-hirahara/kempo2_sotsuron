%plotting the field Ex, Ey, Ez, Bx, By, Bz

if fplot==1
        surf(ex(X2,Y2)'*rne),
        xlabel('X'), ylabel('Y'),zlabel('Ex');
elseif fplot==2   
        surf(ey(X2,Y2)'*rne),
        xlabel('X'), ylabel('Y'),zlabel('Ey');
elseif fplot==3
        surf(ez(X2,Y2)'*rne),
        xlabel('X'), ylabel('Y'),zlabel('Ez');
elseif fplot == 4
        surf((bx(X2,Y2)-bx0)'*rnb),
        xlabel('X'), ylabel('Y'),zlabel('Bx');
elseif fplot == 5
        surf((by(X2,Y2)-by0)'*rnb),
        xlabel('X'), ylabel('Y'),zlabel('By');
elseif fplot == 6
        surf((bz(X2,Y2)-bz0)'*rnb),
        xlabel('X'), ylabel('Y'),zlabel('Bz');
end
if (zmax+zmin) == 0
    zlim([zmin, zmax]);
end

hold off;
         
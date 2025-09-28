function dx = odefun1(t,x)
u = 1;
if(t>1)
    u = 0;
else
    u = 1;
end
dx = zeros(2,1);
dx(1) = x(2);
dx(2) = -2*0.5*10*x(2)-10^2*x(1)+10^2*u;
end
 
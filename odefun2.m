function dx = odefun2(t,x)
b = 1;
r = 10;
A = 0.1;
f = 0.1271;
omega = 2*pi*f;
dx = zeros(2,1);
dx(1) = b*x(2);
dx(2) = x(2)*(x(2)-1)*(1-r*x(2))-x(1)+A/omega*cos(omega*t);
end
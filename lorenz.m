function dx = lorenz(t,x)
dx = zeros(3,1);
alpha = 8/3;
beta = 10;
gamma = 28;
dx(1) = -alpha*x(1)+x(2)*x(3);
dx(2) = -beta*x(2)+beta*x(3);
dx(3) = -x(1)*x(2)+gamma*x(2)-x(3);
end

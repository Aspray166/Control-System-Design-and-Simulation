%% 第二题第一种方法
clc
clear 
tspan = [0 30];
x0 = [0 0 0.001];
[t,x] = ode45('lorenz',tspan,x0);
comet3(x(:,1),x(:,2),x(:,3))
grid
 
%% 第二题第二种方法
clc
clear
% out = sim('task2_2_1.slx');
out = sim('task2_2_2.slx');
yout = zeros(length(out.tout),3);
yout(:,1) = out.yout{1}.Values.Data;
yout(:,2) = out.yout{2}.Values.Data;
yout(:,3) = out.yout{3}.Values.Data;
figure(1)
plot(out.tout,yout)
figure(2)
comet3(yout(:,1),yout(:,2),yout(:,3))
grid
axis([min(yout(:,1)),max(yout(:,1)),min(yout(:,2)),...
    max(yout(:,2)),min(yout(:,3)),max(yout(:,3))])
%% 第三题
load_system('task2_3_1');
find_system('task2_3_1','Type','Block')
P = get_param('task2_3_1/PID Controller','P');
I = get_param('task2_3_1/PID Controller','I');
%% 第四题
%% 第五题
num = 1;
den = [10 1];
G = tf(num,den);
[y,t] = step(G)
plot(t,y)
% syms t s
% f = 1/(10*s^2+1);
% F = ilaplace(f);





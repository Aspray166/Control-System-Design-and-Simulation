%% 实验一
%% 第二题
%求解析解
syms y(x)
eqn = diff(y,x) == 1/x^2-y/x;
cond = y(1) == 1;
ySol(x)= dsolve(eqn,cond)
clc
clear
step = 0.1;
e = 0.001;
% 前项欧拉
Y = zeros(1/step+1,5);
Y(1,1) = 1;
Y(1,2) = 1;
j=1;
for i = 1+step:step:2
    Y(j+1,1) = i;
    Y(j+1,2) = Y(j,2)+step*(1/Y(j,1)^2-Y(j,2)/Y(j,1));
    j = j+1;
end
% 后项欧拉
B = zeros(1/step+1,2);
for i = 1:1/step+1
    B(i,1) = Y(i,2);
end
epsilon = 1;
while epsilon>e
    j = 1;
    B(1,2) = B(1,1);
    for i = 1+step:step:2
        B(j+1,2) = B(j,2)+step*(1/Y(j+1,1)^2-B(j+1,1)/Y(j+1,1));
        j = j+1;
    end
    epsilon = max(abs(B(:,1)-B(:,2)));
    for i = 1:1/step+1
        B(i,1) = B(i,2);
    end
end
for i=1:1/step+1
    Y(i,3) = B(i,2);
end
%梯形公式
T = zeros(1/step+1,2);
for i = 1:1/step+1
    T(i,1) = Y(i,2);
end
epsilon = 1;
while epsilon>e
    j = 1;
    T(1,2) = T(1,1);
    for i = 1+step:step:2
        T(j+1,2) = T(j,2)+0.5*step*((1/Y(j,1)^2-T(j,2)/Y(j,1))+(1/Y(j+1,1)^2-T(j+1,1)/Y(j+1,1)));
        j = j+1;
    end
    epsilon = max(abs(T(:,1)-T(:,2)));
    for i = 1:1/step+1
        T(i,1) = T(i,2);
    end
end
for i=1:1/step+1
    Y(i,4) = T(i,2);
end
%解析解得到的准确值
for i=1:1/step+1
    Y(i,5) = (log(Y(i,1))+1)/Y(i,1);
end
%% 第三题
clear
ksi = 0.5;
omega = 10;
%初值
x = [0;0];y = 0;t0 = 0;tf = 5;  
h = 0.05;%步长
t = t0;%从t0开始迭代，储存时间
M = round(1/h);
N = round((tf-t0)/h)-M;     %迭代步数
r = 1;%单位阶跃输入
num = [omega^2];
den = [1,2*ksi*omega,omega^2];
[A,B,C,~] = tf2ss(num,den);
for i = 1:M
    k1 = A * x+B*r;
    k2 = A * (x+h*k1/2)+B*r;
    k3 = A * (x+h*k2/2)+B*r;
    k4 = A * (x+h*k3)+B*r;
    x = x+h*(k1+2*k2+2*k3+k4)/6;                    %采用四阶龙格库塔法
    y = [y,C*x];                                    %输出值
    t = [t,t(i)+h];
end
r = 0;
for i = M+1:N
    k1 = A * x+B*r;
    k2 = A * (x+h*k1/2)+B*r;
    k3 = A * (x+h*k2/2)+B*r;
    k4 = A * (x+h*k3)+B*r;
    x = x+h*(k1+2*k2+2*k3+k4)/6;                    %采用四阶龙格库塔法
    y = [y,C*x];                                    %输出值
    t = [t,t(i)+h];
end
plot(t,y);
%% 第四题
clear
ksi = 0.5;
omega = 10;
T = 5;
%初值
x = [0;0;0];y = 0;t0 = 0;tf = 50;  
h = 0.05;%步长
t = t0;%从t0开始迭代，储存时间
N = round((tf-t0)/h);     %迭代步数
r = 1;%单位阶跃输入
num = omega^2;
den = conv([1,2*ksi*omega,omega^2],[T,1]);
[A,B,C,D] = tf2ss(num,den);
for i = 1:N
    k1 = A * x+B*r;
    k2 = A * (x+h*k1/2)+B*r;
    k3 = A * (x+h*k2/2)+B*r;
    k4 = A * (x+h*k3)+B*r;
    x = x+h*(k1+2*k2+2*k3+k4)/6;                    %采用四阶龙格库塔法
    y = [y,C*x];                                    %输出值
    t = [t,t(i)+h];
end
plot(t,y);
tr = interp1(y,t',0.9) - interp1(y,t',0.1);%线性插值求上升时间
ts = interp1(y,t',0.95);%调节时间
title(['上升时间' num2str(tr) ' 调节时间' num2str(ts)]);
%% 第五题
clear
tspan = [0 5];
y0 = [0 0];
[t,x] = ode45('odefun1',tspan,y0);
plot(t,x(:,1));
%% 第六题
clear
tspan = [0 100];
y0 = [0.2 0.1];
[t,x] = ode45('odefun2',tspan,y0);
plot(t,x(:,1),'-',t,x(:,2),'--')

%% 第一题
clc
clear
num=10;
den=conv([1,0],[1,1]);
g=tf(num,den);
sys = feedback(g,1)
[num1,den1]=tfdata(sys,'v');
[z,p,k]=tf2zp(num1,den1);
sys2=zpk(z,p,k)
sys3=tf2ss(num1,den1)
[r,p,k]=residue(num1,den1);
[n,d]=residue(r,p,k)
%% 第二题
clc
clear
% 1_1 第一种方法
num = [1 2 3 1];
den = [1 5 2 1 1];
sys1 = tf(num,den);
[z,p,k] = tf2zp(num,den);
ii = find(real(p)>0);
n1 = length(ii);
if(n1>0)
    disp('The Unnstable Poles are:')
    disp(p(ii));
else
    disp('System is Stable');
end
pzmap(num,den);
title('Zero-Pole Map');
%1_2 第二种方法
clc
clear
num = [1 2 3 1];
den = [1 5 2 1 1];
r = roots(den);
ii = find(real(r)>0);
n1 = length(ii);
if(n1>0)
    disp('The Unnstable Poles are:')
    disp(r(ii));
else
    disp('System is Stable');
end
% 2 第一种方法
clc
clear
T = 0.1;
num = [-3 2];
den = [1 -0.2 0.25 0.05];
sys2 = tf(num,den,'Ts',T);
r = roots(den);
ii = find(abs(r)>1);
n1 = length(ii);
if(n1>0)
    disp(['System is Unstable,with' "int2str(n1)" 'unstable pole']);
else
    disp('System is Stable');
end
% 3_1 第一种方法
clc
clear
A = [1 3;5 2];
P = poly(A);
r = roots(P);
ii = find(real(r)>0);
n = length(ii);
if(n>0)
    disp('System is Unstable')
else
    disp('System is Stable');
end
% 3_2 第二种方法
clc
clear
A = [1 3;5 2];Q = eye(size(A));
P = lyap(A,Q);i1 = find(P(1,1)>0);
n1 = length(i1);i2 = find(det(P)>0);
n2 = length(i2);
if(n1>0&&n2>0)
    disp('P>0;Positive Definite;The equilibrium state at the original point is asymptotic stable');
else
    disp('System is Unstable');
end
%% 第三题
clc
clear
T = 0.1;
num = 1;
den = [1 1];
sys = tf(num,den);
[u0,t] = gensig('square',5,30,T);
u = zeros(length(u0),1);
for i = 1:length(u0)-1
    u(i+1) = u(i)+0.5-u0(i+1);
end
umax = max(u);
u = u/umax;
figure(1);
lsim(sys,u,t)
% 附加
u = zeros(1/T+1,1);
m = 0;
for i = 0:T:3
    m = m+1;
    u(m) = i*(i>=0&&i<1)+(i-1)*(i>=1&&i<2)+(i-2)*(i>=2&&i<3);
end
t = (0:T:3)';
figure(2);
lsim(sys,u,t)
%% 第四题
clc
clear
wn = 5;
ksi = [0:0.1:1];
figure(1);
hold on
j = 0;
for i=ksi
    j = j+1;
    num = wn^2;
    den = [1 2*i*wn wn^2];
    [y,~,t] = step(num,den,8);
    plot(t,y);
    if(i==0)
        str{j} = ['等幅振荡'];
        continue;
    end
    [p,m] = find_pole(y);%找到响应曲线中极点及其所对应的位置
    if(~isempty(p))%如果极点的数组不为空，即有超调
        n = p(1)/p(3);%衰减比为第一个峰值比上第二个峰值
        overshoot = (p(1)-1)/1;%超调量为第一个峰值超过稳态值的高度
        tp = t(m(1));%峰值时间为第一个极点对应的时间
        y1 = y(1:m(1));
        t1 = t(1:m(1));%从0开始到第一个峰值划分一个区间
        %响应曲线进入±5%之前一个极点与之后一个极点之间划分一个区间
        p1 = find(p>=1.05,1,'last');
        p2 = find(p<=0.95,1,'last');
        if(isempty(p2))
            if(isempty(p1))
                k = find(y>0.95,1,'first');
                y2 = y(1:k);
                t2 = t(1:k);
            else
                k = p1;
                y2 = y(m(k):m(k+1));
                t2 = t(m(k):m(k+1));
            end
        else
            k = max(p1,p2,'omitnan');
            y2 = y(m(k):m(k+1));
            t2 = t(m(k):m(k+1));
        end
        tr = interp1(y1,t1',1);
        ts = max(interp1(y2,t2',1.05),interp1(y2,t2',0.95));    
        str{j} = ['ksi= ' num2str(i) ' tr=' num2str(tr,3) ' ts=' num2str(ts,3)...
         ' tp=' num2str(tp,3) ' overshoot=' num2str(overshoot*100,3) '% n=' num2str(n,3)];
    else
        k = find(y>0.95,1,'first');
        y1 = y(1:k);
        t1 = t(1:k);
        tr = interp1(y1,t1',0.9) - interp1(y1,t1',0.1);
        ts = interp1(y1,t1',0.95);
        str{j} = ['ksi= ' num2str(i) ' tr=' num2str(tr,3) ' ts=' num2str(ts,3)];    
    end
end
legend(str);
grid on;
hold off
%% 第五题
clc
clear
num = 1;
den = conv(conv([2 1],[1 1]),[0.1 1]);
[r,~] = rlocus(num,den);
figure(1);
plot(r,'x')
axis([-2 1 0.5 3.5])
grid on
[K,poles] = rlocfind(num,den);
figure(2);
[Gm,Pm,Wcg,Wcp] = margin(10*num,den);
margin(10*num,den)
%% 第六题(1)
clc
clear
n1 = 2;
d1 = conv([1 0],conv([0.1 1],[0.3 1]));
sope = tf(n1,d1);
sys = feedback(sope,1);
step(sys)
hold on
%确定期望极点位置
sigma = 0.2;
zeta = ((log(1/sigma))^2/((pi)^2+(log(1/sigma))^2))^(1/2);
wn = 3.5/(zeta*1);
p = [1 2*zeta*wn wn*wn];
r = roots(p);
%求校正器的传递函数
kc = 5;
s_1 = r(1,1);
ngv = polyval(n1,s_1);
dgv = polyval(d1,s_1);
g = ngv/dgv;
zetag = angle(g);
mg = abs(g);
ms = abs(s_1);
zetas = angle(s_1);
tz = (sin(zetas)-kc*mg*sin(zetag-zetas))/(kc*mg*ms*sin(zetag));
tp = -(kc*mg*sin(zetas)+sin(zetag+zetas))/(ms*sin(zetag));
nk=[tz,1];
dk = [tp,1];
Gc = tf(nk,dk)
%校验矫正器
sys = feedback(kc*sope*Gc,1);
step(sys);
grid
legend('校正前','校正后')
%% 第六题（2）
clc
clear
k0 = 3;
n1 = 2;
d1 = conv([1 0],conv([0.1 1],[0.3 1]));
[mag,phase,w] = bode(k0*n1,d1);
figure(1);
margin(mag,phase,w);
hold on
figure(2)
s1 = tf(k0*n1,d1);
sys = feedback(s1,1);
step(sys)
hold on

gama = 45*pi/180;
alfa = (1+sin(gama))/(1-sin(gama));
adb = 20*log10(mag);
am = -10*log10(alfa);
wc = spline(adb,w,am);
T = 1/(wc*sqrt(alfa));
Gc = tf([alfa*T,1],[T,1])

[mag,phase,w] = bode(s1*Gc);
figure(1)
margin(mag,phase,w)
legend('校正前','校正后')
figure(2)
sys = feedback(s1*Gc,1);
step(sys)
legend('校正前','校正后')
%% 第七题
clc
clear
num = 1;
den = conv([1,1],[2,1]);
G = tf(num,den);
Kp = 2;
Ti = [3,6,14,21,25];
for i=1:5
    G1 = tf([Kp,Kp/Ti(i)],[1,0]);
    sys = feedback(G*G1,1);
    step(sys);
    hold on
    str{i} = ['积分时间Ti= ' num2str(Ti(i))];
end
legend(str);
%% 第八题
clc
clear
ess = 0.1;
k = 1/ess;
n1 = k;
d1 = conv([1 0],[1,1]);
[mag,phase,w] = bode(n1,d1);
figure(1);
margin(mag,phase,w);
legend('校正前');
hold on
s1 = tf(n1,d1);

gama = 45*pi/180;
alfa = (1+sin(gama))/(1-sin(gama));
adb = 20*log10(mag);
am = -10*log10(alfa);
wc = spline(adb,w,am);
T = 1/(wc*sqrt(alfa));
Gc = tf([alfa*T,1],[T,1])

[mag,phase,w] = bode(s1*Gc);
figure(2)
margin(mag,phase,w)
legend('校正后')



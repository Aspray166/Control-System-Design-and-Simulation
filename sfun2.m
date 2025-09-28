function [sys,x0,str,ts,simStateCompliance] = sfun2(t,x,u,flag)

num = 1;
den = [10 1];
G = tf(num,den);
[Y,T] = step(G);
switch flag,
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;
  case 1,
    sys=mdlDerivatives(t,x,u);
  case 2,
    sys=mdlUpdate(t,x,u);
  case 3,
    sys=mdlOutputs(t,x,u,Y,T);
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);
  case 9,
    sys=mdlTerminate(t,x,u);
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [0 0];
simStateCompliance = 'UnknownSimState';
function sys=mdlDerivatives(t,x,u)

sys = [];
function sys=mdlUpdate(t,x,u)

sys = [];
function sys=mdlOutputs(t,x,u,Y,T)
N = find(T>=t,1);
sys = Y(N);
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;
function sys=mdlTerminate(t,x,u)

sys = [];

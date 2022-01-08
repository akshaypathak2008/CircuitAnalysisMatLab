f0=1.35e+08;
np=4;
nt=101;
dt=1/(nt*f0);
df=f0/np;
time=0:dt:(np/f0)-dt;
freq=0:df:(nt*f0)-df;
vip=10;
K=1.38e-23;
T=300;
q=1.602e-19;
VT=K*T/q;
IS=1e-015;
VON=1;
VBR=-50;
GBR=10;
G1=1e-03;
G2=1e-03;
Vcc=5;
vi=vip*sin(2*pi*f0*time);
syms x
I=@(x) ((x<=VBR).*GBR.*(x-VBR)+(VBR<x & x<VON).*0+(x>=VON).*IS.*(exp((x-VON)/VT)-1));
x0=0;
opts=optimoptions('fsolve','Algorithm','levenberg-marquardt',...
'FunctionTolerance',1e-12,'StepTolerance',1e-12,'Display','off');
for k=1:length(vi)
    x=fsolve(@KCL,x0,opts);
vd1(k)=x;
x0=x;
end
vo=vi-vd1;
figure
plot(vi,vo)
figure
plot(vd1,I(vd1))
figure
plot(time,vi)
hold on
plot(time,vo)
plot(time,vd1)
function F=KCL(x)
global vi k I
global G1 G2 Vcc
F=G1*(vi(k)-x)+G2*(vi(k)-x-Vcc)-I(x);
end
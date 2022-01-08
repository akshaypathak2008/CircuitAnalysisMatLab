global vp
G1=1e-03;
G2=1e-03;
G3=1e-03;
C1=1e-09;
Vcc=5;
tf=20e-06;
ts=0;
ys=0;
time=0;
cap=0;
Vo=-Vcc;
out=Vo;
c=char('');
syms t v(t)
while ts < tf
vp=G2*Vo/(G1+G2);
EQ=C1*diff(v(t),t) == G3*(Vo-v(t));
[V]=odeToVectorField(EQ);
M=matlabFunction(V,'vars', {'t','Y'});
opts=odeset('Events',@yval);
[T,y,te,ye,ie]=ode45(M,[ts tf],ys,opts);
ts=T(end);
ys=y(end);
time=vertcat(time,T(1:end));
cap=vertcat(cap,y(1:end));
out=vertcat(out,repmat(Vo,length(y),1));
Vo=-Vo;
fprintf(repmat('\b',[1 length(c)]))
c=char([num2str(round(100*ts/tf)),' percent done!']);
fprintf(c)
end
fprintf('\n')
Vi=cap;
Vo=out;
figure
plot(Vi(2:end),Vo(2:end))
figure
plot(time,Vi)
hold on
plot(time,Vo)
function [value,isterminal,direction]=yval(~,y)
global vp
value=y(1)-vp;
isterminal=1;
direction=0;
end
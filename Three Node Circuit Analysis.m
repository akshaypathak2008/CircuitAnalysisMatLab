syms s
syms G C L gm
syms YL YS
syms v1 v2 v3
syms i1 i2 i3
syms ins inG inl
Vvec=[v1 v2 v3];
Nvec=inG;
equ=[i1==YS*(v1-v2)+ins;
i2==YS*(v2-v1)-ins+(1/(s*L))*v2+(G+s*C)*(v2-v3)+inG;...
i3==(G+s*C)*(v3-v2)-inG+gm*v1+YL*v3+inl];
[Nmat,b]=equationsToMatrix(equ,Nvec);
[Ymat,a]=equationsToMatrix(b,Vvec);
YP=subs(Ymat((1:end)~=1,(1:end)~=1),[YS,YL],[0,0]);
NP=-Nmat((1:end)~=1,:);
A=YP;
in=1;
out=2;
A([2,out],:)=A([out,2],:);
A(:,[2,out])=A(:,[out,2]);
A([1,in],:)=A([in,1],:);
A(:,[1,in])=A(:,[in,1]);
B=NP;
B([2,out],:)=B([out,2],:);
B([1,in],:)=B([in,1],:);
[ma,na]=size(A);
[mb,nb]=size(B);
for kk=1:ma-2
for km=1:ma-1
for kn=1:nb
B(km,kn)=B(km,kn)-(B(ma,kn)*A(km,na))/A(ma,na);
end
for kn=1:na-1
A(km,kn)=A(km,kn)-(A(ma,kn)*A(km,na))/A(ma,na);
end
end
ma=ma-1;
na=na-1;
end
IN1=sum(abs(B(1,:)).^2.*Nvec);
IN2=sum(abs(B(2,:)).^2.*Nvec)+inl;
Y11=A(1,1);
Y12=A(1,2);
Y21=A(2,1);
Y22=A(2,2);
disp('KCL equations :')
disp(equ)
disp('two-port Y-parameters :')
disp(A(1:2,1:2))
disp('two-port noise matrix :')
disp(B(1:2,:))
disp('noise equations :')
disp([['ini = ';'ino = '],char(string([IN1;IN2]))])
fprintf('\n')
f1=1e08;
f2=1e10;
f0=1.6e09;
np=4;
nt=101;
dt=1/(nt*f0);
time=0:dt:(np/f0)-dt;
f=logspace(log10(f1),log10(f2),100*(log10(f2)-log10(f1))+1);
G=1e-03;
L=1e-09;
C=10e-12;
gm=1e-03;
s=1i*2*pi*f;
K=1.38e-23;
T=300;
Z0=1;
YL=logspace(log10(0.1/Z0),log10(10/Z0),11).';
YS=(1/50).*ones(size(YL));
inG=4*K*T*G;
ins=4*K*T*real(YS);
inl=4*K*T*real(YL);
vns=4*K*T./real(YS);
vnl=4*K*T./real(YL);
Yin=Y11-Y12.*Y21./(Y22+YL);
NPower=2*YS.*real(Yin)./(abs(YS+Yin)).^2;
inti=abs(eval(IN1)-eval(IN2).*abs(eval(Y12./(Y22+YL))).^2);
vnti=inti./abs(eval(Yin)).^2;
ini=eval(IN1);
vni=eval(IN2).*abs(eval(Y12./(Y11.*(Y22+YL)-Y12.*Y21))).^2;
vnivar=sum(vni.*repmat([0,diff(f)],size(YL)),2);
inivar=sum(ini.*repmat([0,diff(f)],size(YL)),2);
vnsvar=sum(vns.*repmat([0,diff(f)],size(YL)),2);
insvar=sum(ins.*repmat([0,diff(f)],size(YL)),2);
NF=1+(vnivar+inivar./real(YS).^2)./vnsvar;
vnit=sqrt(vnivar).*randn(1,length(time));
init=sqrt(inivar).*randn(1,length(time));
loglog(f,sqrt(vni))
loglog(f,sqrt(ini))
plot(time,vnit)
plot(time,init)
semilogy(1./YL,NF)
mesh(f,1./YL,eval(NPower)*100)
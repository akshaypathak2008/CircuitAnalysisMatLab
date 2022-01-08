% The following MATLAB code lists the complete optimization algorithm in order
% to accurately model the component values of the substrate inductor model using the
% Nelder–Mead algorithm. The MATLAB code first reduces the inductor model to two
% nodes in order to match the number of equations with the measured two-port
% S-parameters data. It finally extracts the model parameters by comparing the measured
% S-parameters data with results generated from the inductor model.

syms G L C Cox Gsi Csi
syms s
syms Z0
Y=[G+s*(C+0.5*Cox) -0.5*s*Cox 0 -G -s*C;...
-0.5*s*Cox 0.5*Gsi+0.5*s*(Csi+Cox) 0 0 0;...
0 0 0.5*Gsi+0.5*s*(Csi+Cox) 0 -0.5*s*Cox;...
-G 0 0 G+1/(s*L) -1/(s*L);...
-s*C 0 -0.5*s*Cox -1/(s*L) 1/(s*L)+s*(C+0.5*Cox)];
in=1;
out=5;
A=Y;
A([2,out],:)=A([out,2],:);
A(:,[2,out])=A(:,[out,2]);
A([1,in],:)=A([in,1],:);
A(:,[1,in])=A(:,[in,1]);
[ma,na]=size(A);
for kk=1:ma-2
for km=1:ma-1
for kn=1:na-1
A(km,kn)=A(km,kn)-(A(ma,kn)*A(km,na))/A(ma,na);
end

end
ma=ma-1;
na=na-1;
end
Y11=A(1,1);
Y12=A(1,2);
Y21=A(2,1);
Y22=A(2,2);
Ydif=(Y11*Y22-Y12*Y21)/(Y11+Y12+Y21+Y22);
deltaS=(1+Z0*Y11).*(1+Z0*Y22)-Z0^2*Y12.*Y21;
S11=((1-Z0*Y11).*(1+Z0*Y22)+Z0^2*Y12.*Y21)./deltaS;
S12=-2*Z0*Y12./deltaS;
S21=-2*Z0*Y21./deltaS;
S22=((1+Z0*Y11).*(1-Z0*Y22)+Z0^2*Y12.*Y21)./deltaS;
Sm=[
1.000e+08 0.135297+0.120255i 0.864700-0.120489i 0.864700-0.120489i 0.135297+0.120255i;...
1.259e+08 0.144849+0.149739i 0.855146-0.150033i 0.855146-0.150033i 0.144849+0.149739i;...
1.585e+08 0.159575+0.185301i 0.840416-0.185671i 0.840416-0.185671i 0.159575+0.185301i;...
1.995e+08 0.181932+0.227148i 0.818055-0.227613i 0.818055-0.227613i 0.181932+0.227148i;...
2.512e+08 0.215094+0.274510i 0.784885-0.275095i 0.784885-0.275095i 0.215094+0.274510i;...
3.162e+08 0.262624+0.324921i 0.737343-0.325658i 0.737343-0.325658i 0.262624+0.324921i;...
3.981e+08 0.327503+0.373535i 0.672445-0.374461i 0.672445-0.374461i 0.327503+0.373535i;...
5.012e+08 0.410405+0.413110i 0.589513-0.414273i 0.589513-0.414273i 0.410405+0.413110i;...
6.310e+08 0.507827+0.435513i 0.492044-0.436971i 0.492044-0.436971i 0.507827+0.435513i;...
7.943e+08 0.611645+0.434773i 0.388152-0.436596i 0.388152-0.436596i 0.611645+0.434773i;...
1.000e+09 0.711326+0.410005i 0.288357-0.412276i 0.288357-0.412276i 0.711326+0.410005i;...
1.259e+09 0.797837+0.365790i 0.201672-0.368604i 0.201672-0.368604i 0.797837+0.365790i;...
1.585e+09 0.866458+0.309667i 0.132793-0.313125i 0.132793-0.313125i 0.866458+0.309667i;...
1.995e+09 0.916904+0.248879i 0.081971-0.253079i 0.081971-0.253079i 0.916904+0.248879i;...
2.512e+09 0.951637+0.188435i 0.046720-0.193454i 0.046720-0.193454i 0.951637+0.188435i;...
3.162e+09 0.974000+0.130809i 0.023680-0.136688i 0.023680-0.136688i 0.974000+0.130809i;...
3.981e+09 0.987049+0.076512i 0.009820-0.083249i 0.009820-0.083249i 0.987049+0.076512i;...
5.012e+09 0.993044+0.024789i 0.002938-0.032356i 0.002938-0.032356i 0.993044+0.024789i;...
6.310e+09 0.993309-0.025832i 0.001798+0.017439i 0.001798+0.017439i 0.993309-0.025832i;...
7.943e+09 0.988191-0.077161i 0.006133+0.067864i 0.006133+0.067864i 0.988191-0.077161i;...
1.000e+10 0.977025-0.131020i 0.016658+0.120617i 0.016658+0.120617i 0.977025-0.131020i];
f=Sm(:,1);
S11m=Sm(:,2);
S12m=Sm(:,3);
S21m=Sm(:,4);
S22m=Sm(:,5);
x=sym('x',[1,6]);
xs=[sym('G'),sym('L'),sym('C'),sym('Cox'),sym('Gsi'),sym('Csi')];
S11x=subs(subs(subs(S11,sym('Z0'),50),xs,x),s,1i*2*pi*f);
S12x=subs(subs(subs(S12,sym('Z0'),50),xs,x),s,1i*2*pi*f);
S21x=subs(subs(subs(S21,sym('Z0'),50),xs,x),s,1i*2*pi*f);
S22x=subs(subs(subs(S22,sym('Z0'),50),xs,x),s,1i*2*pi*f);
fun=@(x)sseval(x,S11x,S12x,S21x,S22x,S11m,S12m,S21m,S22m);
%G,L,C,Cox,Gsi,Csi
%x0=[0.074405,2.5133e-08,2.9749e-14,7.4374e-15,0.00033333,3.4531e-15];
x0=[1,1e-09,1e-12,1e-12,1,1e-12];
options=optimset('Display','iter','TolFun',1e-06,'TolX',1e-06);
bestx=fminsearch(fun,x0,options);
G0=bestx(1);
L0=bestx(2);
C0=bestx(3);
Cox0=bestx(4);
Gsi0=bestx(5);
Csi0=bestx(6);
disp(['G = ',num2str(G0,'%0.3e')])
disp(['L = ',num2str(L0,'%0.3e')])
disp(['C = ',num2str(C0,'%0.3e')])
disp(['Cox = ',num2str(Cox0,'%0.3e')])
disp(['Gsi = ',num2str(Gsi0,'%0.3e')])
disp(['Csi = ',num2str(Csi0,'%0.3e')])
function sse=sseval(x,S11x,S12x,S21x,S22x,S11m,S12m,S21m,S22m)
x1=x(1);
x2=x(2);
x3=x(3);
x4=x(4);
x5=x(5);
x6=x(6);
if 1<0
disp([x1,x2,x3,x4,x5,x6])
end
sse=(sum((real(double(eval(S11x)))-real(S11m)).^2)+...
sum((imag(double(eval(S11x)))-imag(S11m)).^2))*...
(sum((real(double(eval(S12x)))-real(S12m)).^2)+...
sum((imag(double(eval(S12x)))-imag(S12m)).^2))*...
(sum((real(double(eval(S21x)))-real(S21m)).^2)+...
sum((imag(double(eval(S21x)))-imag(S21m)).^2))*...
(sum((real(double(eval(S22x)))-real(S22m)).^2)+...
sum((imag(double(eval(S22x)))-imag(S22m)).^2));
end
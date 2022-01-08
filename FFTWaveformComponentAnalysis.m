% The following MATLAB code calculates the Fast Fourier Transform of a composite
% signal and plots its transient waveform and spectrum.

f0=1e06;
np=10;
nt=101;
dt=1/(nt*f0);
df=f0/np;
time=0:dt:(np/f0)-dt;
freq=0:df:(nt*f0)-df;
N=np*nt;
sig1=7*sin(2*pi*1e06*time);
sig2=5*sin(2*pi*2e06*time);
sig3=3*sin(2*pi*3e06*time);
sig=sig1+sig2+sig3;
SIG1=[];
out=double.empty(0,1);
for k1=1:N
for k2=1:N
out(k2)=sig(k2).*exp(-1i*2*pi*((-N/2)+k2-1)*((-N/2)+k1-1)/N);
end
tmp=[SIG1 sum(out)];
SIG1=tmp;
end
SIG2=fft(sig);
figure
plot(time,sig1)
hold on
plot(time,sig2)
plot(time,sig3)
figure
plot(time,sig)
figure
semilogx(freq(1:round(end/2)),2*abs(SIG1(end-round(end/2)+1:
end))/N)
hold on
semilogx(freq(1:round(end/2)),2*abs(SIG2(1:round(end/2)))/N,'--')

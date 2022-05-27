% %1、利用I型线性相位滤波器设计满足下列指标的FIR高通滤波器。Ωp=0.6πrad, Ωs=0.4πrad，Ap≤ 0.3dB,  As≥40dB
% Wp=0.6*pi;Ws=0.4*pi;Ap=0.3;As=40;
% N=ceil(6.3*pi/(Wp-Ws));
% N=mod(N+1,2)+N;
% M=N-1;
% %hamming窗长度N满足条件
% w=hamming(N)';
% %理想低通截频，计算截频
% Wc=(Wp+Ws)/2;
% k=0:M;
% hd=-(Wc/pi)*sinc(Wc*(k-0.5*M)/pi);
% hd(0.5*M+1)=hd(0.5*M+1)+1;
% h=hd.*w;
% omega=linspace(0,pi,512);
% mag=freqz(h,[1],omega);
% %phasez(h,1)
% plot(omega/pi,20*log10(abs(mag)));grid;%画出增益响应
% omega1=linspace(0,Wp,512);
% h1=freqz(h,[1],omega1);
% omega2=linspace(0,Ws,512);
% h2=freqz(h,[1],omega2);
% fprintf('Ap=%.4f\n',-20*log10(max(abs(h1))));
% fprintf('As=%.4f\n',-20*log10(max(abs(h2))));%计算AS和AP
% a=roundn(-20*log10(max(abs(h1))),-4);
% b=roundn(-20*log10(max(abs(h2))),-4);
% s = sprintf('Ap=%.4e\nAs=%.4e',a,b);
% text(0.5,-20,s);
% % 2、利用I型线性相位滤波器和 Fir1 函数设计满足下列指标的 FIR 高通滤波器。Ωp =0.8πrad, Ωs =0.7πrad Ap ≤ 0.3dB, As ≥ 40dB
% Wp=0.8*pi;Ws=0.7*pi;Ap=0.3;As=40;
% N=ceil(6.2*pi/(Wp-Ws));
% N=mod(N+1,2)+N;
% M=N-1;
% %hamming窗长度N满足条件
% w=hamming(N)';
% %理想低通截频，计算截频
% Wc=(Wp+Ws)/2;
% k=0:M;
% hd=-(Wc/pi)*sinc(Wc*(k-0.5*M)/pi);
% hd(0.5*M+1)=hd(0.5*M+1)+1;
% h=hd.*w;
% omega=linspace(0,pi,512);
% mag=freqz(h,[1],omega);
% %phasez(h,1)
% b=fir1(N,Wc/pi,'high');%不选择窗函数默认使用海明窗
% omega=linspace(0,pi,512);
% mag=freqz(b,[1],omega);
% phasez(b,1)
% plot(omega/pi,20*log10(abs(mag)));grid;%画出增益响应
% omega1=linspace(0,Wp,512);
% h1=freqz(b,[1],omega1);
% omega2=linspace(0,Ws,512);
% h2=freqz(b,[1],omega2);
% fprintf('Ap=%.4f\n',-20*log10(max(abs(h1))));
% fprintf('As=%.4f\n',-20*log10(max(abs(h2))));%计算AS和AP
% a=roundn(-20*log10(max(abs(h1))),-4);
% b=roundn(-20*log10(max(abs(h2))),-4);
% s = sprintf('Ap=%.4e\nAs=%.4e',a,b);
% text(0.2,0,s);
% %3.试用Kaiser窗设计满足下列指标的有2个通带FIR滤波器Ωs1=0.1πrad, Ωp1=0.2πrad, Ωp2=0.4πrad, Ωs2=0.5πrad, Ωs3=0.6πrad, Ωp3=0.7πrad, Ωp4=0.8πrad, Ωs4=0.9πrad, δs=0.008
% f=[0.1 0.2 0.4 0.5 0.6 0.7 0.8 0.9];
% a=[0,1,0,1,0];Rs=0.008;
% dev=Rs*ones(1,length(a));
% [N,Wc,beta,ftype] = kaiserord(f,a,dev);
% h = fir1(N,Wc,ftype,kaiser(N+1,beta));
% [h1,w1]=freqz(h,1);
%  plot(w1/pi,20*log10(abs(h1)));
% xlabel('Normalized frequency');
% ylabel('Gain, db');grid;
% axis([0 1 -80 5]);
%4.利用频率取样法设计Wp=0.6πrad的I型线性相位高通数字滤波器。（分为有过度点和无过度点两种情况）
M=32;Wp=0.6*pi;m=0:M/2;
Wm=2*pi*m./(M+1);
mtr=ceil(Wp*(M+1)/(2*pi));
Ad=double([Wm>=Wp]);
Ad(mtr)=0.28;
Hd=Ad.*exp(-j*0.5*M*Wm);
Hd=[Hd conj(fliplr(Hd(2:end)))];
h=real(ifft(Hd));
w=linspace(0,pi,1000);
H=freqz(h,[1],w);
figure(1)
plot(w/pi,20*log10(abs(H)));
xlabel('Normalized frequency');
ylabel('Gain,db');
axis([0 1 -50 0]);
f1=200;f2=700;f3=800;
fs=2000;
text(0.2,-20,'无过渡点');
grid on;
% 
% M=32;Wp=0.6*pi;m=0:M/2;
% Wm=2*pi*m./(M+1);
% Ad=double([Wm>=Wp]);
% Hd=Ad.*exp(-j*0.5*M*Wm);
% Hd=[Hd conj(fliplr(Hd(2:end)))];
% h=real(ifft(Hd));
% w=linspace(0,pi,1000);
% H=freqz(h,[1],w);
% figure(1)
% plot(w/pi,20*log10(abs(H)));
% xlabel('Normalized frequency');
% ylabel('Gain,db');
% axis([0 1 -50 0]);
% f1=200;f2=700;f3=800;
% fs=2000;
% text(0.2,-20,'有过渡点');
% grid on;

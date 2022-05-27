% n1=0;n2=7;n0=4;
% n=n1:n2;
% N=length(n);
% xn=[(n-n0)>=0];%建立时域信号 
% subplot(3,1,1); 
% stem(n,xn);
% title('x(n)');
% grid on;
% k=0:N-1;
% Xk=fft(xn,N);%用FFT计算信号的DFT 
% subplot(3,1,2); 
% stem(k,abs(Xk));
% title('Xk=DFT((n))');
% grid on;
% xn1=ifft(Xk,N);%用IFFT计算信号的IDFT
% subplot(3,1,3);
% stem(n,xn1); 
% title('x(n)=IDFT(Xk)');
% grid on;
% --------------------------------------------------------------------------
% Fs=10;
% xn=[1,2,3,2,1];N=length(xn);
% D=2*pi*Fs/N; %计算模拟频率分辨率 
% k=floor(-(N-1)/2:(N-1)/2);%频率显示范围对应［－pi,pi］
% X=fftshift(fft(xn,N)); %作FFT运算且移位p
% subplot(1,2,1);plot(k*D,abs(X),'o:'); %横轴化成模拟频率作幅度谱 
% title('幅度频谱');
% xlabel('rad/s');
% subplot(1,2,2);
% plot(k*D,angle(X),'o:');
% title('相位频谱');
% xlabel('rad/s'); 
% --------------------------------------------------------------------------
% Fs=10;N=1000;
% xn=[1,2,3,2,1];Nx=length(xn); 
% xn=[1,2,3,2,1,zeros(1,N-Nx-1)];
% D=2*pi*Fs/N;
% k=floor(-(N-1)/2:(N-1)/2);
% X=fftshift(fft(xn,N));
% subplot(1,2,1);plot(k*D,abs(X)); 
% title('幅度频谱');xlabel('rad/s'); 
% subplot(1,2,2);plot(k*D,angle(X)); 
% title('相位频谱');
% xlabel('rad/s');
% --------------------------------------------------------------------------
% Fs=20;C=[8,16,128]; %输入不同的N值
% for r=0:2;
% N=C(r+1);
% n=0:N-1;
% xn=exp(-0.5*n);%建立x(n) 
% D=2*pi*Fs/N;
% k=floor(-(N-1)/2:(N-1)/2); 
% X=fftshift(fft(xn,N));
% subplot(3,2,2*r+1);
% plot(k*D,abs(X));
% axis([-80,80,0,3]);
% subplot(3,2,2*r+2);
% stairs(k*D,angle(X));
% axis([-80,80,-1,1]);
% end
% %--------------------------------------------------------------------------
% T0=[0.5,0.25,0.125,0.125]; %输入不同的Ts值 
% N0=[256,256,256,2048];%输入不同的N值 
% for r=1:4;
% Ts=T0(r);N=N0(r);%赋Ts和N值 
% n=0:N-1;
% D=2*pi/(Ts*N);%计算模拟频率分辨率
% xa=exp(-0.01*n*Ts).*(sin(2*n*Ts)+sin(2.1*n*Ts)+sin(2.2*n*Ts));
% k=floor(-(N-1)/2:(N-1)/2);
% Xa=Ts*fftshift(fft(xa,N));
% [r,Xa(1)]%输出Xa(1)的数值,供误差计算用 
% subplot(2,2,r);
% plot(k*D,abs(Xa),'b');
% axis([1,3,1.1*min(abs(Xa)),1.1*max(abs(Xa))]);
% end
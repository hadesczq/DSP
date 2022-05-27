% %1、设计一个满足下列指标的BW型、 CB I型、椭圆模拟低通滤波器fp=1kHz, fs=2kHz, Ap≤1dB, As≥40dB
% Wp=2*pi*1000;Ws=2*pi*2000;Ap=1;As=40; % Analog filter specifications
% -----------------------BW型-----------------------------------------------
%  [N,Wc]=buttord(Wp,Ws,Ap,As,'s'); %Computer the order of analog filter
%  [num,den] = butter(N,Wc,'s'); %Compute AF coefficients
% 
% -----------------------CB I型---------------------------------------------
% [N,Wc]=cheb1ord(Wp,Ws,Ap,As,'s'); %Computer the order of analog filter
% [num,den] = cheby1(N,Ap,Wc,'s'); %Compute AF coefficients
% -----------------------椭圆滤波器-----------------------------------------
%  [N,Wc]=ellipord(Wp,Ws,Ap,As,'s'); %Computer the order of analog filter
%  [num,den] = ellip(N,Ap,As,Wc,'s'); %Compute AF coefficients
% --------------------------------------------------------------------------
% disp('Numerator polynomial'); fprintf('%.4e\n',num);
% disp('Denominator polynomial'); fprintf('%.4e\n',den);
% omega=[Wp Ws]; h = freqs(num,den,omega); %Compute Ap and As of AF
% fprintf('Ap= %.4f\n',-20*log10(abs(h(1))));
% fprintf('As= %.4f\n',-20*log10(abs(h(2))));
% a=roundn(-20*log10(abs(h(1))),-4);
% b=roundn(-20*log10(abs(h(2))),-4);
% omega = [0: 200: 3000*2*pi];
% h = freqs(num,den,omega); %Compute the frequency response of the AF
% gain=20*log10(abs(h)); 
% % 创建 figure
% figure1 = figure;
% % 创建 axes
% axes1 = axes('Parent',figure1);
% plot(omega/(2*pi),gain);
% % title({'BW型'});
% %  title({'CB I型'});
%  title({'椭圆型'});
% s = sprintf('N=%d\nAp=%.4e\nAs=%.4e',N,a,b);
% text(2000,-10,s);
% xlabel('Frequency in Hz'); ylabel('Gain in dB');
% %-------------------------------------------------------------------------
% %2、试设计满足下列指标的BW型模拟高通滤波器fp=5kHz, fs=1kHz, Ap≤1dB, As≥40dB。
% fp=5000;fs=1000; Ap=1;As=40; % AF specifications
% wLp=1/(2*pi*fp); wLs=1/(2*pi*fs); % HP specifications to LP’s
% [N,Wc]=buttord(wLp,wLs,Ap,As,'s'); % Computer AF order
% [num,den] = butter(N,Wc,'s'); % Compute AF coefficients
% [numt,dent] = lp2hp(num,den,1); % LP to HP
% omega=[2*pi*fp, 2*pi*fs];
% h = freqs(numt, dent, omega);
% fprintf('Ap= %.4f\n',-20*log10(abs(h(1))));
% fprintf('As= %.4f\n',-20*log10(abs(h(2))));
% a=roundn(-20*log10(abs(h(1))),-4);
% b=roundn(-20*log10(abs(h(2))),-4);
% omega = [0: 200: 6000*2*pi];
% h = freqs(numt, dent, omega); % Compute the frequency response of the AF
% gain=20*log10(abs(h)); plot(omega/(2*pi),gain);
% s = sprintf('Ap=%.4e\nAs=%.4e',a,b);
% text(3000,-40,s);
% xlabel('Frequency in Hz'); ylabel('Gain in dB')
% %-------------------------------------------------------------------------
% 3、试设计满足下列指标的BW型模拟带通滤波器wp1=6 rad/s, wp2=8 rad/s, Ap≤1dB, ws1=4 rad/s, ws2=11 rad/s, As≥32dB
% wp=1; ws=3.3182; Ap=1; As=32; % AF specifications
% w0=sqrt(48); B=2;
% [N,Wc]=buttord(wp,ws,Ap,As,'s'); % Computer low-pass AF order
% [num,den] = butter(N,Wc,'s'); % Compute low-pass AF coefficients
% [numt,dent] = lp2bp(num,den,w0,B); % LP to BP
% omega=[6 8 4 11];
% h = freqs(numt,dent,omega); %Compute Ap and As of the AF
% fprintf('Ap1= %.4f\n',-20*log10(abs(h(1))));
% fprintf('Ap2= %.4f\n',-20*log10(abs(h(2))));
% fprintf('As1= %.4f\n',-20*log10(abs(h(3))));
% fprintf('As2= %.4f\n',-20*log10(abs(h(4))));
% a=roundn(-20*log10(abs(h(1))),-4);
% b=roundn(-20*log10(abs(h(2))),-4);
% c=roundn(-20*log10(abs(h(3))),-4);
% d=roundn(-20*log10(abs(h(4))),-4);
% w=linspace(2,12,1000);
% h=freqs(numt,dent,w);
% plot(w,20*log10(abs(h))) ; grid on;
% s = sprintf('Ap1=%.4e\nAp2=%.4e\nAs1=%.4e\nAs2=%.4e',a,b,c,d);
% text(6,-30,s);
% xlabel('Frequency in rad/s');
% ylabel('Gain in dB')
% %-------------------------------------------------------------------------
% %4、试设计满足下列指标的BW型模拟带阻滤波器wp1=10 rad/s, wp2=30 rad/s, Ap≤1dB, ws1=19 rad/s, ws2=21 rad/s, As≥20dB
% Ap=1;As=20; wp1=10; wp2=30; ws1=19; ws2=21; % BS AF specification
% B=ws2-ws1;w0=sqrt(ws1*ws2);
% wLp1=B*wp1/(w0*w0-wp1*wp1); wLp2=B*wp2/(w0*w0-wp2*wp2);
% wLp=max(abs(wLp1),abs(wLp2)); wLs=1;
% [N,Wc]=buttord(wLp, wLs, Ap, As, 's' ); % Computer low-pass AF order
% [num,den] = butter(N, Wc, 's' ); % Compute low-pass AF coefficients
% [numt,dent]=lp2bs(num,den,w0,B); % LP to BS
% omega=[10 30 19 21];
% h = freqs(numt,dent,omega);
% fprintf('Ap1= %.4f\n',-20*log10(abs(h(1))));
% fprintf('Ap2= %.4f\n',-20*log10(abs(h(2))));
% fprintf('As1= %.4f\n',-20*log10(abs(h(3))));
% fprintf('As2= %.4f\n',-20*log10(abs(h(4))));
% a=roundn(-20*log10(abs(h(1))),-4);
% b=roundn(-20*log10(abs(h(2))),-4);
% c=roundn(-20*log10(abs(h(3))),-4);
% d=roundn(-20*log10(abs(h(4))),-4);
% w=linspace(5,35,1000); h=freqs(numt,dent,w);
% plot(w,20*log10(abs(h)));
% s = sprintf('Ap1=%.4e\nAp2=%.4e\nAs1=%.4e\nAs2=%.4e',a,b,c,d);
% text(6,-30,s);
% xlabel('Frequency in rad/s'); ylabel('Gain in dB')
% %-------------------------------------------------------------------------
% %5、利用BW型模拟低通滤波器设计满足指标Wp=0.2πrad, Ws=0.6πrad, Ap≤2dB,As≥15dB的数字低通滤波器。
% Wp=0.2*pi; Ws=0.6*pi; Ap=2; As=15;%DF BW LP specification
% Fs=1; %Sampling frequency(Hz)
% wp=Wp*Fs; ws=Ws*Fs;%Analog Butterworth specification
% [N,wc]=buttord(wp,ws,Ap,As,'s');%determine the order of AF filter
% [numa,dena]=butter(N,wc,'s');%determine the AF-BW filter
% [numd,dend]=impinvar(numa,dena,Fs);%脉冲响应不变法
% [numd1,dend1]=bilinear(numa,dena,Fs);%双线性变换法
% w=linspace(0,pi,1024);%plot the frequency response
% h=freqz(numd,dend,w);
% h1=freqz(numd1,dend1,w);
% plot(w/pi,20*log10(abs(h)))
% hold on
% plot(w/pi,20*log10(abs(h1)))
% xlabel('frequency in rad');
% ylabel('Gain in dB');% computer Ap As of the designed filter
% w=[Wp Ws];
% h=freqz(numd,dend,w);
% h1=freqz(numd1,dend1,w);
% fprintf('脉冲响应不变法\n');
% fprintf('Ap= %.4f\n',-20*log10( abs(h(1))));
% fprintf('As= %.4f\n',-20*log10( abs(h(2))));
% fprintf('双线性变换法\n');
% fprintf('Ap= %.4f\n',-20*log10( abs(h1(1))));
% fprintf('As= %.4f\n',-20*log10( abs(h1(2))));
% a=roundn(-20*log10(abs(h(1))),-4);
% b=roundn(-20*log10(abs(h(2))),-4);
% c=roundn(-20*log10(abs(h1(1))),-4);
% d=roundn(-20*log10(abs(h1(2))),-4);
% s = sprintf('脉冲响应不变法\nAp1=%.4e\nAp2=%.4e\n双线性变换法\nAs1=%.4e\nAs2=%.4e',a,b,c,d);
% text(0.3,-100,s);
% 
% 

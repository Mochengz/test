clc;
clear all;
close all;
fs=1e8;      %采样频率100MHZ
Tc=8e-6;    %脉冲长度8us
Ts=1/fs;
T=0.5e-3;    %脉冲序列间隔 0.5ms
N=Tc/Ts;
startt=-(N-1)/2;%脉冲起始处
endd=(N-1)/2;%脉冲截止处
fI=10e6;  %脉冲数字起始频率
B=4e6;       %带宽4MHZ
chirp_rate=B/Tc;   %信号调频率5*10^11
t1=startt*(1/fs):1/fs:endd*(1/fs);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t2=0*(1/fs):1/fs:(N-1)*(1/fs); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pulse1=exp(1j*(2*pi*fI*t1+pi*chirp_rate*(t1).^2));
pulse2=exp(1j*(2*pi*fI*t2+pi*chirp_rate*(t2).^2));
figure(1)
plot((0:N-1)*fs/N,abs(fft(pulse1)))
figure(2)
plot((0:N-1)*fs/N,abs(fft(pulse2)))
%% 二维搜索估计调频率和初始频率
            % 峰值搜索
xk=pulse1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xk=pulse1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pp=0.5:0.01:1.5;
u=1:1:N;
PP=zeros(length(u),length(pp));
k=0;
for ii=1:1:length(pp)
         p3 = myfrft(xk,pp(ii));
        %p3 = Fr_FT(xk,pp(ii));
        %p3 =dmyfrft(xk,pp(ii)*pi/2,N,N,Ts);
        PP(:,ii) = p3.*conj(p3);
end
[CX1,CP]=max(max(PP));
%[CX2,CU]=max(max(PP'));
great_p=pp(CP);%最优阶次
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pp=great_p-0.01:0.0005:great_p+0.01;
u=1:1:N;
PP=zeros(length(u),length(pp));
for ii=1:1:length(pp)
         p3 = myfrft(xk,pp(ii));
         %p3 = Fr_FT(xk,pp(ii));
       %p3 =dmyfrft(xk,pp(ii)*pi/2,N,N,Ts);
        PP(:,ii) = p3.*conj(p3);
end
[CX1,CP]=max(max(PP));
[CX2,CU]=max(max(PP'));
great_p=pp(CP);%最优阶次

great_u=u(round(CU));
kd_s=-fs/Tc*cot(pi*great_p/2)%调频率估计值
%u0 = (great_u-(N+1)/2)/sqrt(N);                    %？？？
%f_s=sqrt(fs/Tc)*u0*csc(great_p*pi/2)%中心频率估计值
f_s1=(great_u-N/2)*fs/N*csc(great_p*pi/2)
%% 解线调的方法估计中心频率
figure(3)
pulse1_jie=pulse1.*exp(-1j*pi.*kd_s.*t1.^2);
[a,b]=max(fft(pulse1_jie,512*512));
plot((0:N-1)*fs/N,abs(fft(pulse1_jie)));
f_s_1=b*fs/(512*512)

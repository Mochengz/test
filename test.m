close all
a=0:1:4;%分数阶傅里叶变换阶数

% 生成一个窗函数
% fx=zeros(2048,1);
% fx(512:1548)=1;
% 正弦函数
t=(1:N)/fs;
fx=sin(2*pi*t*50);

%归一化
N =2048;
fs = 2048;
T = N/fs;
S = (T/fs)^(1/2);
delta_x = (T*fs)^(1/2);
delta_t = 1/delta_x;
fx=fx'*S;
u = linspace(-1/2*delta_x,1/2*delta_x,N);

for ai=a
    figure
    F=myfrft(fx',ai);
    plot(u,abs(F))
    title("a="+num2str(ai))
    grid on
    ylim([0,5])
end

% figure
% 
% x_f = linspace(0,fs,N);
% 
% plot(x_f,abs(fft(fx))/N);
% S^(1/2)
% 
% clc
% clear
% fs=100;N=1000;
% n=0:N-1;t=n/fs;
% x=sin(10*t)
% y1=fft(x,N)
% y1 = abs(fftshift(y1)/N)
% plot(n,y1);
% figure
% y=myfrft(x,1)
% plot(n,abs(y))
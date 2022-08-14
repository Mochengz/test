%% 利用FrFT(Ozakats算法)对LFM信号进行检测和参数估计（起始频率和调频率）
%% 只有两级搜索（粗搜索和精细搜索）
clc
clear all
close all

N=1024;
fs=1024;
fc=250;
k=5;
t=(0:N-1)/fs;
sn=exp(1i*pi*k*t.^2+1i*2*pi*fc*t);
snr=10;
sn=awgn(sn,snr);
sn = sn.';
pp=mod(2*acot(-k)/pi,2);
%% 粗搜索
dpc=0.2;       % 粗搜索间隔
s1=2/dpc+1;
Smax=zeros(s1,1);
i=1;
pb=0:dpc:2;    % 要搜索的所有阶次
for p= pb    % 阶次具有周期性0<=p<=2
    S_frft = myfrft(sn, p);
    Smax(i)=max(abs(S_frft));
    i=i+1;
    %figure;plot(fftshift(abs(S_frft)))
end
[index,~]=find(Smax==max(Smax));
pg=pb(index);
%% 精细搜索
dpx=0.002;     % 精细搜索间隔
% s2=1/dpx+1;
% S_frft=zeros(N,s2);
px=pg-dpc:dpx:pg+dpc;
j=1;
for p=px
    S_frft(:,j) = abs(myfrft(sn,p));
    j=j+1;
end
[row,col]=find(S_frft==max(max(S_frft)));  % 最大值所在的行和列
Smax=max(max(S_frft));    % 最大值
pps=px(col);      % 归一化后的匹配阶次


ks=-cot(pps*pi/2);      % 归一化后的调频率估计值
Td=N/fs;
mu=ks*fs/Td;            % 真实的调频率估计值
pzs=mod(2*acot(-mu)/pi,2);      % 真实的匹配阶次

%%%%%%%%%%%%%%%%%% 载频 %%%%%%%%%%%%%%%%%%%%

%f=(0:N-1)*fs/N;   % 频域采样序列
f=(-N/2:N/2-1)*fs/N;
u=f(row);
fh=u*csc(pps*pi/2)     % 中心频率估计值
f0=fh-mu*Td/2          % 起始频率


[X,Y]= meshgrid(px,f);
figure;mesh(X,Y,S_frft)


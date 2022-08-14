%% ����FrFT(Ozakats�㷨)��LFM�źŽ��м��Ͳ������ƣ���ʼƵ�ʺ͵�Ƶ�ʣ�
%% ֻ�������������������;�ϸ������
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
%% ������
dpc=0.2;       % ���������
s1=2/dpc+1;
Smax=zeros(s1,1);
i=1;
pb=0:dpc:2;    % Ҫ���������н״�
for p= pb    % �״ξ���������0<=p<=2
    S_frft = myfrft(sn, p);
    Smax(i)=max(abs(S_frft));
    i=i+1;
    %figure;plot(fftshift(abs(S_frft)))
end
[index,~]=find(Smax==max(Smax));
pg=pb(index);
%% ��ϸ����
dpx=0.002;     % ��ϸ�������
% s2=1/dpx+1;
% S_frft=zeros(N,s2);
px=pg-dpc:dpx:pg+dpc;
j=1;
for p=px
    S_frft(:,j) = abs(myfrft(sn,p));
    j=j+1;
end
[row,col]=find(S_frft==max(max(S_frft)));  % ���ֵ���ڵ��к���
Smax=max(max(S_frft));    % ���ֵ
pps=px(col);      % ��һ�����ƥ��״�


ks=-cot(pps*pi/2);      % ��һ����ĵ�Ƶ�ʹ���ֵ
Td=N/fs;
mu=ks*fs/Td;            % ��ʵ�ĵ�Ƶ�ʹ���ֵ
pzs=mod(2*acot(-mu)/pi,2);      % ��ʵ��ƥ��״�

%%%%%%%%%%%%%%%%%% ��Ƶ %%%%%%%%%%%%%%%%%%%%

%f=(0:N-1)*fs/N;   % Ƶ���������
f=(-N/2:N/2-1)*fs/N;
u=f(row);
fh=u*csc(pps*pi/2)     % ����Ƶ�ʹ���ֵ
f0=fh-mu*Td/2          % ��ʼƵ��


[X,Y]= meshgrid(px,f);
figure;mesh(X,Y,S_frft)


function Faf = myfrft(f, a)
%分数阶傅里叶变换函数
%输入参数f为原始信号，a为阶数
%输出结果为原始信号的a阶傅里叶变换
N = length(f);%总采样点数
shft = rem((0:N-1)+fix(N/2),N)+1;
sN = sqrt(N);
a = mod(a,4);%按周期性将a定义在[0,4]

%特殊情况直接处理
if (a==0), Faf = f; return; end%自身
if (a==2), Faf = flipud(f); return; end%f(-x)
if (a==1), Faf(shft,1) = fft(f(shft))/sN; return; end%f的傅里叶变换
if (a==3), Faf(shft,1) = ifft(f(shft))*sN; return; end%f的逆傅里叶变换

%利用叠加性将阶数变换到0.5 < a < 1.5
if (a>2.0), a = a-2; f = flipud(f); end%a=2是反转
if (a>1.5), a = a-1; f(shft,1) = fft(f(shft))/sN; end%a=1是傅里叶变换
if (a<0.5), a = a+1; f(shft,1) = ifft(f(shft))*sN; end%a=-1是逆傅里叶变换

%开始正式的变换
alpha = a*pi/2;
tana2 = tan(alpha/2);
sina = sin(alpha);

f = [zeros(N-1,1) ; interp(f) ; zeros(N-1,1)];%使用香农插值，拓展为4N
% 对应论文中公式（29）
% 线性调频预调制
chrp = exp(-1i*pi/N*tana2/4*(-2*N+2:2*N-2)'.^2);  % N = delta ^ 2  且已经延拓为4N
f = chrp.*f;
% 线性调频卷积
c = pi/N/sina/4;
Faf = fconv(exp(1i*c*(-(4*N-4):4*N-4)'.^2),f);   % k∈[-2N,2N] l∈[-2N,2N] k-l∈[-4N,4N]  一般情况下，当x(n)及h（n）的"长度"（离散值的个数）分别为N1及N2时，卷积y（n）的长度则为N1+N2-1. 但是fft取了2的整数次方个点数 补零后去零
Faf = Faf(4*N-3:8*N-7)*sqrt(c/pi); %取卷积结果的中间1/3总长度的数据 总长为12N  sqrt(c/pi)为最后式子的分母 2*delta*sin(fai)^1/2
% 线性调频后调制
Faf = chrp.*Faf;
% 乘以最前面的A_Phi项
Faf = exp(-1i*(1-a)*pi/4)*Faf(N:2:end-N+1); %从2N点数抽取为N个点数输出
end

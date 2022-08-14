function z = fconv(x,y)%利用fft快速计算卷积
N = length([x(:);y(:)])-1;%计算最大点数
P = 2^nextpow2(N);%补零
z = ifft( fft(x,P) .* fft(y,P));%频域相乘，时域卷积
z = z(1:N);%去零
end
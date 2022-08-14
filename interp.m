function xint=interp(x)%香农插值
% sinc interpolation
N = length(x);
y = zeros(2*N-1,1);
y(1:2:2*N-1) = x;
xint = fconv(y(1:2*N-1), sinc([-(2*N-3):(2*N-3)]'/2));%计算卷积
xint = xint(2*N-2:end-2*N+3);
end

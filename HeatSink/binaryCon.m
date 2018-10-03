function [B, dB] = binaryCon(x, n)

B = 0;
dB = zeros(n, 1);
for i = 1:n
    B = B + sin(pi*x(i))/n;
    dB(i) = pi*cos(pi*x(i))/n;
end
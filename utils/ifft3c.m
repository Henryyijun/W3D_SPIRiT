function res = ifft3c(x)

S = size(x);
fctr = S(1)*S(2)*S(3);

x = reshape(x,S(1),S(2),prod(S(3:end)));

res = zeros(size(x));
res = sqrt(fctr)*fftshift(ifftn(ifftshift(x)));

res = reshape(res,S);




function res = fft2c(x)

S = size(x);
fctr = S(1)*S(2);

x = reshape(x,S(1),S(2),prod(S(3:end)));

res = zeros(size(x));
global gpu_enable;
if gpu_enable
    res = gpuArray(res);
end
for n=1:size(x,3)
	res(:,:,n) = 1/sqrt(fctr)*fftshift(fft2(ifftshift(x(:,:,n))));
end

res = reshape(res,S);




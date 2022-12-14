clear
clc
c = randn(3) + randn(3)*1i
cr = real(c); % a
ci = imag(c); % a

ifc = ifft2(c);

% ifcr = ifft2(cr);
% ifci = ifft2(ci);
% ifcrr = real(ifcr);
% ifcri = imag(ifcr);
% ifcir = real(ifci);
% ifcii = imag(ifci);
% yr = ifcrr-ifcii;
% yi = ifcri+ifcir;

ft_yr = fft2(real(ifc));
ft_yi = fft2(imag(ifc));

x_r = real(ft_yr) - imag(ft_yi);
x_i = imag(ft_yr) + real(ft_yi);
c2 = complex(x_r, x_i)
% ifcr+ifci*1i
% issame(ifc, ifcr+ifci*1i)
% 
% a = real(fft2(ifcrr)) - imag(fft2(ifcri));
% b = real(fft2(ifcir)) - imag(fft2(ifcii));
% fft2(ifc);
% result = a+b*1i;
% norm(fft2(ifc)-result)
% issame(fft2(ifc), result);


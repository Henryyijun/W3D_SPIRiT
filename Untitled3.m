clear
clc
c = randn(4) + randn(4)*1i

Bc = B(c)
BtBc = Bt(Bc)

norm(abs(c-BtBc))

function  [res] = B(x) 
  res = IFFT2_3D_N(x);
end

function [res] = Bt(x)
  re_x = real(x);
  im_x = imag(x);

  t1 = FFT2_3D_N(re_x);
  t2 = FFT2_3D_N(im_x);
  a = real(t1) - imag(t2);
  b = imag(t1) + real(t2);
  res = complex(a, b);
end
   
function [u] = SPIRiT_PD3O_Real2(g, mask, lev, lambda, max_iter, ker, ker_tra)

%%

un_mask = ~mask;

% z=WF^{-1}g
[z]= B(g, 1, lev);

A = @(x) VidHaarDec3S(IFFT2_3D_N(un_mask .* x), lev);
At = @(x) un_mask .* FFT2_3D_N(VidHaarRec3S(x, lev));
c = size(g,3);
%%
Lip =  Est_Lip_Ker_Mat_C(g, ker, [size(ker, 1), size(ker, 2)]);
gamma = 1.99 / (9*(1+Lip)^2);% (1+Lip)^2;
delta = 0.99 / gamma;

uk = g;
sk = A(zeros(size(uk)));
WinSize = 3; 
AveKer = ones(WinSize)/(WinSize*WinSize);
Thr = zeros(size(g,1),size(g,2),size(g,3),lev*4);

%%
global gpu_enable
if gpu_enable == true
  uk = gpuArray(uk);
  g = gpuArray(g);
  sk = gpuArray(sk);
  AveKer = gpuArray(AveKer);
  Thr = gpuArray(Thr);
  z = gpuArray(z);
  un_mask = gpuArray(un_mask);
  ker = gpuArray(ker);
  ker_tra = gpuArray(ker_tra);
end
for k = 1:max_iter
%   Qu_g = uk;
  Qu_g = un_mask .* uk + g;
  C_I_Qu_g = Kernel_Rec_ks_C_I_Pro(Qu_g, ker, 1);
  gradient_f = Kernel_Rec_ks_C_I_Pro(C_I_Qu_g, ker_tra, un_mask(:,:,1));
  uk_gamma_gradient_f = uk - gamma * gradient_f;

  gamma_At_sk_uk_gamma_gradient_f = gamma*Bt(sk, un_mask, lev) - uk_gamma_gradient_f;
  [tp1] = B(gamma_At_sk_uk_gamma_gradient_f, un_mask, lev);
  xk = sk - delta*tp1;
  delta_xk_z = xk/delta + z;
 
  if k < 40 && mod(k,3) == 1
    y = VidHaarDec3S(IFFT2_3D_N(un_mask .* uk + g), lev);
    for l = 1:lev
      Sigma = imfilter(abs(y(:,:,:,[2,3,4,5]+5*(l-1))), AveKer,'replicate');
      Thr(:,:,:,[1,2,3,4]+4*(l-1)) =  lambda(1).*8.^(l-1)./Sigma;
      % Thr(:,:,:,[1,2,3,4]+4*(l-1)) = Sigma; 
    end
  end
  
  for l = 1:lev
     sk(:,:,:,[2,3,4,5]+5*(l-1)) = xk(:,:,:,[2,3,4,5]+5*(l-1)) - ...
       delta * wthresh(delta_xk_z(:,:,:,[2,3,4,5]+5*(l-1)), 's', Thr(:,:,:,[1,2,3,4]+4*(l-1))/delta ) + delta * z(:,:,:,[2,3,4,5]+5*(l-1));
  end
  uk = uk_gamma_gradient_f - gamma*Bt(sk, un_mask, lev);
  
  figure(22)
  imshow(sos(IFFT2_3D_N(un_mask .*uk + g)), []);
  fprintf("At iteration %d\n", k);
end

u = sos(IFFT2_3D_N(un_mask .*uk + g));
uk = gather(uk);
u = gather(u);

function  [res] = B(x, mask, lev) 
  res = VidHaarDec3S(IFFT2_3D_N(mask .* x), lev);
end 

function [res] = Bt(x, mask, lev)
  re_x = real(x);
  im_x = imag(x);

  t1 = mask .* FFT2_3D_N(VidHaarRec3S(re_x, lev));
  t2 = mask .* FFT2_3D_N(VidHaarRec3S(im_x, lev));
  a = real(t1) - imag(t2);
  b = imag(t1) + real(t2);

  res = complex(a, b);
end
    

end

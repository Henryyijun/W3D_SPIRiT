function [u] = SPIRiT_PD3O_Real(g, mask, lev, lambda, max_iter, ker, ker_tra)

%%

un_mask = ~mask;

% z=WF^{-1}g
[z1, z2]= B(g, 1, lev); 

A = @(x) VidHaarDec3S(IFFT2_3D_N(un_mask .* x), lev);
At = @(x) un_mask .* FFT2_3D_N(VidHaarRec3S(x, lev));
%%
Lip =  Est_Lip_Ker_Mat_C(g, ker, [size(ker, 1), size(ker, 2)]);
gamma = 1.99 / ((1+Lip)^2);% (1+Lip)^2;
delta = 0.5;
uk = g;
sk = A(zeros(size(uk)));
sk1 = sk;
sk2 = sk;
WinSize = 3; 
AveKer = ones(WinSize)/(WinSize*WinSize);
Thr1 = zeros(size(g,1),size(g,2),size(g,3),lev*4);
Thr2 = zeros(size(g,1),size(g,2),size(g,3),lev*4);
%%

for k = 1:max_iter
  Qu_g = un_mask .* uk + g;
  C_I_Qu_g = Kernel_Rec_ks_C_I_Pro(Qu_g, ker, 1);
  gradient_f = Kernel_Rec_ks_C_I_Pro(C_I_Qu_g, ker_tra, un_mask(:,:,1));
  uk_gamma_gradient_f = uk - gamma * gradient_f;

  gamma_At_sk_uk_gamma_gradient_f = gamma*Bt(sk1, sk2, un_mask, lev) - uk_gamma_gradient_f;
  [tp1, tp2] = B(gamma_At_sk_uk_gamma_gradient_f, un_mask, lev);
  xk1 = sk1 - delta*tp1;
  xk2 = sk2 - delta*tp2;
  
  y1 = VidHaarDec3S(IFFT2_3D_N(un_mask .* real(uk) + real(g)), lev);
  y2 = VidHaarDec3S(IFFT2_3D_N(un_mask .* imag(uk) + imag(g)), lev);
  delta_xk_z1 = xk1/delta + (z1);
  delta_xk_z2 = xk2/delta + (z2);
    
%   if k < 40 && mod(k,8) ==1
%     for l = 1:lev
%       Sigma = imfilter(abs(y1(:,:,:,[2,3,4,5]+5*(l-1))), AveKer,'replicate');
%       Sigma =  lambda(1)*8^(l-1)/Sigma;
%       Thr1(:,:,:,[1,2,3,4]+4*(l-1)) = Sigma; 
%       
%       Sigma = imfilter(abs(y2(:,:,:,[2,3,4,5]+5*(l-1))), AveKer,'replicate');
%       Sigma =  lambda(1)*8^(l-1)/Sigma;
%       Thr2(:,:,:,[1,2,3,4]+4*(l-1)) = Sigma; 
%     end
%   end
%   
%   for l = 1:lev
%      sk1(:,:,:,[2,3,4,5]+5*(l-1)) = xk1(:,:,:,[2,3,4,5]+5*(l-1)) - ...
%        delta * wthresh(delta_xk_z1(:,:,:,[2,3,4,5]+5*(l-1)), 's', Thr1(:,:,:,[1,2,3,4]+4*(l-1))/delta ) + delta * z1(:,:,:,[2,3,4,5]+5*(l-1));
%      
%      sk2(:,:,:,[2,3,4,5]+5*(l-1)) = xk2(:,:,:,[2,3,4,5]+5*(l-1)) - ...
%        delta * wthresh(delta_xk_z2(:,:,:,[2,3,4,5]+5*(l-1)), 's', Thr2(:,:,:,[1,2,3,4]+4*(l-1))/delta ) + delta * z2(:,:,:,[2,3,4,5]+5*(l-1));
%   end
  sk1 = xk1 - delta*wthresh(delta_xk_z1, 's', lambda(1)) + z1;
  sk2 = xk2 - delta*wthresh(delta_xk_z2, 's', lambda(1)) + z2;
  uk = uk_gamma_gradient_f - gamma*Bt(sk1, sk2, un_mask, lev);
  
  figure(22)
  imshow(sos(IFFT2_3D_N(un_mask .*uk + g)), []);
  fprintf("At iteration %d\n", k);
end

u = sos(IFFT2_3D_N(un_mask .*uk + g));


function  [res1, res2] = B(x, mask, lev) 
%   re_x = real(x); 
%   im_x = imag(x);
%   res1 = VidHaarDec3S(IFFT2_3D_N(mask .* re_x), lev); % real 
%   res2 = VidHaarDec3S(IFFT2_3D_N(mask .* im_x), lev); % imag
res1 = VidHaarDec3S(IFFT2_3D_N(mask .* re_x), lev);
end

function [res] = Bt(x1, x2, mask, lev)
%   re_x1 = real(x1);
%   im_x1 = imag(x1);
%   re_x2 = real(x2);
%   im_x2 = imag(x2);
  
  t1 = real(mask .* FFT2_3D_N(VidHaarRec3S(x1, lev)));
  t2 = real(mask .* FFT2_3D_N(VidHaarRec3S(x2, lev)));
%   t1 = (mask .* FFT2_3D_N(VidHaarRec3S(x1, lev)));
%   t2 = (mask .* FFT2_3D_N(VidHaarRec3S(x2, lev)));
%   t3 = mask .* FFT2_3D_N(VidHaarRec3S(re_x2, lev));
%   t4 = mask .* FFT2_3D_N(VidHaarRec3S(im_x2, lev));
%   res1 = real(t1) - imag(t2);
%   res2 = real(t3) - imag(t4);
%   
  res = t1 + t2*1i;
end
    

end

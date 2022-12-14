function [u] = PD3O(g, mask, lev, lambda, max_iter, ker, ker_tra)


%%
artifact_image = IFFT2_3D_N(g);
un_mask = ~mask;

% z=WF^{-1}g
z = VidHaarDec3S(artifact_image, lev);

% A = WF^{-1}Q
A = @(x) VidHaarDec3S(IFFT2_3D_N(un_mask .* x), lev);
At = @(x) un_mask .* FFT2_3D_N(VidHaarRec3S(x, lev));


%%
Lip =  Est_Lip_Ker_Mat_C(g, ker, [size(ker, 1), size(ker, 2)]);
gamma = 1.99/(15);
delta = 0.99 / gamma;

%%
uk = g;
sk = A(zeros(size(uk)));
WinSize = 3; 
AveKer = ones(WinSize)/(WinSize*WinSize);
Thr = zeros(size(g,1),size(g,2),size(g,3),lev*4);
for k = 1:max_iter
  Qu_g = un_mask .* uk + g;
  C_I_Qu_g = Kernel_Rec_ks_C_I_Pro(Qu_g, ker, 1);
  gradient_f = Kernel_Rec_ks_C_I_Pro(C_I_Qu_g, ker_tra, un_mask(:,:,1));
  
  uk_gamma_gradient_f = uk - gamma * gradient_f;
  
  gamma_At_sk_uk_gamma_gradient_f = gamma*At(sk) - uk_gamma_gradient_f;
  xk = sk - delta * A(gamma_At_sk_uk_gamma_gradient_f);
  
  y = VidHaarDec3S(IFFT2_3D_N(un_mask .* uk + g), lev);
  delta_xk_z = xk/delta + z;
  
    
  if k < 40 && mod(k,8) ==1
    for l = 1:lev
      Sigma = imfilter(abs(y(:,:,:,[2,3,4,5]+5*(l-1))), AveKer,'replicate');
      Sigma =  lambda(1)*8^(l-1)/Sigma;
      Thr(:,:,:,[1,2,3,4]+4*(l-1)) = Sigma; 
    end
  end
  abs_delta_xk_z = abs(delta_xk_z);
  for l = 1:lev
     sk(:,:,:,[2,3,4,5]+5*(l-1)) = xk(:,:,:,[2,3,4,5]+5*(l-1)) - ...
       delta  .* delta_xk_z(:,:,:,[2,3,4,5]+5*(l-1))./abs_delta_xk_z(:,:,:,[2,3,4,5]+5*(l-1)) ...
       .* (wthresh(delta_xk_z(:,:,:,[2,3,4,5]+5*(l-1)), 's', Thr(:,:,:,[1,2,3,4]+4*(l-1))/delta ) -... 
        delta * z(:,:,:,[2,3,4,5]+5*(l-1)));
  end
  
  uk = uk_gamma_gradient_f - gamma*At(sk);
  figure(22)
  imshow(SoS(ifft2_3D(un_mask .*uk + g)), []);
end

u = SoS(ifft2_3D(un_mask .*uk + g));
end


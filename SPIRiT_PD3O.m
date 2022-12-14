function [u] = SPIRiT_PD3O(g, mask, lev, lambda, max_iter, ker, ker_tra)

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
gamma = 1.99/(9*(1+Lip)^2);
delta = 0.99 / gamma;

%%
uk = g;
sk = A(zeros(size(uk)));
WinSize = 3; 
AveKer = ones(WinSize)/(WinSize*WinSize);
Thr = zeros(size(sk));

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
  Qu_g = un_mask .* uk + g;
  C_I_Qu_g = Kernel_Rec_ks_C_I_Pro(Qu_g, ker, 1);
  gradient_f = Kernel_Rec_ks_C_I_Pro(C_I_Qu_g, ker_tra, un_mask(:,:,1));
  
  uk_gamma_gradient_f = uk - gamma * gradient_f;
  
  gamma_At_sk_uk_gamma_gradient_f = gamma*At(sk) - uk_gamma_gradient_f;
  xk = sk - delta * A(gamma_At_sk_uk_gamma_gradient_f);
  
  % y = VidHaarDec3S(IFFT2_3D_N(un_mask .* uk + g), lev);
  delta_xk_z = xk + delta*z;

  if k < 40 && mod(k,3) == 1
    y = VidHaarDec3S(IFFT2_3D_N(un_mask .* uk + g), lev);
    for l = 1:lev
      Sigma = imfilter(abs(y(:,:,:,[2,3,4,5]+5*(l-1))), AveKer,'replicate');
      Sigma =  lambda(1)*8^(l-1)/Sigma;
      Thr(:,:,:,[2,3,4,5]+5*(l-1)) = Sigma; 
    end
  end
  
  sk = delta_xk_z - wthresh(delta_xk_z, 's', Thr);
    
  % if k < 40 && mod(k,3) ==1
  %   for l = 1:lev
  %     Sigma = imfilter(abs(y(:,:,:,[2,3,4,5]+5*(l-1))), AveKer,'replicate');
  %     Sigma =  lambda(1)*8^(l-1)/Sigma;
  %     Thr(:,:,:,[1,2,3,4]+4*(l-1)) = Sigma; 
  %   end
  % end
  
  % for l = 1:lev
  %    sk(:,:,:,[2,3,4,5]+5*(l-1)) = xk(:,:,:,[2,3,4,5]+5*(l-1)) - ...
  %      delta * wthresh(delta_xk_z(:,:,:,[2,3,4,5]+5*(l-1)), 's', Thr(:,:,:,[1,2,3,4]+4*(l-1))/delta ) + delta * z(:,:,:,[2,3,4,5]+5*(l-1));
  % end
  
  uk = uk_gamma_gradient_f - gamma*At(sk);
  figure(22)
  imshow(sos(IFFT2_3D_N(un_mask .*uk + g)), []);
  fprintf("At iteration %d\n", k);
end

u = sos(IFFT2_3D_N(un_mask .*uk + g));
end


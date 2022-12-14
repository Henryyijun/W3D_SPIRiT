function [ u ] = Spirit_3D_S( g,Mask,Lev,lambda,Rho,Iter_Out,Iter_In,...
    Ker, Ker_Tra)
%
%

% inital
ImDown = ifft2_3D(g);
Mask_Un = ~Mask;
W_IFFT_g = VidHaarDec3S(ImDown, Lev);
y = zeros(size(W_IFFT_g));
z = y;

Thr = zeros(size(g,1),size(g,2),size(g,3),Lev*4);

global gpu_enable
if gpu_enable == true
    ImDown = gpuArray(ImDown);
    Mask_Un = gpuArray(Mask_Un);
    Mask = gpuArray(Mask);
    W_IFFT_g = gpuArray(W_IFFT_g);
    y = gpuArray(y);
    z = gpuArray(z);
    Ker = gpuArray(Ker);
    Ker_Tra = gpuArray(Ker_Tra);
    Thr = gpuArray(Thr);
    g = gpuArray(g);
end

A_g = Kernel_Rec_ks_C_I_Pro(g, Ker,1);
Qo_At_A_g = Kernel_Rec_ks_C_I_Pro(A_g, Ker_Tra,Mask_Un(:,:,1));
u = g; 

for k = 1:Iter_Out
    disp(['Iter:',num2str(k)]);
    
%     if(k<6)
%         Iter_In = 2; 
%     else
%         Iter_In = 3; 
%     end
    u = Update_u(Mask,u,y,z,Lev,Qo_At_A_g,Rho,Ker, Ker_Tra, Iter_In);
     
    Qo_u = u.*Mask_Un;
    Im = ifft2_3D(Qo_u);
    Gu = VidHaarDec3S(Im, Lev);

  [ y,Thr ] = Update_y_Itr(Gu,z,Rho,lambda, Lev, W_IFFT_g,k,Thr) ;  
    
    z = Update_z(Gu,y,z,Rho, W_IFFT_g);

    % figure(198), imshow(SoS(ifft2_3D(Qo_u+g)), []);
    
end

end

function [ u ] = Update_u(Mask,u,y,z,Lev,Qo_At_A_g,Rho,Ker, Ker_Tra, Iter_In)
Mask_Un = ~Mask;
% FFT_Wt_z = fft2_3D(VidHaarRec3(z, Lev));
% FFT_Wt_y = fft2_3D(VidHaarRec3(y, Lev));
 FFT_Wt_z_rho_Wt_y = fft2_3D(VidHaarRec3S(z-Rho*y, Lev));

% A_g = Kernel_Rec_ks_C_I_Pro(g, Ker,1);
% Qo_At_A_g = Kernel_Rec_ks_C_I_Pro(A_g, Ker_Tra,Mask_Un(:,:,1));
target_b = Qo_At_A_g + FFT_Wt_z_rho_Wt_y.*Mask_Un;
% target_b = Qo_At_A_g + (FFT_Wt_z - Rho.*FFT_Wt_y).*Mask_Un;
% target_b = -target_b(:);
% CG
u = CG(@aprod, -target_b, u, Iter_In, Mask, Mask_Un, Rho,Ker, Ker_Tra);
% u = reshape(u, size(g));
end



function [ y,Thr ] = Update_y_Itr(Gu,z,Rho,lambda, Lev, c,Itr,Thr)
y = Gu+c+z./Rho;

WinSize = 3; 
AveKer = ones(WinSize)/(WinSize*WinSize);
global gpu_enable;
if gpu_enable == true
    AveKer = gpuArray(AveKer);
end
for Lev_Num = 1:Lev
    if(Itr<10&mod(Itr,2)==1)
     Sigma = imfilter(abs(y(:,:,:,[2,3,4,5]+5*(Lev_Num-1))), AveKer,'replicate');
     Sigma =  lambda(1)*8^(Lev_Num-1)/Rho/Sigma;
     Thr(:,:,:,[1,2,3,4]+4*(Lev_Num-1)) =Sigma; 
    end
       y(:,:,:,[2,3,4,5]+5*(Lev_Num-1)) = wthresh(y(:,:,:,[2,3,4,5]+5*(Lev_Num-1)), 's', Thr(:,:,:,[1,2,3,4]+4*(Lev_Num-1)) );
      
    
end

end




function [ z ] = Update_z(Gu,y,z,Rho, c)

z = z + Rho.*(Gu+c-y);

end


function [res] = aprod(x,Mask, Mask_Un, Rho,Ker, Ker_Tra)
% x = reshape(x, size(Mask));
Qo_u = x.*Mask_Un;
A_Qo_u = Kernel_Rec_ks_C_I_Pro(Qo_u, Ker, 1);
Qo_At_A_Qo_u = Kernel_Rec_ks_C_I_Pro(A_Qo_u, Ker_Tra, Mask_Un(:,:,1));
res = Qo_At_A_Qo_u + Rho.*Qo_u;
% res = res(:);
end

function [x0] = CG(Afun,b,x0,Iter,Mask, Mask_Un, Rho, Ker, Ker_Tra)
r0 = b - Afun(x0,Mask, Mask_Un, Rho, Ker, Ker_Tra);
p0 = r0;
r0t_r0 = r0(:)'*r0(:);

for k = 1:Iter
    Ap0 = Afun(p0,Mask, Mask_Un, Rho,Ker, Ker_Tra);
    a0 = r0t_r0/(p0(:)'*Ap0(:));
    x0 = x0 + a0*p0;
    r0 = r0 - a0*Ap0;
    r0t_r0_k = r0(:)'*r0(:);
    
    b0 = r0t_r0_k/r0t_r0;
    p0 = r0 + b0*p0;
    r0t_r0 = r0t_r0_k;
    
end
end

clear
close all
%% Read the data
addpath utils;
load data\05_t2_tse_tra_512_s33_3mm_13.mat
% ksFull = ksfull; % k-space data

% Energy normalization
ImData = ifft2_3D(ksfull);
% figure(98);
% imshow(abs(sqrt(sum(ImData.^2,3))),[]);
% ImData = ImageNormalization(ImData);
ImData = NormalizedCoeByEnergy(ImData);
ksFull = fft2_3D(ImData);
size_ksFull = size(ksFull);


%% Compute the Sensitivity
center_size = 20;
center_width = 256;
KsCenter = crop(ksFull, center_size, center_width, size(ksFull, 3));
SenCoi = Sensitivity_Compute(KsCenter, size_ksFull);
SenCoi_SOS = sqrt(sum(abs(SenCoi).^(2),3));

%% Downsample + noise
[ImgRow, ImgCol, ImgDim] = size(ksFull);
% Mask = imread('./data/CartesianMask512_0_33.png');  % Qp sampling matrix of k-space data
Mask = imread('./data/CartesianMask512_0_20_3.png');  % Qp sampling matrix of k-space data
Mask = imrotate(Mask, 90);
% Iden = ones(size(Mask));
Mask_x = ~Mask; % the unsampled matrix Qq

ImgTruth = SoS(ifft2_3D(ksFull));
figure('Name','underlying image');
imshow(abs(ImgTruth),[])

ksData = zeros(size(ksFull));
for i = 1: size(ksFull, 3)
    ksData(:, :, i) = Mask .* ksFull(:, :, i);
end
ImgData = ifft2_3D(ksData); % downsampled data image domain
% show the aliased image and coil images
Clipsos = sqrt(sum(abs(ImgData).^(2),3));
% figure('Name','Aliased image');
% imshow(abs(Clipsos),[])

%% Find the missing line index in ACS
Target_center_size = 512;
Mask_center = crop(Mask, Target_center_size, center_width);
Missing = find(Mask(:,1) == 0);
fprintf(['Size of missing line is.' num2str(size(Missing, 1)) '.\n']);
Fill_line_num = 30; % How many lines fill once

%% Generate the data
Level = 2; % decompose Level
g = ksData; % g is the sampled k-sapce data

%% Iteration
% Initial
uk = Clipsos; % u0 = u_sos
Q0_M_u_g = repmat(Mask_x,1,1,size(SenCoi, 4)) .* M_oper(uk, SenCoi) + g;
vk = VidHaarDec3S(ifft2_3D(Q0_M_u_g), Level);
beta_k = zeros(size(vk)); % bk = W3D*Ft(Q0_M_u_g) - v0
max_Itr = 40; % the maximum number of iterations
rho = 1;
lambda = 0.025;
Sense_lambda = 0.01;
Sigmma = zeros([size(vk),Level]);
Sigmma_update_num = 5;% how many iterations to update Sigmma
crop_num = center_size;
Sense_update_num = 8; % how many iterations to update Sensitivity
Sense_update_time = 0;
% make Mask_x into 3D
Mask_x = repmat(Mask_x, 1,1,size(SenCoi, 3));
Thr_Par = zeros([size(SenCoi),8]);
WinSize = 3;
AveKer = ones(WinSize)/(WinSize*WinSize);
% loop
for k = 1: max_Itr
    % 1.update u
    W3Dt_v_b = VidHaarRec3S(vk - beta_k / rho, Level);
    FD_W3D_v_b = fft2_3D(W3Dt_v_b);
    Q0t_FD_W3D_v_b = Mask_x .* FD_W3D_v_b;
    Mt_Q0t_FD_W3D_v_b_g_r = M_Trans_oper(real(rho * Q0t_FD_W3D_v_b + g), SenCoi);
    Mt_Q0t_FD_W3D_v_b_g_i = M_Trans_oper(imag(rho * Q0t_FD_W3D_v_b + g), SenCoi);
    uk = real(Mt_Q0t_FD_W3D_v_b_g_r) - imag(Mt_Q0t_FD_W3D_v_b_g_i); % real - image
    
    % Update Sense
    if ((mod(k, Sense_update_num) == 0) && (k >= Sense_update_num))
        Sense_update_time = Sense_update_time + 1;
        fprintf('Update Sensitivity.\n');
        % A. Create the k-space location that needs to be filled
        % fill pattern 1 -- Fill a few lines at a time
        kspace_fill_Position = GetPosition(Missing, Sense_update_time*Fill_line_num, [ImgRow, ImgCol]);
        % fill pattern 2 -- Fill all the missing area
%         kspace_fill_Position = GetPosition(Missing, size(Missing, 1), [ImgRow, ImgCol]);
        % B. Create Mask_Sense
        % type 1 -- Not fixed window
        Mask_Sense = Get_Mask_Sense(size(SenCoi), center_size + Sense_update_time*Fill_line_num, center_width);
        % type 2 -- fix window
%         Mask_Sense = Get_Mask_Sense(size(SenCoi), ImgRow/2, center_width);
        % type 3 -- full window
%         Mask_Sense = Get_Mask_Sense(size(SenCoi), ImgRow, ImgCol);
        % C. update the Sensitivity
        if(Sense_update_time <= 100)
            SenCoi = Sense_3D_model(uk, SenCoi, g, Mask_x, Sense_lambda);
        else
            [SenCoi,~] = Sensitivity_Update(uk, SenCoi, g, kspace_fill_Position, Mask_Sense);
        end
%        max_SenCoi = max(max(max(abs(SenCoi))))
%        min_SenCoi = min(min(min(abs(SenCoi))))
%         figure, imshow(abs(SenCoi(:,:,1)),[])
        figure(99)
        imshow(abs(uk),[])
    end

    % 2.update v
    Q0_M_u = Mask_x .* M_oper(uk, SenCoi);
    W3D_FDt_Q0_M_u_g = VidHaarDec3S(ifft2_3D(Q0_M_u + g), Level);
    vk = W3D_FDt_Q0_M_u_g + beta_k / rho; % update low-pass too
    if(k<15)
        if(mod(k,4)==1)
            for Lev = 1:Level
                Sigma_Loc = imfilter(abs(vk(:,:,:,(2:5)+5*(Lev-1))),AveKer,'circular');
                Thr_Par(:,:,:,(1:4)+4*(Lev-1)) = 8^(Lev-1)*lambda ./ Sigma_Loc;
            end
            fprintf('Update e.\n');
        end  
    end
    for Lev = 1:Level
        vk(:,:,:,(2:5)+5*(Lev-1)) = wthresh(vk(:,:,:,(2:5)+5*(Lev-1)), 's', Thr_Par(:,:,:,(1:4)+4*(Lev-1)));
    end

    % 3.update b
    beta_k = beta_k + rho * ( W3D_FDt_Q0_M_u_g - vk);
    fprintf(['Ieration is: ' num2str(abs(k)) '.\n']);
end
figure(98);
imshow(abs(uk), [], 'border','tight');

%%-----------------------------------------------------------------------%%
% imwrite(abs(uk/max(uk(:))),'.\Result\Fill and Window\Fill_Line_Win_Fix.png')
% imwrite(abs(SenCoi(:,:,1)),'.\Result\Fill and Window\Sensitivity\Sense_Fill_Line_Win_Fix.png');

%         save(['.\Result\Sense\mat\lambda_' num2str(lambda) '.mat'],'SenCoi')
%         save(['.\Result\Sense\u\Itr_' num2str(k) '.mat'],'uk')
%         save(['.\Result\Sense\g.mat'],'g');

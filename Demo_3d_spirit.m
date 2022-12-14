
clear
% -----------------load full k-space image---------------
Slice  = 17;%29 ;
FilName = ['./data/05_t2_tse_tra_512_s33_3mm_',num2str(Slice) '.mat'];
load(FilName);
           
% ------------process---------------------
ksfull = NormalizedCoeByEnergy(ksfull);
KsData = ksfull(:,:,:); %2:4:end
%%  Parameter settings 
R = [4];
ACS_Line = [48];
Mask_Type = 'unif'; %unif  random
Lev = 2;
Rho = 1;
% lambda = [0.014, 0.02];
lambda = [0.0520, 0.0053];
Iter_Out = 50;
Iter_In = 5;        


%% Preparation
% ---------------------get mask
switch Mask_Type
    case 'unif'
        % acs: acs data     AcsLoc: index of acs rows
        [~,AcsLoc] = Acs_idx(KsData,ACS_Line);
        % mask: 1 represents exit data
        [Mask,~] = FindMask(size(KsData),R,AcsLoc);
         SamRate = round(sum(Mask(:))/size(Mask,1)/size(Mask,2)/size(Mask,3)*100);  %Keep the smapling rate 
    case 'random'
      Mask = imread('./data/CartesianMask512_0_15_2.png');  %Get a sampling model of k-space data
 
        % get size of calibration area from mask
        %  SamRate = round(sum(Mask(:))/size(Mask,1)/size(Mask,2)*100);  %Keep the smapling rate 
         Mask = repmat(Mask', [1 1 size(KsData,3)]);
         Mask = logical(Mask);
        [CalibSize, ~] = getCalibSize_1D_Edt( (Mask(:,:,1)));
        ACS_Line = CalibSize(1);
    otherwise
        error('Please input "unif" or "random"');
end
                          
Mask_Un  = ~Mask;              % Unsample positions
% -----------------get spirit kernel
Ker_Size= [9 9];
[Ker,Ker_Tra]=Kernel_Estimation(KsData,Ker_Size,ACS_Line); %(:,100:412,:) (240:272,240:272,:)
Lip_C = Est_Lip_Ker_Mat_C(KsData, Ker, Ker_Size);
g= zeros(size(KsData));
g(Mask) = KsData(Mask);

%%
% ImData = ifft2_3D(KsData);
% ImData = NormalizedCoeByEnergy(ImData);
% g1 = fft2_3D(ImData);
g1 = KsData .* Mask;
% g1 = KsData .* Mask;
% g1 = NormalizedCoeByEnergy(g1);
Img = IFFT2_3D_N(g1);
ImgSoS = SoS(Img);
figure(98);
imshow(ImgSoS,[]);
lambda1 = 0.1;
global gpu_enable;
gpu_enable = true;
disp("begin");
tic
u1 = SPIRiT_PD3O(g1,Mask,Lev,lambda1,Iter_Out,Ker, Ker_Tra);
toc;
figure(900)
imshow(u1,[]);

res1 = u1 ./ max(u1(:));
imwrite(res1, "./result/PD3O.png");
% 
% 
% ImData = ifft2_3D(KsData);
% ImData = NormalizedCoeByEnergy(ImData);
% g2 = fft2_3D(ImData);
g2 = KsData .* Mask;
lambda2 = 0.55;
tic;
[u2] = Spirit_3D_S( g2,Mask,Lev,lambda2,Rho,Iter_Out,Iter_In,...
    Ker, Ker_Tra);
toc
res2 = sos(ifft2_3D(u2));
res2 = res2 ./ max(res2(:));
imwrite(res2, "./result/ADMM.png");
figure(901)
imshow(sos(ifft2_3D(u2)),[]);

% img_3d = SoS(IFFT2_3D_N(u.*Mask_Un+g));
% img_3d = img_3d/max(img_3d(:));
% imwrite(img_3d,Rec_Name);
% 
% figure(4);
% imshow(img_3d, []);




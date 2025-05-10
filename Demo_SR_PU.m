 clc
  clear all;
  close all;
warning('off');
sf    =   4;
addpath(genpath('data\'));
addpath(genpath('LFCTNN-main\'));
for data_num=1
clear img
switch data_num    
    case 1
         S=imread('.\data\original_rosis.tif');
         img=S(1:256,1:256,11:end);
        name= 'original_rosis';
    case 2
         load 'WDC.mat'
         img=WDC(:,:,1:103);
         name='DCmall';
     case 3
        load 'RemoteImage.mat'
        img=Hw;
        name='RemoteImage';
end  
img=double(img);
O_Img = img/max(img(:));
S= O_Img ;
[M N L]=size(S);
sz=[M N];
%%
load('.\data\R.mat');
F=R;
F=F(:,1:end-10);
for band = 1:size(F,1)
        div = sum(F(band,:));
        for i = 1:size(F,2)
            F(band,i) = F(band,i)/div;
        end
end
%%
Comparison_LFCTNN_FC = 1;
Comparison_LFCTNN_FL = 1;
%%
s0=1;
downsampling_scale = sf;
% psf        =    ones(sf)/(sf^2);
psf        =    fspecial('gaussian',7,2);
par.fft_B      =    psf2otf(psf,sz);
par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);
[M,N,L] = size(S);
S_bar = hyperConvert2D(S);
hyper= par.H(S_bar);
SNRh=30;
sigma = sqrt(sum(hyper(:).^2)/(10^(SNRh/10))/numel(hyper));
rng(10,'twister')
hyper = hyper+ sigma*randn(size(hyper));
HSI=hyperConvert3D(hyper,M/downsampling_scale, N/downsampling_scale );
rng(10,'twister')
Y = F*S_bar;
SNRm=35;
sigmam = sqrt(sum(Y(:).^2)/(10^(SNRm/10))/numel(Y));
Y = Y+ sigmam*randn(size(Y));
MSI=hyperConvert3D(Y,M,N);
case1 = 'uni';
case2 = num2str(SNRh);
case3 =  num2str(SNRm);
%%
if Comparison_LFCTNN_FC == 1;
ex14='LFCTNN_FC';
cd 'LFCTNN-main'
addpath(genpath(cd));
aa = 0;
para.K=400;
para.eta=[3e-2, 5e-3, 3e-2]; %% for LFCTNN-FC
para.p=10;
para.patchsize = 6;
t0=clock;
Z14 =  Gernalized_LFCTNN_FC(HSI,MSI,F, par.fft_B,sf,S,para,aa);
t15=etime(clock,t0);    
cd ..
[psnr14,rmse14, ergas14, sam14, uiqi14,ssim14,DD14,CC14] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z14)), 0, 1.0/sf);
disp(['PSNR_',ex14,'_',name ' = ',num2str(psnr14)]);
disp(['SSIM_',ex14,'_',name ' = ',num2str(ssim14)]);
disp(['ergas_',ex14,'_',name ' = ',num2str(ergas14)]);
disp(['sam_',ex14,'_',name ' = ',num2str(sam14)]);
disp(['uiqi_',ex14,'_',name ' = ',num2str(uiqi14)]);
disp(['DD_',ex14,'_',name ' = ',num2str(DD14)]);
disp(['CC_',ex14,'_',name ' = ',num2str(CC14)]);
% save(['result_' name '_' ex14 '_' case1 '_' case2 '_' case3  '.mat'],'psnr14','ssim14','ergas14','sam14','uiqi14','DD14','CC14','Z14');
end

%% 15
if Comparison_LFCTNN_FL == 1;
ex15='LFCTNN_FL';
cd 'LFCTNN-main'
addpath(genpath(cd));
aa = 0;
para.K=400;
para.eta=[3e-2, 1e-2, 3e-2]; %% for LFCTNN-FL
para.p=10;
para.patchsize = 6;
t0=clock;
Z15 =  Gernalized_LFCTNN_FL(HSI,MSI,F, par.fft_B,sf,S,para,aa);
t15=etime(clock,t0);    
cd ..
[psnr15,rmse15, ergas15, sam15, uiqi15,ssim15,DD15,CC15] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z15)), 0, 1.0/sf);
disp(['PSNR_',ex15,'_',name ' = ',num2str(psnr15)]);
disp(['SSIM_',ex15,'_',name ' = ',num2str(ssim15)]);
disp(['ergas_',ex15,'_',name ' = ',num2str(ergas15)]);
disp(['sam_',ex15,'_',name ' = ',num2str(sam15)]);
disp(['uiqi_',ex15,'_',name ' = ',num2str(uiqi15)]);
disp(['DD_',ex15,'_',name ' = ',num2str(DD15)]);
disp(['CC_',ex15,'_',name ' = ',num2str(CC15)]);
% save(['result_' name '_' ex15 '_' case1 '_' case2 '_' case3  '.mat'],'psnr15','ssim15','ergas15','sam15','uiqi15','DD15','CC15','Z15');
    end
end  
clc
clear all;
close all;
warning('off');
sf              =   16;
sf_str = num2str(sf);
addpath(genpath('data\'));
addpath(genpath(cd));
pathstr = '.\data\cave';
dirname = fullfile(pathstr, '*.mat');
imglist = dir(dirname);
index = zeros(length(imglist), 5);
for yy = 1:length(imglist)
    img = load(fullfile(pathstr, imglist(yy).name));
    img = img.label;
    name='cave';
    ii_str = num2str(yy); % ½« ii ×ª»»Îª×Ö·û´®
    O_ImgUint8=double(img);
    O_Img = O_ImgUint8/max(O_ImgUint8(:));
    S= O_Img ;
    [M N L]=size(S);
    sz=[M N];
%%
F=create_F();
%%
Comparison_GTNN =1;
s0=1;
downsampling_scale = sf;
psf        =    fspecial('gaussian',7,2);
par.fft_B      =    psf2otf(psf,sz);
par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);
par.P=create_F();
[M,N,L] = size(S);
S_bar = hyperConvert2D(S);
hyper= par.H(S_bar);
SNRh=20;
sigma = sqrt(sum(hyper(:).^2)/(10^(SNRh/10))/numel(hyper));
rng(10,'twister')
hyper = hyper+ sigma*randn(size(hyper));
HSI=hyperConvert3D(hyper,M/downsampling_scale, N/downsampling_scale );
rng(10,'twister')
Y = F*S_bar;
SNRm=20;
sigmam = sqrt(sum(Y(:).^2)/(10^(SNRm/10))/numel(Y));
Y = Y+ sigmam*randn(size(Y));
MSI=hyperConvert3D(Y,M,N);


%% 15
if Comparison_GTNN== 1;
 ex15='LFCTNN';
cd 'LFCTNN-main'
addpath(genpath(cd));
aa = 0;
for jj = [1] %%  2 for LFCTNN-FL;  1 for LFCTNN-FC
jj_str = num2str(jj); 
para=[];
para.K=400;

 for e1=[1] %% 3 for LFCTNN-FL;  1 for LFCTNN-FC
para.eta=[e1,e1,e1];
para.p=10;
iii = 1;
t0=clock;
 [Z15,ZZZZZ] =  Gernalized_LFCTNN_CAVE(HSI,MSI,F, par.fft_B,sf,S,para,jj,aa,iii);
  t15=etime(clock,t0)    
cd ..
  [psnr15,rmse15, ergas15, sam15, uiqi15,ssim15,DD15,CC15] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z15)), 0, 1.0/sf);
disp(['PSNR_',ex15,'_',name ' = ',num2str(psnr15)]);
disp(['SSIM_',ex15,'_',name ' = ',num2str(ssim15)]);
disp(['ergas_',ex15,'_',name ' = ',num2str(ergas15)]);
disp(['sam_',ex15,'_',name ' = ',num2str(sam15)]);
disp(['uiqi_',ex15,'_',name ' = ',num2str(uiqi15)]);
disp(['DD_',ex15,'_',name ' = ',num2str(DD15)]);
disp(['CC_',ex15,'_',name ' = ',num2str(CC15)]);
save(['result_' name '_' ex15 '_' ii_str '_' sf_str '_' jj_str '.mat'],'psnr15','ssim15','ergas15','sam15','uiqi15','DD15','CC15','Z15');
aa = aa +1;
mpsnr15(yy)=psnr15;
 end
        end
    end 
end
% end
% % mpsnr15(yy)=psnr15;
% end
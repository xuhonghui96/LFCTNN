function  [Z,ZZZZZ] =  Gernalized_TSVD_Subpace_FUS( HSI, MSI,R,FBm,sf,S,para,ii,aaa,iii)
K=para.K;
eta=para.eta(1);
eta2=para.eta(2);
eta3=para.eta(3);
p=para.p;
mu=1e-3;
mu0=1e-3;
patchsize=8;
overlap=4;

HSI3=Unfold(HSI,size(HSI),3);
[D,~]=svds(HSI3,p);
RD=R*D;
L1=size(D,2);

[~,~,D2]=svds(HSI3,1);
D2=reshape(D2,[size(HSI,1),size(HSI,2)]);
D2=D2/(max(D2(:))-min(D2(:)))-min(D2(:))/(max(D2(:))-min(D2(:)));
D2=imresize(D2,sf);
path='G:\2024ÊµÑématlab\para\';
filename = [path 'cave-ctnn4-yy2-3_10.xls']; 

bparams.block_sz = [patchsize, patchsize];
bparams.overlap_sz=[overlap overlap];
[nr, nc,~]=size(MSI);
L=size(HSI,3);
num1=(nr-patchsize)/(patchsize-overlap)+1;
num2=(nc-patchsize)/(patchsize-overlap)+1;
 bparams.block_num=[num1 num2]; 
 sizeag=size(MSI,3);
MSIG(:,:,1:sizeag)=MSI;
MSIG(:,:,sizeag+1)=D2;
predenoised_blocks = ExtractBlocks(MSIG, bparams);
Y2=Unfold(predenoised_blocks,size(predenoised_blocks),4);
if K==1
    aa=ones(num1*num2,1);
else
  [aa ]=fkmeans(Y2,K);
end

HSI_int=zeros(nr,nc,L);
HSI_int(1:sf:end,1:sf:end,:)=HSI;
FBmC  = conj(FBm);
FBs  = repmat(FBm,[1 1 L]);
FBs1  = repmat(FBm,[1 1 L1]);
FBCs=repmat(FBmC,[1 1 L]);
FBCs1=repmat(FBmC,[1 1 L1]);
HHH=ifft2((fft2(HSI_int).*FBCs));
HHH1=hyperConvert2D(HHH);




%% iteration

MSI3=Unfold(MSI,size(MSI),3);

n_dr=nr/sf;
n_dc=nc/sf;

HR_load1=imresize(HSI, sf,'bicubic');

V1=D'*hyperConvert2D(HR_load1);
V2=V1;
V3=V1;



G1=zeros(size(V2));
G2=G1;
G3=G1;
DTD=D'*D;
CCC=DTD\(RD'*MSI3+D'*HHH1);
 C1=DTD\(RD'*RD+3*mu*eye(size(D,2))); 
 [Q,Lambda]=eig(C1);
Lambda=reshape(diag(Lambda),[1 1 L1]);
InvLbd=1./repmat(Lambda,[ sf*n_dr  sf*n_dc 1]);
B2Sum=PPlus(abs(FBs1).^2./( sf^2),n_dr,n_dc);
InvDI=1./(B2Sum(1:n_dr,1:n_dc,:)+repmat(Lambda,[n_dr n_dc 1]));



for i=1:30
%     HR_HSI3=mu*(V1+G1/(mu)+V2+G2/(mu)+V3+G3/(mu));
    HR_HSI3=mu*(V1+G1/(2*mu)+V2+G2/(2*mu)+V3+G3/(2*mu));
C3=CCC+DTD\HR_HSI3;
C30=fft2(reshape((Q\C3)',[nr nc L1   ])).*InvLbd;
temp  = PPlus_s(C30/( sf^2).*FBs1,n_dr,n_dc);
invQUF = C30-repmat(temp.*InvDI,[ sf  sf 1]).*FBCs1; % The operation: C5bar- temp*(\lambda_j d Im+\Sum_i=1^d Di^2)^{-1}Dv^H)
VXF    = Q*reshape(invQUF,[nc*nc L1])';
A = reshape(real(ifft2(reshape(VXF',[nr nc L1   ]))),[nc*nc L1])'; 
%   [ZE1] = Sylvester(C1,psfY.B, sf,n_dr,n_dc,C3);  
    
Zt=hyperConvert3D(D*A,nr, nc );


%   rmse2(i)=getrmse(double(im2uint8(S)),double(im2uint8(Zt))) 
psnr(i)=csnr(double(im2uint8(S)),double(im2uint8(Zt)),0,0)
s=strcat('A',num2str(i+aaa*40));
%% spatial similarties

  

 B1=hyperConvert3D(A-G1/(2*mu),nr, nc );
 predenoised_blocks1 = ExtractBlocks(B1, bparams);
   Z_block1=zeros( bparams.block_sz(1), bparams.block_sz(2),L1, bparams.block_num(1)* bparams.block_num(2));
 
   B2=hyperConvert3D(A-G2/(2*mu),nr, nc );
 predenoised_blocks2 = ExtractBlocks(B2, bparams);
 Z_block2=zeros( bparams.block_sz(1), bparams.block_sz(2),L1, bparams.block_num(1)* bparams.block_num(2));
 
  B3=hyperConvert3D(A-G3/(2*mu),nr, nc );
 predenoised_blocks3 = ExtractBlocks(B3, bparams);
  Z_block3=zeros( bparams.block_sz(1), bparams.block_sz(2),L1, bparams.block_num(1)* bparams.block_num(2));


for mn=1:max(aa)
    gg=find(aa==mn);
 XES=predenoised_blocks1(:,:,:,gg);  
   [a, b, c, d ]=size(XES);
    XES = reshape(XES,[a*b c d]);
%       U=sqrt(d)*dct(eye(d));
switch ii
    case 1
      V1=Log_prox_tnn_FC(XES, eta/2/mu,iii);%_2  
    case 2
      V1=Log_prox_tnn_FL(XES, eta/2/mu,iii);%_2  
end


V1=reshape(V1,[a b c d]); 
Z_block1(:,:,:,gg)=V1;
 



   XES=predenoised_blocks2(:,:,:,gg);
   [a, b, c, d ]=size(XES);
    XES = reshape(XES,[a*b c d]);
    XES=permute(XES,[1,3,2]);
switch ii
    case 1
      V2=Log_prox_tnn_FC(XES, eta2/2/mu,iii);%_2  
    case 2
      V2=Log_prox_tnn_FL(XES, eta2/2/mu,iii);%_2  
end
       V2=permute(V2,[1,3,2]);
V2=reshape(V2,[a b c d]); 
  Z_block2(:,:,:,gg)=V2;
 
    XES=predenoised_blocks3(:,:,:,gg);
   [a, b, c, d ]=size(XES);
    XES = reshape(XES,[a*b c d]);
    XES=permute(XES,[2,3,1]);
switch ii
    case 1
      V3=Log_prox_tnn_FC(XES, eta3/2/mu,iii);%_2  
    case 2
      V3=Log_prox_tnn_FL(XES, eta3/2/mu,iii);%_2  
end
       V3=permute(V3,[3,1,2]);
V3=reshape(V3,[a b c d]); 
  Z_block3(:,:,:,gg)=V3;
end
    
V1= JointBlocks(Z_block1, bparams);
V1=hyperConvert2D(V1);
G1=G1+2*mu*(V1-A);

V2= JointBlocks(Z_block2, bparams);
V2=hyperConvert2D(V2);
G2=G2+2*mu*(V2-A);

V3= JointBlocks(Z_block3, bparams);
V3=hyperConvert2D(V3);
G3=G3+2*mu*(V3-A);
Z=hyperConvert3D(D*A,nr, nc );
ZZ{i}=Z;
if i>1
    ZZZ1 = ZZ{i};
    ZZZ2 = ZZ{i-1};
    ZZZZZ(i)=norm(ZZZ1(:) - ZZZ2(:)) / norm(ZZZ2(:));
end
end

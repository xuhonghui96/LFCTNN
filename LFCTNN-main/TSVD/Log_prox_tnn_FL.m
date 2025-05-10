function [X,tnn,trank] = Log_prox_tnn(Y,rho)

[n1,n2,n3] = size(Y);
X = zeros(n1,n2,n3);

Y = fft(Y,[],3);
if n3 == 1
    [UU,~,~] = svd(Y,'econ');
else 
    YO = Unfold(Y,[n1 n2 n3],3);
    [mm,nn]=size(YO);
    if mm > nn
       [UU,~,~] = svd(YO);
    else
       [UU,~,~] = svd(YO,'econ');
    end
    O = tenmat(Y,[3]); %square norm
    SS = O.data;
    Y = UU'*SS;
    Y = tensor(tenmat(Y, O.rdims, O.cdims, O.tsize));
    Y = Y.data;
end


tnn = 0;
trank = 0;

for i = 1 : n3
   X(:,:,i) = WNNM( Y(:,:,i), rho, 1);
end
tnn = tnn/n3;


if n3==1
    X = X;
else
    O = tenmat(X,[3]); %square norm
    SS = O.data;
    Y = UU*SS;
    Y = tensor(tenmat(Y, O.rdims, O.cdims, O.tsize));
    X = Y.data;
end
X = ifft(X,[],3);

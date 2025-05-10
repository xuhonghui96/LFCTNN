function [X,tnn,trank] = Log_prox_tnn(Y,rho)

[n1,n2,n3] = size(Y);
X = zeros(n1,n2,n3);
Y = fft(Y,[],3);
transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
Y = lineartransform(Y,transform);
tnn = 0;
trank = 0;

for i = 1 : n3
   X(:,:,i) = WNNM( Y(:,:,i), rho,1);
end
tnn = tnn/n3;

X = inverselineartransform(X,transform);
X = ifft(X,[],3);

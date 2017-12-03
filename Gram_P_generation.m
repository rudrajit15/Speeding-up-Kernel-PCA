function [K,prob] = Gram_P_generation(n,m,d,gamma)
%images = loadMNISTImages('train-images-idx3-ubyte');
S = rand(n,d);
D = zeros(d,d);
for i = 1:d
    D(i,i) = 1-(i-1)/d;
end
Z = rand(m,m);
[U,~,~] = svd(Z);
U = U(:,1:d)';
N = rand(n,m);
%A = S*D*U+N/gamma;
A = S*D*U;
K = A'*A;
denom = trace(K'*K);
prob = zeros(m,1);
for i = 1:m
    prob(i) = (norm(K(:,i)))^2/denom;
end
function [eig_val,eig_vec,eig_val2,eig_vec2] = Nystrom2(G,c_id,p,k)
n = size(G,1);
c = size(c_id,1);
S = zeros(n,c);
D = zeros(c,c);
t0=1:c;
%p1 = p(c_id(t0));
%t1 = sqrt(1./(c*p(c_id(t0))));
%D1 = diag(t1);
%D = zeros(c,c);
D((t0-1)*c+t0) = sqrt(1./(c*p(c_id(t0))));
%for t = 1:c
%    D(t,t) = 1/sqrt(c*p(c_id(t)));
%    S(c_id(t),t) = 1;
%end
%disp(trace(D1-D));
S((t0-1)*n+c_id') = 1;

%t = 1:c;
%v = 1./sqrt(c*p(c_id(t)));
%D = diag(v);
%S((t-1)*n+c_id(t)) = 1;
C = G*S*D;
W = D*S'*C;

%[Uw,Sw,Vw] = svd(W);
%for i = 1:k
%        cu = Uw(:,i);
%        cv = Vw(:,i);
%        if(sum(cu.*cv) < 0)
%            Sw(i,i) = -1/Sw(i,i);
%        else
%            Sw(i,i) = 1/Sw(i,i);
%        end
%end

%Uw = Uw(:,1:k);
%Sw = Sw(1:k,1:k);

[Uw,Sw] = eigs(W,k);
%for i = 1:k
%   Sw(i,i) = 1/Sw(i,i);
%end
tem = find(Sw);
Sw(tem) = 1./Sw(tem);


W_plus = Uw*Sw*Uw';

t = C*Uw;
[Uc,Sc,Vc] = svd(t,'econ');
%Uc = Uc(:,1:k);
%Sc = Sc(1:k,1:k);
%Vc = Vc(:,1:k);
Dc = Sc*Vc'*Sw*Vc*Sc;

[Uf,Df] = eig(Dc);
[~,I] = sort(diag(Df),'descend');
Df = Df(:, I);
Uf = Uf(:, I);

eig_vec = Uc*Uf;
eig_val = zeros(k,1);
for i = 1:k
    eig_val(i) = Df(i,i);
end


G_k = C*W_plus*C';
%[eig_vec2,eig_val2] = eig(G_k);
[eig_vec2,eig_val2] = eig(G);
[eig_val2,I] = sort(diag(eig_val2),'descend');
eig_vec2 = eig_vec2(:, I);
%eig_val2 = eig_val2(1:k);
%eig_vec2 = eig_vec2(1:k,1:k);

eig_val_diff = zeros(k,1);
for i = 1:k
    eig_val_diff(i) = eig_val2(i) - eig_val(i);
end
eig_val_diff

eig_vec_diff = zeros(k,1);
for q = 1:k
    temp = abs(eig_vec2(:,q)) - abs(eig_vec(:,q));
    %eig_vec_diff(q,1) = sqrt((temp'*temp)/length(eig_vec2(:,q)));
    eig_vec_diff(q,1) = sqrt((temp'*temp));
end

%cross_correl = G_e_vec'*Nys_e_vec;
%imshow(cross_correl);
eig_vec_diff


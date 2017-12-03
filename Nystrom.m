function [eig_val,eig_vec] = Nystrom(G,c_id,p,k)
n = size(G,1);
c = size(c_id,1);
S = zeros(n,c);
D = zeros(c,c);
%for t = 1:c
%    D(t,t) = 1/sqrt(c*p(c_id(t)));
%    S(c_id(t),t) = 1;
%end
%t=1:c;
%D((t-1)*c+t) = 1/(sqrt(c*p(c_id(t))));
%S((t-1)*n+c_id(t)) = 1;
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


C = G*S*D;
W = D*S'*C;
[Uw,Sw] = eigs(W,k);
%for i = 1:k
%    if(Sw(i,i) ~= 0)
%        Sw(i,i) = 1/Sw(i,i);
%    end
%end
%tem = find(Sw);
%Sw(tem) = 1/Sw(tem);
tem = find(Sw);
Sw(tem) = 1./Sw(tem);

t = C*Uw;
[Uc,Sc,Vc] = svd(t,'econ');
Dc = Sc*Vc'*Sw*Vc*Sc;

[Uf,Df] = eig(Dc);
%[Df,I] = sort(diag(Df),'descend');
%Uf = Uf(:,I);

eig_vec = Uc*Uf;
eig_val = Df;
%function [eig_val,eig_vec,eig_val2,eig_vec2] = main(K,prob)
function [eig_val,eig_vec] = main(K,prob)
n = 200;
m = 6500;
d = 50;
gamma = 2;
epsilon = 1.5;
k = 50;
c = floor(64*k/(epsilon^4));
%c = 4000;
%[K,prob] = Gram_P_generation(n,m,d,gamma);
tic
[c_ind] = non_uniform_sampling(prob,c);
%[eig_val,eig_vec,eig_val2,eig_vec2] = Nystrom2(K,c_ind,prob,d);
[eig_val,eig_vec] = Nystrom(K,c_ind,prob,d);
%[eig_val,eig_vec] = eigs(K,d);
toc
end

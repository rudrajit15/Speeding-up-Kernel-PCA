function [c_ind] = non_uniform_sampling(probs,c)
n = numel(probs);
cdf = zeros(n,1);
cdf(1,1) = probs(1);
for i = 2:n
    cdf(i,1) = cdf(i-1,1) + probs(i);
end
%cdf1 = zeros(n,1);
%cdf = cumsum(probs);
%i0=2:n;
%cdf1(i0) = cdf1(i0-1)+cdf1(i0);
%disp(sum(cdf-cdf1));
%cdf(n,1)
c_ind = zeros(c,1);
sequence = rand(c,1);
for t = 1:c
     [~,ind] = min(abs(cdf-sequence(t,1)));
     if(cdf(ind) < sequence(t,1))
         ind = ind + 1;
     end
     c_ind(t,1) = ind;
end

    


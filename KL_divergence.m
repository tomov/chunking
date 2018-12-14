function surprise = KL_divergence(P, Q)
% Compute the KL divergence from Q to P i.e. D_KL(P || Q)
% = sum over i, 
% assumes each column corresponds to a different value of the RV
% if given matrices, treats each row as a separate pair of P, Q
%
% INPUT:
% P = posterior distributions, one per row
% Q = prior distributions, one per row
%
% OUTPUT:
% surprise = D_KLs, one per row
%

assert(isequal(size(P), size(Q)));

if ndims(P) == 3
    P = reshape(P, size(P,1)*size(P,2), size(P,3))';
    Q = reshape(Q, size(Q,1)*size(Q,2), size(Q,3))';
end

%{
assert(abs(mean(sum(P,2) - 1)) < 1e-6); 
assert(abs(mean(sum(Q,2) - 1)) < 1e-6);
assert(min(min(P)) >= 0);
assert(max(max(P)) <= 1);
assert(min(min(Q)) >= 0);
assert(max(max(Q)) <= 1);
%}

h = P .* (log2(P) - log2(Q)); 
h(isnan(h)) = 0; % lim_{x->0} x log(x) = 0
%assert(sum(sum(isinf(h))) == 0); % shouldn't be any inf's left

surprise = sum(h, 2);

%assert(sum(isinf(surprise)) == 0); % sanity
%assert(sum(isnan(surprise)) == 0); % sanity

end


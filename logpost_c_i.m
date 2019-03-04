% P(H|D) for updates of c_i
% i.e. with new c's up to c_i, the candidate c_i, then old c's after (and old rest of H)
%
function logp = logpost_c_i(c_i, i, H, D, h)
    H.c(i) = c_i;

    cnt = get_H_cnt(H, D);
    H.cnt = cnt;
    disp('H after update in matlab:');
    H
    logp = logpost(H, D, h);
end


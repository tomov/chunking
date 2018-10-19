function H = populate_H(H, D)

% fill up the missing deets in H 

H.cnt = get_H_cnt(H, D);
H.E = get_H_E(H, D);
H.b = get_H_b(H, D);

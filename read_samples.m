load count_expirement_1.mat

% 4 is connected to 5, 6, 7

n = nsubjects;
c1 = count_4_1;
c2 = count_4_7;

p = 2 * binocdf(min(c1,c2), n, 0.5);
y = binoinv([0.025 0.975], n, 0.5);
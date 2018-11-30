clear all;

sem = @(x) std(x) / sqrt(length(x));

N = 30; % participants
h.alpha = 5;
nwalks = 18; % how many random walks or hamiltonians for each subject (based on paper)

D = init_D_from_txt('schapiro.txt');

%hamils = hamiltonians(D);
load('schapiro_hamils.mat', 'hamils'); % pregenerated hamiltonian paths

% community transitions
comm_trans = [
3 6;
6 3;
4 11;
11 4;
10 12;
12 10];

for s = 1:N % for each simulated subject
    [H, P] = sample(D, h, 10);
    H_all{s} = H;
    P_all{s} = P;

    [~,I] = max(P); % MAP H
    H = H(I);

    paths = random_walks(nwalks);
end

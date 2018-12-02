% simulate experiment 1 from schapiro 2013

clear all;

rng default;

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
    fprintf('subject %d\n', s);

    [H, P] = sample(D, h, 10);
    H_all{s} = H;
    P_all{s} = P;

    [~,I] = max(P); % MAP H
    H = H(I);

    paths = random_walks(D, nwalks, D.G.N-1); % random walks
    i = randsample(length(hamils), nwalks);
    paths = [paths; hamils(i)]; % hamiltonian paths

    comm_press = 0;
    comm_nopress = 0;
    other_press = 0;
    other_nopress = 0;
    comm_press_hamil = 0;
    comm_nopress_hamil = 0;
    other_press_hamil = 0;
    other_nopress_hamil = 0;
    for i = 1:length(paths) % for each path
        path = paths{i};
        assert(length(path) == D.G.N);
        for j = 1:length(path) - 1 % for each transition
            u = path(j);
            v = path(j+1);
            if ismember([u v], comm_trans, 'rows') 
                % this is a community transition
                if H.c(u) ~= H.c(v)
                    % ...and also a cross-cluster transition according to H
                    comm_press = comm_press + 1;
                    if i > nwalks
                        % separate counts for hamiltonians
                        comm_press_hamil = comm_press_hamil + 1;
                    end
                else
                    comm_nopress = comm_nopress + 1;
                    if i > nwalks
                        % separate counts for hamiltonians
                        comm_nopress_hamil = comm_nopress_hamil + 1;
                    end
                end
            else
                % this is NOT a community transition
                if H.c(u) ~= H.c(v)
                    % ...but it is a cross-cluster transition according to H
                    other_press = other_press + 1;
                    if i > nwalks
                        % separate counts for hamiltonians
                        other_press_hamil = other_press_hamil + 1;
                    end
                else
                    other_nopress = other_nopress + 1;
                    if i > nwalks
                        % separate counts for hamiltonians
                        other_nopress_hamil = other_nopress_hamil + 1;
                    end
                end
            end
        end
    end

    comm_p(s) = comm_press / (comm_press + comm_nopress);
    other_p(s) = other_press / (other_press + other_nopress);
    comm_p_hamil(s) = comm_press_hamil / (comm_press_hamil + comm_nopress_hamil);
    other_p_hamil(s) = other_press_hamil / (other_press_hamil + other_nopress_hamil);
end

save('schapiro.mat');


load('schapiro.mat');

figure;

m = [mean(comm_p) mean(other_p);
     mean(comm_p_hamil) mean(other_p_hamil)];
se = [sem(comm_p) sem(other_p);
      sem(comm_p_hamil) sem(other_p_hamil)];

h = bar(m);
xs = [h(1).XData + h(1).XOffset, ...
      h(2).XData + h(2).XOffset];
hold on;
m = m(:); se = se(:); 
errorbar(xs, m, se, 'linestyle', 'none');
hold off;
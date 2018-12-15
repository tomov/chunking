
clear all;

rng default;

sem = @(x) std(x) / sqrt(length(x));

h = init_hyperparams();


%D = init_D_from_txt('solway4.txt');
D = init_D_from_txt('solway4.txt');
D.G.E(:) = 0; % erase edges; start empty

M = 100; % # particles
nsamples = 10; % # rejuvination steps

filename = sprintf('KL_rdms_M=%d_nsamples=%d_solway4.mat', M, nsamples);

for i = 1:M
    H(i) = init_H(D, h);
    P(i) = logpost(H(i), D, h);
end
P = P - logsumexp(P);

tic
for k = 1:size(D.G.edges,1)

    k

    % read new edge
    u = D.G.edges(k,1);
    v = D.G.edges(k,2);
    D.G.E(u,v) = 1;
    D.G.E(v,u) = 1;

    % update posteriors
    Q = P;
    for i = 1:M
        P(i) = logpost(H(i), D, h);
    end
    P = P - logsumexp(P);

    % KL
    KL(k) = KL_divergence(exp(P), exp(Q));
    fprintf('KL(%d) = %f\n', k, KL(k));

    % MCMC rejuvination
    for i = 1:M
        disp(i);
        [samples, post] = sample(D, h, nsamples, 1, 1, H(i));
        H(i) = samples(end);
        P(i) = post(end);
    end
    P = P - logsumexp(P);

    P_all(k,:) = P;
end

toc

RDM = squareRDMs(pdist(P_all, 'cosine'));

save(filename);

%D = init_D_from_txt('hourglass.txt');
%D = init_D_from_txt('starwars.txt');
D = init_D_from_txt('three.txt');

% GP kernel
% see Kondor (2002) Diffusion kernels on graphs
beta = 1;
H = -laplacian(graph(D.G.E)); % kernel generator
K = expm(beta * H); % kernel

% set up GP
meanfun = {@meanDiscrete, D.G.N};
covfun = {@covDiscrete, D.G.N};
likfun = @likGauss;

% covariance hyperparams = K but slightly modified
L = chol(K);
L(1:D.G.N+1:end) = log(diag(L));
covhyp = L(triu(true(D.G.N)));

%x = [3]'; % state 3
%y = [30]'; % ...delivered 30
%xs = [1 4]'; % how much will states 1 and 4 deliver?

%x = [5]';
%y = [30]';
%xs = [1 7 10]';

x = [5]';
y = [30]';
xs = [1 6 12]';


% GP hyperparams
% expect reward 15 from each state
hyp = struct('mean', ones(1, D.G.N) * 15, 'cov', covhyp, 'lik', -1);

% GP
[mu s2] = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y, xs);
mu

function H = MAP_H(D, h)

% TODO use EM to find actual max marginal posterior P(c | D) = sum over other H params P(c, other H params | D)

rng default; % for reproducibility

[samples, post] = sample(D, h);

[~,i] = max(post);
H = samples(i);


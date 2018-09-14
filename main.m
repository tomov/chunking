h.alpha = 1.5;

D = init_D();

H = init_H(D, h);

rng default;

[samples, post] = sample(D, h);

[~,i] = max(post);
H = samples(i);

H

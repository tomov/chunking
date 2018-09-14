h.alpha = 1.5;

D = init_D();

H = init_H(D, h);

rng default;

samples = sample(D, h);

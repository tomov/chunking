h.alpha = 1.5;

D = init_D_from_txt('hourglass.txt');

%D(2) = init_D_from_csv('/Users/momchil/Dropbox/Research/chunking/occluder3D_data/100931.csv');

H = init_H(D, h);

rng default;

[samples, post] = sample(D, h);

[~,i] = max(post);
H = samples(i);

H

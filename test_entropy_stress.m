D = init_D_from_txt('active.txt');
h = init_hyperparams;

H = sample_c(D, h, 1000);
for i = 1:100
    H = sample_c(D, h, 10, 1, 1, H(end));
end

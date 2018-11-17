

clear all;
rng default;

h.alpha = 1.5;
h.var_theta = 10;
h.var_mu = 10;

D(1) = init_D_from_txt('hourglass_rewards.txt');
for i = 1:length(D)
    tic 
    [samples, post] = sample(D(i), h, 1000);
    for j = 1:length(samples)
        H(i,j) = samples(j);
        P(i,j) = post(j);
    end
    toc
end

save demo4.mat

% load demo1.mat;

figure;

k = 5;
for i = 1:length(D)
    post = P(i,:);
    [~,I] = maxk(post, k);
    for j = 1:k
        subplot(length(D),k, (i-1)*k+j);
        plot_H(H(i,I(j)), D(i));
        if j == ceil(k/2);
            %ylabel(D(i).name);
            title(D(i).name);
        end
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
    end
end

%figure;
%        subplot(length(D),5, (i-1)*5+j);
%        plot_H(H(i,j), D(i));
%        if j == 1
%            ylabel(D(i).name);
%        end

%D(2) = init_D_from_csv('/Users/momchil/Dropbox/Research/chunking/occluder3D_data/100931.csv');


%H = MAP_H(D, h);
%plot_H(H, D);
%title(D.name);

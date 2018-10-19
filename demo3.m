% data from Rohan's chunking experiment

clear all;

h.alpha = 1.5;

D = init_Ds_from_data_Rohan('occluder3D_data');

D = D(8:9);

for i = 1:length(D)
    [samples, post] = sample(D(i), h, 1000);
    for j = 1:length(samples)
        H(i,j) = samples(j);
        P(i,j) = post(j);
    end
end

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

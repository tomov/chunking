clear all;
% {
h.alpha = 1.5;

D(1) = init_D_from_txt('hourglass.txt');

for i = 1:1
    tic 
    [D(i), samples, post] = sample_graph_update(D(i), h);
    for j = 1:length(samples)
        H(i,j) = samples(j);
        P(i,j) = post(j);
    end
    toc
end

save demo_update_graph.mat
% }

load demo_update_graph.mat;

figure;

k = 1;
for i = 1:length(D)
    post = P(i,:);
    [~,I] = maxk(post, k);
    for j = 1:k
        subplot(length(D),k, (i-1)*k+j);
        plot_H(H(i,I(j)), D(i));
        if j == ceil(k/2)
            %ylabel(D(i).name);
            title(D(i).name);
        end
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
    end
end
clear all;
close all;

%{
init_all_plots;

sem = @(x) std(x) / sqrt(length(x));

for i = 1:length(pl)
    for j = 1:length(pl(i).dirnames)
        if ~isnan(pl(i).m(j))
            continue;
        end

        D = init_Ds_from_data(pl(i).dirnames{j});
        pl(i).D{j} = D;

        % TODO dedupe w/ demo5
        clear H;
        clear P;
        for k = 1:length(D)
            k
            [samples, post] = sample(D(k), h, 1000);
            for l = 1:length(samples)
                H(k,l) = samples(l);
                P(k,l) = post(l);
            end
        end
        pl(i).H{j} = H;
        pl(i).H{j} = P;

        s = pl(i).starts(j);
        g = pl(i).goals(j);
        clear move;
        for k = 1:length(D)
            [~,I] = maxk(P(k,:), 1);
            [path, hpath] = hbfs(s, g, H(k,I(1)), D(k));
            move(k) = path(2);
        end

        % TODO dedupe w/ analyze_all_data.m
        c1 = sum(move == pl(i).nexts(j)); % count 1
        c2 = sum(move ~= pl(i).nexts(j)); % count 2
        n = sum(which);
        p = 2 * binocdf(min(c1,c2), n, 0.5);

        y = binoinv([0.025 0.975], n, 0.5);
        pl(i).ci(j) = (y(2) - y(1)) / 2;
        pl(i).n(j) = n;
        pl(i).m(j) = c1;
        pl(i).p(j) = p;
    end
end

save('model_all_data.mat');
%}

load('model_all_data.mat');

plot_all_data;

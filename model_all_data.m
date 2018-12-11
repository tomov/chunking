clear all;
close all;

rng default;

init_all_plots;

h.alpha = 2;
nsamples = 40;
filename = sprintf('model_all_data_%dsamples_MAP_%dalpha.mat', nsamples, h.alpha);

sem = @(x) std(x) / sqrt(length(x));

for i = 1:length(pl)
    for j = 1:length(pl(i).dirnames)
        if ~isnan(pl(i).m(j))
            continue;
        end

        fprintf('Modeling %d,%d: %s, %s\n', i, j, pl(i).title, pl(i).xticklabels{j});
        D = init_Ds_from_data(pl(i).dirnames{j});
        pl(i).D{j} = D;

        % TODO dedupe w/ demo5
        clear H;
        clear P;
        for k = 1:length(D)
            fprintf('      subject %d\n', k);
            tic
            [samples, post] = sample(D(k), h, nsamples);
            for l = 1:length(samples)
                H(k,l) = samples(l);
                P(k,l) = post(l);
            end
            toc
        end
        pl(i).H{j} = H;
        pl(i).P{j} = P;

        s = pl(i).starts(j);
        g = pl(i).goals(j);
        clear move;
        for k = 1:length(D)
            [~,I] = maxk(P(k,:), 1); % MAP H
            %I = length(P(k,:)); % last H
            [path, hpath] = hbfs(s, g, H(k,I(1)), D(k));
            move(k) = path(2);
        end

        % TODO dedupe w/ analyze_all_data.m
        c1 = sum(move == pl(i).nexts(j)); % count 1
        c2 = sum(move ~= pl(i).nexts(j)); % count 2
        n = length(move);
        switch pl(i).tests(j)
            case 1 % right-tailed
                p = 1 - binocdf(c1, n, 0.5);
            case 2 % left-tailed
                p = binocdf(c1, n, 0.5);
            case 3 % two-tailed
                p = 2 * binocdf(min(c1,c2), n, 0.5);
            otherwise
                assert(false);
        end

        y = binoinv([0.025 0.975], n, 0.5);
        pl(i).ci(j) = (y(2) - y(1)) / 2;
        pl(i).n(j) = n;
        pl(i).m(j) = c1;
        pl(i).p(j) = p;
    end
end

save(filename);


plot_all_data;

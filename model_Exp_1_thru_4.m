clear all;
close all;

rng default;

init_all_plots;

h = init_hyperparams;
nsamples = 10000;
take_map = false;

%filename = sprintf('model_all_data_%dsamples_MAP_%dalpha.mat', nsamples, h.alpha);
if take_map
    filename = sprintf('model_Exp_1_thru_4_samples=%d_alpha=%.4f_MAP.mat', nsamples, h.alpha);
else
    filename = sprintf('model_Exp_1_thru_4_samples=%d_alpha=%.4f_last.mat', nsamples, h.alpha);
end
filename

sem = @(x) std(x) / sqrt(length(x));

for i = 1:length(pl) - 1 % TODO not counting mines10_map; find better way (we need it for analyzing the behavioral data...)
    for j = 1:length(pl(i).dirnames)

        if ~isnan(pl(i).m(j))
            continue;
        end
        if contains(pl(i).title, 'mines') % we model it separately in mines10.m
            continue;
        end

        fprintf('Inferring H %d,%d: %s, %s\n', i, j, pl(i).title, pl(i).xticklabels{j});
        tic

        [D, filenames] = init_Ds_from_data(pl(i).dirnames{j});
        pl(i).D{j} = D;

        % TODO dedupe w/ demo5
        for subj = 1:length(D)
            fprintf('      subject %d\n', subj);

            if take_map
                [H, P] = sample_c(D(subj), h, nsamples);
                [~,I] = max(P);
                H = H(I);
            else
                [H, P] = sample_c(D(subj), h, 1, nsamples);
            end
            pl(i).H{j}(subj) = H;
        end

        toc

    end
end

filename
save(filename, '-v7.3');

for i = 1:length(pl) - 1 % TODO not counting mines10_map; find better way (we need it for analyzing the behavioral data...)
    for j = 1:length(pl(i).dirnames)

        if ~isnan(pl(i).m(j))
            continue;
        end
        if contains(pl(i).title, 'mines') % we model it separately
            continue;
        end

        fprintf('Running HBFS %d,%d: %s, %s\n', i, j, pl(i).title, pl(i).xticklabels{j});
        tic

        s = pl(i).starts(j);
        g = pl(i).goals(j);
        clear move;
        for subj = 1:length(D)
            H = pl(i).H{j}(subj);
            [path, hpath] = hbfs(s, g, H, D(subj));
            move(subj) = path(2);
        end

        % TODO dedupe w/ analyze_all_data.m
        c1 = sum(move == pl(i).nexts(j)); % count 1
        c2 = sum(move ~= pl(i).nexts(j)); % count 2
        n = length(move);
        switch pl(i).tests(j)
            case 1 % right-tailed
                p = 1 - binocdf(c1, n, 0.5);
                assert(false);
            case 2 % left-tailed
                p = binocdf(c1, n, 0.5);
                assert(false);
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

        toc
    end
end

filename
save(filename, '-v7.3');


plot_all_data;

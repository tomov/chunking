
figure;

col = [];
for i = 1:length(pl)
    subplot(length(pl), 1, i);

    m = pl(i).m ./ pl(i).n;
    ci = pl(i).ci ./ pl(i).n;
    clear se;
    for j = 1:length(pl(i).m)
        se(j) = std([ones(1,pl(i).m(j)) zeros(1,pl(i).n(j) - pl(i).m(j))]) / sqrt(pl(i).n(j));
    end
  
    hold on;
    for j = 1:length(pl(i).m)
        h = bar(j, m(j));
        if pl(i).fake(j)
            set(h, 'facealpha', 0.5);
        end
        if isempty(col)
            col = get(h, 'facecolor');
        else
            set(h, 'facecolor', col);
        end

        text(j, 1, sprintf('p = %.3f\nN = %d', pl(i).p(j), pl(i).n(j)));
    end
    errorbar(m, se, 'linestyle', 'none', 'color', 'black');
    %errorbar(m, ci, 'linestyle', 'none', 'color', 'black');
    plot([-5 length(pl(i).m) + 5], [0.5 0.5], '--', 'color', [0.8 0.8 0.8]);
    hold off;

    title(pl(i).title);

    xticks(pl(i).xticks);
    xticklabels(pl(i).xticklabels);
    xtickangle(45);
    xlim([0 length(pl(i).m) + 1]);

    yticks(pl(i).yticks);
    yticklabels(pl(i).yticklabels);
    ylim([-0.1 1.1]);
    ylabel(pl(i).ylabel);

    % stats
    %
    fprintf('\n ---- %s -----\n\n', pl(i).title);

    for j = 1:length(pl(i).m)
        fprintf('m = %.3f +- %.3f, binomial two-sided n = %d, p = %d\n', m(j), se(j), pl(i).n(j), pl(i).p(j));
    end
    [tbl, chi2stat, pval] = chi2(pl(i).m, pl(i).n);
    tbl
    N = sum(pl(i).n);
    df = (2 - 1) * (length(pl(i).m) - 1);
    fprintf('test of independence (crosstab) for all: chi2(%d, %d) = %.3f, p = %.4f\n', df, N, chi2stat, pval);

    for j = 1:length(pl(i).m)
        for k = j+1:length(pl(i).m)
            fprintf('   for %s vs. %s:\n', pl(i).xticklabels{j}, pl(i).xticklabels{k});
            [tbl, chi2stat, pval] = chi2(pl(i).m([j k]), pl(i).n([j k]));
            N = sum(pl(i).n([j k]));
            df = (2 - 1) * (2 - 1);
            fprintf('              chi2(%d, %d) = %.3f, p = %.4f\n', df, N, chi2stat, pval);
        end
    end

end


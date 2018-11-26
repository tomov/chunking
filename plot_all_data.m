
figure;

col = [];
for i = 1:length(pl)
    subplot(length(pl), 1, i);

    m = pl(i).m ./ pl(i).n;
    ci = pl(i).ci ./ pl(i).n;
  
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
    errorbar(m, ci, 'linestyle', 'none', 'color', 'black');
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
end

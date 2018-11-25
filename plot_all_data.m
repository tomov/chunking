
figure;

for i = 1:length(pl)
    subplot(length(pl), 1, i);

    m = pl(i).m ./ pl(i).n;
    ci = pl(i).ci ./ pl(i).n;
    
    bar(m);
    hold on;
    errorbar(m, ci, 'linestyle', 'none');
    plot([-5 length(pl(i).m) + 5], [0.5 0.5], '--', 'color', [0.6 0.6 0.6]);
    hold off;

    title(pl(i).title);

    xticks(pl(i).xticks);
    xticklabels(pl(i).xticklabels);
    xtickangle(45);
    xlim([0 length(pl(i).m) + 1]);

    yticks(pl(i).yticks);
    yticklabels(pl(i).yticklabels);
    ylim([-0.1 1.1]);
end

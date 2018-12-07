h.alpha = 1.5;

D(1) = init_D_from_txt('hourglass1.txt');

nsubjects = 40;

count_4_1 = 0;
count_4_7 = 0;
both = 0;
for i = 1:nsubjects
    [~, samples, ~] = sample_graph_update(D(1), h);
    H = samples(1);
    if (H.c(4) == H.c(1))
        count_4_1 = count_4_1 + 1;
    end
    if (H.c(4) == H.c(7))
        count_4_7 = count_4_7 + 1;
    end
    if (H.c(4) == H.c(7)) && (H.c(4) == H.c(1))
       both = both + 1; 
    end
end

save count_expirement_1.mat
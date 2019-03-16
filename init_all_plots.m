pl(1).m = [NaN 50];
pl(1).ci = [NaN 10];
pl(1).n = [NaN 100];
pl(1).p = [NaN 0.7];
pl(1).fake = [0 1];
pl(1).nrows = [81 81];
pl(1).title = 'subway 10, full view';
pl(1).ylabel = ['P(right)'];
pl(1).xticklabels = {'chunks', 'no chunks'};
pl(1).xticks = [1 2];
pl(1).yticklabels = {'0', '0.5', '1'};
pl(1).yticks = [0 0.5 1];
pl(1).dirnames = {'exp/results/subway10_map/', ''};
pl(1).starts = [6 6];
pl(1).goals = [1 1];
pl(1).nexts = [5 5];
pl(1).tests = [3 3]; % 1 = right tailed, 2 = left tailed, 3 = two-tailed



pl(2).m = [NaN NaN NaN NaN];
pl(2).ci = [NaN NaN NaN NaN];
pl(2).n = [NaN NaN NaN NaN];
pl(2).p = [NaN NaN NaN NaN];
pl(2).fake = [0 0 0 0];
pl(2).nrows = [81 21 81 81];
pl(2).title = 'subway 9, full view';
pl(2).ylabel = ['P(bad)'];
pl(2).xticklabels = {'bad chunks', 'no chunks (1)', 'no chunks (2)', 'good chunks'};
pl(2).xticks = [1 2 3 4];
pl(2).yticklabels = {'0 (good)', '0.5', '1 (bad)'};
pl(2).yticks = [0 0.5 1];
pl(2).dirnames = {'exp/results/subway9_map/', 'exp/results/subway9_map_control_1/', 'exp/results/subway9_map_control_2/', 'exp/results/subway9_map_goodchunks/'};
pl(2).starts = [6 6 6 6];
pl(2).goals = [1 1 1 1];
pl(2).nexts = [5 5 5 5];
pl(2).tests = [3 3 3 3]; % 1 = right tailed, 2 = left tailed, 3 = two-tailed



pl(3).m = [NaN 50];
pl(3).ci = [NaN 10];
pl(3).n = [NaN 100];
pl(3).p = [NaN 0.8];
pl(3).fake = [0 1];
pl(3).nrows = [83 81];
pl(3).title = 'subway 10';
pl(3).ylabel = ['P(right)'];
pl(3).xticklabels = {'chunks', 'no chunks'};
pl(3).xticks = [1 2];
pl(3).yticklabels = {'0', '0.5', '1'};
pl(3).yticks = [0 0.5 1];
pl(3).dirnames = {'exp/results/subway10_repro/', ''};
pl(3).starts = [6 6];
pl(3).goals = [1 1];
pl(3).nexts = [5 5];
pl(3).tests = [3 3]; % 1 = right tailed, 2 = left tailed, 3 = two-tailed


pl(4).m = [NaN NaN NaN];
pl(4).ci = [NaN NaN NaN];
pl(4).n = [NaN NaN NaN];
pl(4).p = [NaN NaN NaN];
pl(4).fake = [0 0 0];
pl(4).nrows = [81 81 81];
pl(4).title = 'subway 9';
pl(4).ylabel = ['P(bad)'];
pl(4).xticklabels = {'bad chunks', 'no chunks (2)', 'good chunks'};
pl(4).xticks = [1 2 3 4];
pl(4).yticklabels = {'0 (good)', '0.5', '1 (bad)'};
pl(4).yticks = [0 0.5 1];
pl(4).dirnames = {'exp/results/subway9/', 'exp/results/subway9_control/', 'exp/results/subway9_goodchunks/'};
pl(4).starts = [6 6 6];
pl(4).goals = [1 1 1];
pl(4).nexts = [5 5 5];
pl(4).tests = [3 3 3]; % 1 = right tailed, 2 = left tailed, 3 = two-tailed

% see mines10.m
pl(5).m = [NaN 50];
pl(5).ci = [NaN 10];
pl(5).n = [NaN 100];
pl(5).p = [NaN 0.7];
pl(5).fake = [0 1];
pl(5).nrows = [101 81];
pl(5).title = 'mines 10, full view';
pl(5).ylabel = ['P(right)'];
pl(5).xticklabels = {'chunks', 'no chunks'};
pl(5).xticks = [1 2];
pl(5).yticklabels = {'0', '0.5', '1'};
pl(5).yticks = [0 0.5 1];
pl(5).dirnames = {'exp/results/mines10_map/', ''};
pl(5).starts = [6 6];
pl(5).goals = [1 1];
pl(5).nexts = [5 5];
pl(5).tests = [3 3]; % 1 = right tailed, 2 = left tailed, 3 = two-tailed





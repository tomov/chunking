function [h, x, y] = plot_solway4_graph(H, D)

r1 = 4;
r2 = 2;
r3 = 1;
cx0 = cos([0:2]*2*pi/3 + pi/2) * r1;
cy0 = sin([0:2]*2*pi/3 + pi/2) * r1;
cx1 = cos([0:2]*2*pi/3 + pi/2) * r2;
cy1 = sin([0:2]*2*pi/3 + pi/2) * r2;
cx2 = cos([0:2]*2*pi/3 + pi/2) * r3;
cy2 = sin([0:2]*2*pi/3 + pi/2) * r3;

x = zeros(1, D.G.N);
y = zeros(1, D.G.N);
for i = 1:D.G.N
    i0 = floor((i - 1) / 9) + 1;
    i1 = floor(mod(i - 1, 9) / 3) + 1;
    i2 = mod(i - 1, 3) + 1;
    x(i) = x(i) + cx0(i0) + cx1(i1) + cx2(i2);
    y(i) = y(i) + cy0(i0) + cy1(i1) + cy2(i2);
end

H.c = [1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3];


h = plot_H(H, D);
set(h, 'XData', x);
set(h, 'YData', y);
labelnode(h, 1:D.G.N, '');
for i = 1:D.G.N
    highlight(h, i, 'nodecolor', [0.1 0.1 0.1], 'markersize', 12);
end
hold on;

h = plot_H(H, D);
set(h, 'XData', x);
set(h, 'YData', y);
labelnode(h, 1:D.G.N, '');
for i = 1:D.G.N
    if H.c(i) == 1
        highlight(h, i, 'nodecolor', [80 150 114]/255);
    elseif H.c(i) == 3
        highlight(h, i, 'nodecolor', [16 118 188]/255);
    else
        highlight(h, i, 'nodecolor', [246 135 31]/255);
    end

    for j = i+1:D.G.N
        if D.G.E(i,j)
            highlight(h, i, j, 'edgecolor', [0.1 0.1 0.1], 'linewidth', 2);
        end
    end
end

hold off;

set(gca, 'xtick', []);
set(gca, 'ytick', []);
xlim([min(x) - 1, max(x) + 1]);
ylim([min(y) - 1, max(y) + 1]);

clear all;

rng default;

sem = @(x) std(x) / sqrt(length(x));

N = 30; % participants
h.alpha = 5;

%{
D = init_D_from_txt('schapiro.txt');

for i = 1:10
    [H, P] = sample(D, h, 1000);
    H_all{i} = H;
    P_all{i} = P;

    [~,I] = max(P);
    H_map(i) = H(I);
end

save('schapiro_viz.mat');
%}

load('schapiro_viz.mat');

r1 = 6;
r2 = 3;
cx = cos([0:2]*2*pi/3 + pi/2) * r1;
cy = sin([0:2]*2*pi/3 + pi/2) * r1;
x(1:5) = cos([0:4]*2*pi/5 + pi/2) * r2 + cx(1);
y(1:5) = sin([0:4]*2*pi/5 + pi/2) * r2 + cy(1);
x(6:10) = cos([0:4]*2*pi/5 + pi/2) * r2 + cx(2);
y(6:10) = sin([0:4]*2*pi/5 + pi/2) * r2 + cy(2);
x(11:15) = cos([0:4]*2*pi/5 + pi/2) * r2 + cx(3);
y(11:15) = sin([0:4]*2*pi/5 + pi/2) * r2 + cy(3);


figure;

for i = 1:length(H_map)
    subplot(1, length(H_map), i);
    h = plot_H(H_map(i), D);
    set(h, 'XData', x);
    set(h, 'YData', y);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
end

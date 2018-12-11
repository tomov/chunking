clear all;



figure('pos', [100 100 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

p = zeros(13,50);
x = zeros(size(p));
m = zeros(1,size(p,2));
b = zeros(1,size(p,2));

t = 5;
p(12,t+1) = 1;
p(11,t+2) = 1;
p(13,t+2) = 1;
p(6,t+3) = 1;
p(5,t+4) = 1;
p(4,t+5) = 1;

x(6,1:t+10) = 1;
x(5,t+11:t+15) = 1;
x(4,t+16:t+20) = 1;

b(t+9) = 3;
%b(t+9) = 1;
m(t+10) = 1;

b(t+14) = 1;
m(t+15) = 1;

t = 20;
p(12,t+1) = 1;
p(11,t+2) = 1;
p(13,t+2) = 1;
p(4,t+3) = 1;
p(3,t+4) = 1;
p(2,t+5) = 1;
p(1,t+6) = 1;

x(4,t+1:t+11) = 1;
x(3,t+12:t+16) = 1;
x(2,t+17:t+21) = 1;
x(1,t+22:t+28) = 1;

b(t+10) = 3;
%b(t+10) = 1;
m(t+11) = 1;

b(t+15) = 1;
m(t+16) = 1;

b(t+20) = 1;
m(t+21) = 1;

Str = zeros(1,length(b)*3);
Str(3:3:end) = b;
M1 = zeros(1,length(m)*3);
M1(1:3:end) = m;


h = subplot(4,1,1);
pos = get(h, 'position');
pos(1) = pos(1) * -0.15;
pos(2) = pos(2) * 0.96;
pos(3) = pos(3) * 1.35;
pos(4) = pos(4) * 1.35;
subplot(4,1, 1, 'position', pos);

PICpng = imread('circuitH.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  





subplot(8,1,3);

plot(M1);
ylabel('M1');
ylim([-0.5 1.5]);
set(gca,'ytick', []);
set(gca,'xtick', []);


subplot(8,1,4);

plot(Str);
ylabel('Str');
ylim([-1 4]);
set(gca,'ytick', []);
set(gca,'xtick', []);


h = subplot(2,1,2);
pos = get(h, 'position');
pos(1) = pos(1) * 1.0;
pos(2) = pos(2) * 1.6;
pos(3) = pos(3) * 1.0;
pos(4) = pos(4) * 1.0;
subplot(2,1, 2, 'position', pos);

I = 1 - min(x * 0.5 + p, 1);
C(:,:,1) = I;
C(:,:,2) = I;
C(:,:,3) = I;
imagesc(C);
yticks(1:13);
set(gca, 'TickLength', [0 0]);
ylabel('Hippocampus');
xlabel('time');



ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.08, 0.95, 'A', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.41, 0.95, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.79, 0.95, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.08, 0.71, 'D', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.08, 0.60, 'E', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.08, 0.50, 'F', 'FontSize', lettersize, 'FontWeight', 'bold');





% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
print('figures/neuro.pdf', '-dpdf');

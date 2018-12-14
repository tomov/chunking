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
Str(Str == 1) = 0;
M1 = zeros(1,length(m)*3);
M1(1:3:end) = m;


%
% within-trial
%




% circuit subway 10


h = subplot(3,2,1);
pos = get(h, 'position');
pos(1) = pos(1) * 1.00;
pos(2) = pos(2) * 1.00;
pos(3) = pos(3) * 1.00;
pos(4) = pos(4) * 1.00;
subplot(3,2,1, 'position', pos);

PICpng = imread('circuit_subway10.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  


title('Within trial', 'fontsize', fontsize);


% hippocampus

h = subplot(3,2,3);
pos = get(h, 'position');
pos(1) = pos(1) * 1.0;
pos(2) = pos(2) * 1.1;
pos(3) = pos(3) * 1.0;
pos(4) = pos(4) * 1.0;
subplot(3,2,3, 'position', pos);

I = 1 - min(x * 0.5 + p, 1);
C(:,:,1) = I;
C(:,:,2) = I;
C(:,:,3) = I;
imagesc(C);
yticks(1:13);
xticks([2 (find(b > 0) + 2)]);
xticklabels([6 5 4 3 2 1]);
set(gca, 'TickLength', [0 0]);
ylabel('Hippocampus');
xlabel('state');



% Str and M1

subplot(6,2,9);

plot(Str);
ylabel('Str');
ylim([-1 4]);
set(gca,'ytick', []);
set(gca,'xtick', []);

subplot(6,2,11);

plot(M1);
ylabel('M1');
ylim([-0.5 1.5]);
set(gca,'ytick', []);
set(gca,'xtick', []);

xlabel('time within trial');





%
% across-trial
%


% circuit example

h = subplot(3,2,2);
pos = get(h, 'position');
pos(1) = pos(1) * 1.00;
pos(2) = pos(2) * 1.00;
pos(3) = pos(3) * 1.00;
pos(4) = pos(4) * 1.00;
subplot(3,2,2, 'position', pos);

PICpng = imread('ex_circuit.png');
[rows columns numberOfColorChannels] = size(PICpng);
imshow(PICpng, 'InitialMagnification', 'fit');  

title('Across trials', 'fontsize', fontsize);


% KL in PPC

load('KL_rdms.mat');

KL(1) = 0;
subplot(6,2,6);
plot(KL);
xlim([1 length(KL)]);
ylim([min(KL)-50 max(KL) + 50]);
set(gca,'ytick', []);
set(gca,'xtick', []);
xlabel('trial number');
ylabel('PPC');



% RDM

subplot(2,2,4);
image(scale01(rankTransform_equalsStayEqual(RDM,1)),'CDataMapping','scaled')
colormap(gca, RDMcolormap)
colorbar;
xlabel('trial number');
ylabel('trial number');



ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.07, 0.95, 'A', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.07, 0.69, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.07, 0.39, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.52, 0.95, 'D', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.52, 0.69, 'E', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.52, 0.50, 'F', 'FontSize', lettersize, 'FontWeight', 'bold');
%text(0.41, 0.95, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
%text(0.79, 0.95, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');
%text(0.08, 0.71, 'D', 'FontSize', lettersize, 'FontWeight', 'bold');
%text(0.08, 0.60, 'E', 'FontSize', lettersize, 'FontWeight', 'bold');
%text(0.08, 0.50, 'F', 'FontSize', lettersize, 'FontWeight', 'bold');





% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
print('figures/neuro.pdf', '-dpdf');

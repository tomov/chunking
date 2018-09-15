clear all

h.alpha = 1.5;

D(1) = init_D_from_txt('hourglass.txt');
D(2) = init_D_from_txt('solway1.txt');
D(3) = init_D_from_txt('solway2.txt');
D(4) = init_D_from_txt('schapiro.txt');

tic

for i = 1:length(D)
    [samples, post] = sample(D(i), h);
    for j = 1:length(samples)
        H(i,j) = samples(j);
        P(i,j) = post(j);
    end
end

toc

%figure;
%        subplot(length(D),5, (i-1)*5+j);
%        plot_H(H(i,j), D(i));
%        if j == 1
%            ylabel(D(i).name);
%        end

%D(2) = init_D_from_csv('/Users/momchil/Dropbox/Research/chunking/occluder3D_data/100931.csv');


%H = MAP_H(D, h);
%plot_H(H, D);
%title(D.name);

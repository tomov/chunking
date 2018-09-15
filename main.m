h.alpha = 1.5;

D(1) = init_D_from_txt('hourglass.txt');
D(2) = init_D_from_txt('solway1.txt');
D(3) = init_D_from_txt('solway2.txt');
D(4) = init_D_from_txt('schapiro.txt');

figure;
for i = 1:length(D)
    [samples, post] = sample(D(i), h);
    for j = 1:5
        H(i,j) = samples(end+1-j);

        subplot(length(D),5, (i-1)*5+j);
        plot_H(H(i,j), D(i));
        if j == 1
            ylabel(D(i).name);
        end
    end
end


%D(2) = init_D_from_csv('/Users/momchil/Dropbox/Research/chunking/occluder3D_data/100931.csv');


%H = MAP_H(D, h);
%plot_H(H, D);
%title(D.name);

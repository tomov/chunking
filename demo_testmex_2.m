% with realistic data

load('mines10_alpha=2_nsamples=100.mat')

H = map_H{1};
D.G.edges = [];
for i = 1:D.G.N
    for j = i+1:D.G.N
        if D.G.E(i,j)
            D.G.edges = [D.G.edges; i j];
        end
    end
end

Hout = testmex(D, h, 100, 1, 1, H)

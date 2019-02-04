load('demo_testmex.mat');

Hout = testmex(D, h, 12)

% or
Hout = testmex(D, h, 100, 1, 1, H)

Hout = testmex(D, h, 100, 1, 1, Hout)

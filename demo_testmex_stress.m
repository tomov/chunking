% checks for memory leaks, etc
demo_testmex

for i = 1:1000
    Hout = testmex(D, h, 100, 1, 1, Hout)
end

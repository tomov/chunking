models = {'ISL', 'PF1', 'PF2', 'PF3', 'PF4'};

files = {'mat/isl_MH_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat', ...
         'mat/pf1_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat', ...
         'mat/pf2_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat', ...
         'mat/pf3_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat', ...
         'mat/pf4_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat'};

LME = [];
LME_probes = [];

for i = 1:numel(files);

    load(files{i}, 'lme', 'lme_probes');
    LME = [LME lme];
    LME_probes = [LME_probes lme_probes];

end


[alpha,exp_r,xp,pxp,bor] = bms(LME);

{models}
pxp
bor


[alpha,exp_r,xp,pxp,bor] = bms(LME_probes);

{models}
pxp
bor


save('mat/isl_bms.mat');

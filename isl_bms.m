addpath('/n/home04/mtomov13/libs/mfit');

models = {'ISL', 'PF1', 'PF2', 'PF3', 'PF4'};

files = {'mat/isl_MH_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat', ...
         'mat/pf1_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat', ...
         'mat/pf2_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat', ...
         'mat/pf3_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat', ...
         'mat/pf4_alpha=1.0000_nsamples=1000_div_eps=0.6000_last_np=1000.mat'};
[~, name] = system('hostname');

LME = [];
LME_probes = [];

for i = 1:numel(files);

    filename = files{i};
    if ~isempty(strfind(name, 'omchil')) || ~isempty(strfind(name, 'dhcp-'))
        % local
        filename
    else
        % cannon
        filename = fullfile(getenv('MY_SCRATCH'), 'chunking', filename);
        filename
    end

    load(filename, 'lme', 'lme_probes');

    %size(lme)

    LME = [LME lme(1:104,:)];
    LME_probes = [LME_probes lme_probes(1:104,:)];

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

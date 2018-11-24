function [p_1, p_2] = predict(D, h, nsamples, burnin, lag, m, tau, ids)
    %
    % D, h, nsmaples, burnin, and lag are the parameters for the sample fn
    % m: # of times to sample H
    % tau: temperature for Softmax function
    % ids: vector of length 2 where ids[1] is the node ID for the first
    %      neigbhor and ids[2] is the node ID for the second neighbor.
    %

    % Take m samples of H
    for i = 1:m
        sample(D, h, nsamples, burnin, lag); % TODO fill in params
    end
    
    % Expected value of mu_1 given observed reward 
    % (equal to sum over m of theta_c1 in each sample of H)
    mu_1_exp = 0;
    for i = 1:m
        curSampleH = sample(D, h, nsamples, burnin, lag); 
        mu_1_exp = mu_1_exp + curSampleH.theta_c1;
    end
    
    % Expected value of mu_2 given observed reward 
    % (equal to sum over m of theta_c2 in each sample of H)
    mu_2_exp = 0;
    for i = 1:m
        curSampleH = sample(D, h, nsamples, burnin, lag); 
        mu_2_exp = mu_2_exp + curSampleH.theta_c2;
    end
    
    % Pass through softmax function
    % P(choose 1 | E[mu_1], E[mu_2]) = exp(E[mu_1]/tau)/(exp(E[mu_1]/tau) +
    %   exp(E[mu_3]/tau)); draw from this multinomial distribution (??)
    prob = exp(mu_1_exp/tau)/((mu_2_exp/tau) + exp(mu_3_exp/tau));
    
    % Return p_1 and p_2
    p_1 = None;
    p_2 = None;
end

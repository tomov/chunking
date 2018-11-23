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
        sample(); % TODO fill in params
    end
    
    % Expected value of mu_1 given observed reward 
    % (equal to sum over m of theta_c1 in each sample of H)
    
    % Expected value of mu_2 given observed reward 
    % (equal to sum over m of theta_c2 in each sample of H)
    
    % Pass through softmax function
    % P(choose 1 | E[mu_1], E[mu_2]) = exp(E[mu_1]/tau)/(exp(E[mu_1]/tau) +
    %   exp(E[mu_3]/tau)); draw from this multinomial distribution (??)
    
    % Return p_1 and p_2
    p_1 = None;
    p_2 = None;
end

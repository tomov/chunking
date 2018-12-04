function [p, mu, H, post] = predict(D, h, M, burnin, lag, tau, ids)
    %
    % D, h, nsmaples, burnin, and lag are the parameters for the sample fn
    % M: # of times to sample H
    % tau: temperature for Softmax function
    % ids: vector of length 2 where ids[1] is the node ID for the first
    %      neigbhor and ids[2] is the node ID for the second neighbor.
    %

    % Take m samples of H
    % w is probabilities (normalized weights)
    [H, post] = sample(D, [], h, M, burnin, lag);
    %w = exp(w)/sum(exp(w));
    
    
    % Expected value of mu_1 given observed reward 
    % (equal to sum over m of theta_c1 in each sample of H)
    [~, max_index] = max(post);
    H_max = H(max_index);
    mu = zeros(D.G.N);
    for i = 1:D.G.N
%         mu(i) = 0;
%         for m = 1:M
%             mu(i) = mu(i) + H(m).theta(H(m).c(i));
%         end
%         mu(i) = mu(i)/M;
        mu(i) =  H_max.theta(H_max.c(i));
    end 
    
    % Expected value of mu_2 given observed reward 
    % (equal to sum over m of theta_c2 in each sample of H)
    
    % Pass through softmax function
    % P(choose 1 | E[mu_1], E[mu_2]) = exp(E[mu_1]/tau)/(exp(E[mu_1]/tau) +
    %   exp(E[mu_3]/tau)); draw from this multinomial distribution (??)
    p = exp(mu(ids)/tau)/(sum(exp(mu(ids)/tau)));
    
    % Return p_1 and p_2
    %p_1 = None;
    %p_2 = None;
end

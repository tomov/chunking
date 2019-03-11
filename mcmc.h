// MCMC sampler for hierarchies & helper functions
//

#include "printmex.h"
#include "datastructs.h"

#include <cmath>


// P(H|D) for updates of c_i
// i.e. with new c's up to c_i, the candidate c_i, then old c's after (and old rest of H)
//
// TODO optimize the crap out of it
double logpost_c_i(int c_i, int i, Hierarchy& H, const Data &D, const Hyperparams &h)
{
    double old_c_i = H.c[i];
    H.c[i] = c_i;
    double logP = H.LogPost(D, h);
    H.c[i] = old_c_i;
    return logP;
}

// proposal PMF for c_i
// inspired by Algorithm 5 from Neal 1998: MCMC for DP mixtures
// TODO use algorithm 7 instead -- this one just adds everyone to a single component and it sucks
//
// notice proposal distr q(x'|x) = q(x') i.e. it's independent of previous cluster assignment => implements Independent Metropolis-Hastings
// TODO optimize the crap out of it
std::vector<double> propP_c_i(int i, const Hierarchy& H, const Data &D, const Hyperparams &h)
{
    std::vector<double> cnt(H.cnt.begin(), H.cnt.end()); // TODO optimize -- no need to copy, could be done in O(1) 
    assertThis(H.c[i] - 1 < cnt.size());
    assertThis(cnt[H.c[i] - 1] > 0);
    cnt[H.c[i] - 1]--; // careful with off-by-one!

    DEBUG_PRINT("cnt = [");
    for (int k = 0; k < cnt.size(); k++)
    {
        DEBUG_PRINT("%lf ", cnt[k]);
    }
    DEBUG_PRINT("]\n");

    std::vector<int> z;
    for (int k = 0; k < cnt.size(); k++)
    {
        if (cnt[k] == 0)
        {
            z.push_back(k); // reuse empty bins -- notice this is legit b/c we're not reusing parameters; empty bin = new cluster; in fact, we have to take care of that in case c_i_old was the only one by itelf
        }
    }

    if (z.empty())
    {
        cnt.push_back(h.alpha);
    }
    else
    {
        for (int i = 0; i < z.size(); i++)
        {
            //cnt[z[i]] = h.alpha; // <-- this is what sample.m used to do 
            cnt[z[i]] = h.alpha / z.size(); // notice all the empty bins have equal probability = alpha, but that's fine b/c it doesn't matter which one we use as the new cluster; we just have to make sure their total probability is not too high, otherwise effective alpha is greater
        }
    }

    double sum = 0;
    for (int k = 0; k < cnt.size(); k++)
    {
        sum += cnt[k];
    }

    std::vector<double> P(cnt.size());
    for (int k = 0; k < cnt.size(); k++)
    {
        P[k] = cnt[k] / sum;
    }
    return P;
}

// draw proposals for c_i
//
int proprnd_c_i(int /*c_i_old*/, int i, const Hierarchy& H, const Data &D, const Hyperparams &h)
{
    std::vector<double> P = propP_c_i(i, H, D, h);
    int c_i_new = CatRnd(P) + 1; // careful with off-by-one!

    // TODO bridges
    return c_i_new;
}

// proposal distribution for c_i; note that it doesn't really depend on c_i_old
//
double logprop_c_i(int c_i_new, int /*c_i_old*/, int i, const Hierarchy& H, const Data &D, const Hyperparams &h)
{
    std::vector<double> P = propP_c_i(i, H, D, h);
    double logP = log(P[c_i_new - 1]); // off by one!!
    return logP;
}



// P(H|D) for updates of p
//
double logpost_p(double p, Hierarchy& H, const Data &D, const Hyperparams &h)
{
    double old_p = H.p;
    H.p = p;
    double logP = H.LogPost(D, h);
    H.p = old_p;
    return logP;
}

// P(H|D) for updates of q
//
double logpost_q(double q, Hierarchy& H, const Data &D, const Hyperparams &h)
{
    double old_q = H.q;
    H.q = q;
    double logP = H.LogPost(D, h);
    H.q = old_q;
    return logP;
}

// P(H|D) for updates of tp
//
double logpost_tp(double tp, Hierarchy& H, const Data &D, const Hyperparams &h)
{
    double old_tp = H.tp;
    H.tp = tp;
    double logP = H.LogPost(D, h);
    H.tp = old_tp;
    return logP;
}

// P(H|D) for updates of hp
//
double logpost_hp(double hp, Hierarchy& H, const Data &D, const Hyperparams &h)
{
    double old_hp = H.hp;
    H.hp = hp;
    double logP = H.LogPost(D, h);
    H.hp = old_hp;
    return logP;
}



// draw proposals for p; random walk, bounded between 0 and 1
// TODO rm useless params; from below too
//
double proprnd_p(double p_old, const Hierarchy &H, const Data &D, const Hyperparams &h)
{
    while (true) // TODO can use universality of uniform inverse CDF thingy
    {
        double p_new = NormRnd(p_old, 0.1); // TODO const TODO adaptive
        if (p_new <= 1 && p_new >= 0)
        {
            return p_new; // keep params within bounds

        }
    }
}

// proposal PDF for p
// accounts for truncating that keeps params within bounds 
//
double logprop_p(double p_new, double p_old, const Hierarchy &H, const Data &D, const Hyperparams &h)
{
    double Z = NormCDF(1, p_old, 0.1) - NormCDF(0, p_old, 0.1); // TODO consts TODO adaptive
    double logP = log(NormPDF(p_new, p_old, 0.1)) - log(Z);
    return logP;
}


// P(H|D) for updates of theta
//
double logpost_theta(double theta_k, int k, Hierarchy& H, const Data &D, const Hyperparams &h)
{
    double theta_old = H.theta[k];
    H.theta[k] = theta_k;
    double logP = H.LogPost(D, h);
    H.theta[k] = theta_old;
    return logP;
}

// P(H|D) for updates of mu
// TODO optimize all the logposts....
//
double logpost_mu(double mu_i, int i, Hierarchy& H, const Data &D, const Hyperparams &h)
{
    double mu_old = H.mu[i];
    H.mu[i] = mu_i;
    double logP = H.LogPost(D, h);
    H.mu[i] = mu_old;
    return logP;
}


// proposals for mu, theta; random walk
//
double proprnd_unbounded(double p_old, const Hierarchy &H, const Data &D, const Hyperparams &h)
{
    double p_new = NormRnd(p_old, 1); // TODO const TODO adaptive
    return p_new;
}

// unbounded proposal f'n
// 
double logprop_unbounded(double p_new, double p_old, const Hierarchy &H, const Data &D, const Hyperparams &h)
{
    double logP = log(NormPDF(p_new, p_old, 1)); // TODO const TODO adaptive
    return logP;
}




// logpost_new = f(x')
// logpost_old = f(x)
// logprop_new = q(x'|x)
// logprop_old = q(x|x')
//
bool MetropolisHastingsFlip(double logpost_new, double logpost_old, double logprop_new, double logprop_old)
{
    double logA = std::min(log(1), (logpost_new - logprop_new) - (logpost_old - logprop_old)); // note logs
    double A = exp(logA);
    assertThis(A >= 0 - EPS && A <= 1 + EPS, "A >= 0 && A <= 1");

    double U = UnifRnd();

    return (U < A);
}




std::vector<Hierarchy*> 
sample(const Data &D, const Hyperparams &h, const int nsamples, const int burnin, const int lag, Hierarchy &H, std::vector<double> /*out*/ &post)
{
    post.clear();
    std::vector<Hierarchy*> samples;

    // Roberts & Rosenthal (2009)
    // Metropolis-within-Gibbs
    //
    // as a reminder, Metropolis-Hastings acceptance probability:
    // A = min(1, [f(x') / q(x'|x)] / [f(x) / q(x|x')]
    // where f = target distr, q = proposal distr, x' = proposal
    //
    for (int n = 0; n < nsamples * lag + burnin; n++)
    {
        // update clusters (Neal 1998: MCMC for DP mixtures)
        //
        for (int i = 0; i < D.G.N; i++)
        {
            int c_i_new = proprnd_c_i(H.c[i], i, H, D, h);
            int c_i_old = H.c[i];

            double logpost_new = H.LogPost_c_i(c_i_new, i, D, h); // f(x') TODO can probs speed up, but connectivity messes things up
            double logpost_old = H.LogPost(D, h); // f(x)

            double logprop_new = logprop_c_i(c_i_new, c_i_old, i, H, D, h); // q(x'|x)
            double logprop_old = logprop_c_i(c_i_old, c_i_new, i, H, D, h); // q(x|x')

            if (MetropolisHastingsFlip(logpost_new, logpost_old, logprop_new, logprop_old))
            { 
                // accept
                // TODO notice if it's a new cluster, we redraw a new theta that is different from the one
                // we used to compute the logpost_new; that's probably okay, but not clear from Neal (1998)
                H.Update_c_i(c_i_new, i, D, h);
            }
            else
            {
                // reject
                assertThis(H.c[i] == c_i_old);
            }
        }

        // TODO everywhere do not recompute the entire logpost; only local change
        // BUT do sanity check first

        // update H.p 
        //
        {
            double p_new = proprnd_p(H.p, H, D, h);
            double p_old = H.p;

            double logpost_new = logpost_p(p_new, H, D, h);
            double logpost_old = H.LogPost(D, h); // TODO calculate once only

            double logprop_new = logprop_p(p_new, p_old, H, D, h);
            double logprop_old = logprop_p(p_old, p_new, H, D, h);

            if (MetropolisHastingsFlip(logpost_new, logpost_old, logprop_new, logprop_old))
            { 
                // accept
                H.p = p_new;
            }
            else
            {
                // reject
                assertThis(fabs(H.p - p_old) < EPS);
            }
        }

        // update H.q TODO refactor
        //
        {
            double q_new = proprnd_p(H.q, H, D, h);
            double q_old = H.q;

            double logpost_new = logpost_q(q_new, H, D, h);
            double logpost_old = H.LogPost(D, h); // TODO calculate once only

            double logprop_new = logprop_p(q_new, q_old, H, D, h);
            double logprop_old = logprop_p(q_old, q_new, H, D, h);

            if (MetropolisHastingsFlip(logpost_new, logpost_old, logprop_new, logprop_old))
            { 
                // accept
                H.q = q_new;
            }
            else
            {
                // reject
                assertThis(fabs(H.q - q_old) < EPS);
            }
        }

        // update H.hp TODO refactor
        //
        {
            double hp_new = proprnd_p(H.hp, H, D, h);
            double hp_old = H.hp;

            double logpost_new = logpost_hp(hp_new, H, D, h);
            double logpost_old = H.LogPost(D, h); // TODO calculate once only

            double logprop_new = logprop_p(hp_new, hp_old, H, D, h);
            double logprop_old = logprop_p(hp_old, hp_new, H, D, h);

            if (MetropolisHastingsFlip(logpost_new, logpost_old, logprop_new, logprop_old))
            { 
                // accept
                H.hp = hp_new;
            }
            else
            {
                // reject
                assertThis(fabs(H.hp - hp_old) < EPS);
            }
        }

        // update H.tp TODO refactor
        //
        {
            double tp_new = proprnd_p(H.tp, H, D, h);
            double tp_old = H.tp;

            double logpost_new = logpost_tp(tp_new, H, D, h);
            double logpost_old = H.LogPost(D, h); // TODO calculate once only

            double logprop_new = logprop_p(tp_new, tp_old, H, D, h);
            double logprop_old = logprop_p(tp_old, tp_new, H, D, h);

            if (MetropolisHastingsFlip(logpost_new, logpost_old, logprop_new, logprop_old))
            { 
                // accept
                H.tp = tp_new;
            }
            else
            {
                // reject
                assertThis(fabs(H.tp - tp_old) < EPS);
            }
        }

        // TODO do
        // update thetas
        //
        /*
        for (int k = 0; k < H.theta.size(); k++)
        {
            if (H.cnt[k] > 0)
            {
                double theta_k_new = proprnd_unbounded(H.theta[k], H, D, h);
                double theta_k_old = H.theta[k];

                double logpost_new = logpost_theta(theta_k_new, k, H, D, h);
                double logpost_old = H.LogPost(D, h); // TODO optim

                double logprop_new = logprop_unbounded(theta_k_new, theta_k_old, H, D, h);
                double logprop_old = logprop_unbounded(theta_k_old, theta_k_new, H, D, h);

                if (MetropolisHastingsFlip(logpost_new, logpost_old, logprop_new, logprop_old))
                { 
                    // accept
                    H.theta[k] = theta_k_new;
                }
                else
                {
                    // reject
                    assertThis(fabs(H.theta[k] - theta_k_old) < EPS);
                }
            }
        }

        // update mus
        //
        for (int i = 0; i < H.N; i++)
        {
            double mu_i_new = proprnd_unbounded(H.mu[i], H, D, h);
            double mu_i_old = H.mu[i];

            double logpost_new = logpost_mu(mu_i_new, i, H, D, h);
            double logpost_old = H.LogPost(D, h); // TODO optim

            double logprop_new = logprop_unbounded(mu_i_new, mu_i_old, H, D, h);
            double logprop_old = logprop_unbounded(mu_i_old, mu_i_new, H, D, h);

            if (MetropolisHastingsFlip(logpost_new, logpost_old, logprop_new, logprop_old))
            { 
                // accept
                H.mu[i] = mu_i_new;
            }
            else
            {
                // reject
                assertThis(fabs(H.mu[i] - mu_i_old) < EPS);
            }
        }
        */

        // TODO bridges

        if (n >= burnin && n % lag == 0)
        {
            samples.push_back(new Hierarchy(H));
            post.push_back(H.LogPost(D, h)); // TODO optim
        }
    }

    return samples;
}


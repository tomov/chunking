// Hierarchy, Data, and Hyperparams classes
//
#include "printmex.h"
#include <string>
#include <memory>
#include <random>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <boost/math/distributions.hpp>

// TODO separate into .h and .cpp

// see https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
// for random number generation
std::random_device rd;  //Will be used to obtain a seed for the random number engine
//std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::mt19937 gen(0); // for reproducibility

// random draw U ~ Unif(a, b)
// see https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
// 
double UnifRnd(double a = 0, double b = 1)
{
    std::uniform_real_distribution<double> dis(a, b);
    double U = dis(gen);
    DEBUG_PRINT("unif = %lf\n", U);
    return U;
}

// random draw X ~ Beta(alpha, beta)
// Uses universality of the normal, i.e. F(X) ~ Unif(0,1), where X ~ CDF F
// see https://stackoverflow.com/questions/4181403/generate-random-number-based-on-beta-distribution-using-boost
//
double BetaRnd(double alpha = 1, double beta = 1)
{
    boost::math::beta_distribution<double> dist(alpha, beta);

    double U = UnifRnd();
    double X = quantile(dist, U);
    return X;
}

double BetaPDF(double x, double alpha = 1, double beta = 1)
{
    boost::math::beta_distribution<double> dist(alpha, beta);
    return boost::math::pdf(dist, x);
}


// random draw X ~ N(mu, sigma)
// see http://www.cplusplus.com/reference/random/normal_distribution/
//
double NormRnd(double mu, double sigma)
{
    std::normal_distribution<double> dist(mu, sigma);
    double X = dist(gen);
    return X;
}

double NormPDF(double x, double mu, double sigma)
{
    boost::math::normal_distribution<double> dist(mu, sigma);
    return boost::math::pdf(dist, x);
}

// random draw X ~ Cat(p), where p is a vector of (unnormalized) categorical probabilities
// see http://www.cplusplus.com/reference/random/discrete_distribution/
//
int CatRnd(const std::vector<double> &p)
{
    std::discrete_distribution<int> dist(p.begin(), p.end());
    int X = dist(gen);
    return X;
}


using namespace matlab::mex;
using namespace matlab::data;

class Data
{
    public:
        Data(StructArray const matlabStructArrayD);
        ~Data();

        struct Edge
        {
            int u, v;
            Edge(int _u, int _v) : u(_u), v(_v) {}
        };

        struct Task
        {
            int s, g;
            Task(int _s, int _g) : s(_s), g(_g) {}
        };
        
        struct Graph
        {
            int **E;
            int N;
            std::vector<Edge> edges;
            std::vector<int> *adj;
        };

        //std::string name;
        Graph G;
        std::vector<Task> tasks;
        std::vector<double> *rewards;
};


// TODO pass entryIndex and use instead of [0]
Data::Data(StructArray const matlabStructArrayD)
{
    // const TypedArray<char*> _name = matlabStructArrayD[0]["name"];

    // convert G
    //
    const StructArray matlabStructArrayG = matlabStructArrayD[0]["G"];
    const TypedArray<double> _N = matlabStructArrayG[0]["N"];
    const TypedArray<double> _E = matlabStructArrayG[0]["E"];
    const TypedArray<double> _edges = matlabStructArrayG[0]["edges"];

    G.N = (int)_N[0];
	DEBUG_PRINT("G.N = %d\n", G.N);

	DEBUG_PRINT("G.E = \n");
    G.E = new int*[G.N];
    G.adj = new std::vector<int>[G.N];
    for (int i = 0; i < G.N; i++)
    {
        G.E[i] = new int[G.N];
        for (int j = 0; j < G.N; j++)
        {
            G.E[i][j] = (int)_E[i][j];
            G.adj[i].push_back(j);
            DEBUG_PRINT("%d ", G.E[i][j]);
        }
        DEBUG_PRINT("\n");
    }

    for (int i = 0; i < _edges.getDimensions()[0]; i++)
    {
        int u = _edges[i][0];
        int v = _edges[i][1];
        G.edges.push_back(Edge(u, v));
        DEBUG_PRINT("G.edge %d %d\n", u, v);
    }

    // convert tasks  
    //
    const StructArray matlabStructArrayTasks = matlabStructArrayD[0]["tasks"];
    const TypedArray<double> _s = matlabStructArrayTasks[0]["s"];
    const TypedArray<double> _g = matlabStructArrayTasks[0]["g"];

    for (int i = 0; i < _s.getNumberOfElements(); i++)
    {
        int s = _s[i];
        int g = _g[i];
        tasks.push_back(Task(s, g));
        DEBUG_PRINT("task %d %d\n", s, g);
    }

    // convert rewards
    // 
    const CellArray matlabStructArrayRewards = matlabStructArrayD[0]["r"];
   
    rewards = new std::vector<double>[G.N];
    for (int i = 0; i < G.N; i++)
    {
        const TypedArray<double> _r = matlabStructArrayRewards[i];
        DEBUG_PRINT("r{%d} = [", i);
        for (int j = 0; j < _r.getNumberOfElements(); j++)
        {
            double r = _r[j];
            rewards[i].push_back(r);
            DEBUG_PRINT("%lf ", r);
        }
        DEBUG_PRINT("]\n");
    }
}

Data::~Data()
{
    for (int i = 0; i < G.N; i++)
    {
        delete [] G.E[i];
    }
    delete [] G.E;
    delete [] G.adj;
    delete [] rewards;
}


class Hyperparams
{
    public:
        Hyperparams(StructArray const matlabStructArrayHyperparams);

        double alpha;
        double std_theta;
        double theta_mean;
        double std_mu;
        double std_r;
};

Hyperparams::Hyperparams(StructArray const matlabStructArrayHyperparams)
{
    const TypedArray<double> _alpha = matlabStructArrayHyperparams[0]["alpha"];
    alpha = _alpha[0];

    const TypedArray<double> _std_theta = matlabStructArrayHyperparams[0]["std_theta"];
    std_theta = _std_theta[0];

    const TypedArray<double> _theta_mean = matlabStructArrayHyperparams[0]["theta_mean"];
    theta_mean = _theta_mean[0];

    const TypedArray<double> _std_mu = matlabStructArrayHyperparams[0]["std_mu"];
    std_mu = _std_mu[0];

    const TypedArray<double> _std_r = matlabStructArrayHyperparams[0]["std_r"];
    std_r = _std_r[0];

    DEBUG_PRINT("h = %lf %lf %lf %lf %lf\n", alpha, std_theta, theta_mean, std_mu, std_r);
}



class Hierarchy
{
    public:
        Hierarchy(int _N);
        Hierarchy(const Hierarchy &H); // copy constructor
        ~Hierarchy();

        void InitFromMATLAB(StructArray const matlabStructArrayH);
        void InitFromPrior(const Data &D, const Hyperparams &h);
        void PopulateCnt();

        void Print();

        double LogPrior(const Data &D, const Hyperparams &h) const;
        double LogLik(const Data &D, const Hyperparams &h) const;
        double LogPost(const Data &D, const Hyperparams &h) const;

        int N;
        int *c;
        std::vector<int> cnt;
        double p;
        double q;
        double hp; // p'
        double tp; // p''

        std::vector<double> theta;
        double *mu;
};

// populate cluster cnt counts from cluster assignments c
//
void Hierarchy::PopulateCnt()
{
    int K = *std::max_element(this->c, this->c + N); // # of clusters
    this->cnt.clear();
    this->cnt.resize(K);
    for (int i = 0; i < this->N; i++)
    {
        this->cnt[this->c[i] - 1]++;
    }
    DEBUG_PRINT("H.cnt = [");
    for (int k = 0; k < K; k++)
    {
        DEBUG_PRINT("%d ", this->cnt[k]);
    }
    DEBUG_PRINT("]\n");
}

// TODO pass entryIndex and use instead of [0]
void Hierarchy::InitFromMATLAB(StructArray const matlabStructArrayH)
{
    const TypedArray<double> _c = matlabStructArrayH[0]["c"];
    for (int i = 0; i < _c.getNumberOfElements(); i++)
    {
        this->c[i] = _c[i];
    }

    const TypedArray<double> _p = matlabStructArrayH[0]["p"];
    this->p = _p[0];

    const TypedArray<double> _q = matlabStructArrayH[0]["q"];
    this->q = _q[0];

    const TypedArray<double> _tp = matlabStructArrayH[0]["tp"];
    this->tp = _tp[0];

    const TypedArray<double> _hp = matlabStructArrayH[0]["hp"];
    this->hp = _hp[0];


    const TypedArray<double> _theta = matlabStructArrayH[0]["theta"];
    this->theta.clear();
    for (int k = 0; k < _theta.getNumberOfElements(); k++)
    {
        this->theta.push_back(_theta[k]);
    }
   
    const TypedArray<double> _mu = matlabStructArrayH[0]["mu"];
    for (int i = 0; i < _mu.getNumberOfElements(); i++)
    {
        this->mu[i] = _mu[i];
    }

    this->PopulateCnt();
}

void Hierarchy::Print()
{
    DEBUG_PRINT("H.c = [");
    for (int i = 0; i < this->N; i++)
    {
        DEBUG_PRINT("%d ", this->c[i]);
    }
    DEBUG_PRINT("]\n");

    DEBUG_PRINT("H.cnt = [");
    for (int k = 0; k < this->cnt.size(); k++)
    {
        DEBUG_PRINT("%d ", this->cnt[k]);
    }
    DEBUG_PRINT("]\n");

    DEBUG_PRINT("H.p q tp hp = [%lf %lf %lf %lf]\n", this->p, this->q, this->tp, this->hp);

    DEBUG_PRINT("H.theta = [");
    for (int k = 0; k < this->theta.size(); k++)
    {
        DEBUG_PRINT("%lf ", this->theta[k]);
    }
    DEBUG_PRINT("]\n");

    DEBUG_PRINT("H.mu = [");
    for (int i = 0; i < this->N; i++)
    {
        DEBUG_PRINT("%lf ", this->mu[i]);
    }
    DEBUG_PRINT("]\n");
}

// transpiled from init_H.m
// careful with off-by-ones
//
void Hierarchy::InitFromPrior(const Data &D, const Hyperparams &h)
{
    this->c[0] = 1;
    this->cnt.clear();
    this->cnt.push_back(1);
    for (int i = 1; i < D.G.N; i++)
    {
        std::vector<double> p(this->cnt.begin(), this->cnt.end()); // TODO optimize -- no need to copy, could be done in O(1)
        p.push_back(h.alpha);
        int c_new = CatRnd(p) + 1; // careful with off-by-one
        if (c_new - 1 >= this->cnt.size())
        {
            this->cnt.push_back(1);
        }
        else
        {
            this->cnt[c_new - 1]++;
        }
        this->c[i] = c_new;
    }
    std::random_shuffle(this->c, this->c + this->N);


    this->p = BetaRnd(1, 1);
    this->q = BetaRnd(1, 1);
    this->tp = BetaRnd(1, 1);
    this->hp = BetaRnd(1, 1);

    int K = this->cnt.size();
    this->theta.clear();
    for (int k = 0; k < K; k++)
    {
        this->theta.push_back(NormRnd(h.theta_mean, h.std_theta));
    }

    for (int i = 0; i < D.G.N; i++)
    {
        this->mu[i] = NormRnd(this->theta[this->c[i] - 1], h.std_mu);
    }
}


Hierarchy::Hierarchy(int _N)
{
    N = _N;
    c = new int[N];
    mu = new double[N];
}


Hierarchy::Hierarchy(const Hierarchy &H)
{
    N = H.N;
    cnt = H.cnt;
    theta = H.theta;
    p = H.p;
    q = H.q;
    hp = H.hp;
    tp = H.tp;

    c = new int[N];
    for (int i = 0; i < this->N; i++)
    {
        c[i] = H.c[i];
    }

    mu = new double[N];
    for (int i = 0; i < this->N; i++)
    {
        mu[i] = H.mu[i];
    }
}



Hierarchy::~Hierarchy()
{
    DEBUG_PRINT("deleting c for %p\n", this);
    delete [] c;
    DEBUG_PRINT("deleting mu for %p\n", this);
    delete [] mu;
    DEBUG_PRINT("deleted all! for %p\n", this);
}

double Hierarchy::LogPrior(const Data &D, const Hyperparams &h) const
{
    assert(D.G.N == this->N);

    double logP = 0;

    // cluster assignments
    //
    std::vector<int> cnt(this->cnt.size()); // temporary count
    cnt[this->c[0] - 1] = 1;
    for (int i = 1; i < this->N; i++)
    {
        int c = this->c[i];
        assert(c <= this->cnt.size());
        if (cnt[c - 1] == 0)
        {
            logP += log(h.alpha) - log(i + h.alpha);
        }
        else
        {
            logP += log(cnt[c - 1]) - log(i + h.alpha);
        }
        cnt[c - 1]++;
    }

    // TODO optimize by having beta for object
    // TODO or marginalize over them
    // TODO actually jk just rm them
    logP += log(BetaPDF(this->p, 1, 1)) + log(BetaPDF(this->q, 1, 1)) + log(BetaPDF(this->tp, 1, 1)) + log(BetaPDF(this->hp, 1, 1));// TODO const

    // cluster rewards
    //
    assert(this->cnt.size() == this->theta.size());
    for (int k = 0; k < this->theta.size(); k++)
    {
        // TODO optimize with norm dist objects for each k for H
        logP += log(NormPDF(this->theta[k], h.theta_mean, h.std_theta));
    }

    // state rewards
    //
    for (int i = 0; i < this->N; i++)
    {
        logP += log(NormPDF(this->mu[i], this->theta[this->c[i] - 1], h.std_mu));
    }

    // prevent -Infs = impossible events; equivalent to using a Gaussian + uniform mixture
    //
    if (isinf(logP))
    {
        logP = 1e-100;
    }

    return logP;
}

double Hierarchy::LogLik(const Data &D, const Hyperparams &h) const
{
    assert(D.G.N == this->N);

    double logP = 0;

    // connectivity
    //
    for (int i = 0; i < D.G.N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (this->c[i] == this->c[j])
            {
                if (D.G.E[i][j])
                {
                    logP += log(this->p);
                }
                else 
                {
                    logP += log(1 - this->p);
                }
            }
            else
            {
                if (D.G.E[i][j])
                {
                    logP += log(this->p * this->q);
                }
                else 
                {
                    logP += log(1 - this->p * this->q);
                }
            }
        }
    }

	DEBUG_PRINT("at 1 -> %.6lf\n", logP);

    // transitive closures
    // TODO optimize the crap out of this
    //
    int A[this->N][this->N];
    for (int i = 0; i < D.G.N; i++)
    {
        for (int j = 0; j < D.G.N; j++)
        {
            if (this->c[i] == this->c[j])
            {
                A[i][j] = D.G.E[i][j];
            }
        }
    }
    for (int k = 0; k < D.G.N; k++) 
    {
        for (int i = 0; i < D.G.N; i++)
        {
            if (this->c[i] != this->c[k])
            {
                continue;
            }
            for (int j = 0; j < D.G.N; j++)
            {
                if (this->c[i] != this->c[j])
                {
                    continue;
                }
                A[i][j] = A[i][j] || (A[i][k] && A[k][j]);
            }
        }
    }

    // penalize disconnected chunks
    //
    for (int i = 0; i < D.G.N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (this->c[i] == this->c[j] && !A[i][j])
            {
                logP -= 100;
            }
        }
    }

	DEBUG_PRINT("at 2 -> %.6lf\n", logP);

    // get_H_E
    int K = this->cnt.size();
    int E[K][K];
    memset(E, 0, sizeof(E));
    for (int i = 0; i < this->N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (this->c[i] != this->c[j] && D.G.E[i][j])
            {
                DEBUG_PRINT("hierarchical edge (%d %d)\n", this->c[i], this->c[j]);
                E[this->c[i] - 1][this->c[j] - 1] = 1;
                E[this->c[j] - 1][this->c[i] - 1] = 1;
            }
        }
    }

    // (hierarchical) edges
    // TODO sample them too, or marginalize over them
    //
    DEBUG_PRINT(" this->hp = %.4lf\n", this->hp);
    for (int k = 0; k < K; k++)
    {
        if (this->cnt[k] == 0)
        {
            continue;
        }
        for (int l = 0; l < k; l++)
        {
            if (this->cnt[l] == 0)
            {
                continue;
            }
            if (E[k][l])
            {
                logP += log(this->hp);
                DEBUG_PRINT(" %.4lf for edge (%d %d)\n", log(this->hp), k + 1, l + 1);
            }
            else
            {
                logP += log(1 - this->hp);
                DEBUG_PRINT(" %.4lf for edge -(%d %d)\n", log(1 - this->hp), k + 1, l + 1);
            }
        }
    }

	DEBUG_PRINT("at 3 -> %.6lf\n", logP);

    // bridges
    //
    for (int k = 0; k < K; k++)
    {
        if (this->cnt[k] == 0)
        {
            continue;
        }
        for (int l = 0; l < k; l++)
        {
            if (this->cnt[l] == 0)
            {
                continue;
            }
            if (E[k][l])
            {
                logP += log(1) - log(this->cnt[k]) - log(this->cnt[l]);
                logP -= log(this->p * this->q); // b/c the bridge is always there, but we penalized / overcounted the corresponding edge when accounting for the connectivity of G
            }
        }
    }

	DEBUG_PRINT("at 4 -> %.6lf\n", logP);

    // tasks
    //
    for (int i = 0; i < D.tasks.size(); i++)
    {
        int s = D.tasks[i].s;
        logP += log(1) - log(this->N);

        int g = D.tasks[i].g;
        if (this->c[s] == this->c[g])
        {
            logP += log(1);
        }
        else
        {
            logP += log(this->tp);
        }
        double denom = this->cnt[this->c[s] - 1] + (this->N - this->cnt[this->c[s] - 1]) * this->tp;
        logP -= log(denom);
    }

	DEBUG_PRINT("at 5 -> %.6lf\n", logP);

    // rewards
    //
    for (int i = 0; i < this->N; i++)
    {
        for (int o = 0; o < D.rewards[i].size(); o++)
        {
            logP += log(NormPDF(D.rewards[i][o], this->mu[i], h.std_r));
        }
    }

	DEBUG_PRINT("at 6 -> %.6lf\n", logP);

    return logP;
}


double Hierarchy::LogPost(const Data &D, const Hyperparams &h) const
{
    double logP = this->LogPrior(D, h) + this->LogLik(D, h);
    return logP;
}

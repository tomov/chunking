// Hierarchy, Data, and Hyperparams classes
//

#ifndef DATA_STRUCTS_H 
#define DATA_STRUCTS_H

#include <string>
#include <memory>
#include <random>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <boost/math/distributions.hpp>

#include "helpermex.h"
#include "printmex.h"

// TODO separate into .h and .cpp


// see https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
// for random number generation
std::random_device rd;  //Will be used to obtain a seed for the random number engine
//std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::mt19937 gen(0); // for reproducibility

const double EPS = 1e-9;

// random draw U ~ Unif(a, b)
// see https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
// 
double UnifRnd(double a = 0, double b = 1)
{
    std::uniform_real_distribution<double> dis(a, b);
    double U = dis(gen);
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

double NormCDF(double x, double mu, double sigma)
{
    boost::math::normal_distribution<double> dist(mu, sigma);
    return boost::math::cdf(dist, x);
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
        // types -- see https://www.mathworks.com/help/matlab/apiref/matlab.data.arraytype.html
        static const std::vector<std::string> fieldNames;
        static const std::vector<ArrayType> fieldTypes;

        static void check(StructArray const &matlabStructArrayD, std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr);

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
            static const std::vector<std::string> fieldNames; // D.G
            static const std::vector<ArrayType> fieldTypes;

            int **E;
            int **hidden_E;
            int N;
            std::vector<Edge> edges;
            std::vector<Edge> hidden_edges;
            std::vector<int> *adj;
        };

        //std::string name;
        Graph G;
        std::vector<Task> tasks;
        std::vector<double> *rewards;
};

const std::vector<std::string> Data::fieldNames = {"name", "G", "tasks", "r"};
const std::vector<ArrayType> Data::fieldTypes = {ArrayType::CHAR, ArrayType::STRUCT, ArrayType::STRUCT, ArrayType::CELL};

const std::vector<std::string> Data::Graph::fieldNames = {"N", "E", "edges", "hidden_E", "hidden_edges"}; // D.G
const std::vector<ArrayType> Data::Graph::fieldTypes = {ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE};

void Data::check(StructArray const &matlabStructArrayD, std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr)
{
    checkStructureElements(matlabStructArrayD, "D", Data::fieldNames, Data::fieldTypes, matlabPtr);

    size_t total_num_of_elements = matlabStructArrayD.getNumberOfElements();
    for (size_t i=0; i<total_num_of_elements; i++) 
    {
        // check D.G
        const StructArray structFieldG = matlabStructArrayD[i]["G"];
        checkStructureElements(structFieldG, "D.G", Data::Graph::fieldNames, Data::Graph::fieldTypes, matlabPtr);

        const TypedArray<double> _N = structFieldG[0]["N"];
        const TypedArray<double> _E = structFieldG[0]["E"];
        int N = (int)_N[0];
        if (_E.getNumberOfElements() != N * N)
        {
            displayError("D.G.E must have D.G.N^2 elements.", matlabPtr);
        }

        // check D.tasks
        const StructArray matlabStructArrayTasks = matlabStructArrayD[i]["tasks"];
        const TypedArray<double> _s = matlabStructArrayTasks[0]["s"];
        const TypedArray<double> _g = matlabStructArrayTasks[0]["g"];
        if (_s.getNumberOfElements() != _g.getNumberOfElements())
        {
            displayError("D.tasks.s and D.tasks.g must have the same number of elements.", matlabPtr);
        }

        // check D.r
        const CellArray matlabStructArrayRewards = matlabStructArrayD[i]["r"];
        if (matlabStructArrayRewards.getNumberOfElements() != N)
        {
            displayError("D.r should have D.N elements", matlabPtr);
        }
    }
}

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
    const TypedArray<double> _hidden_E = matlabStructArrayG[0]["hidden_E"];
    const TypedArray<double> _hidden_edges = matlabStructArrayG[0]["hidden_edges"];

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

	DEBUG_PRINT("G.hidden_E = \n");
    G.hidden_E = new int*[G.N];
    G.adj = new std::vector<int>[G.N];
    for (int i = 0; i < G.N; i++)
    {
        G.hidden_E[i] = new int[G.N];
        for (int j = 0; j < G.N; j++)
        {
            G.hidden_E[i][j] = (int)_hidden_E[i][j];
            // important for connectivity that hidden edges are not present in E
            ASSERT(!G.hidden_E[i][j] || !G.E[i][j], "!G.hidden_E[i][j] || !G.E[i][j]");
            DEBUG_PRINT("%d ", G.hidden_E[i][j]);
        }
        DEBUG_PRINT("\n");
    }

    for (int i = 0; i < _hidden_edges.getDimensions()[0]; i++)
    {
        int u = _hidden_edges[i][0];
        int v = _hidden_edges[i][1];
        G.hidden_edges.push_back(Edge(u, v));
        DEBUG_PRINT("G.hidden_edge %d %d\n", u, v);
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
        delete [] G.hidden_E[i];
    }
    delete [] G.E;
    delete [] G.adj;
    delete [] rewards;
}


class Hyperparams
{
    public:
        static const std::vector<std::string> fieldNames; // h
        static const std::vector<ArrayType> fieldTypes;

        static void check(StructArray const &matlabStructArrayHyperparams, std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr);

        Hyperparams(StructArray const matlabStructArrayHyperparams);

        double alpha;
        double std_theta;
        double theta_mean;
        double std_mu;
        double std_r;
        double eps;
};

const std::vector<std::string> Hyperparams::fieldNames = {"alpha", "std_theta", "theta_mean", "std_mu", "std_r", "eps"}; // h
const std::vector<ArrayType> Hyperparams::fieldTypes = {ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE};

void Hyperparams::check(StructArray const &matlabStructArrayHyperparams, std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr)
{
    checkStructureElements(matlabStructArrayHyperparams, "h", Hyperparams::fieldNames, Hyperparams::fieldTypes, matlabPtr);
}

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

    const TypedArray<double> _eps = matlabStructArrayHyperparams[0]["eps"];
    eps = _eps[0];

    DEBUG_PRINT("h = %lf %lf %lf %lf %lf %lf\n", alpha, std_theta, theta_mean, std_mu, std_r, eps);
}



class Hierarchy
{
    public:
        // types -- see https://www.mathworks.com/help/matlab/apiref/matlab.data.arraytype.html
        static const std::vector<std::string> fieldNames; // H
        static const std::vector<ArrayType> fieldTypes;

        static void check(StructArray const &matlabStructArrayH, const Data &D, std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr);

        Hierarchy(int _N);
        Hierarchy(const Hierarchy &H); // copy constructor
        ~Hierarchy();

        bool Equals(const Hierarchy& H) const;
        void InitFromMATLAB(StructArray const matlabStructArrayH);
        void InitFromPrior(const Data &D, const Hyperparams &h);
        void PopulateCnt();

        void Print() const;

        double LogPrior(const Data &D, const Hyperparams &h) const;
        double LogLik(const Data &D, const Hyperparams &h) const;
        double LogPost(const Data &D, const Hyperparams &h) const;

        double LogPost_c_i(int c_i_new, int i, const Data &D, const Hyperparams &h);
        void Update_c_i(int c_i_new, int i, const Data &D, const Hyperparams &h);
        void Update_c_i(int c_i_new, int i, const Data &D, const Hyperparams &h, int /*out*/ &c_i_old, double /*out*/ &theta_old, std::vector<int> /*out*/ &E_old);
        void Undo_c_i(int c_i_new, int i, const Data &D, const Hyperparams &h, int c_i_old, double theta_old, const std::vector<int> &E_old);
        
        void Sanity(const Data &D, const Hyperparams &h);
        void Sanity();

        int N;
        int **E;
        int *c;
        std::vector<int> cnt;
        double p;
        double q;
        double hp; // p'
        double tp; // p''

        std::vector<double> theta;
        double *mu;
};

const std::vector<std::string> Hierarchy::fieldNames = {"c", "p", "q", "tp", "hp", "theta", "mu", "E"}; // H
const std::vector<ArrayType> Hierarchy::fieldTypes = {ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE};

void Hierarchy::check(StructArray const &matlabStructArrayH, const Data &D, std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr)
{
    checkStructureElements(matlabStructArrayH, "H", Hierarchy::fieldNames, Hierarchy::fieldTypes, matlabPtr);

    const TypedArray<double> _c = matlabStructArrayH[0]["c"];
    if (_c.getNumberOfElements() != D.G.N)
    {
        displayError("H.c should have D.G.N elements", matlabPtr);
    }

    const TypedArray<double> _theta = matlabStructArrayH[0]["theta"];
    // TODO this is wrong; needs to be max(c)
    // also, wtf try to pass Hout as input argument -> Busy
    //if (_theta.getNumberOfElements() != D.G.N)
    //{
    //    displayError("H.theta should have D.G.N elements");
    //}

    const TypedArray<double> _mu = matlabStructArrayH[0]["mu"];
    if (_mu.getNumberOfElements() != D.G.N)
    {
        displayError("H.mu should have D.G.N elements", matlabPtr);
    }

    const TypedArray<double> _E = matlabStructArrayH[0]["E"];
    if (_E.getNumberOfElements() != D.G.N * D.G.N)
    {
        displayError("H.E must have D.G.N^2 elements.", matlabPtr);
    }
}

// P(H|D) for updates of c_i
// i.e. with new c's up to c_i, the candidate c_i, then old c's after (and old rest of H)
//
double Hierarchy::LogPost_c_i(int c_i_new, int i, const Data &D, const Hyperparams &h) // notice it's not const b/c it temporarily changes H for efficiency
{
#if (DEBUG)
    Hierarchy H(*this);
#endif

    int c_i_old;
    double theta_old;
    std::vector<int> E_old;
    this->Update_c_i(c_i_new, i, D, h, c_i_old, theta_old, E_old);

    double logP = this->LogPost(D, h); // TODO much more efficiently, maybe

    this->Undo_c_i(c_i_new, i, D, h, c_i_old, theta_old, E_old);

#if (DEBUG)
    ASSERT(this->Equals(H), "this->Equals(H)"); // TODO rm in prod
#endif
    this->Sanity(D, h);

    return logP;
}


// set c[i] = c_i_new and update stuff accordingly
void Hierarchy::Update_c_i(int c_i_new, int i, const Data &D, const Hyperparams &h)
{
    int c_i_old;
    double theta; // dummies
    std::vector<int> E_old;
    Update_c_i(c_i_new, i, D, h, c_i_old, theta, E_old);
}

// set c[i] = c_i_new and update stuff accordingly
// also returns stuff that can redo the op
// careful with off-by-one everywhere
// 
void Hierarchy::Update_c_i(int c_i_new, int i, const Data &D, const Hyperparams &h, int /*out*/ &c_i_old, double /*out*/ &theta_old, std::vector<int> /*out*/ &E_old)
{
    c_i_old = this->c[i];
    this->c[i] = c_i_new;

    ASSERT(c_i_old - 1 >= 0, "c_i_old - 1 >= 0, Update_c_i");
    ASSERT(c_i_old - 1 < this->cnt.size(), "c_i_old - 1 < this->cnt.size(), Update_c_i");
    ASSERT(this->cnt[c_i_old - 1] > 0, "this->cnt[c_i_old - 1] > 0, Update_c_i");
    this->cnt[c_i_old - 1]--;

    // removed a singleton cluster TODO make theta logic similar
    E_old.clear();
    if (this->cnt[c_i_old - 1] == 0)
    {
        for (int k = 0; k < D.G.N; k++)
        {
            ASSERT(this->E[c_i_old - 1][k] == this->E[k][c_i_old - 1], "this->E[c_i_old - 1][k] == this->E[k][c_i_old - 1]");
            E_old.push_back(this->E[c_i_old - 1][k]);
            this->E[c_i_old - 1][k] = -1;
            this->E[k][c_i_old - 1] = -1;
        }
    }

    ASSERT(c_i_new - 1 >= 0, "c_i_new - 1 >= 0, Update_c_i");
    ASSERT(c_i_new - 1 <= this->cnt.size(), "c_i_new - 1 <= this->cnt.size(), Update_c_i"); // note we allow equality

    // creating new cluster
    if (c_i_new - 1 == this->cnt.size())
    {
        this->cnt.push_back(0);
        this->theta.push_back(nan(""));
    }

    this->cnt[c_i_new - 1]++;

    theta_old = this->theta[c_i_new - 1];

    if (this->cnt[c_i_new - 1] == 1)
    {
        // created a new cluster, so create a new theta -- note that we could be reusing an old one, so we take care to save the theta in case we need to undo this, e.g. when computing the MCMC updates
        //this->theta[c_i_new - 1] = NormRnd(h.theta_mean, h.std_theta); // InitFromPrior 
        this->theta[c_i_new - 1] = this->mu[i]; // sort-of empirical prior TODO is this legit? 
        ASSERT(this->cnt.size() <= D.G.N, "this->cnt.size() <= D.G.N");
        for (int k = 0; k < this->cnt.size(); k++)
        {
            if (this->cnt[k] > 0)
            {
                this->E[c_i_new - 1][k] = UnifRnd() < this->p;
                this->E[k][c_i_new - 1] = this->E[c_i_new - 1][k];
            }
        }
        this->E[c_i_new - 1][c_i_new - 1] = 0;
    }

    DEBUG_PRINT(" update_c_i -- c_i_new = %d, i = %d; c_i_old = %d\n", c_i_new, i, c_i_old);
    this->Print();
    this->Sanity(D, h); // TODO DEBUG only
}

// undo Update_c_i; notice we go in reverse order
// careful with off-by-one everywhere
//
void Hierarchy::Undo_c_i(int c_i_new, int i, const Data &D, const Hyperparams &h, int c_i_old, double theta_old, const std::vector<int> &E_old)
{
    ASSERT(c_i_new - 1 >= 0, "c_i_new - 1 >= 0, Undo_c_i");
    ASSERT(c_i_new - 1 < this->cnt.size(), "c_i_new - 1 < this->cnt.size(), Undo_c_i"); // notice strict < here
    if (this->cnt[c_i_new - 1] == 1)
    {
        this->theta[c_i_new - 1] = theta_old;

        ASSERT(c_i_new - 1 < D.G.N, "c_i_new - 1 < D.G.N");
        for (int k = 0; k < D.G.N; k++)
        {
            this->E[k][c_i_new - 1] = -1;
            this->E[c_i_new - 1][k] = -1;
        }
    }

    ASSERT(this->cnt[c_i_new - 1] > 0, "this->cnt[c_i_new - 1] > 0, Undo_c_i");
    this->cnt[c_i_new - 1]--;

    if (c_i_new == this->cnt.size() && isnan(theta_old))
    {
        // we added a new cluster -> remove it
        this->cnt.pop_back();
        this->theta.pop_back();
    }

    ASSERT(c_i_old - 1 >= 0, "c_i_old - 1 >= 0, Undo_c_i");
    ASSERT(c_i_old - 1 < this->cnt.size(), "c_i_old - 1 < this->cnt.size(), Undo_c_i");
    if (this->cnt[c_i_old - 1] == 0)
    {
        ASSERT(E_old.size() == D.G.N, "E_old.size() == D.G.N");
        for (int k = 0; k < D.G.N; k++)
        {
            this->E[c_i_old - 1][k] = E_old[k];
            this->E[k][c_i_old - 1] = E_old[k];
        }
    }

    this->cnt[c_i_old - 1]++;

    this->c[i] = c_i_old;

    DEBUG_PRINT(" undo_c_i -- c_i_new = %d, i = %d; c_i_old = %d\n", c_i_new, i, c_i_old);
    this->Print();

    this->Sanity(D, h); // TODO DEBUG only
}


void Hierarchy::Sanity(const Data &D, const Hyperparams &h)
{
#if DEBUG == 1
    ASSERT(D.G.N == this->N, "D.G.N == this->N, Sanity");
    this->Sanity();
#endif
}

void Hierarchy::Sanity()
{
#if DEBUG == 1
    ASSERT(this->cnt.size() == this->theta.size(), "this->cnt.size() == this->theta.size(), Sanity");

    int K = *std::max_element(this->c, this->c + N); // # of clusters
    ASSERT(this->cnt.size() >= K, "this->cnt.size() >= K, Sanity");

    std::vector<int> cnt(this->cnt.size());
    for (int i = 0; i < this->N; i++)
    {
        cnt[this->c[i] - 1]++;
    }

    int sum = 0;
    for (int k = 0; k < K; k++)
    {
        sum += this->cnt[k];
        ASSERT(this->cnt[k] == cnt[k], "this->cnt[k] == cnt[k], Sanity");
        ASSERT(!isnan(this->theta[k]), "!isnan(this->theta[k])");
    }
    ASSERT(sum == this->N, "sum == this->N, Sanity");

    for (int i = 0; i < this->N; i++)
    {
        ASSERT(this->c[i] - 1 >= 0, "this->c[i] - 1 >= 0, Sanity");
        ASSERT(this->c[i] - 1 < this->theta.size(), "this->c[i] - 1 < this->theta.size(), Sanity");
        ASSERT(this->c[i] - 1 < this->cnt.size(), "this->c[i] - 1 < this->cnt.size(), Sanity");
    }

    for (int k = 0; k < this->N; k++)
    {
        for (int l = 0; l <= k; l++)
        {
            ASSERT(this->E[k][l] == this->E[l][k], "this->E[k][l] == this->E[l][k]");
            if (k < this->cnt.size() && this->cnt[k] > 0 && l < this->cnt.size() && this->cnt[l] > 0)
            {
                ASSERT(this->E[k][l] != -1, "this->E[k][l] != -1, Sanity");
            }
        }
    }
#endif
}


bool Hierarchy::Equals(const Hierarchy& H) const
{
    if (this->N != H.N) return false;
    if (fabs(this->p - H.p) > EPS) return false;
    if (fabs(this->q - H.q) > EPS) return false;
    if (fabs(this->tp - H.tp) > EPS) return false;
    if (fabs(this->hp - H.hp) > EPS) return false;

    for (int i = 0; i < this->N; i++)
    {
        if (this->c[i] != H.c[i]) return false;
        if (fabs(this->mu[i] - H.mu[i]) > EPS) return false;
    }

    for (int k = 0; k < this->cnt.size(); k++)
    {
        if (this->cnt[k] != H.cnt[k]) return false;
        if (fabs(this->theta[k] - H.theta[k]) > EPS) return false;
    }

    for (int k = 0; k < this->N; k++)
    {
        for (int l = 0; l < this->N; l++)
        {
            if (this->E[k][l] != H.E[k][l])
            {
                return false;
            }
        }
    }

    return true;
}


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
    // pad with zeros in case there were empty clusters; this occurs
    // e.g. if we return H to MATLAB after sampling, which erases this->cnt
    //
    // EDIT: DON'T bad -- might screw up the prior b/c it will count a bunch of useless thetas
    //while (this->cnt.size() < this->theta.size())
    //{
    //    this->cnt.push_back(0);
    //}
    // remove thetas for empty clusters; happens b/c of MATLAB
    // agni
    if (this->theta.size() > K)
    {
        this->theta.resize(K);
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

    const TypedArray<double> _E = matlabStructArrayH[0]["E"];
    for (int k = 0; k < this->N; k++)
    {
        for (int l = 0; l < this->N; l++)
        {
            this->E[k][l] = _E[k][l];
        }
    }

    this->PopulateCnt();

    this->Sanity();
}

void Hierarchy::Print() const
{
#if (DEBUG)
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

    DEBUG_PRINT("H.E = \n");
    for (int k = 0; k < this->N; k++)
    {
        for (int l = 0; l < this->N; l++)
        {
            DEBUG_PRINT("%d ", this->E[k][l]);
        }
        DEBUG_PRINT("\n");
    }
#endif
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
        ASSERT(this->c[i] - 1 >= 0, "this->c[i] - 1 >= 0, InitFromPrior");
        ASSERT(this->c[i] - 1 < this->theta.size(), "this->c[i] - 1 < this->theta.size(), InitFromPrior");
        this->mu[i] = NormRnd(this->theta[this->c[i] - 1], h.std_mu);
    }

    for (int k = 0; k < D.G.N; k++)
    {
        for (int l = 0; l < k; l++)
        {
            if (k >= this->cnt.size() || l >= this->cnt.size() || this->cnt[k] == 0 || this->cnt[l] == 0)
            {
                this->E[k][l] = -1;
            }
            else
            {
                this->E[k][l] = UnifRnd() < this->p;
            }
            this->E[l][k] = this->E[k][l];
        }
        if (k >= this->cnt.size() || this->cnt[k] == 0)
        {
            this->E[k][k] = -1;
        }
        else
        {
            this->E[k][k] = 0;
        }
    }

    this->Sanity(D, h);
}


Hierarchy::Hierarchy(int _N)
{
    N = _N;
    c = new int[N];
    mu = new double[N];
    E = new int*[N];
    for (int k = 0; k < N; k++)
    {
        E[k] = new int[N];
    }
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

    E = new int*[N];
    for (int k = 0; k < N; k++)
    {
        E[k] = new int[N];
        for (int l = 0; l < N; l++)
        {
            E[k][l] = H.E[k][l];
        }
    }
}



Hierarchy::~Hierarchy()
{
    delete [] c;
    delete [] mu;
    for (int k = 0; k < N; k++)
    {
        delete [] E[k];
    }
    delete [] E;
}

double Hierarchy::LogPrior(const Data &D, const Hyperparams &h) const
{
    ASSERT(D.G.N == this->N, "D.G.N == this->N, LogPrior");

    double logP = 0;

    // cluster assignments
    //
    std::vector<int> cnt(this->cnt.size()); // temporary count
    ASSERT(this->c[0] - 1 < cnt.size(), "this->c[0] - 1 < cnt.size()");
    cnt[this->c[0] - 1] = 1;
    for (int i = 1; i < this->N; i++)
    {
        int c = this->c[i];
        ASSERT(c - 1 < this->cnt.size(), "c - 1 < this->cnt.size()");
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
    // TODO restore when hyperprior is introduced; right now this is just 0
    // logP += log(BetaPDF(this->p, 1, 1)) + log(BetaPDF(this->q, 1, 1)) + log(BetaPDF(this->tp, 1, 1)) + log(BetaPDF(this->hp, 1, 1)); // TODO const


    // hierarchical edges
    // TODO same problem as thetas -- more clusters just additionally reduces the prior???
    //
    int K = this->cnt.size();
    ASSERT(K <= D.G.N, "K <= D.G.N");
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
            ASSERT(this->E[k][l] != -1, "this->E[k][l] != -1, LogPrior");
            if (this->E[k][l])
            {
                logP += log(this->hp);
            }
            else
            {
                logP += log(1 - this->hp);
            }
        }
    }


    // cluster rewards
    //
    // TODO discuss w/ Sam 
    /*
    ASSERT(this->cnt.size() == this->theta.size(), "this->cnt.size() == this->theta.size()");
    for (int k = 0; k < this->theta.size(); k++)
    {
        // don't do for empty clusters TODO think about it more carefully TODO might be a bug in MATLAB version
        //if (this->cnt[k] > 0) // agni <------------- OMG IT'S THIS!!!!! that screws up mines10 ...
        {
            // TODO optimize with norm dist objects for each k for H
            //DEBUG_PRINT("theta k [%d] = %.4lf\n", k, this->theta[k]);
            double p = log(NormPDF(this->theta[k], h.theta_mean, h.std_theta));
            if (isinf(p))
            {
                // prevent -Infs = impossible events; equivalent to using a
                // Gaussian + uniform mixture
                // do it in a "soft" way so MCMC can recover one by one
                //
                logP += 1e-100;
            }
            else
            {
                logP += p;
            }
            //logP += p; // agni

            DEBUG_PRINT("N(theta[%d]): logp += %e = %e\n", k + 1, p, logP);
        }
    }
    */

    // state rewards
    //
    for (int i = 0; i < this->N; i++)
    {
        ASSERT(this->c[i] - 1 >= 0, "this->c[i] - 1 >= 0");
        ASSERT(this->c[i] - 1 < this->theta.size(), "this->c[i] - 1 < this->theta.size()");
        //DEBUG_PRINT("mu i [%d] = %.4lf\n", i, this->mu[i]);
        double p = log(NormPDF(this->mu[i], this->theta[this->c[i] - 1], h.std_mu));

        if (isinf(p))
        {
            // prevent -Infs = impossible events; equivalent to using a
            // Gaussian + uniform mixture
            // do it in a "soft" way so MCMC can recover one by one
            //
            logP += 1e-100;
        }
        else
        {
            logP += p;
        }

        DEBUG_PRINT("N(mu[%d], theta[%d]): logp += %e = %e\n", i+1, this->c[i], p, logP);
    }

    ASSERT(!isinf(logP), "!isinf(logP) in LogPrior");

    return logP;
}

double Hierarchy::LogLik(const Data &D, const Hyperparams &h) const
{
    ASSERT(D.G.N == this->N, "D.G.N == this->N, LogLik");

    double logP = 0;

    // connectivity
    //
    for (int i = 0; i < D.G.N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (D.G.hidden_E[i][j])
            {
                // hidden edge -> don't count
                continue;
            }
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
                logP += -100;
            }
        }
    }

    // get_H_E and H_hidden_E 
    int K = this->cnt.size();
    int E[K][K];
    memset(E, 0, sizeof(E));
    int hidden_E[K][K];
    memset(hidden_E, 0, sizeof(hidden_E));
    for (int i = 0; i < this->N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (this->c[i] != this->c[j] && D.G.E[i][j])
            {
                E[this->c[i] - 1][this->c[j] - 1]++;
                E[this->c[j] - 1][this->c[i] - 1]++;
            }
            if (this->c[i] != this->c[j] && D.G.hidden_E[i][j])
            {
                hidden_E[this->c[i] - 1][this->c[j] - 1]++;
                hidden_E[this->c[j] - 1][this->c[i] - 1]++;
            }
        }
    }

    // MOVED TO PRIOR -- sampled properly now TODO remove competely
    /*
    // (hierarchical) edges
    // TODO sample them too, or marginalize over them
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
                logP += log(this->hp);
            }
            else
            {
                logP += log(1 - this->hp);
            }
        }
    }
    */

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
            if (this->E[k][l]) // there's a hierarchical edge between clusters k and l
            {
                if (E[k][l]) // there are actual edges between them => one of them is a bridge
                {
                    logP += log(1) - log(this->cnt[k]) - log(this->cnt[l]);
                    logP += - log(this->p * this->q); // b/c the bridge is always there, but we penalized / overcounted the corresponding edge when accounting for the connectivity of G
                }
                else if (hidden_E[k][l]) // there are no visible edges but there are hidden/unobserved edges = potential edges
                {
                    logP += log(1) - log(this->cnt[k]) - log(this->cnt[l]); // count the bridge but don't adjust, b/c we don't penalize hidden edges
                }
                else  // no visible nor hidden edges between clusters => bad; penalize a lot
                {
                    logP += -100;
                }
            }
        }
    }


    // tasks
    //
    for (int i = 0; i < D.tasks.size(); i++)
    {
        int s = D.tasks[i].s;
        ASSERT(s > 0 && s <= this->N, "s > 0 && s <= this->N");
        logP += log(1) - log(this->N);

        int g = D.tasks[i].g;
        ASSERT(g > 0 && g <= this->N, "g > 0 && g <= this->N");
        if (this->c[s - 1] == this->c[g - 1]) // off by one!
        {
            logP += log(1);
        }
        else
        {
            logP += log(this->tp);
        }
        double denom = this->cnt[this->c[s - 1] - 1] * 1 + (this->N - this->cnt[this->c[s - 1] - 1]) * this->tp; // careful with double off-by-ones
        logP -= log(denom);
    }

    // rewards
    //
    for (int i = 0; i < this->N; i++)
    {
        for (int o = 0; o < D.rewards[i].size(); o++)
        {
            // Pr(r = x | rest of H)
            double p = log(NormPDF(D.rewards[i][o], this->mu[i], h.std_r));
            if (isinf(p))
            {
                // prevent -Infs = impossible events; equivalent to using a
                // Gaussian + uniform mixture
                // do it in a "soft" way so MCMC can recover one by one
                //
                logP += 1e-100;
            }
            else
            {
                logP += p;
            }
        }
    }

    ASSERT(!isinf(logP), "!isinf(logP) in LogLik");

    return logP;
}


double Hierarchy::LogPost(const Data &D, const Hyperparams &h) const
{
    double logP = this->LogPrior(D, h) + this->LogLik(D, h);
    return logP;
}

#endif

// compile with:
//
// mex testmex.cpp printmex.cpp -I/usr/local/boost-1.64.0/include/
//
// see https://stackoverflow.com/questions/16127060/what-is-the-default-location-for-boost-library-when-installed-using-macport-on-m
// and https://www.mathworks.com/matlabcentral/answers/7955-using-boost-libraries-with-mex-function-in-matlab

/* ========================================================================
 *  copy of 
 * phonebook.cpp
 * example for illustrating how to manipulate structure.
 *
 * takes a (MxN) structure matrix which has first field as 
 * character array(name), and second field as scalar double (phone number). 
 * This function returns a new structure (1x1)containing following fields: 
 * for character array input, it will be (MxN) cell array; 
 * and for numeric double (noncomplex, scalar) input, it will be (MxN)
 * cell array where each field is numeric array of type double.
 *
 * Build : from MATLAB
 *         >> mex phonebook.cpp
 * Usage with example : from MATLAB
 *         >> friends(1).name = 'Jordan Robert';
 *         >> friends(1).phone = 3386;
 *         >> friends(2).name = 'Mary Smith';
 *         >> friends(2).phone = 3912;
 *         >> friends(3).name = 'Stacy Flora';
 *         >> friends(3).phone = 3238;
 *         >> friends(4).name = 'Harry Alpert';
 *         >> friends(4).phone = 3077;
 *         >> phonebook(friends)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2017 The MathWorks, Inc.
 *=======================================================================*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "printmex.h"
#include <string>
#include <memory>
#include <random>
#include <algorithm>
#include <boost/math/distributions.hpp>

#define DEBUG 1
#define DEBUG_PRINT(args ...) if (DEBUG) printThis(args)

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
    boost::math::beta_distribution<> dist(alpha, beta);

    double U = UnifRnd();
    double X = quantile(dist, U);
    return X;
}

// random draw X ~ N(mu, sigma)
// see http://www.cplusplus.com/reference/random/normal_distribution/
//
double NormRnd(double mu, double sigma)
{
    std::normal_distribution<double> dis(mu, sigma);
    double X = dis(gen);
    return X;
}

// random draw X ~ Cat(p), where p is a vector of (unnormalized) categorical probabilities
// see http://www.cplusplus.com/reference/random/discrete_distribution/
//
int CatRnd(const std::vector<double> &p)
{
    std::discrete_distribution<int> dis(p.begin(), p.end());
    int X = dis(gen);
    return X;
}


using namespace matlab::mex;
using namespace matlab::data;

// TODO put in .h file
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
    for (int i = 0; i < G.N; i++)
    {
        G.E[i] = new int[G.N];
        for (int j = 0; j < G.N; j++)
        {
            G.E[i][j] = (int)_E[i][j];
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
    delete [] rewards;
}


// TODO move to .h file
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
        ~Hierarchy();

        void InitFromMATLAB(StructArray const matlabStructArrayH);
        void InitFromPrior(const Data &D, const Hyperparams &h);

        void PopulateCnt();

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
    for (int i = 0; i < K; i++)
    {
        DEBUG_PRINT("%d ", this->cnt[i]);
    }
    DEBUG_PRINT("]\n");
}

// TODO pass entryIndex and use instead of [0]
void Hierarchy::InitFromMATLAB(StructArray const matlabStructArrayH)
{
    const TypedArray<double> _c = matlabStructArrayH[0]["c"];
    DEBUG_PRINT("H.c = [");
    for (int i = 0; i < _c.getNumberOfElements(); i++)
    {
        this->c[i] = _c[i];
        DEBUG_PRINT("%d ", this->c[i]);
    }
    DEBUG_PRINT("]\n");

    const TypedArray<double> _p = matlabStructArrayH[0]["p"];
    this->p = _p[0];

    const TypedArray<double> _q = matlabStructArrayH[0]["q"];
    this->q = _q[0];

    const TypedArray<double> _tp = matlabStructArrayH[0]["tp"];
    this->tp = _tp[0];

    const TypedArray<double> _hp = matlabStructArrayH[0]["hp"];
    this->hp = _hp[0];

    DEBUG_PRINT("H.p q tp hp = [%lf %lf %lf %lf]\n", this->p, this->q, this->tp, this->hp);

    const TypedArray<double> _theta = matlabStructArrayH[0]["theta"];
    DEBUG_PRINT("H.theta = [");
    this->theta.clear();
    for (int k = 0; k < _theta.getNumberOfElements(); k++)
    {
        this->theta.push_back(_theta[k]);
        DEBUG_PRINT("%lf ", this->theta[k]);
    }
    DEBUG_PRINT("]\n");
   
    const TypedArray<double> _mu = matlabStructArrayH[0]["mu"];
    DEBUG_PRINT("H.mu = [");
    for (int i = 0; i < _mu.getNumberOfElements(); i++)
    {
        this->mu[i] = _mu[i];
        DEBUG_PRINT("%lf ", this->mu[i]);
    }
    DEBUG_PRINT("]\n");

    this->PopulateCnt();
}

// transpiled from init_H.m
// careful with off-by-ones
//
void Hierarchy::InitFromPrior(const Data &D, const Hyperparams &h)
{
    this->c[0] = 1;
    this->cnt.clear();
    this->cnt.push_back(1);
    for (int i = 0; i < D.G.N; i++)
    {
        std::vector<double> p(this->cnt.begin(), this->cnt.end()); // TODO optimize -- no need to copy, could be done in O(1)
        p.push_back(h.alpha);
        int c_new = CatRnd(p);
        if (c_new >= this->cnt.size())
        {
            this->cnt.push_back(1);
        }
        else
        {
            this->cnt[c_new - 1]++;
        }
    }
    // TODO randomize?
    //

    this->p = BetaRnd(1, 1);
    this->q = BetaRnd(1, 1);
    this->tp = BetaRnd(1, 1);
    this->hp = BetaRnd(1, 1);

    int K = this->cnt.size();
    this->theta.clear();
    DEBUG_PRINT("H.theta = [");
    for (int k = 0; k < K; k++)
    {
        this->theta.push_back(NormRnd(h.theta_mean, h.std_theta));
        DEBUG_PRINT("%lf ", this->theta[k]);
    }
    DEBUG_PRINT("]\n");

    DEBUG_PRINT("H.mu = [");
    for (int i = 0; i < D.G.N; i++)
    {
        this->mu[i] = NormRnd(this->theta[this->c[i]], h.std_mu);
        DEBUG_PRINT("%lf ", this->mu[i]);
    }
    DEBUG_PRINT("]\n");
}

Hierarchy::Hierarchy(int _N)
{
    N = _N;
    c = new int[N];
    mu = new double[N];
}


Hierarchy::~Hierarchy()
{
    delete [] c;
    delete [] mu;
}


class MexFunction : public Function {
private:
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;

    // types -- see https://www.mathworks.com/help/matlab/apiref/matlab.data.arraytype.html
    const std::vector<std::string> fieldNamesH = {"c", "p", "q", "tp", "hp", "theta", "mu"}; // H
    const std::vector<ArrayType> fieldTypesH = {ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE};

    const std::vector<std::string> fieldNamesD = {"name", "G", "tasks", "r"}; // D
    const std::vector<ArrayType> fieldTypesD = {ArrayType::CHAR, ArrayType::STRUCT, ArrayType::STRUCT, ArrayType::CELL};

    const std::vector<std::string> fieldNamesG = {"N", "E", "edges"}; // D.G
    const std::vector<ArrayType> fieldTypesG = {ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE};

    const std::vector<std::string> fieldNamesh = {"alpha", "std_theta", "theta_mean", "std_mu", "std_r"}; // h
    const std::vector<ArrayType> fieldTypesh = {ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE, ArrayType::DOUBLE};

public:
  /* Constructor for the class. */
  MexFunction()
  {
    matlabPtr = getEngine();
  }
  
  void displayError(std::string errorMessage)
  {
    ArrayFactory factory;
    matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
            0, std::vector<Array>({
      factory.createScalar(errorMessage) }));
  }
  
  
  /* This is the gateway routine for the MEX-file. */
  void 
  operator()(ArgumentList outputs, ArgumentList inputs) { 
    
    checkArguments (outputs,inputs);

    // check D
    StructArray const matlabStructArrayD = inputs[0];
    checkStructureElements(matlabStructArrayD, "D", fieldNamesD, fieldTypesD);

    size_t total_num_of_elements = matlabStructArrayD.getNumberOfElements();
    for (size_t i=0; i<total_num_of_elements; i++) 
    {
        // check D.G
        const StructArray structFieldG = matlabStructArrayD[i]["G"];
        checkStructureElements(structFieldG, "D.G", fieldNamesG, fieldTypesG);

        const TypedArray<double> _N = structFieldG[0]["N"];
        const TypedArray<double> _E = structFieldG[0]["E"];
        int N = (int)_N[0];
        if (_E.getNumberOfElements() != N * N)
        {
            displayError("D.G.E must have D.G.N^2 elements.");
        }

        // check D.tasks
        const StructArray matlabStructArrayTasks = matlabStructArrayD[i]["tasks"];
        const TypedArray<double> _s = matlabStructArrayTasks[0]["s"];
        const TypedArray<double> _g = matlabStructArrayTasks[0]["g"];
        if (_s.getNumberOfElements() != _g.getNumberOfElements())
        {
            displayError("D.tasks.s and D.tasks.g must have the same number of elements.");
        }

        // check D.r
        const CellArray matlabStructArrayRewards = matlabStructArrayD[i]["r"];
        if (matlabStructArrayRewards.getNumberOfElements() != N)
        {
            displayError("D.r should have D.N elements");
        }
    }

    // init D
    Data D(matlabStructArrayD);

    // check h
    StructArray const matlabStructArrayHyperparams = inputs[1];
    checkStructureElements(matlabStructArrayHyperparams, "h", fieldNamesh, fieldTypesh);

    // init h
    Hyperparams h(matlabStructArrayHyperparams);

    // check nsamples
    int nsamples = 10000; // default
    if (inputs.size() > 2)
    {
        const TypedArray<double> _nsamples = inputs[2];
        nsamples = _nsamples[0];
    }
    DEBUG_PRINT("nsamples = %d\n", nsamples);

    // check burnin
    int burnin = 1; // default
    if (inputs.size() > 3)
    {
        const TypedArray<double> _burnin = inputs[3];
        burnin = _burnin[0];
    }
    DEBUG_PRINT("burnin = %d\n", burnin);

    // check lag
    int lag = 1; // default
    if (inputs.size() > 4)
    {
        const TypedArray<double> _lag = inputs[4];
        lag = _lag[0];
    }
    DEBUG_PRINT("lag = %d\n", lag);

    // check H
    Hierarchy H(D.G.N);
    if (inputs.size() > 5)
    {
        StructArray const matlabStructArrayH = inputs[5];
        checkStructureElements(matlabStructArrayH, "H", fieldNamesH, fieldTypesH);

        const TypedArray<double> _c = matlabStructArrayH[0]["c"];
        if (_c.getNumberOfElements() != D.G.N)
        {
            displayError("H.c should have D.G.N elements");
        }

        const TypedArray<double> _theta = matlabStructArrayH[0]["theta"];
        if (_theta.getNumberOfElements() != D.G.N)
        {
            displayError("H.theta should have D.G.N elements");
        }

        const TypedArray<double> _mu = matlabStructArrayH[0]["mu"];
        if (_mu.getNumberOfElements() != D.G.N)
        {
            displayError("H.mu should have D.G.N elements");
        }

        H.InitFromMATLAB(matlabStructArrayH);
    }
    else
    {
        H.InitFromPrior(D, h);
    }

    // read up on https://www.mathworks.com/help/matlab/apiref/matlab.data.arrayfactory.html?searchHighlight=createarray&s_tid=doc_srchtitle#bvn7dve-1
    ArrayFactory factory;   

    // see https://www.mathworks.com/help/matlab/matlab_external/create-struct-arrays-1.html
    StructArray resultH = factory.createStructArray({ 1,1 }, MexFunction::fieldNamesH ); // dims, fieldNames

    resultH[0]["c"] = factory.createArray<double>({1, (size_t)H.N}, (const double*)H.c, (const double*)(H.c + H.N));
    DEBUG_PRINT("c\n");
    resultH[0]["p"] = factory.createScalar<double>(H.p);
    DEBUG_PRINT("p\n");
    resultH[0]["q"] = factory.createScalar<double>(H.q);
    DEBUG_PRINT("q\n");
    resultH[0]["tp"] = factory.createScalar<double>(H.tp);
    DEBUG_PRINT("tp\n");
    resultH[0]["hp"] = factory.createScalar<double>(H.hp);
    DEBUG_PRINT("hp\n");
    double *theta = new double[H.theta.size()];
    DEBUG_PRINT("pre theta 1\n");
    for (int i = 0; i < H.theta.size(); i++)
    {
        theta[i] = H.theta[i];
        DEBUG_PRINT("pre theta    s %d\n", i);
    }
    DEBUG_PRINT("preee theta");
    resultH[0]["theta"] = factory.createArray<double>({1, H.theta.size()}, (const double*)theta, (const double*)(theta + H.theta.size()));
    DEBUG_PRINT("thetaaa\n");
    delete [] theta;
    DEBUG_PRINT("theta\n");
    resultH[0]["mu"] = factory.createArray<double>({1, (size_t)H.N}, (const double*)H.mu, (const double*)(H.mu + H.N));
    DEBUG_PRINT("mu\n");

    /*
    size_t n = 10; // # of H's
   
    StructArray resultH = factory.createStructArray({ 1,n }, MexFunction::fieldNamesH ); // dims, fieldNames
    for (size_t i = 0; i < n; i++)
    {
        // see https://www.mathworks.com/help/matlab/apiref/matlab.data.arrayfactory.html#bvmdqqr-1
        // and https://www.mathworks.com/help/matlab/apiref/matlab.data.arrayfactory.html
        resultH[i]["c"] = factory.createArray<double>({1, 6}, {1, 2, 3, 4, 5, 6});
        resultH[i]["p"] = factory.createScalar<double>(0.1);
        resultH[i]["q"] = factory.createScalar<double>(0.1);
        resultH[i]["tp"] = factory.createScalar<double>(0.1);
        resultH[i]["hp"] = factory.createScalar<double>(0.1);
        resultH[i]["theta"] = factory.createArray<double>({1, 4}, {0.1, 35.3, 34.5, 12.33});
        resultH[i]["mu"] = factory.createArray<double>({1, 4}, {1000.1, 9935.3, 9934.5, 9912.33});
    }
    */
    
    outputs[0] = resultH;
  }


  // check fields for any structure
  //
  void checkStructureElements(StructArray const & matlabStructArray, const std::string & name, const std::vector<std::string> & expectedFieldNames, const std::vector<ArrayType> & expectedFieldTypes)
  {
      std::ostringstream stream;
      size_t nfields = matlabStructArray.getNumberOfFields();
      auto fields = matlabStructArray.getFieldNames();
      size_t total_num_of_elements = matlabStructArray.getNumberOfElements();
      std::vector<std::string> fieldNames(fields.begin(), fields.end());

      char err[100];

      /* Produce error if structure has wrong number of fields. */
      if(nfields != expectedFieldNames.size())
      {
          sprintf(err, "Struct %s must contain %lu fields.", name.c_str(), expectedFieldNames.size());
          displayError(err);
      }

      for (size_t i = 0; i < expectedFieldNames.size(); i++)
      {
          auto it = find(fieldNames.begin(), fieldNames.end(), expectedFieldNames[i]);

          /* Produce error if field is missing. */
          if (it == fieldNames.end())
          {
              sprintf(err, "Struct %s must contain field '%s'.", name.c_str(), expectedFieldNames[i].c_str());
              displayError(err);
          }
          else
          {
              for (size_t entryIndex=0; entryIndex<total_num_of_elements; entryIndex++) 
              {
                  const Array structField = matlabStructArray[entryIndex][expectedFieldNames[i]];

                  /* Produce error if name field in structure is empty. */
                  if (structField.isEmpty()) 
                  {
                      sprintf(err, "Struct %s has empty field %s on index %zu.", name.c_str(), expectedFieldNames[i].c_str(), entryIndex);
                      displayError(err);
                  }

                  /* Produce error if name is not a valid character array. */
                  if (structField.getType() != expectedFieldTypes[i])
                  {
                      sprintf(err, "Struct %s has field %s on index %zu has invalid type", name.c_str(), expectedFieldNames[i].c_str(), entryIndex);
                      displayError(err);
                  }
              }
          }
      }
  }
  
 
  // check function arguments
  // 
  void checkArguments(ArgumentList outputs, ArgumentList inputs) {
      if (inputs.size() < 2)
      {
        displayError("Specify at least D and h as input arguments.");
      }
      if (inputs.size() > 6)
      {
        displayError("Too many input arguments.");
      }
      if (outputs.size() > 2)
      {
        displayError("Too many outputs specified.");
      }

      if (inputs[0].getType() != ArrayType::STRUCT)
      {
        displayError("D must be a structure.");
      }
      if (inputs[1].getType() != ArrayType::STRUCT)
      {
        displayError("h must be a structure.");
      }
      if (inputs.size() > 2 && inputs[2].getType() != ArrayType::DOUBLE)
      {
        displayError("nsamples must be a number.");
      }
      if (inputs.size() > 3 && inputs[3].getType() != ArrayType::DOUBLE)
      {
        displayError("burnin must be a number.");
      }
      if (inputs.size() > 4 && inputs[4].getType() != ArrayType::DOUBLE)
      {
        displayError("lag must be a number.");
      }
      if (inputs.size() > 5 && inputs[5].getType() != ArrayType::STRUCT)
      {
        displayError("H must be a number.");
      }
  }
};


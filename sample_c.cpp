// compile with:
//
// mex sample_c.cpp printmex.cpp -I/usr/local/boost-1.64.0/include/
// 
// for debugging (mxAssert):
// mex -g sample_c.cpp printmex.cpp -I/usr/local/boost-1.64.0/include/
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

    std::vector<int> z;
    double sum = 0;
    for (int k = 0; k < cnt.size(); k++)
    {
        sum += cnt[k];
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
            //cnt[z[i]] = h.alpha; // <-- this is what sample.m does; but it's wrong; but that works better b/c alpha is too small...
            // this PLUS not deleting from H.cnt and H.theta -> so there's more chance of adding new clusters... 

            cnt[z[i]] = h.alpha / z.size(); // notice all the empty bins have equal probability = alpha, but that's fine b/c it doesn't matter which one we use as the new cluster; we just have to make sure their total probability is not too high, otherwise effective alpha is greater
        }
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
    double logP = log(P[c_i_new]);
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
                H.Update_c_i(c_i_new, i, D, h);
            }
            else
            {
                // reject
                assertThis(H.c[i] == c_i_old);
            }
        }

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

        // update thetas
        //
        for (int k = 0; k < H.theta.size(); k++)
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

        // TODO bridges

        samples.push_back(new Hierarchy(H));
        post.push_back(H.LogPost(D, h)); // TODO optim
    }

    return samples;
}




using namespace matlab::mex;
using namespace matlab::data;


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
    std::srand(0); // for reproducibility
    
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
        // TODO this is wrong; needs to be max(c)
        // also, wtf try to pass Hout as input argument -> Busy
        //if (_theta.getNumberOfElements() != D.G.N)
        //{
        //    displayError("H.theta should have D.G.N elements");
        //}

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

    H.Print();


    std::vector<double> post;
    std::vector<Hierarchy*> samples = sample(D, h, nsamples, burnin, lag, H, post);
    //sample(D, h, nsamples, burnin, lag, H);


    // read up on https://www.mathworks.com/help/matlab/apiref/matlab.data.arrayfactory.html?searchHighlight=createarray&s_tid=doc_srchtitle#bvn7dve-1
    ArrayFactory factory;   

    // see https://www.mathworks.com/help/matlab/matlab_external/create-struct-arrays-1.html

    /*
    StructArray resultH = factory.createStructArray({ 1, 1 }, MexFunction::fieldNamesH ); // dims, fieldNames
    std::vector<double> c(H.c, H.c + H.N);
    resultH[0]["c"] = factory.createArray<std::vector<double>::iterator, double>({1, (size_t)H.N}, c.begin(), c.end());
    resultH[0]["p"] = factory.createScalar<double>(H.p);
    resultH[0]["q"] = factory.createScalar<double>(H.q);
    resultH[0]["tp"] = factory.createScalar<double>(H.tp);
    resultH[0]["hp"] = factory.createScalar<double>(H.hp);
    resultH[0]["theta"] = factory.createArray<std::vector<double>::iterator, double>({1, H.theta.size()}, H.theta.begin(), H.theta.end());
    resultH[0]["mu"] = factory.createArray<double>({1, (size_t)H.N}, (const double*)H.mu, (const double*)(H.mu + H.N));
    */


    // return H = the samples
    //
    StructArray resultH = factory.createStructArray({ 1, samples.size() }, MexFunction::fieldNamesH ); // dims, fieldNames
    for (int i = 0; i < samples.size(); i++)
    {
        Hierarchy *sample = samples[i];

        //resultH[i]["c"] = factory.createArray<int>({1, (size_t)sample->N}, (const int*)sample->c, (const int*)(sample->c + sample->N)); // double is the default in MATLAB; having int here introduces complications...
        std::vector<double> c(sample->c, sample->c + sample->N);
        resultH[i]["c"] = factory.createArray<std::vector<double>::iterator, double>({1, (size_t)sample->N}, c.begin(), c.end());
        resultH[i]["p"] = factory.createScalar<double>(sample->p);
        resultH[i]["q"] = factory.createScalar<double>(sample->q);
        resultH[i]["tp"] = factory.createScalar<double>(sample->tp);
        resultH[i]["hp"] = factory.createScalar<double>(sample->hp);
        resultH[i]["theta"] = factory.createArray<std::vector<double>::iterator, double>({1, sample->theta.size()}, sample->theta.begin(), sample->theta.end());
        resultH[i]["mu"] = factory.createArray<double>({1, (size_t)sample->N}, (const double*)sample->mu, (const double*)(sample->mu + sample->N));

        // IMPORTANT!! free memory; otherwise memory leak => MATLAB hangs
        delete sample;
    }

    outputs[0] = resultH;

    // return post = the log posterior
    //
    assertThis(samples.size() == post.size(), "samples.size() == post.size()");
    Array resultP = factory.createArray<std::vector<double>::iterator, double>({1, post.size()}, post.begin(), post.end());

    outputs[1] = resultP;
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
                      sprintf(err, "Struct %s field %s on index %zu has invalid type", name.c_str(), expectedFieldNames[i].c_str(), entryIndex);
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
        displayError("H must be a structure.");
      }
  }
};


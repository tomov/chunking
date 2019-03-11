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

#include "mcmc.h"
// important -- needs to come after mcmc.h; otherwise boost libraries complain
#include "helpermex.h"


using namespace matlab::mex;
using namespace matlab::data;


class MexFunction : public Function {
private:
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;

public:
  /* Constructor for the class. */
  MexFunction()
  {
    matlabPtr = getEngine();
  }
  
  void displayError(std::string errorMessage)
  {
      ::displayError(errorMessage, matlabPtr);
  }

  /* This is the gateway routine for the MEX-file. */
  void 
  operator()(ArgumentList outputs, ArgumentList inputs) { 
    
    checkArguments (outputs,inputs);

    // D
    StructArray const matlabStructArrayD = inputs[0];
    Data::check(matlabStructArrayD, matlabPtr);
    Data D(matlabStructArrayD);

    // h
    StructArray const matlabStructArrayHyperparams = inputs[1];
    Hyperparams::check(matlabStructArrayHyperparams, matlabPtr);
    Hyperparams h(matlabStructArrayHyperparams);

    // nsamples
    int nsamples = 10000; // default
    if (inputs.size() > 2)
    {
        const TypedArray<double> _nsamples = inputs[2];
        nsamples = _nsamples[0];
    }
    DEBUG_PRINT("nsamples = %d\n", nsamples);

    // burnin
    int burnin = 1; // default
    if (inputs.size() > 3)
    {
        const TypedArray<double> _burnin = inputs[3];
        burnin = _burnin[0];
    }
    DEBUG_PRINT("burnin = %d\n", burnin);

    // lag
    int lag = 1; // default
    if (inputs.size() > 4)
    {
        const TypedArray<double> _lag = inputs[4];
        lag = _lag[0];
    }
    DEBUG_PRINT("lag = %d\n", lag);

    // H
    Hierarchy H(D.G.N);
    if (inputs.size() > 5)
    {
        StructArray const matlabStructArrayH = inputs[5];
        Hierarchy::check(matlabStructArrayH, D, matlabPtr);
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
    StructArray resultH = factory.createStructArray({ 1, samples.size() }, Hierarchy::fieldNames ); // dims, fieldNames
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

        double E[sample->N * sample->N];
        for (int k = 0; k < sample->N; k++)
        {
            for (int l = 0; l < sample->N; l++)
            {
                E[k + sample->N * l] = sample->E[k][l]; // TODO orientation?
            }
        }
        resultH[i]["E"] = factory.createArray<double>({(size_t)sample->N, (size_t)sample->N}, (const double*)E, (const double*)(E + sample->N * sample->N));

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


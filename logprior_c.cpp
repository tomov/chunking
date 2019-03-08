// compile with:
//
// mex logprior_c.cpp printmex.cpp -I/usr/local/boost-1.64.0/include/
// TODO dedupe with loglik_c and sample_c
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
// needs to come last I think
#include "helpermex.h"





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
    std::srand(0); // for reproducibility
    
    checkArguments (outputs,inputs);

    // D
    StructArray const matlabStructArrayD = inputs[1];
    Data::check(matlabStructArrayD, matlabPtr);
    Data D(matlabStructArrayD);

    // h
    StructArray const matlabStructArrayHyperparams = inputs[2];
    Hyperparams::check(matlabStructArrayHyperparams, matlabPtr);
    Hyperparams h(matlabStructArrayHyperparams);


    // H
    Hierarchy H(D.G.N);
    StructArray const matlabStructArrayH = inputs[0];
    Hierarchy::check(matlabStructArrayH, D, matlabPtr);
    H.InitFromMATLAB(matlabStructArrayH);

    H.Print();


    // compute Logprior
    //
    double logp = H.LogPrior(D, h);


    // read up on https://www.mathworks.com/help/matlab/apiref/matlab.data.arrayfactory.html?searchHighlight=createarray&s_tid=doc_srchtitle#bvn7dve-1
    ArrayFactory factory;   

    Array result = factory.createScalar<double>(logp);

    outputs[0] = result;
  }


  // check function arguments
  // 
  void checkArguments(ArgumentList outputs, ArgumentList inputs) {
      if (inputs.size() < 3)
      {
        displayError("Specify H, D and h as input arguments.");
      }
      if (inputs.size() > 3)
      {
        displayError("Too many input arguments.");
      }
      if (outputs.size() > 1)
      {
        displayError("Too many outputs specified.");
      }

      if (inputs[0].getType() != ArrayType::STRUCT)
      {
        displayError("H must be a structure.");
      }
      if (inputs[1].getType() != ArrayType::STRUCT)
      {
        displayError("D must be a structure.");
      }
      if (inputs[2].getType() != ArrayType::STRUCT)
      {
        displayError("h must be a structure.");
      }
  }
};


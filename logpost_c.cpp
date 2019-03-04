// compile with:
//
// mex logpost_c.cpp printmex.cpp -I/usr/local/boost-1.64.0/include/
// TODO dedupe with logprior_c and sample_c and loglik
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
#define DEBUG 0 // <-- uncomment to print debug stuff; must be right after printmex.h

#include "datastructs.h"




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

    // check D TODO dedupe with sample_c.cpp
    StructArray const matlabStructArrayD = inputs[1];
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
    StructArray const matlabStructArrayHyperparams = inputs[2];
    checkStructureElements(matlabStructArrayHyperparams, "h", fieldNamesh, fieldTypesh);

    // init h
    Hyperparams h(matlabStructArrayHyperparams);

    // init & check H TODO dedupe with sample_c.cpp
    Hierarchy H(D.G.N);

    StructArray const matlabStructArrayH = inputs[0];
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

    H.Print();



    // compute Logpost
    //
    double logp = H.LogPost(D, h);


    // read up on https://www.mathworks.com/help/matlab/apiref/matlab.data.arrayfactory.html?searchHighlight=createarray&s_tid=doc_srchtitle#bvn7dve-1
    ArrayFactory factory;   

    Array result = factory.createScalar<double>(logp);

    outputs[0] = result;
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


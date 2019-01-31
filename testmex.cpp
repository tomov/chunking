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
#include<string>
#include<memory>

using namespace matlab::mex;
using namespace matlab::data;

// TODO put in .h file
class Data
{
    public:
        Data(StructArray const &matlabStructArrayD);
        ~Data();

        void initFromTxt(std::string filename);

    private:

        struct Edge
        {
            int u, v;
        };

        struct Task
        {
            int s, g;
        };
        
        struct Graph
        {
            int **E;
            int N;
            std::vector<Edge> edges;
        };

        std::string name;
        Graph G;
        std::vector<Task> tasks;
        std::vector<double> *rewards;
};


Data::Data(StructArray const &matlabStructArrayD)
{
    const TypedArray<double> _N = matlabStructArrayD[0]["N"];
    G.N = (int)_N[0];

	printThis("sheeeeiiit\n");

    G.E = new int*[G.N];
    for (int i = 0; i < G.N; i++)
    {
        G.E[i] = new int[G.N];
    }
}

Data::~Data()
{
    for (int i = 0; i < G.N; i++)
    {
        delete [] G.E[i];
    }
    delete [] G.E;
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
  
  /* Helper function to print output string on MATLAB command prompt. */
  void displayOnMATLAB(std::ostringstream stream)
  {
    ArrayFactory factory;
    matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("fprintf"),0, std::vector<Array>
            ({ factory.createScalar(stream.str())}));
  }
  
  /* Helper function to generate an error message from given string,
   * and display it over MATLAB command prompt.
   */
  void displayError(std::string errorMessage)
  {
    ArrayFactory factory;
    matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
            0, std::vector<Array>({
      factory.createScalar(errorMessage) }));
  }
  
  /* Helper function to information about an empty field in the structure. */
  void emptyFieldInformation(std::string fieldName, size_t index)
  {
    std::ostringstream stream;
    stream<<"Field: "<<std::string(fieldName)<<" of the element at index: "
          <<index+1<<" is empty."<<std::endl;
    displayOnMATLAB(std::move(stream));
  }
  
  /* Helper function to information about an invalid field in the structure. */
  void invalidFieldInformation(std::string fieldName, size_t index)
  {
    std::ostringstream stream;
    stream<<"Field: "<<std::string(fieldName)<<" of the element at index: "
          <<index+1<<" contains wrong value."<<std::endl;
    displayOnMATLAB(std::move(stream));
  }
  
  
  /* This is the gateway routine for the MEX-file. */
  void 
  operator()(ArgumentList outputs, ArgumentList inputs) { 
    
    checkArguments (outputs,inputs);

    // check D
    StructArray const matlabStructArrayD = inputs[0];
    checkStructureElements(matlabStructArrayD, "D", fieldNamesD, fieldTypesD);

    // check D.G
    size_t total_num_of_elements = matlabStructArrayD.getNumberOfElements();
    for (size_t i=0; i<total_num_of_elements; i++) 
    {
        const StructArray structFieldG = matlabStructArrayD[i]["G"];
        checkStructureElements(structFieldG, "D.G", fieldNamesG, fieldTypesG);
    }

    printThis("slkdjflks %d %s\n", 34, "sdf");

    // check h
    StructArray const matlabStructArrayh = inputs[1];
    checkStructureElements(matlabStructArrayh, "h", fieldNamesh, fieldTypesh);

    // check nsamples

    ArrayFactory factory;   
    /*

    checkStructureElements(matlabStructArray);
    auto fields = matlabStructArray.getFieldNames();
    size_t total_num_of_elements = matlabStructArray.getNumberOfElements();   
    CellArray phoneNumberStringArray = 
            factory.createCellArray({ 1,total_num_of_elements });
    CellArray nameStringArray = 
            factory.createCellArray({ 1,total_num_of_elements });
    std::vector<std::string> fieldNames(fields.begin(), fields.end());
    */

    size_t n = 10; // # of H's
   
    // see https://www.mathworks.com/help/matlab/matlab_external/create-struct-arrays-1.html
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
  
  /* Make sure that the passed structure has valid data. */
  // TODO rm 
  void checkStructureElements(StructArray const & matlabStructArray)
  {
    std::ostringstream stream;
    size_t nfields = matlabStructArray.getNumberOfFields();
    auto fields = matlabStructArray.getFieldNames();
    size_t total_num_of_elements = matlabStructArray.getNumberOfElements();
    std::vector<std::string> fieldNames(fields.begin(), fields.end());
    
    /* Produce error if structure has more than 2 fields. */
    if(nfields != 2) {
      displayError("Struct must consist of 2 entries."
                   "(First: char array, Second: numeric double scalar).");
    }
    
    /* Walk through each structure element. */
    for (size_t entryIndex=0; entryIndex<total_num_of_elements; entryIndex++) {
      const Array structField1 = 
              matlabStructArray[entryIndex][fieldNames[0]];
      const Array structField2 = 
              matlabStructArray[entryIndex][fieldNames[1]];
      
      /* Produce error if name field in structure is empty. */
      if (structField1.isEmpty()) {
        emptyFieldInformation(fieldNames[0],entryIndex);
        displayError("Empty fields are not allowed in this program." 
                     "This field must contain character array.");
      }
      
      /* Produce error if phone number field in structure is empty. */
      if(structField2.isEmpty()) {
        emptyFieldInformation(fieldNames[1],entryIndex);
        displayError("Empty fields are not allowed in this program." 
                     "This field must contain numeric double scalar.");
      }
      
      /* Produce error if name is not a valid character array. */
      if(structField1.getType()!= ArrayType::CHAR) {
        invalidFieldInformation(fieldNames[0],entryIndex);
        displayError("This field must contain character array.");
      }
      /* Produce error if phone number is not a valid double scalar. */
      if (structField2.getType() != ArrayType::DOUBLE
          || structField2.getNumberOfElements() != 1) {
        invalidFieldInformation(fieldNames[1],entryIndex);
        displayError("This field must contain numeric double scalar.");
      }
    }
  }
  
  /* This function makes sure that user has provided structure as input,
   * and is not expecting more than one output in results.
   */
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


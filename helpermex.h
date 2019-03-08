#ifndef HELPER_MEX_H
#define HELPER_MEX_H

#include "mex.hpp"

using namespace matlab::mex;
using namespace matlab::data;


void displayError(std::string errorMessage, std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr)
{
  ArrayFactory factory;
  matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
          0, std::vector<Array>({
    factory.createScalar(errorMessage) }));
}
  
// check fields for any structure
//
void checkStructureElements(StructArray const & matlabStructArray, const std::string & name, const std::vector<std::string> & expectedFieldNames, const std::vector<ArrayType> & expectedFieldTypes, std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr)
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
        displayError(err, matlabPtr);
    }

    for (size_t i = 0; i < expectedFieldNames.size(); i++)
    {
        auto it = find(fieldNames.begin(), fieldNames.end(), expectedFieldNames[i]);

        /* Produce error if field is missing. */
        if (it == fieldNames.end())
        {
            sprintf(err, "Struct %s must contain field '%s'.", name.c_str(), expectedFieldNames[i].c_str());
            displayError(err, matlabPtr);
        }
        else
        {
            for (size_t entryIndex=0; entryIndex<total_num_of_elements; entryIndex++) 
            {
                const Array structField = matlabStructArray[entryIndex][expectedFieldNames[i]];

                /* Produce error if name field in structure is empty. */
                // jk no need to; sometimes we expect it; also we check the counts in the check() functions in datastructs.h
                /*
                if (structField.isEmpty()) 
                {
                    sprintf(err, "Struct %s has empty field %s on index %zu.", name.c_str(), expectedFieldNames[i].c_str(), entryIndex);
                    displayError(err, matlabPtr);
                }
                */

                /* Produce error if name is not a valid character array. */
                if (structField.getType() != expectedFieldTypes[i])
                {
                    sprintf(err, "Struct %s field %s on index %zu has invalid type", name.c_str(), expectedFieldNames[i].c_str(), entryIndex);
                    displayError(err, matlabPtr);
                }
            }
        }
    }
}
  
 

#endif

#include "printmex.h"
#include "mex.h"
#include "stdarg.h"

// helper thingy to printf to matlab; cannot do it directly in cpp file b/c including mex.h (the C version) and mex.hpp (the C++ version) is incompatible. We need former for mexPrintf, latter for C++ stuff
//

void printThis( const char* format, ... ) {
    va_list args;
    va_start( args, format );
    mexPrintf( format, args );
    va_end( args );
}

#ifndef PRINT_MEX_H
#define PRINT_MEX_H

// TODO figure out where to define DEBUG 
#define DEBUG 0
#define DEBUG_PRINT(args ...) if (DEBUG) printThis(args)
#define ASSERT(args ...) if (DEBUG) assertThis(args)

void printThis(const char* format, ...);
void assertThis(int expr, const char* msg = "");

#endif

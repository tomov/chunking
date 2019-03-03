#ifndef PRINT_MEX_H
#define PRINT_MEX_H

#define DEBUG 1
#define DEBUG_PRINT(args ...) if (DEBUG) printThis(args)

void printThis(const char* format, ...);
void assertThis(int expr, const char* msg = "");

#endif

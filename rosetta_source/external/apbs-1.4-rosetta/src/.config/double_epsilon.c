#include <stdio.h>

int main()
{
    double doubleEpsilon = 1.0;
    while(1.0 + doubleEpsilon / 2.0 != 1.0)
        doubleEpsilon /= 2.0;
    printf( "%e", doubleEpsilon );
    return 0;
}

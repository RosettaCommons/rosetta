#include <stdio.h>

int main()
{
    float floatEpsilon = 1.0;
    while(1.0 + floatEpsilon / 2.0 != 1.0)
        floatEpsilon /= 2.0;
    printf( "%e", floatEpsilon );
    return 0;
}

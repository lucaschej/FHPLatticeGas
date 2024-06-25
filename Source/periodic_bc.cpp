// periodic_bc.cpp
#include "periodic_bc.h"

//Computes periodic boundary conditions.
int periodic_bc(int k, int maxk)
{
    if (k < 0)
    {
        return maxk-1;
    }
    else if(k == maxk)
    {
        return 0;
    }
    else
    {
        return k;
    }
}

//Periodic conditions for bit-by-bit checking
int periodic_bc(int b, int maxb, int& x)
{
    if (b < 0)
    {
        x = periodic_bc(x-1, XMAX); //Goes backward on x
        return maxb-1; //Get last bit
    }
    else if(b == maxb)
    {
        x = periodic_bc(x+1, XMAX); //Goes forward on x
        return 0; //Get first bit
    }
    else
    {
        return b;
    }
}
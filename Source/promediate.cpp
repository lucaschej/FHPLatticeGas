// promediate.cpp
#include "promediate.h"

//Promediate over measured values in the main loop to get stable data
void promediate(double mean_data[NMEANDATA], vector<double> ni_to_write[7], vector<double> u_to_write[2], double& rho_eq, double u_eq[2], double Ni_eq[7], int fraction_its, int n_its)
{
    int i,j;

    //Promediate in the equilibrium
    if (i >= fraction_its * n_its)
    {
        for (j=0; j < 7; j++)
        {
            Ni_eq[j] += mean_data[j];
        }
        rho_eq += mean_data[7];
        u_eq[0] += mean_data[8];
        u_eq[1] += mean_data[9];
    }

    //Division by the number of particles will be done after the main loop, to save time

    //Add the ni and u measured values to do a graph
    for (j=0; j < 7; j++)
    {
        ni_to_write[j].push_back(mean_data[j]);
    }
    u_to_write[0].push_back(mean_data[8]);
    u_to_write[1].push_back(mean_data[9]);

    return;
}
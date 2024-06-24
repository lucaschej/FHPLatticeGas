// get_results.cpp
#include "get_results.h"

//Writes the values on the variables given by reference.
void get_results_FHPI(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re)
{
    d = rho / 6.0;
    cs = 1.0/ sqrt(2.0);
    gd = ((1.0 - 2.0 * d) / (1.0 - d)) / 2.0;
    nud = (1.0 / (d *  pow(1.0 - d, 3.0))) / 12.0 - 1.0/8.0;
    etad = 0.0;
    re = cs * gd / nud;

    return;
}

void get_results_FHPII(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re)
{
    d = rho / 7.0;
    cs = sqrt(3.0 / 7.0);
    gd = 7.0*(1.0-2.0*d)/(12.0 * (1.0 - d));
    nud = 1.0 / (28.0 * d * pow(1.0 - d, 3.0) * (1.0 - 4.0 * d / 7.0)) - 1.0/8.0;
    etad = 1.0 / (98.0 * d * pow(1.0 - d, 4.0)) - 1.0/28.0;
    re = cs * gd / nud;

    return;
}

void get_results_FHPIII(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re)
{
    d = rho / 7.0;
    cs = sqrt(3.0 / 7.0);
    gd = 7.0*(1.0-2.0*d)/(12.0 * (1.0 - d));
    nud = 1.0 / (28.0 * d * (1.0 - d) * (1.0 - 8.0 * d * (1.0 - d)/ 7.0)) - 1.0/8.0;
    etad = 1.0 / (98.0 * d * (1.0 - d) * (1.0 - 2.0 * d * (1.0 - d))) - 1.0/28.0;
    re = cs * gd / nud;

    return;
}
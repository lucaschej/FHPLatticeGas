#include "config.h"
#include "initialize.h"
#include "collisionFHP_I.h"
#include "collisionFHP_II.h"
#include "collisionFHP_III.h"
#include "propagation.h"
#include "measure.h"

#include<iostream>
#include<cstdlib>
#include<random>
#include<ctime>
#include<fstream>
#include<vector>
#include<string>
#include<chrono>

using namespace std;

//Inits the variables
int initialize(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> >& nsbit, mt19937 generator,  double prob);

//Get the results for g(d), nu(d), etc for the models
void get_results_FHPI(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re);
void get_results_FHPII(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re);
void get_results_FHPIII(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re);

//Get the k-th bit of n
int bit_at(uint64_t n, int k);

//Write all the data
void writeData(vector< vector<uint64_t> > cell[7], string filename);
void writeGrid();
void writeData(vector<double> data_to_write[7], int N, string filename);
void writeResults(double rho, string filename);

//Measures the most important values of the distribution
//void measure(vector< vector<uint64_t> > cell[7], double ci[7][2], vector< vector<double> > Ni[7], vector< vector<double> >&  rho, vector< vector<double> >  u[2], double mean_data[NMEANDATA]);
void promediate(double mean_data[NMEANDATA], vector<double> ni_to_write[7], vector<double> u_to_write[2], double& rho_eq, double u_eq[2], double Ni_eq[7], int fraction_its, int n_its);

//Gives the equilibrium values for the mean occupation numbers
void equilibrium_ni(double rho, double ci[7][2], double u[2], double ni[7], double nieq[7]);


int main(void)
{
    int i,j;
    double k;
    //Six directions + rest particles
    //  a   b     |   0   1
    //   \ /      |    \ /
    //f --g-- c   |  5--6--2
    //   / \      |    / \
    //  e   d     |   4   3
    vector< vector<uint64_t> > cell[7];
    vector< vector<uint64_t> > result_cell[7];
    vector< vector<uint64_t> > nsbit (XMAX, vector<uint64_t>(YMAX));

    //Velocities for the grid above
    double ci[7][2] =
    {
      {-1.0/2.0, sqrt(3.0)/2.0},
      {1.0/2.0, sqrt(3.0)/2.0},
      {1.0, 0.0},
      {1.0/2.0, -sqrt(3.0)/2.0},
      {-1.0/2.0, -sqrt(3.0)/2.0},
      {-1.0, 0.0}
    };

    //Variables to get the numerical data
    double N_particles; //Number of particles

    vector< vector<double> > Ni[7]; //Mean occupation numbers
    double Ni_the[7]; //Theoretical mean occupation numbers in equilibrium

    vector< vector<double> > rho (XMAX*64, vector<double>(YMAX)); //Density
    //vector< vector<double> > rho (YMAX, vector<double>(YMAX)); //Density
    vector< vector<double> > u[2]; //Velocity
    double mean_data[NMEANDATA]; //Mean of Ni (0-6), rho (7), u (8-9) in all the space
    vector<double> ni_to_write[7]; //Used to write the Ni in every iteration
    vector<double> u_to_write[2]; //Used to write u in every iteration

    double rho_eq; //Values promediated over all the iterations...
    double u_eq[2];
    double Ni_eq[7];

    //Vector init for all the variables
    for (i=0; i < 7; i++)
    {
        cell[i] = vector< vector<uint64_t> >(XMAX, vector<uint64_t>(YMAX));
        result_cell[i] = vector< vector<uint64_t> >(XMAX, vector<uint64_t>(YMAX));
        Ni[i] = vector< vector<double> >(XMAX*64, vector<double>(YMAX));
        //Ni[i] = vector< vector<double> >(YMAX, vector<double>(YMAX));
    }
    for (i=0; i < 2; i++)
    {
        u[i] = vector< vector<double> >(XMAX*64, vector<double>(YMAX));
        //u[i] = vector< vector<double> >(XMAX, vector<double>(XMAX));
    }


    //Init of the random generator and distribution
    mt19937 gen (chrono::high_resolution_clock::now().time_since_epoch().count());

    k = 0.2;

    //for (k=0.01; k < 0.5; k += 0.02)
    //{
        //Init the grid and get all the particles
        N_particles = initialize(cell, result_cell, nsbit, gen, k);

        //Init values for promediating for every iteration
        rho_eq = 0.0;
        u_eq[0] = 0.0;
        u_eq[1] = 0.0;
        for (i=0; i <7;i++)
        {
            Ni_eq[i] = 0.0;
        }

        //Do iterations
        for (i=0; i <= n_its; i++)
        {
            collisionFHP_I(cell, result_cell, nsbit, gen); //Collision
            propagation(result_cell,cell,nsbit); //Propagation

            measure(cell, ci, Ni, rho, u, mean_data); //Measure
            promediate(mean_data, ni_to_write, u_to_write, rho_eq, u_eq, Ni_eq, fraction_its, n_its); //Promediate values to obtain _eq ones

            //To do animations
            writeData(cell,"data"+to_string(i)+".txt");

        }

        //Write ni in every iteration
        //writeData(ni_to_write, 7, "ni");

        //Write u in every iteration
        //writeData(u_to_write, 2, "speed");


        //Finish the mean values
        for (j=0; j < 7; j++)
        {
            Ni_eq[j] /= (1.0-fraction_its) * n_its;
        }
        rho_eq /= (1.0-fraction_its) * n_its;
        u_eq[0] /= (1.0-fraction_its) * n_its;
        u_eq[1] /= (1.0-fraction_its) * n_its;


        cout << rho_eq << endl;

        //Calculate equilibrium values for the numerical data obtained
        equilibrium_ni(rho_eq, ci, u_eq, Ni_eq, Ni_the);

        //Append the results to a file.
        //writeResults(rho_eq, "results3.txt");
    //}

    return 0;
}

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

void writeGrid()
{
    ofstream output;

    double L = 1.0; //side of the grid
    int i,j,k; //counters
    double xpos; //Position to write.

    output.open("grid.txt");
    for (i=0; i < XMAX; i++)
    {
        for (j=0; j < YMAX; j++)
        {
            for (k=0; k < 64; k++)
            {
                xpos = k + 64 * i; //Get every position inside the bit array
                output << 0.5 * j * L + xpos * L << " " << j * L * 0.8660254038 << endl; //Write the output
            }
        }
    }
    output.close();

    return;
}

void writeData(vector< vector<uint64_t> > cell[7], string filename)
{
    ofstream output;
    uint64_t particle; //Is there a particle in?
    double L = 1.0; //side of the grid
    int i,j,k; //counters
    int b; //Every bit of particle
    double xpos; //Position to write.

    output.open(filename);
    for (i=0; i < XMAX; i++)
    {
        for (j=0; j < YMAX; j++)
        {
            //Particle = 1 if there's at least one particle
            particle = cell[0][i][j]|cell[1][i][j]|cell[2][i][j]|cell[3][i][j]|cell[4][i][j]|cell[5][i][j]|cell[6][i][j];

            //for (k=0; k < 64; k++)
            for (k=0; k < 1; k++)
            {
                b = bit_at(particle, k);
                if (b == 1)
                {
                    xpos = k + 1 * i; //Get every position inside the bit array
                    //xpos = k + 64 * i; //Get every position inside the bit array
                    output << 0.5 * j * L + xpos * L << " " << j * L * 0.8660254038 << endl; //Write the output
                }
            }
        }
    }
    output.close();


    return;
}

//Write the first N data stored in the vector. With N=7, Ni can be written, and
//with N=2, also u
void writeData(vector<double> data_to_write[7], int N, string filename)
{
    ofstream output;
    int i,j;


    for (i=0; i < N; i++)
    {
        output.open(filename + to_string(i) + ".txt");
        for (j=0; j < data_to_write[i].size(); j++)
        {
            output << j << " " << data_to_write[i][j] << endl;
        }
        output.close();
    }


    return;
}

void writeResults(double rho, string filename)
{
    double d, cs, gd, nud, etad, re;
    ofstream file;

    get_results_FHPIII(rho, d, cs, gd, nud, etad, re);

    file.open(filename, ios_base::app);
    file << d << " " << cs << " " << gd << " " << nud << " " << etad << " " << re << endl;
    file.close();

    return;
}

//Returns the k-th bit of n. What it does is to put the k-th as first,
//then check the k-th & 1 (because 1<<63 is 1000...000) and we can obtain then
//0000...000 or 1000...0000 so we do again >>63 to get 1 or 0.
int bit_at(uint64_t n, int k)
{
    return (((n<<k) & ((uint64_t(1)<<63))) >>63);
}

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

void equilibrium_ni(double rho, double ci[7][2], double u[2], double ni[7], double nieq[7])
{
    int i;
    double Grho = (6.0-2.0*rho)/(3.0*(6-rho));
    double cdotu; //Scalar product ci * u
    for (i=0; i < 7; i++)
    {
        cdotu = ci[i][0]*u[0]+ci[i][1]*u[1];
        nieq[i] = rho/6.0 + rho*(cdotu)/3.0 + rho*Grho*((cdotu*cdotu) - (u[0]*u[0]+u[1]*u[1])/2.0);
    }
    return;
}

// include.cpp
#include "measure.h"
#include "periodic_bc.h"

//Get the k-th bit of n
int bit_at(uint64_t n, int k);

//Measures the mean occupation numbers, density and velocity
void measure(vector< vector<uint64_t> > cell[7], double ci[7][2], vector< vector<double> > Ni[7], vector< vector<double> >&  rho, vector< vector<double> >  u[2], double mean_data[NMEANDATA])
{
    int i,j,k,b;
    int bnext, bprev, ynext, yprev;
    int neiga, neigb, neigc, neigd, neige, neigf; //Following the directions of the FHP grid
    double sum;
    vector< vector<double> > flux[2];

    for (i=0; i < 2; i++)
    {
        flux[i] = vector< vector<double> >(XMAX*64, vector<double>(YMAX));
        //flux[i] = vector< vector<double> >(YMAX, vector<double>(YMAX));
    }

    //Init sums for the mean data
    for (i=0; i < NMEANDATA; i++)
    {
            mean_data[i] = 0.0;
    }
    //Make measures for every (x,y) in space
    for (i=0; i < XMAX; i++)
    {
        for (j=0; j < YMAX; j++)
        {

            for (b=0; b < 64; b++)
            {
                rho[b + i*64][j] = 0.0; //Initalize the rho in this  (x,y)
                flux[0][b + i*64][j] = 0.0;
                flux[1][b + i*64][j] = 0.0; //The same for flux vector

                bnext = periodic_bc(b+1, 64, i);
                bprev = periodic_bc(b-1, 64, i);
                ynext = periodic_bc(j+1, YMAX);
                yprev = periodic_bc(j-1, YMAX);

                //Calculate Ni for every direction:
                for (k=0; k < 7; k++)
                {
                    neiga = bit_at(cell[k][i][ynext],bprev);
                    neigb = bit_at(cell[k][i][ynext],b);
                    neigc = bit_at(cell[k][i][j],bnext);
                    neigd = bit_at(cell[k][i][yprev],bnext);
                    neige = bit_at(cell[k][i][yprev],b);
                    neigf = bit_at(cell[k][i][j],bprev);

                    Ni[k][b + i*64][j] = (1.0*(neiga+neigb+neigc+neigd+neige+neigf))/6.0; //Obtain the Ni at x,y
                    rho[b + i*64][j]  += Ni[k][b + i*64][j]; //Add the result to rho

                    mean_data[k] += Ni[k][b + i*64][j]; //Computes sum of all Ni for the entire space

                    flux[0][b + i*64][j] += ci[k][0] * Ni[k][b + i*64][j];
                    flux[1][b + i*64][j] += ci[k][1] * Ni[k][b + i*64][j]; //The same for the vector flux

                }
                //Get the velocities from the data calculated, if density not equals 0.0
                u[0][b + i*64][j] = rho[b + i*64][j] != 0 ? flux[0][b + i*64][j] / rho[b + i*64][j] : 0.0;
                u[1][b + i*64][j] = rho[b + i*64][j] != 0 ? flux[1][b + i*64][j] / rho[b + i*64][j] : 0.0;



                //Sum for the entire space for rho and u
                mean_data[7] += rho[b + i*64][j];
                mean_data[8] += u[0][b + i*64][j];
                mean_data[9] += u[1][b + i*64][j];


            }


        }
    }

    //Finish the mean dividing by the number of cells
    for (i=0; i < NMEANDATA; i++)
    {
            //mean_data[i] /= (XMAX*64)*YMAX;
            mean_data[i] /= (YMAX)*YMAX;
    }

    return;
}
// propagation.cpp
#include "propagation.h"
#include "periodic_bc.h"

//Propagation of particles. Result is exported again to cell from the result_cell
//computed in collisions.
void propagation(vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > cell[7], vector< vector<uint64_t> > nsbit)
{
    int x,y, xnext, xprev, ynext, yprev;

    for (x=0; x < XMAX; x++)
    {
        for (y=0; y < YMAX; y++)
        {
            //Use periodic boundary conditions
            xnext = periodic_bc(x+1, XMAX);
            ynext = periodic_bc(y+1, YMAX);
            xprev = periodic_bc(x-1, XMAX);
            yprev = periodic_bc(y-1, YMAX);

            //Propagator:
            cell[0][x][ynext] = (result_cell[0][x][y]<<1)^(result_cell[0][xnext][y]>>63);
            cell[1][x][ynext] = result_cell[1][x][y];
            cell[2][xnext][y] = (result_cell[2][x][y]<<63)^(result_cell[2][xnext][y]>>1);
            cell[3][x][yprev] = (result_cell[3][x][y]>>1)^(result_cell[3][xprev][y]<<63);
            cell[4][x][yprev] = result_cell[4][x][y];
            cell[5][xprev][y] = (result_cell[5][x][y]>>63)^(result_cell[5][xprev][y]<<1);
            cell[6][x][y] = result_cell[6][x][y]; //Rest particles don't move

        }
    }
}

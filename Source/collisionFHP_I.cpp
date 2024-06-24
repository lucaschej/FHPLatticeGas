// collisionFHP_I.cpp
#include "collisionFHP_I.h"

void collisionFHP_I(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > nsbit, mt19937 &generator)
{
    int x,y;

    uint64_t rnd, no_rnd, nsb; //Random, negate of random, and non-solid bit

    uint64_t a,b,c,d,e,f; //Single cells:
    uint64_t tb_ad, tb_be, tb_cf; //Two-body possible collisions
    uint64_t triple; //Three body collision
    uint64_t cha, chb, chc, chd, che, chf; //Change in each cell

    uniform_int_distribution<uint64_t> dist(uint64_t(0), ~uint64_t(0)); //from 000...0000 to 111...1111

    for (x=0; x < XMAX; x++)
    {
        for (y=0; y < YMAX; y++)
        {

            //Shorter names
            a = cell[0][x][y];
            b = cell[1][x][y];
            c = cell[2][x][y];
            d = cell[3][x][y];
            e = cell[4][x][y];
            f = cell[5][x][y];

            //If there're particles in i and i+3 but not in any other cell, then
            //we have a two body collision
            tb_ad = (a&d)&(~(b|c|e|f));
            tb_be = (b&e)&(~(a|c|d|f));
            tb_cf = (c&f)&(~(a|b|d|e));


            //If we don't have any contiguous occupied or non-occupied sites, then
            //we have a three body collision
            triple = (a^b)&(b^c)&(c^d)&(d^e)&(e^f);

            //Get the random bit and wall bit
            rnd = dist(generator);
            no_rnd = ~rnd; //Save time!
            nsb = nsbit[x][y];


            //Change the configuration using the collisions.
            cha = (tb_ad|triple|(rnd&tb_be)|(no_rnd&tb_cf)&nsb);
            chd = (tb_ad|triple|(rnd&tb_be)|(no_rnd&tb_cf)&nsb);
            chb = (tb_be|triple|(rnd&tb_cf)|(no_rnd&tb_ad)&nsb);
            che = (tb_be|triple|(rnd&tb_cf)|(no_rnd&tb_ad)&nsb);
            chc = (tb_cf|triple|(rnd&tb_ad)|(no_rnd&tb_be)&nsb);
            chf = (tb_cf|triple|(rnd&tb_ad)|(no_rnd&tb_be)&nsb);

            //Update the configuration with the change
            result_cell[0][x][y] = a ^ cha;
            result_cell[1][x][y] = b ^ chb;
            result_cell[2][x][y] = c ^ chc;
            result_cell[3][x][y] = d ^ chd;
            result_cell[4][x][y] = e ^ che;
            result_cell[5][x][y] = f ^ chf;
        }
    }

    return;
}

// collisionFHP_III.cpp
#include "collisionFHP_III.h"

void collisionFHP_III(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > nsbit, mt19937 generator)
{
    int x,y;

    uint64_t rnd, no_rnd, nsb; //Random, negate of random, and non-solid bit

    uint64_t a,b,c,d,e,f,r; //Single cells:
    uint64_t ad, be, cf; //Auxiliar cells;
    uint64_t tb_ad, tb_be, tb_cf, tb_ra, tb_rb, tb_rc, tb_rd, tb_re, tb_rf; //Two-body possible collisions
    uint64_t fb_adbe, fb_adcf, fb_becf; //Possible four-body collisions
    uint64_t ra, rb, rc, rd, re, rf; //Collision of particles giving a rest one
    uint64_t spect_ad_b, spect_ad_c, spect_ad_e, spect_ad_f;
    uint64_t spect_be_a, spect_be_c, spect_be_d, spect_be_f;
    uint64_t spect_cf_a, spect_cf_b, spect_cf_d, spect_cf_e; //Two body head-on collision with spectator
    uint64_t spectator_change_ad, spectator_change_be, spectator_change_cf; //Changes due to spectators
    uint64_t fb_change_ad, fb_change_be, fb_change_cf; //Changes due to 4-body collisions
    uint64_t triple; //Three body collision
    uint64_t cha, chb, chc, chd, che, chf, chr; //Change in each cell

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
            r = cell[6][x][y];

            ad = a&d;
            be = b&e;
            cf = c&f;

            //Get the random bit and wall bit
            rnd = dist(generator);
            no_rnd = ~rnd; //Save time!
            nsb = nsbit[x][y];

            //If there're particles in i and i+3 but not in any other cell, then
            //we have a two body collision
            tb_ad = ad&(~(b|c|e|f));
            tb_be = be&(~(a|c|d|f));
            tb_cf = cf&(~(a|b|d|e));

            //With rest:
            tb_ra = (r&a&~(b|c|d|e|f));
            tb_rb = (r&b&~(a|c|d|e|f));
            tb_rc = (r&c&~(a|b|d|e|f));
            tb_rd = (r&d&~(a|b|c|e|f));
            tb_re = (r&e&~(a|b|c|d|f));
            tb_rf = (r&f&~(a|b|c|d|e));

            //Handle the four body collisions
            fb_adbe = ad&be&(~(c|f));
            fb_adcf = ad&cf&(~(b|e));
            fb_becf = be&cf&(~(a|d));

            //Changes to do due to this collisions. When the ad is present, it changes with p=0.5. In becf it always appears:
            fb_change_ad = (rnd&fb_adbe)|(no_rnd&fb_adcf)|fb_becf;
            fb_change_be = (rnd&fb_becf)|(no_rnd&fb_adbe)|fb_adcf;
            fb_change_cf = (rnd&fb_adcf)|(no_rnd&fb_becf)|fb_adbe;

            //Three body collisions with spectator. The array-form is for legibility. Each two body collision can have 4-subcases.
            spect_ad_b = ad&b&(~(c|e|f));
            spect_ad_c = ad&c&(~(b|e|f));
            spect_ad_e = ad&e&(~(b|c|f));
            spect_ad_f = ad&f&(~(b|c|e));
            spect_be_a = be&a&(~(c|d|f));
            spect_be_c = be&c&(~(a|d|f));
            spect_be_d = be&d&(~(a|c|f));
            spect_be_f = be&f&(~(a|c|d));
            spect_cf_a = cf&a&(~(b|d|e));
            spect_cf_b = cf&b&(~(a|d|e));
            spect_cf_d = cf&d&(~(a|b|e));
            spect_cf_e = cf&e&(~(a|b|d));

            //Changes to do due to this collisions. The "a" will change if we have a head-on collision in ad, or if we have a head on collision
            //in be or cf which leads to a rotation to ad
            spectator_change_ad = (spect_ad_b|spect_ad_e|spect_ad_c|spect_ad_f)|(spect_be_c|spect_be_f)|(spect_cf_b|spect_cf_e);
            spectator_change_be = (spect_be_a|spect_be_c|spect_be_d|spect_be_f)|(spect_ad_c|spect_ad_f)|(spect_cf_a|spect_cf_d);
            spectator_change_cf = (spect_cf_a|spect_cf_b|spect_cf_d|spect_cf_e)|(spect_ad_b|spect_ad_e)|(spect_be_a|spect_be_d);

            //If we don't have any contiguous occupied or non-occupied sites, then
            //we have a three body collision
            triple = (a^b)&(b^c)&(c^d)&(d^e)&(e^f);

            //Collision of a two particles at i and i+2 which leads to a rest and a moving one:
            ra = (f&b&~(r|a|c|d|e));
            rb = (a&c&~(r|b|d|e|f));
            rc = (b&d&~(r|a|c|e|f));
            rd = (c&e&~(r|a|b|d|f));
            re = (d&f&~(r|a|b|c|e));
            rf = (e&a&~(r|b|c|d|f));


            //Change the configuration using the collisions.
            cha = (tb_ad|triple|(rnd&tb_be)|(no_rnd&tb_cf)|tb_ra|tb_rb|tb_rf|ra|rb|rf|spectator_change_ad|fb_change_ad&nsb);
            chd = (tb_ad|triple|(rnd&tb_be)|(no_rnd&tb_cf)|tb_rd|tb_rc|tb_re|rd|rc|re|spectator_change_ad|fb_change_ad&nsb);
            chb = (tb_be|triple|(rnd&tb_cf)|(no_rnd&tb_ad)|tb_rb|tb_ra|tb_rc|rb|ra|rc|spectator_change_be|fb_change_be&nsb);
            che = (tb_be|triple|(rnd&tb_cf)|(no_rnd&tb_ad)|tb_re|tb_rd|tb_rf|re|rd|rf|spectator_change_be|fb_change_be&nsb);
            chc = (tb_cf|triple|(rnd&tb_ad)|(no_rnd&tb_be)|tb_rc|tb_rb|tb_rd|rc|rb|rd|spectator_change_cf|fb_change_cf&nsb);
            chf = (tb_cf|triple|(rnd&tb_ad)|(no_rnd&tb_be)|tb_rf|tb_ra|tb_re|rf|ra|re|spectator_change_cf|fb_change_cf&nsb);
            chr = (tb_ra|tb_rb|tb_rc|tb_rd|tb_re|tb_rf|ra|rb|rc|rd|re|rf);

            //Update the configuration with the change
            result_cell[0][x][y] = a ^ cha;
            result_cell[1][x][y] = b ^ chb;
            result_cell[2][x][y] = c ^ chc;
            result_cell[3][x][y] = d ^ chd;
            result_cell[4][x][y] = e ^ che;
            result_cell[5][x][y] = f ^ chf;
            result_cell[6][x][y] = r ^ chr;
        }
    }
}
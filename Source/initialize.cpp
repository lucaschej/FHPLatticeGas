// initialize.cpp
#include "config.h"
#include <vector>
#include <random>
#include <cstdint>
#include <cmath>

using namespace std;

// Inicializa las celdas con velocidades aleatorias
int initialize(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> >& nsbit, mt19937 generator, double prob)
{
    int i, j, k, b; // Contadores
    int n_particles;
    uint64_t aux_bit; // Variables auxiliares
    uint64_t bit_to_add;
    uniform_real_distribution<double> rbit(0, 1); // Se usa para agregar bits aleatorios

    n_particles = 0; // Contador de partículas
    for (i = 0; i < XMAX; i++)
    {
        for (j = 0; j < YMAX; j++)
        {
            for (k = 0; k < 7; k++)
            {
                uint64_t aux_bit = 0; // Inicializar
                // Obtener un número aleatorio con un p% de unos
                for (b = 0; b < 64; b++)
                {
                    bit_to_add = rbit(generator) <= prob ? uint64_t(1) : uint64_t(0);
                    aux_bit = bit_to_add ^ (aux_bit << 1); // Agregar el uno o cero
                }
                cell[k][i][j] = aux_bit; // Agregar a la celda
                n_particles += __builtin_popcount(cell[k][i][j]); // Contar número de unos
                result_cell[k][i][j] = 0;
            }

            nsbit[i][j] = ~uint64_t(0); // No sólido = No pared -> 1 en el fluido
        }
    }
    return n_particles;
}

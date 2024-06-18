#!/bin/bash

# Ejecutar el script Allclean
./Allclean.sh

# Compilar los archivos .cpp
g++ -o lattice_gas lattice_gas.cpp initialize.cpp

# Ejecutar el ejecutable
./lattice_gas

# Ejecutar el script de Python
python3 plotter.py


#!/bin/bash

# Nombre del script
script_name=$(basename "$0")

# Carpeta en la que se encuentra el script
script_dir=$(dirname "$0")

# Archivos protegidos que no se deben borrar
protected_files=(
    "lattice_gas.cpp"
    "config.h"
    "initialize.h"
    "initialize.cpp"
    "Allrun.sh"
    "Allclean.sh"
    "plotter.py"
    "collisionFHP_I.cpp"
    "collisionFHP_I.h"
    "collisionFHP_II.cpp"
    "collisionFHP_II.h"
    "collisionFHP_III.cpp"
    "collisionFHP_III.h"
    "propagation.cpp"
    "propagation.h"
    "measure.h"
    "measure.cpp"
    "get_results.h"
    "get_results.cpp"
    "promediate.cpp"
    "promediate.h"
)

# Navegar a la carpeta donde se encuentra el script
cd "$script_dir"

# Borrar todos los archivos excepto el script y los archivos protegidos
for file in *; do
    if [[ "$file" != "$script_name" && ! " ${protected_files[@]} " =~ " $file " ]]; then
        rm -f "$file"
    fi
done

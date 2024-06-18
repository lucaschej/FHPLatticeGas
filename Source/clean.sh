#!/bin/bash

# Nombre del script
script_name=$(basename "$0")

# Carpeta en la que se encuentra el script
script_dir=$(dirname "$0")

# Nombre del archivo que no se debe borrar
protected_file_1="lattice_gas.cpp"
protected_file_2="config.h"
protected_file_3="initialize.h"
protected_file_4="initialize.cpp"
# Navegar a la carpeta donde se encuentra el script
cd "$script_dir"

# Borrar todos los archivos excepto el script y el archivo protegido
for file in *; do
  if [[ "$file" != "$script_name" && "$file" != "$protected_file_1" && "$file" != "$protected_file_2" && "$file" != "$protected_file_3" && "$file" != "$protected_file_4" ]]; then
    rm -f "$file"
  fi
done


#!/bin/bash

# Ejecutar el script Allclean
./Allclean.sh

# Obtener todos los archivos .cpp en el directorio actual
cpp_files=$(find . -maxdepth 1 -name "*.cpp")

# Compilar los archivos .cpp encontrados
g++ -o lattice_gas ${cpp_files}

# Verificar si la compilación tuvo éxito
if [ $? -eq 0 ]; then
    echo "Compilación exitosa. Ejecutando el programa..."
    # Ejecutar el ejecutable
    ./lattice_gas

    # Verificar si la ejecución tuvo éxito
    if [ $? -eq 0 ]; then
        echo "Programa ejecutado correctamente."
        # Ejecutar el script de Python
        python3 plotter.py
    else
        echo "Error al ejecutar el programa lattice_gas."
    fi
else
    echo "Error durante la compilación."
fi


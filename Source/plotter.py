import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# Obtener la ruta actual
current_path = os.getcwd()

# Nombres de los archivos
files = [f'data{i}.txt' for i in range(10)]

# Parámetros de la malla extendidos
n = 40
m = 25
dx = 1
dy = np.sqrt(3) / 2

# Inicializar listas de coordenadas de la malla
x_coords = []
y_coords = []

# Generar las coordenadas de los nodos
for j in range(m):
    for i in range(n):
        x = i * dx + (j % 2) * (dx / 2)
        y = j * dy
        x_coords.append(x)
        y_coords.append(y)

# Filtrar y graficar los datos de cada archivo
for i, file_name in enumerate(files):
    # Leer el archivo
    file_path = os.path.join(current_path, file_name)
    data = pd.read_csv(file_path, sep=' ', header=None, names=['x', 'y'])

    # Crear la figura y los ejes
    fig, ax = plt.subplots(figsize=(15, 10))

    # Graficar la malla hexagonal
    ax.scatter(x_coords, y_coords, edgecolor='lightgray', facecolor='none', s=20, linewidths=0.5)
    for j in range(m):
        for k in range(n):  # Usar k como índice para evitar conflicto con el índice i
            x = k * dx + (j % 2) * (dx / 2)
            y = j * dy
            # Vecinos posibles
            neighbors = [
                (x + dx, y),
                (x - dx, y),
                (x + dx / 2, y + dy),
                (x - dx / 2, y + dy),
                (x + dx / 2, y - dy),
                (x - dx / 2, y - dy)
            ]
            for neighbor in neighbors:
                if (0 <= neighbor[0] <= (n - 1) * dx + (m % 2) * (dx / 2) and
                    0 <= neighbor[1] <= (m - 1) * dy):
                    ax.plot([x, neighbor[0]], [y, neighbor[1]], color='lightgray', linewidth=0.5)

    # Graficar los datos
    ax.scatter(data['x'], data['y'], color='blue', s=50)  # Cambiar color a azul y aumentar tamaño de los puntos
    ax.set_title(f'Coordenadas desde {file_name}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_xlim(-1, n)
    ax.set_ylim(-1, (m - 1) * dy + dy)
    ax.grid(False)

    # Guardar la imagen en el directorio actual
    save_path = os.path.join(current_path, f'x_{i+1}.png')
    plt.savefig(save_path)
    plt.close()  # Cerrar la figura para liberar memoria

print("Las imágenes han sido generadas y guardadas con éxito.")


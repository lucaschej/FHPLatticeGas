import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
animation.writers.list()  # Verificar los escritores disponibles
import pylab as pylab
from tqdm import tqdm

##########################################################################

rows = 300
columns = 300
n_particles = 300

def initialize():
    vectorized_object = np.vectorize(lambda obj: node(False))
    config = np.zeros((rows, columns)).astype(object) 
    config[:] = vectorized_object(config[:])
    
    #adding the walls:
    vectorized_object_wall = np.vectorize(lambda obj: node(True))
    
    config[0, :] = vectorized_object_wall(config[0, :])
    config[rows-1, :] = vectorized_object_wall(config[rows-1, :])
    config[:, 0] = vectorized_object_wall(config[:, 0])
    config[:, columns-1] = vectorized_object_wall(config[:, columns-1])
            
     
    return config
# a function to observe and visualize the grid:
def observe(config):
    vectorized_x = np.vectorize(lambda obj: obj.occupied)
    view_grid = vectorized_x(config)
    
    plt.figure(figsize = (10, 8))
    pylab.cla()
    pylab.imshow(view_grid, vmin = 0, vmax = 1, cmap = pylab.cm.binary)

##########################################################################

class node: 
    def __init__(self, wall):
        self.ni_s = [0, 0, 0, 0, 0, 0]
        self.wall = wall
        if not wall:
            self.occupied = 0
        else:
            self.occupied = 1

##########################################################################

def get_neighbors(coordinates, config):
    row, column = coordinates
    neighbors = []
    neighbors.append([row, (column+1)%columns]) #right
    neighbors.append([(row-1)%rows, (column+1)%columns]) # top right
    neighbors.append([(row-1)%rows, (column-1)%columns]) #top left
    neighbors.append([row, (column-1)%columns]) #left
    neighbors.append([(row+1)%rows, (column-1)%columns]) #bottom left
    neighbors.append([(row+1)%rows, (column+1)%columns]) # bottom right
    
    return neighbors

##########################################################################

#getting all possible coordinates and choosing N radnomly:
i_coords, j_coords = np.meshgrid(range(rows), range(columns), indexing='ij')
coordinate_grid = np.array([i_coords, j_coords])
all_coordinates = []
for i in range(1, rows-1):
    for j in range(1, columns-1):
        all_coordinates.append(list(coordinate_grid[:, i, j]))
        
        
def fill_random(N, config, d = False):
    samples = random.sample(all_coordinates, N)
    for row,column in samples:
        config[row, column].occupied = 1
        
        our_neighbors = get_neighbors([row, column], config)
        while True:
            if not d:
                direction = random.choice([i for i in range(6)])
            else:
                direction = d-1
            n_row, n_col = our_neighbors[direction]
            if not config[n_row,n_col].wall:
                if config[n_row,n_col].occupied == 0:
                    config[n_row,n_col].ni_s[direction]=1
                    break
     
    return config

##########################################################################

my_grid = initialize()
my_grid = fill_random(n_particles, my_grid)
observe(my_grid)

##########################################################################

def D_i(n, i):
    """
    The head-on collision function
    """
    return n[i%6]*n[(i+3)%6]*(1-n[(i+1)%6])*(1-n[(i+2)%6])*(1-n[(i+4)%6])
def T_i(n, i):
    """
    The three-body collision function
    """
    return n[i]*n[(i+2)%6]*n[(i+4)%6]*(1-n[(i+1)%6])*(1-n[(i+3)%6])*(1-n[(i+5)%6])
def collision_factor(coords, config, i):
    """
    The collision factor Ω
    """
    q = random.choice([0, 1])
    n = config[coords[0], coords[1]].ni_s
    coll = -D_i(n, i)+q*D_i(n, (i-1)%6) + (1-q)*D_i(n, (i+1)%6) - T_i(n, i)+T_i(n, (i+3)%6)
    return coll

##########################################################################

def update(config):
    next_config = initialize()
    
    for row in range(rows):
        for column in range(columns): # for each cell in the grid
            if not config[row, column].wall: # if it's not a wall
                 #get all of its neighbors    
                our_neighbors = get_neighbors([row, column], config)
                
                
                bounce_direc = [] # make this list to store if it's going to bounce off a wall
                
                for i in range(6):# for each neighbor
              # start by letting the particle that is coming along i (from i+3) get in
                    if config[row, column].ni_s[i] >0:
                    # then update so that the coming particle can enter the cell
                        next_config[row, column].occupied = 1
                        th_row, th_col = our_neighbors[(i+3)%6] 
                       # the cell that sent that particle
                        if not next_config[th_row, th_col].wall: 
                            next_config[th_row, th_col].occupied = 0         # update it to 0 since the particle moved already
                        next_config[row, column].ni_s[i] = 0
                    #then for each neghibor, check if something will pop inside them from this cell and change their n_i
                    other_row, other_col = our_neighbors[i]
                    #calculate its value using the FHP equation
                    next_config[other_row,other_col].ni_s[i] = config[row, column].ni_s[i]+\
                                collision_factor([row, column], config, i)
                    
                    # check if any particle is hitting the wall:
                    if next_config[other_row,other_col].wall and next_config[other_row,other_col].ni_s[i]>0:
                        next_config[other_row,other_col].ni_s[i] = 0
                        bounce_direc.append((i+3)%6) #store their bounce off directions
                        
                # if any particles were hiting the wall, send them back:        
                if bounce_direc:
                    next_config[row, column].ni_s[bounce_direc[0]] = 1
    return next_config

##########################################################################

def turn_state(config):
    """
    A function to visualize the states
    """
    view_grid = np.copy(config)
    for i in range(len(config)):
        for j in range(len(config)):
            if type(config[i,j]) is not int:
                view_grid[i, j] = config[i,j].occupied
    view_grid = view_grid.astype(int)
    return view_grid

def build_animation(states, filename):
    fig = plt.figure()
    plt.axis('off')

    ims = []
    for state in states:
        im = plt.imshow(state, animated=True, cmap='binary')
        ims.append([im])

    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
    ani.save(filename, writer='pillow')

# Main code
my_grid = initialize()
my_grid = fill_random(n_particles, my_grid)
observe(my_grid)

states = []
states.append(turn_state(my_grid))
next_grid = update(my_grid)
states.append(turn_state(next_grid))
states.append(turn_state(next_grid))

for i in tqdm(range(200)):
    next_grid = update(next_grid)
    states.append(turn_state(next_grid))

# Cambiar el nombre del archivo a .gif
build_animation(states, "video.gif")


##########################################################################
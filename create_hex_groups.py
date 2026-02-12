import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
#from matplotlib.colormaps import get_cmap
from pylmgc90 import pre
import sys
import random
import json


def rotation_matrix(phi):
    return np.array([
        [np.cos(phi), -np.sin(phi)],
        [np.sin(phi),  np.cos(phi)]
    ])


def create_hexagon(center, radius, phi):
    angles = np.linspace(0, 2*np.pi, 7)[:-1]
    angles += phi
    x = center[0] + radius * np.cos(angles)
    y = center[1] + radius * np.sin(angles)
    return np.column_stack((x, y))

def generate_hexagon_groups(nrows=10, ncols=12, R=1.0, dR=0, disk_R_min=3, disk_R_max=10, phi=0):
    #hexagons = []
    #disks = []
    hexagons_pos = []
    for i in range(nrows):
        for j in range(ncols):
            pos = np.array([1+i*3/2, (1+i%2+2*j)*np.sqrt(3)/2])*R
            hexagons_pos.append(pos)
            #verts = create_hexagon(pos, R, phi)
            #hexagons.append(verts)
    
    Lx = (0.5 + nrows*3/2)*R
    Ly = (1 + 2*ncols) * np.sqrt(3)/2 * R
    
    nb_particles = int(10*(Lx/disk_R_min)*(Ly/disk_R_min))
    radii = pre.granulo_Random(nb_particles, disk_R_min, disk_R_max, seed)

    nb_laid_particles, coors, radii = pre.depositInBox2D(radii,2*Lx,2*Ly)
    coors = coors - np.array([Lx/2, Ly/2])
        
    # store disks that are within the box [(R,R), (Lx-R, Ly-R)]
    disks_pos = []
    disks_rad = []
    for i in range(nb_laid_particles):
        center = coors[i]
        r = radii[i]
        if center[0] + r >= R and center[0] - r <= Lx-R and center[1] + r >= R and center[1] - r <= Ly-R:
            disks_pos.append(center)
            disks_rad.append(r)
    
    plot_hex_disk('hexagon_circle_all.png', hexagons_pos, R, dR, disks_pos, disks_rad, Lx, Ly, phi)
    
    hexagons_pos = np.array(hexagons_pos)
    hexagons_mask = np.zeros(len(hexagons_pos)).astype(bool)
    
    hexagons_pos_grouped = []
    for i in range(len(disks_pos)):
        pos = disks_pos[i]
        rad = disks_rad[i] + dR
        dists = np.linalg.norm(hexagons_pos-pos, axis=1)
        mask = dists<=rad
        hexagons_pos_grouped.append(hexagons_pos[mask & ~hexagons_mask])
        hexagons_mask = hexagons_mask | mask
    
    plot_hexagon_groups('hexagons_grouped.png', hexagons_pos_grouped, disks_pos, disks_rad, Lx, Ly, phi, dR)
    
    N_init = len(hexagons_pos)
    N_final = sum([len(hexagons) for hexagons in hexagons_pos_grouped])
    
    porosity = (N_init - N_final)/N_init
    print('N_init:', N_init)
    print('N_final:', N_final)
    print('porosity:', porosity)
    
    # prepare data to save as a json file
    data = {
            'parameters': {
                'radius': R,
                'phi': phi,
                'num_grs': len(hexagons_pos_grouped),
                'Lx': Lx,
                'Ly': Ly
            },
            'groups': [],
            'hexagons': []
        }
    hex_id = 1 # Hexagon id starts from 1 in lmgc90
    for i, positions in enumerate(hexagons_pos_grouped):
        N = len(positions)
        IDs = [hex_id + i for i in range(N)]
        data['groups'].append(IDs)
        hex_id = hex_id + N
        for j, pos in enumerate(positions):
            data['hexagons'].append({
                'id': IDs[j],
                'center': tuple(pos),
                'gr_id' : i
            })
    
    # Save to file
    filename = 'hexagons_data_grouped.json'
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)
        
    
def plot_hex_disk(filename, hexagons_pos, R, dR, disks_pos, disks_rad, Lx, Ly, phi=0):
    
    fig, ax = plt.subplots()
    for pos in hexagons_pos:
        verts = create_hexagon(pos, R, phi)
        hex_patch = Polygon(verts, closed=True,
                            edgecolor='r', fill=False)
        ax.add_patch(hex_patch)
    
    # Draw disks
    for i in range(len(disks_pos)):
        circle = plt.Circle(disks_pos[i], disks_rad[i],
                            edgecolor='k', fill=False, linestyle='-')
        ax.add_patch(circle)
        
        circle = plt.Circle(disks_pos[i], disks_rad[i]+dR,
                            edgecolor='k', fill=False, linestyle='--')
        ax.add_patch(circle)

    ax.set_aspect('equal')
    ax.set_xlim(-Lx*0.05, Lx*1.05)
    ax.set_ylim(-Ly*0.05, Ly*1.05)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()

def plot_hexagon_groups(filename, hexagons_pos_grouped, disks_pos, disks_rad, Lx, Ly, phi=0, dR=0):
    fig, ax = plt.subplots()
    num_groups = len(hexagons_pos_grouped)
    
    colors = []
    
    for i in range(num_groups):
        color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
        colors.append(color)
        grouped_pos = hexagons_pos_grouped[i]
        for pos in grouped_pos: # grouped_pos is an array
            verts = create_hexagon(pos, R, phi)
            hex_patch = Polygon(verts, closed=True,
                                edgecolor=color, fill=False)
            ax.add_patch(hex_patch)
    
    # Draw disks
    for i in range(len(disks_pos)):
        circle = plt.Circle(disks_pos[i], disks_rad[i]+dR,
                            edgecolor='k', fill=False, linestyle='--')
        ax.add_patch(circle)

    ax.set_aspect('equal')
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    fig, ax = plt.subplots()
    for i in range(num_groups):
        grouped_pos = hexagons_pos_grouped[i]
        for pos in grouped_pos:
            verts = create_hexagon(pos, R, phi)
            hex_patch = Polygon(verts, closed=True,
                                edgecolor='k', fill=True, facecolor=colors[i])
            ax.add_patch(hex_patch)
    ax.set_aspect('equal')
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    plt.savefig(filename[:-4]+'_without_circles.png', dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    
# Example
nrows = 10
ncols = 8
R = 1e-1
dR = 0.5*R
disk_R_min = 2*R 
disk_R_max = 5*R 
phi = 0

if '--norand' in sys.argv:
  seed = 1
else:
  seed = None

generate_hexagon_groups(nrows, ncols, R, dR, disk_R_min, disk_R_max, phi)

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches #
import sys, random, json
from pylmgc90 import chipy, pre
from time import time


def get_vertices(pos, radius, n=6, phi=0):
    dtheta = 2*np.pi/n
    verts = []
    for i in range(n):
        v = (round(pos[0] + radius*np.cos(phi + i*dtheta), 5), round(pos[1] + radius*np.sin(phi + i*dtheta), 5))
        verts.append(v)
    return verts
    
def imposedFy(t):
    # Set time evolution of the force on the top wall 
    # [k<10: Fy=0; 10<k<110: Fy=-f*(k-10)/100; k>110: Fy=-f]
    #
    if t <10*dt :
      return 0.
    elif 10*dt <= t and t <110*dt :
      return Fy*(t-(10*dt))/(100*dt)
    else :
      return Fy

def get_mecaMailx_nodes_of_edge(pos_u, pos_v, idBody):
    """
    Returns the nodes that are inside the triangle made by vertices pos_u, pos_v and the centre of the hexagon.
    Uses Barycentric Coordinate System to determine if point P(px, py) is inside triangle ABC.
    """
    nb_nodes = chipy.mecaMAILx_GetNbNodes(idBody)
    if nb_nodes < 1:
        print('nb_nodes < 1')
        breakpoint()
        sys.exit()
    pos_nodes = np.zeros((nb_nodes, 2))   # 2D
    for idNode in range(1, nb_nodes + 1):
        pos_nodes[idNode - 1, :] = chipy.mecaMAILx_GetNodeCooref(idBody, idNode)
        
    cx, cy = pos_nodes[:,0].mean(), pos_nodes[:,1].mean()
    
    ax, ay = pos_u
    bx, by = pos_v
    
    # For numerical precision use pos_u = pos_u + epsilon,  pos_v = pos_v + epsilon
    # epsilon = 5% of centre to vertex distance
    ax = ax + (ax-cx)*0.05
    ay = ay + (ay-cy)*0.05
    bx = bx + (bx-cx)*0.05
    by = by + (by-cy)*0.05
    
    # Vectors from A to other vertices
    v0x, v0y = cx - ax, cy - ay 
    v1x, v1y = bx - ax, by - ay 
    # Compute dot products
    dot00 = v0x * v0x + v0y * v0y 
    dot01 = v0x * v1x + v0y * v1y 
    dot11 = v1x * v1x + v1y * v1y 
    # Compute the denominator
    denom = (dot00 * dot11 - dot01 * dot01) 
    inv_denom = 1.0 / denom #
    
    node_IDs = [] # starts from 0 ## mind it ##
    for i in range(nb_nodes):
        px, py = pos_nodes[i]
        
        # Vector from A to the node
        v2x, v2y = px - ax, py - ay
        
        # Compute dot products
        dot02 = v0x * v2x + v0y * v2y
        dot12 = v1x * v2x + v1y * v2y
    
        # Compute barycentric coordinates (u, v)
        u = (dot11 * dot02 - dot01 * dot12) * inv_denom
        v = (dot00 * dot12 - dot01 * dot02) * inv_denom
    
        # Check if the node is in the triangle
        if (u >= 0) and (v >= 0) and (u + v <= 1):
            node_IDs.append(i)
        #
    return np.array(node_IDs)

def calculate_probability_weights_nodes(idBody):
    nb = chipy.mecaMAILx_GetNbNodes(idBody)
    node_pos = []
    for i in range(1, nb+1):
        pos = chipy.mecaMAILx_GetNodeCooref(idBody, i)
        node_pos.append(pos)
    node_pos = np.array(node_pos)
    centre = (node_pos[:,0].mean(), node_pos[:,1].mean())
    node_pos = node_pos - centre # weight is proportional to the distance from the centre.
    w = node_pos[:,0]**2 + node_pos[:,1]**2
    #w = np.sqrt(w)
    w[w==w.min()] = w[w!=w.min()].min() # to eliminate 0 weights
    w = w/w.sum()
    return w

def select_nodes(idBody, nb_select, p, neighborsToNodes, boundaryNodes):
    # p -> weights
    #
    probable_nodes = [] 
    node_ids = np.arange(1, chipy.mecaMAILx_GetNbNodes(idBody)+1)
    for idNode in node_ids[display_T[idBody-1] >= Tmax]:
        probable_nodes.extend(neighborsToNodes[idNode])
    probable_nodes.extend(boundaryNodes)
    probable_nodes = np.unique_values(probable_nodes)
    p = p[probable_nodes-1]
    p = p/p.sum()
    selected_nodes = np.random.choice(probable_nodes, nb_select, replace=False, p=p)
    return selected_nodes
    
def run_DEM_steps():
    #
    chipy.utilities_logMes('COMPUTE Fext')
    chipy.ComputeFext()
    #
    chipy.utilities_logMes('COMPUTE Fint')
    chipy.ComputeBulk()
    #
    chipy.utilities_logMes('ASSEMB RHS')
    chipy.AssembleMechanicalRHS()
    #
    chipy.utilities_logMes('COMPUTE Free Vlocy')
    chipy.ComputeFreeVelocity()
    #
    chipy.utilities_logMes('SELECT PROX TACTORS')
    chipy.SelectProxTactors()
    #
    chipy.utilities_logMes('RESOLUTION' )
    chipy.RecupRloc()
    #
    chipy.ExSolver(solver_type, norm, tol, relax, gs_it1, gs_it2)
    chipy.UpdateTactBehav()
    #
    chipy.StockRloc()
    #
    chipy.utilities_logMes('COMPUTE DOF, FIELDS, etc.')
    chipy.ComputeDof()
    #
    chipy.utilities_logMes('UPDATE DOF, FIELDS')
    chipy.UpdateStep()
    #
    chipy.utilities_logMes('WRITE OUT')
    chipy.WriteOut(freq_write)
    #
    chipy.utilities_logMes('VISU & POSTPRO')
    chipy.WriteDisplayFiles(freq=freq_display,T=('mecafe','node',display_T))
    chipy.WritePostproFiles()

def update_T():
    pass

def update_edge_type(chosen_edges):
    #        
    # change type if all mecaNodes have reacted
    for e in chosen_edges:
        #
        # If all the FEM nodes associated to the edge has fully reacted change the edge type to type 0. 
        # It will not be selected again.
        
        # neighbor hexagons
        ID_nn = list(G.edges[e]["nn_hex"].keys())
        
        # temps of the nodes that belong the edge
        T_nn1 =  display_T[ID_nn[0]-1][G.edges[e]["nn_hex"][ID_nn[0]]] # from hexagon1 
        # temporary solution for edge with 1 polygon *** ID_nn[-1] ***
        T_nn2 =  display_T[ID_nn[-1]-1][G.edges[e]["nn_hex"][ID_nn[-1]]] # from hexagon2 
        
        # update the type
        if all(T_nn1 >= Tmax) and all(T_nn2 >= Tmax): 
            edges_of_type[G.edges[e]["type"]].remove(e)
            G.edges[e]["type"] = 0
            edges_of_type[0].append(e)
        
    # for type 4 edges update the types of its neighbours
    chosen_edges = [e for e in chosen_edges if G.edges[e]["type"]==4]
    for edge in chosen_edges:
        u, v = edge
        nn_edges = list(G.edges(u)) + list(G.edges(v))
        nn_edges.remove((u,v))
        nn_edges.remove((v,u))
        nn_edges = [tuple(sorted(e)) for e in nn_edges]
        for e in nn_edges:
            if G.edges[e]["type"] == 1:
                # update the type of single edge
                if G.edges[e]["edge_gr"] == None:
                    G.edges[e]["type"] = 4
                    edges_of_type[4].append(e)
                    edges_of_type[1].remove(e)
                
                # if the edge is in a connected chain of edges update type for all of them
                else: 
                    edge_gr = group_of_edges[G.edges[e]["edge_gr"]]
                    for e2 in edge_gr:
                        G.edges[e2]["type"] = 4
                        edges_of_type[4].append(e2)
                        edges_of_type[1].remove(e2) #!! ValueError: list.remove(x): x not in list
    #

def select_mecaNodes(chosen_edges, NselectMecaNodes):
    IdBodies = []
    mecaNodes = []
    for e in chosen_edges:
        G.edges[e]["nselected"] += 1
        for ID in G.edges[e]["nn_hex"].keys():
            nodes = G.edges[e]["nn_hex"][ID] # FEM nodes
            T = display_T[ID-1][nodes]
            
            # skip if fully reacted
            if all(T >= Tmax): continue
            #
            IdBodies.append(ID)
            #
            nodes = nodes[T<Tmax]  
            
            if len(nodes) <= NselectMecaNodes:
                mecaNodes.append(nodes)
            else:
                w = w_nodes_allPolyg[ID-1][nodes] 
                w = w/w.sum()
                mecaNodes.append(np.random.choice(nodes, NselectMecaNodes, p=w, replace=False))
                #
    return IdBodies, mecaNodes


##===================================================##

## Prepare data for networx

filename = 'hexagons_data_grouped.json'
with open(filename, 'r') as foo:
    data = json.load(foo)

radius = data['parameters']['radius']
Lx = data['parameters']['Lx']
Ly = data['parameters']['Ly']
phi = data['parameters']['phi']
hexagons_info =  data['hexagons']
groups = data['groups']
num_grs = data['parameters']['num_grs']


hexagons_vert_data = []
hexagons_pos = []
group_ids_of_hexagons = []
for hexagon in hexagons_info:
    if hexagon['gr_id'] == -1: continue
    hexagons_pos.append(hexagon['center'])
    hexagons_vert_data.append(get_vertices(hexagon['center'], radius, n=6, phi=phi))
    group_ids_of_hexagons.append(hexagon['gr_id'])

nodes_pos_nx = []
edges_all_hex = []

for hex_vertices in hexagons_vert_data:
    for vert in hex_vertices:
        if vert not in nodes_pos_nx:
            nodes_pos_nx.append(vert)

for hex_vertices in hexagons_vert_data:
    num_verts = len(hex_vertices)
    edges = []
    n0 = nodes_pos_nx.index(hex_vertices[0])
    for i  in range(1,num_verts+1):
        n1 = nodes_pos_nx.index(hex_vertices[i%num_verts])
        edges.append(tuple(sorted([n0,n1])))
        n0 = n1
    edges_all_hex.append(edges)

# Networx
G = nx.Graph()

# 1. Build the network from vertex connections
G.add_nodes_from(range(len(nodes_pos_nx)))
for edges in edges_all_hex:
    G.add_edges_from(edges) 

# 1.1 store neighbor hexagons data of the edges
for e in G.edges:
    nn_hexagons = [i+1 for i,es in enumerate(edges_all_hex) if e in es]
    G.edges[e]["nn_hex"] = { ID:[] for ID in nn_hexagons} # store hexagon IDs and the FEM nodes correspond to this edge
    G.edges[e]["edge_gr"] = None # number of group if member of a egde chain
    G.edges[e]["nselected"] = 0 # number of times the edge has been selected
        

# store edges by their types as keys
# Create a dictionary where each value defaults to an empty list
edges_of_type = {0:[], 1:[], 2:[], 3:[], 4:[]}

for e in G.edges:
    nn_hex = list(G.edges[e]["nn_hex"].keys())
    if len(nn_hex)==1:
        G.edges[e]["type"] = 3
        edges_of_type[3].append(e)
    else:
        g1 = group_ids_of_hexagons[nn_hex[0]-1]
        g2 = group_ids_of_hexagons[nn_hex[1]-1]
        if g1 == g2:
            G.edges[e]["type"] = 1
            edges_of_type[1].append(e)
        else:
            G.edges[e]["type"] = 2
            edges_of_type[2].append(e)

#-------------------------------------------------------------#
# Make the outer edges (near the boundary) of type 3 to type 1 
#-------------------------------------------------------------#

x_min = radius
x_max = Lx - radius
y_min = radius
y_max = Ly - radius

for e in G.edges:
    if G.edges[e]["type"] == 3:
        u, v = e
        pos_u = nodes_pos_nx[u]
        if pos_u[0] < x_min or pos_u[0] > x_max or pos_u[1] < y_min or pos_u[1] > y_max:
            G.edges[e]["type"] = 1
            edges_of_type[1].append(e)
            edges_of_type[3].remove(e)

#-------------------------------------------------------------#

# change the type of the edges, that are connected to type 2 or type 3 edges, to type 4
for edge in edges_of_type[2] + edges_of_type[3]:
    u, v = edge
    nn_edges = list(G.edges(u)) + list(G.edges(v))
    nn_edges.remove((u,v))
    nn_edges.remove((v,u))
    nn_edges = [tuple(sorted(e)) for e in nn_edges]
    for e in nn_edges:
        if G.edges[e]["type"] == 1:
            G.edges[e]["type"] = 4
            edges_of_type[4].append(e)
            edges_of_type[1].remove(e)

# group of edges ## broken edges
group_of_edges = dict()
indices = np.random.choice(np.arange(len(edges_of_type[1])), 2)
for ID in indices:
    u, v = edges_of_type[1][ID]
    nn_edges = list(G.edges(u)) + list(G.edges(v))
    nn_edges.remove((v,u))
    nn_edges = [tuple(sorted(e)) for e in nn_edges if G.edges[sorted(e)]["type"] == 1]
    group_of_edges[ID] = nn_edges
    for e in nn_edges:
        G.edges[e]["edge_gr"] = ID  

#'''
#------------------ Visualize the network ----------------------#
# edge colors
colors = ['r','gray','tab:green','cornflowerblue','k']
edge_colors =[]
for e in G.edges:
    edge_colors.append(colors[G.edges[e]["type"]])
edge_type = ['-' if G.edges[e]["edge_gr"]==None else '--' for e in G.edges]

# 4. Visualization using the actual coordinates as positions
plt.figure(figsize=(7.5, 6.2))
nx.draw(G, nodes_pos_nx, with_labels=False, node_size=50, edge_color=edge_colors, style=edge_type, width=3)

ax = plt.gca()
gr_color = ["#{:06x}".format(random.randint(0, 0xFFFFFF)) for i in range(num_grs)]

for i, pos in enumerate(hexagons_pos):
    ax.text(pos[0], pos[1], i+1, color=gr_color[group_ids_of_hexagons[i]])

type_1_edge = mpatches.Patch(color=colors[1], lw=0, label='Edge type = 1')
type_2_edge = mpatches.Patch(color=colors[2], lw=0, label='Edge type = 2')
type_3_edge = mpatches.Patch(color=colors[3], lw=0, label='Edge type = 3')
type_4_edge = mpatches.Patch(color=colors[4], lw=0, label='Edge type = 4')
plt.legend(handles=[type_1_edge, type_2_edge, type_3_edge, type_4_edge])


plt.axis('equal')
plt.savefig('hex_network.png', dpi=300, bbox_inches='tight')
#plt.show()
plt.close()
#------------------------------------------------------#
#'''

## DEM environment :: declare variables
#
dim = 2                   # space dimension
mhyp = 1                  # modeling hypothesis (1 = plain strain, 2 = plain stress, 3 = axi-symmetry)
dt = 1.e-3                # time evolution parameters
theta = 0.55              # theta integrator parameter
deformable = 1            # deformable  yes=1, no=0
Rloc_tol = 5.e-2          # interaction parameters

# nlgs parameters
tol = 1e-4
relax = 1.0
norm = 'Quad '
gs_it1 = 50
gs_it2 = 50
#solver_type='Stored_Delassus_Loops         '
solver_type='Stored_Delassus_Loops'

T0 = 0.                   # initial temp
Tmax = 5000.              # max temp
nb_steps_T_incr = 100     # steps to increase T from T0 to Tmax
nb_loops_max = 200        # max while loops 
nb_select_mecaNodes = 4   # select FEM nodes at each simulation step
#nb_select_edges = 10      # select edges (networkx) at each simulation step
reaction_rate = 0.2       # 20% of active edges
Fy = -1e6                 # Force on the top wall

freq_write   = 1500       # write parameter
freq_display = nb_steps_T_incr         # display parameters

# Set time evolution of the force on the top wall
nb_steps_Fy = nb_loops_max * nb_steps_T_incr + 3*nb_steps_T_incr
tf = nb_steps_Fy * dt
pre.writeEvolution(f=imposedFy, instants=np.linspace(0., tf, nb_steps_Fy) ,path='DATBOX/', name='Fy.dat')

t0 = time()

## DEM :: Initialize

chipy.Initialize()       # Initializing
chipy.checkDirectories() # checking/creating mandatory subfolders

chipy.utilities_DisableLogMes() # logMes

chipy.SetDimension(dim,mhyp) # Set space dimension
#
chipy.utilities_logMes('INIT TIME STEPPING')
chipy.TimeEvolution_SetTimeStep(dt)
chipy.Integrator_InitTheta(theta)
#
chipy.ReadDatbox(deformable) # read and load
#
# open display & postpro
chipy.utilities_logMes('DISPLAY & WRITE')
chipy.OpenDisplayFiles()
chipy.OpenPostproFiles()

chipy.PT2Dx_SetDisplayRadius(0.002)
chipy.CLALp_SetNonSymmetricDetection()
#
# set initial temp
display_T = []
for i in range(1,chipy.mecaMAILx_GetNbMecaMAILx()+1):
    T=T0*np.ones(chipy.mecaMAILx_GetNbNodes(i))
    chipy.mecaMAILx_SetScalarFieldByNode(i,1,T)
    display_T.append(T)

# since constant compute elementary mass matrices once
chipy.utilities_logMes('COMPUTE MASS')
chipy.ComputeMass()

# since constant compute elementary stiffness matrices once
chipy.utilities_logMes('COMPUTE STIFFNESS')
chipy.ComputeBulk()

# since constant compute iteration matrix once
chipy.utilities_logMes('ASSEMB KT')
chipy.AssembleMechanicalLHS()

## Link mecaMailx nodes to the Networx edges
# 
for e in G.edges:
    node_pos1 = nodes_pos_nx[e[0]]
    node_pos2 = nodes_pos_nx[e[1]]
    for ID in G.edges[e]["nn_hex"].keys():
        mecaNodes = get_mecaMailx_nodes_of_edge(node_pos1, node_pos2, ID)
        G.edges[e]["nn_hex"][ID] = mecaNodes

#
## Simulation part
#

# Initialize loading
for k in range(3*nb_steps_T_incr): # 3*nb_steps_T_incr
    if k%10 == 0: print('Loading step:', k)
    chipy.utilities_logMes('INCREMENT STEP')
    chipy.IncrementStep()
    #
    run_DEM_steps()
#

# Assign Gamma-distributed weights to the edges
scale2 =  1.0
scale3 =  1.0
scale4 =  1.0
shape2 = 5.0
shape3 = 10.0
shape4 = 2.0

T_incrmnt = (Tmax - T0)/nb_steps_T_incr
wall_pos = []
Nb = chipy.mecaMAILx_GetNbMecaMAILx()   # Number of polygons
w_nodes_allPolyg = [calculate_probability_weights_nodes(idBody) for idBody in range(1,Nb+1)] # probability weights for nodes

# number of FEM nodes that correspod to an edge
ID_nn_ = list(G.edges[edges_of_type[3][0]]["nn_hex"].keys())
nb_nodes_per_edge = len(G.edges[edges_of_type[3][0]]["nn_hex"][ID_nn_[0]])
# The number of times an edge is needed to be selected to ensure that all the FEM nodes has reacted
nselect_edge_max = np.ceil(nb_nodes_per_edge/nb_select_mecaNodes) 

for loop_id in range(nb_loops_max):
    print('loop id:', loop_id)
    
    edges234 = edges_of_type[2] + edges_of_type[3] + edges_of_type[4]
    
    # find accessible edge length (surface area in 3D)
    nselected_234 = [G.edges[e]["nselected"] for e in edges234]
    nselected_234 = np.array(nselected_234)
    edge_length = (nselect_edge_max - nselected_234)/nselect_edge_max
    edge_length = edge_length.sum()
    
    nb_select_edges = int(reaction_rate * edge_length)
    
    #if len(edges234)==0:
    #    break
    #if len(edges234) <= nb_select_edges:
    #    chosen_edges = edges234
    #else:
    #    weights2 = np.random.gamma(shape2, scale2, len(edges_of_type[2]))
    #    weights3 = np.random.gamma(shape3, scale3, len(edges_of_type[3]))
    #    weights4 = np.random.gamma(shape4, scale4, len(edges_of_type[4]))
    #    weights234 = np.concatenate([weights2, weights3, weights4])
    #    
    #    # Choose edges randomly
    #    prob = weights234 / np.sum(weights234)
    #    chosen_idx = np.random.choice(len(edges234), nb_select_edges, p=prob) # select nb_select_edges edges
    #    chosen_edges = [edges234[idx] for idx in chosen_idx]
    
    if nb_select_edges < 1: 
        break
        
    weights2 = np.ones(len(edges_of_type[2])) * shape2
    weights3 = np.ones(len(edges_of_type[3])) * shape3
    weights4 = np.ones(len(edges_of_type[4])) * shape4
    weights234 = np.concatenate([weights2, weights3, weights4])
    #    
    # Choose edges randomly
    prob = weights234 / np.sum(weights234)
    chosen_idx = np.random.choice(len(edges234), nb_select_edges, p=prob) # select nb_select_edges edges
    chosen_edges = [edges234[idx] for idx in chosen_idx]
    
    # select mecaNodes to increase T
    IdBodies, mecaNodes = select_mecaNodes(chosen_edges, nb_select_mecaNodes)
    
    for j in range(nb_steps_T_incr):
        #
        chipy.utilities_logMes('INCREMENT STEP')
        chipy.IncrementStep()
        #
        # Increase T
        for n, ID in enumerate(IdBodies):
            display_T[ID-1][mecaNodes[n]] += T_incrmnt
            display_T[ID-1][display_T[ID-1]>Tmax] = Tmax
            chipy.mecaMAILx_SetScalarFieldByNode(ID,1,display_T[ID-1])
            #
        run_DEM_steps()
        
        #store top wall position
        wall_pos.append([chipy.TimeEvolution_GetTime(), chipy.JONCx_GetCoor(4)[1]])
    
    # update edge type (type 1 of the neighbors of selected type 4 becomes type 4), fully reacted edges become type 0 
    update_edge_type(chosen_edges)
    
    # save top wall position with time
    if loop_id%10 == 0:
        np.savetxt('POSTPRO/top_wall_position.txt', np.array(wall_pos), header='time   postion')


# close display & postpro
#
chipy.CloseDisplayFiles()
chipy.ClosePostproFiles()

# this is the end
chipy.Finalize()

print('='*50)
print('total loop:',loop_id, 'last edges:', chosen_edges)
print('='*50)

# save top wall position with time
np.savetxt('POSTPRO/top_wall_position.txt', np.array(wall_pos), header='time   postion')


t = time() - t0
h = int(t/3600)
m = int((t - h*3600)/60)
s = t - h*3600 - m*60
print('Time needed to simulate: %ihr %im %is' %(h,m,s))


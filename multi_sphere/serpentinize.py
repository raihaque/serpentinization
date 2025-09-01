import numpy as np
# importing chipy module
from pylmgc90 import chipy
import pickle
import sys
from time import time

t0 = time()

# Initializing
chipy.Initialize()

# checking/creating mandatory subfolders
chipy.checkDirectories()

# logMes
# chipy.utilities_DisableLogMes()

#
# defining some variables
#

# space dimension
dim = 2

# modeling hypothesis ( 1 = plain strain, 2 = plain stress, 3 = axi-symmetry)
mhyp = 1

# time evolution parameters
dt = 1e-5
nb_steps = 100 # sphere expansion step
#radius = 1e-2
rad_incr_percentage = 25 # 10%
N_select = 10 # depends on reaction rate

# theta integrator parameter
theta = 0.5

# interaction parameters
Rloc_tol = 5.e-2

# nlgs parameters
tol = 1e-4
relax = 1.0
norm = 'Quad '
gs_it1 = 500
gs_it2 = 20
solver_type='Stored_Delassus_Loops         '

# write parameter
freq_write   = 100

# display parameters
freq_display = 25

#
# read and load
#

# Set space dimension
chipy.SetDimension(dim,mhyp)
#
chipy.utilities_logMes('INIT TIME STEPPING')
chipy.TimeEvolution_SetTimeStep(dt)
chipy.Integrator_InitTheta(theta)
#
chipy.ReadDatbox(deformable=False)

#
# open display & postpro
#

chipy.utilities_logMes('DISPLAY & WRITE')
chipy.OpenDisplayFiles()
chipy.OpenPostproFiles()

#
# simulation part ...
#

# ... calls a simulation time loop
# since constant compute elementary mass once
chipy.utilities_logMes('COMPUTE MASS')
chipy.ComputeMass()

def pressurize(n):
  chipy.utilities_logMes('WRITE OUT')
  chipy.WriteOut(freq_write)

  chipy.WriteDisplayFiles(freq_display)
  chipy.WritePostproFiles()
  for k in range(n):
    chipy.IncrementStep()
    
    chipy.ComputeFext()
    chipy.ComputeBulk()
    chipy.ComputeFreeVelocity()
    
    chipy.SelectProxTactors() 
    chipy.RecupRloc(Rloc_tol) 
    chipy.ExSolver(solver_type, norm, tol, relax, gs_it1, gs_it2)
    chipy.UpdateTactBehav()
    chipy.StockRloc() 
    chipy.ComputeDof()
    chipy.UpdateStep() 

    chipy.utilities_logMes('WRITE OUT')
    chipy.WriteOut(freq_write)

    chipy.WriteDisplayFiles(freq_display)
    chipy.WritePostproFiles()

def expand_disks(IDs, X, V, n):
  
  for k in range(n):
    chipy.IncrementStep()
    
    chipy.ComputeFext()
    chipy.ComputeBulk()
    chipy.ComputeFreeVelocity()
    
    for ID in IDs:
      chipy.DISKx_SetVdilation(int(ID), V[ID-1])
    
    chipy.SelectProxTactors() 

    chipy.RecupRloc(Rloc_tol) 

    chipy.ExSolver(solver_type, norm, tol, relax, gs_it1, gs_it2)
    chipy.UpdateTactBehav()

    chipy.StockRloc() 

    chipy.ComputeDof()

    chipy.UpdateStep() 
    
    dx = X*k/n
    for ID in IDs:
      chipy.DISKx_SetXdilation(int(ID), dx[ID-1])

    chipy.utilities_logMes('WRITE OUT')
    chipy.WriteOut(freq_write)

    chipy.WriteDisplayFiles(freq_display)
    chipy.WritePostproFiles()
    
def update_probable_IDs(probable_IDs, selected_IDs, contacts, serpentinized):
    probable_IDs = probable_IDs - set(selected_IDs)
    new_IDs = []
    for ID in selected_IDs:
        new_IDs = new_IDs + contacts[ID] # contact is a dictionary
        serpentinized[ID-1] = True # fortran indexing starts from 1
    new_IDs = set(new_IDs)
    for ID in new_IDs:
        if serpentinized[ID-1]:
            new_IDs = new_IDs - {ID}
    probable_IDs.update(new_IDs)
    
    return probable_IDs, serpentinized

def get_neighbor_info(Nb):
    print('find neighbor info..')
    positions = []
    for bid in range(1, Nb+1):      # IDs start at 1
        X = chipy.RBDY2_GetBodyVector('Coor_', bid)   # x, y, theta
        #y = chipy.RBDY2_GetBodyVector('X', bid, 1)   # y-coordinate
        positions.append(X[:2])
    
    radii = []
    for bid in range(1, Nb+1):
        r = chipy.DISKx_GetRadius(bid)
        radii.append(r)
    
    positions = np.array(positions)
    radii = np.array(radii)
    top_pos = positions[:,1] + radii
    top_pos = top_pos.max()
    r_max = radii.max()
    
    in_contact = {i:[] for i in range(1,Nb+1)}
    for i in range(1, Nb+1):
        for j in range(i+1, Nb+1):
            dsq = sum((np.array(positions[i-1]) - np.array(positions[j-1]))**2)
            if dsq <= 1.01*(radii[i-1] + radii[j-1])**2: # epsilon = 1% of radius
                in_contact[i].append(j)
                in_contact[j].append(i)
    
    top_layer = []
    for i in range(1, Nb+1):
        if positions[i-1][1] >= top_pos - r_max:
            top_layer.append(i)
    return in_contact, top_layer, radii

# press without expansion
pressurize(60)

#with open('./DATBOX/neighbors.pkl','rb') as foo:
#    in_contact, toplayer, radii = pickle.load(foo)

Nb = chipy.DISKx_GetNbDISKx()
in_contact, toplayer, radii = get_neighbor_info(Nb)

rad_increment = radii * rad_incr_percentage/100
rad_increment_rate = rad_increment/(nb_steps*dt)



#'''
serpentinized = [False for i in range(Nb)]
probable_IDs = set(toplayer)

loop_id = 0
while not all(serpentinized):
    loop_id += 1 
    if loop_id > 100:
        print('Warning! loop_id > 50. Exiting while loop.')
        print('loop:',loop_id, 'ids:', probable_IDs)
        break 
        
    if len(probable_IDs) > N_select:
        selected_IDs = np.random.choice(list(probable_IDs), N_select)
    else:
        selected_IDs = list(probable_IDs)
    
    if len(selected_IDs)==0:
        print('Warning! len(selected_IDs) = 0')
        print('loop:',loop_id, 'ids:', probable_IDs)
        break
        
    print('-'*50)
    print('loop:',loop_id, 'ids:', probable_IDs)
    print('-'*50)
    
    expand_disks(IDs=selected_IDs, X=rad_increment, V=rad_increment_rate, n=nb_steps)
    probable_IDs, serpentinized = update_probable_IDs(probable_IDs, selected_IDs, in_contact, serpentinized)
    
print('='*50)
print('total loop:',loop_id, 'last ids:', probable_IDs)
print('='*50)
#'''


#
# close display & postpro
#
chipy.CloseDisplayFiles()
chipy.ClosePostproFiles()

# this is the end
chipy.Finalize()

t = time() - t0
h = int(t/3600)
m = int((t - h*3600)/60)
s = t - h*3600 - m*60
print('time needed to simulate: %ihr %im %is' %(h,m,s))

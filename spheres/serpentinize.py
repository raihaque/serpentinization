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
nb_steps = 200
radius = 1e-2
rad_incr = 0.1e-2 # 10%

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
freq_write   = 20

# display parameters
freq_display = 20

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

def expand_disks(IDs, X, n, dt):
  V = X/(n*dt)
  
  for k in range(n):
    chipy.utilities_logMes('INCREMENT STEP')
    chipy.IncrementStep()
    
    chipy.utilities_logMes('COMPUTE Fext')
    chipy.ComputeFext()
    chipy.utilities_logMes('COMPUTE Fint')
    chipy.ComputeBulk()
    chipy.utilities_logMes('COMPUTE Free Vlocy')
    chipy.ComputeFreeVelocity()
    
    for ID in IDs:
      chipy.DISKx_SetVdilation(int(ID), V)
    
    chipy.utilities_logMes('SELECT PROX TACTORS')
    chipy.SelectProxTactors() 

    chipy.utilities_logMes('RESOLUTION' )
    chipy.RecupRloc(Rloc_tol) 

    chipy.ExSolver(solver_type, norm, tol, relax, gs_it1, gs_it2)
    chipy.UpdateTactBehav()

    chipy.StockRloc() 

    chipy.utilities_logMes('COMPUTE DOF, FIELDS, etc.')
    chipy.ComputeDof()

    chipy.utilities_logMes('UPDATE DOF, FIELDS')
    chipy.UpdateStep() 
    
    dx = X*k/n
    for ID in IDs:
      chipy.DISKx_SetXdilation(int(ID), dx)

    chipy.utilities_logMes('WRITE OUT')
    chipy.WriteOut(freq_write)

    chipy.utilities_logMes('VISU & POSTPRO')
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

with open('./DATBOX/neighbors.pkl','rb') as foo:
    in_contact, toplayer = pickle.load(foo)


Nb = chipy.DISKx_GetNbDISKx() 
#print('number of disks = ', Nb)
#print('in_contact:')
#print(in_contact)
#sys.exit()
N_select = 4
serpentinized = [False for i in range(Nb)]
probable_IDs = set(toplayer)

loop_id = 0
while not all(serpentinized):
    loop_id += 1 
    if loop_id > 50:
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
    
    expand_disks(IDs=selected_IDs, X=rad_incr, n=nb_steps, dt=dt)
    probable_IDs, serpentinized = update_probable_IDs(probable_IDs, selected_IDs, in_contact, serpentinized)
    
print('='*50)
print('total loop:',loop_id, 'last ids:', probable_IDs)
print('='*50)

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

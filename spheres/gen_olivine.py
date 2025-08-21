import os
import math
import numpy as np
from pylmgc90 import pre
#import h5py
import pickle

# create folders
if not os.path.isdir('./DATBOX'):
  os.mkdir('./DATBOX')

dim = 2

nb_particles = 100

radius = 1e-2
wall_width = 1e-2
Fy = -10e9

mat = pre.material(name='PDURx', materialType='RIGID', density=3000)
mod = pre.model(name='rigid', physics='MECAx', element='Rxx2D', dimension=dim)

# generate the hexagons
bodies = pre.avatars()

Nx = round(np.sqrt(nb_particles*np.sqrt(3)/2))
Ny = round(np.sqrt(nb_particles*2/np.sqrt(3)))
print('tot particles: {} * {} = {}'.format(Nx, Ny, Nx * Ny))
print('max pos:', np.array([1+Nx*2, (1+Ny*np.sqrt(3)/2)])*radius)
positions = [] # position of spheres
for j in range(Ny):
    for i in range(Nx):
        pos = np.array([1+i*2+j%2, 1+j*np.sqrt(3)])*radius
        body = pre.rigidDisk(r=radius, center=pos, model=mod, material=mat, color='BLUEx')
        bodies.addAvatar(body)
        positions.append(pos)

lx = (1 + 2*Nx)*radius
ly = lx
top_pos = (2 + (Ny-1) * np.sqrt(3)) * radius

# define walls

mat_wall    = pre.material(name='bndry', materialType='RIGID', density=10000)

left   = pre.rigidJonc(axe1=ly, axe2=wall_width, center=[-wall_width, 0.5*ly], model=mod, material=mat_wall, color='WALLx')
left.imposeDrivenDof(component=[1,2,3], dofty='vlocy')
left.rotate(description='axis', alpha=math.pi/2., axis=[0., 0., 1.], center=[-wall_width, 0.5*ly])
bodies.addAvatar(left)


right  = pre.rigidJonc(axe1=ly, axe2=wall_width, center=[lx+wall_width, 0.5*ly], model=mod, material=mat_wall, color='WALLx')
right.imposeDrivenDof(component=[1,2,3], dofty='vlocy')
right.rotate(description='axis', alpha=-math.pi/2., axis=[0., 0., 1.], center=[lx+wall_width, 0.5*ly])
bodies.addAvatar(right)


bottom = pre.rigidJonc(axe1=lx, axe2=wall_width, center=[0.5*lx, -wall_width], model=mod, material=mat_wall, color='WALLx')
bottom.imposeDrivenDof(component=[1,2,3], dofty='vlocy')
bodies.addAvatar(bottom)

top = pre.rigidJonc(axe1=lx, axe2=wall_width, center=[0.5*lx, top_pos + wall_width], model=mod, material=mat_wall, color='WALLx')
top.imposeDrivenDof(component=[1,3], dofty='vlocy')
top.imposeDrivenDof(component=2, dofty='force',description='predefined',ct=Fy)
# evolution : giving a file containing t,f(t)
#top.imposeDrivenDof(component=2, dofty='force',description='evolution',evolutionFile='Fy.dat')

bodies.addAvatar(top)

mats = pre.materials()
mats.addMaterial(mat, mat_wall)
mods = pre.models()
mods.addModel(mod)
svs   = pre.see_tables()
tacts = pre.tact_behavs()

#'''
# ---- GLUEING BETWEEN PARTICLES WITH IQS_MOHR_DS_CLB ----
# Define 'glued' contact behavior
ld_glue = pre.tact_behav(
    name='mhr_c',
    law='IQS_MOHR_DS_CLB',
    cohn=1.5e9,      # normal cohesion (glue strength in normal direction) #si units
    coht=2.0e8,      # tangential cohesion (glue strength in tangential direction)
    stfr=0.6,       # static friction coefficient after glue break
    dyfr=0.4        # dynamic friction coefficient after glue break
)
tacts += ld_glue

# Assign this behavior to particle-particle (POLYG-POLYG) contacts
sv_glue = pre.see_table(
    CorpsCandidat='RBDY2', candidat='DISKx', colorCandidat='BLUEx',
    CorpsAntagoniste='RBDY2', antagoniste='DISKx', colorAntagoniste='BLUEx',
    behav=ld_glue, alert=0.01
)
svs += sv_glue #'''

'''
# (no glue or cohesion)
# interaction definition:
ldkdk = pre.tact_behav(name='iqsc0', law='IQS_CLB', fric=0.3)
tacts+= ldkdk

svdkdk = pre.see_table(CorpsCandidat='RBDY2', candidat='DISKx', colorCandidat='BLUEx',
                       CorpsAntagoniste='RBDY2', antagoniste='DISKx', colorAntagoniste='BLUEx',
                       behav=ldkdk, alert=.01)
svs+=svdkdk 
#'''


# Polygon - wall contact

ldkjc = pre.tact_behav(name='iqsc1', law='IQS_CLB', fric=0.4)
tacts+= ldkjc

svdkjc = pre.see_table(CorpsCandidat='RBDY2', candidat='DISKx', colorCandidat='BLUEx',
                       CorpsAntagoniste='RBDY2', antagoniste='JONCx', colorAntagoniste='WALLx',
                       behav=ldkjc, alert=0.01)
svs+=svdkjc


post = pre.postpro_commands()


# writing files
pre.writeDatbox(dim, mats, mods, bodies, tacts, svs, post=post, gravy=[0., -9.8, 0.])


#-------- store the force on the top wall --------#

dt = 1e-5
nb_steps = 1500
t1 = nb_steps * dt

def imposedFy(t):
    return Fy

pre.writeEvolution(f=imposedFy, instants=np.linspace(0., t1, nb_steps) ,path='DATBOX/', name='Fy.dat')

#try:
#  pre.visuAvatars(bodies)
#except:
#  pass


# find neighbours and store the information
N = len(positions)
in_contact = {i:[] for i in range(1,N+1)}
for i in range(1, N+1):
    for j in range(i+1, N+1):
        dsq = sum((positions[i-1] - positions[j-1])**2)
        if dsq <= 4.04*radius**2: # epsilon = 1% of radius
            in_contact[i].append(j)
            in_contact[j].append(i)

top_layer = [Nx*(Ny-1) + i for i in range(1,Nx+1)]

#with h5py.File('./DATBOX/neighbors.h5', 'w') as f:
#    #dset = f.create_dataset('in_contact', data=in_contact)
#    for k, v in in_contact.items():
#        f.create_dataset(k, data=v)


with open('./DATBOX/neighbors.pkl', 'wb') as foo:
    pickle.dump((in_contact, top_layer), foo)
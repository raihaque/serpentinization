import os
import math
import numpy as np
from pylmgc90 import pre
import json
import gmsh

def getmesh(center, nb_vertices, radius, phi=0) :
  
  # mesh dimension
  lcmin = radius/3
  lcmax = radius/3

  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 1)
  gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lcmin)
  gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lcmax)
  gmsh.option.setNumber("Mesh.Algorithm", 2)

  gmsh.model.add('Grain')
  gmsh.model.setCurrent('Grain')

  # points of the fill contour
  lp=[]
  dtt = 2*math.pi/nb_vertices
  for i in range(nb_vertices):
    lp.append(gmsh.model.geo.addPoint(center[0]+radius*math.cos(phi + i*dtt), center[1]+radius*math.sin(phi + i*dtt), 0, lcmax))

  # lines of the contour
  lc=[]            
  for i in range(nb_vertices):
     lc.append(gmsh.model.geo.addLine(lp[i], lp[(i+1)%nb_vertices]))

  # surface creation
  contour=gmsh.model.geo.addCurveLoop(lc)

  surface=gmsh.model.geo.addPlaneSurface([contour])

  gmsh.model.geo.synchronize()

  # names of the components, limits and contacts
  fill = gmsh.model.addPhysicalGroup(2,[surface])
  gmsh.model.setPhysicalName(2, fill, 'fill')

  side = gmsh.model.addPhysicalGroup(1,lc)
  gmsh.model.setPhysicalName(1, side, 'side')

  # mesh generation and export
  gmsh.model.mesh.generate(2)
  gmsh.write('fill.msh')
  gmsh.finalize()

  # fill mesh reading
  return pre.readMesh('fill.msh', dim)
  
# create folders
if not os.path.isdir('./DATBOX'):
  os.mkdir('./DATBOX')

filename = 'hexagons_data_grouped.json'
with open(filename, 'r') as foo:
    data = json.load(foo)

radius = data['parameters']['radius']
lx = data['parameters']['Lx']
ly = lx
top_pos = data['parameters']['Ly']
phi = data['parameters']['phi']
hexagon_data =  data['hexagons']
groups = data['groups']
num_grs = data['parameters']['num_grs']

if len(groups) > 100000:
    print('Warning! Number of groups is more than 100000.')
    sys.exit()

dim = 2

#radius     = 1e-1   # particle size ## read from json file
wall_width = radius  #
Fy         =-1e6     # Force on the top wall
density_p  = 3000.   # particle density
density_w  = 10000.  # wall density
g = -9.8             # gravity
alert = radius       # used for contact search in pre.see_table
young_modulus = 10e7 #

cn = 2e8             # Normal stiffness (bonds between same group)
ct = 1e8             # Tangential stiffness (bonds between same group)
smax = 1e7           # Maximum cohesive stress (bonds between same group)
w = 6e-6             # Critical separation (bonds between same group)

cn2 = 2e7            # Normal stiffness (bonds between different groups)
ct2 = 1e7            # Tangential stiffness (bonds between different groups)
smax2 = 1e7          # Maximum cohesive stress (bonds between different groups)
w2 = 6e-6            # Critical separation (bonds between different groups)
stfr = 0.7           # static friction 
dyfr = 0.6           # dynamic friction coefficient after glue break


mats = pre.materials()

# define walls
mat_wall = pre.material(name='bndry', materialType='RIGID', density=density_w)

matD = pre.material(name='stone', materialType='ELAS_DILA', elas='standard',
                    young=young_modulus, nu=0.3, anisotropy='isotropic',
                    density=density_p,dilatation=1e-5,T_ref_meca=20.)

mats.addMaterial(matD,mat_wall)

mods = pre.models()

modR = pre.model(name='rigid', physics='MECAx', element='Rxx2D', dimension=dim)


modD = pre.model(name='defox', physics='MECAx',
                 element='T3xxx', dimension=dim,
                 external_model='MatL_', kinematic='small', material='elasd', 
                 anisotropy='iso__', mass_storage='coher',external_fields=['TEMPERATURE'])

mods.addModel(modR, modD)

# generate the hexagons
bodies = pre.avatars()

for hexagon in hexagon_data:
    pos = hexagon['center']
    mesh = getmesh(center=pos, nb_vertices=6, radius=radius, phi=phi)
    color = '{}'.format(hexagon['gr_id']).rjust(5,'x')
    body = pre.buildMeshedAvatar(mesh=mesh, model=modD, material=matD)
    body.addContactors(group='side', shape='ALpxx', color=color, reverse=True)
    body.addContactors(group='side', shape='CLxxx', color=color, reverse=True, weights=[0.01,0.99])                
                                 
    bodies.addAvatar(body)

left   = pre.rigidJonc(axe1=ly, axe2=wall_width, center=[-wall_width, 0.5*ly], model=modR, material=mat_wall, color='WALLx')
left.imposeDrivenDof(component=[1,2,3], dofty='vlocy')
left.rotate(description='axis', alpha=math.pi/2., axis=[0., 0., 1.], center=[-wall_width, 0.5*ly])
bodies.addAvatar(left)

#'''
right  = pre.rigidJonc(axe1=ly, axe2=wall_width, center=[lx+wall_width, 0.5*ly], model=modR, material=mat_wall, color='WALLx')
right.addContactors(shape='PT2Dx', color='REDxx', shift=[0., 0.])
right.imposeDrivenDof(component=[2,3], dofty='vlocy')
#right.imposeDrivenDof(component=[1,2,3], dofty='vlocy')
right.rotate(description='axis', alpha=-math.pi/2., axis=[0., 0., 1.], center=[lx+wall_width, 0.5*ly])
bodies.addAvatar(right)

right2=pre.avatar( dimension=dim)
right2.addBulk( pre.rigid2d() )
right2.addNode( pre.node( coor=np.array([lx+3*wall_width, 0.5*ly]),number=1) )
right2.defineGroups()
right2.defineModel(model=modR)
right2.defineMaterial(material=mat_wall)
right2.addContactors(shape='PT2Dx', color='REDxx', shift=[0., 0.])
right2.imposeDrivenDof(component=[1,2,3], dofty='vlocy')
bodies.addAvatar(right2)
#'''
#
'''
right   = pre.rigidJonc(axe1=ly, axe2=wall_width, center=[lx+wall_width, 0.5*ly], model=modR, material=mat_wall, color='WALLx')
right.imposeDrivenDof(component=[1,2,3], dofty='vlocy')
right.rotate(description='axis', alpha=math.pi/2., axis=[0., 0., 1.], center=[lx+wall_width, 0.5*ly])
bodies.addAvatar(right)'''

bottom = pre.rigidJonc(axe1=lx, axe2=wall_width, center=[0.5*lx, -wall_width], model=modR, material=mat_wall, color='WALLx')
bottom.imposeDrivenDof(component=[1,2,3], dofty='vlocy')
bodies.addAvatar(bottom)

top = pre.rigidJonc(axe1=lx, axe2=wall_width, center=[0.5*lx, top_pos + wall_width], model=modR, material=mat_wall, color='WALLx')
top.imposeDrivenDof(component=[1,3], dofty='vlocy')
#top.imposeDrivenDof(component=2, dofty='force', description='predefined', ct=Fy)
# evolution : giving a file containing t,f(t)
top.imposeDrivenDof(component=2, dofty='force',description='evolution',evolutionFile='Fy.dat')

bodies.addAvatar(top)

svs   = pre.see_tables()
tacts = pre.tact_behavs()

# interaction definition:
colors = ['{}'.format(i).rjust(5,'x') for i in range(len(groups))]

#lgg = pre.tact_behav(name='iqsc0', law='GAP_SGR_CLB', fric=dyfr)

# contact law between same group
lgg = pre.tact_behav(name='mpczm', law='MP3_CZM', cn=cn, ct=ct, smax=smax, w=w, stfr=stfr, dyfr=dyfr)
tacts+= lgg

for color in colors:
    svgg = pre.see_table(CorpsCandidat='MAILx', candidat='CLxxx', colorCandidat=color,
                         CorpsAntagoniste='MAILx', antagoniste='ALpxx', colorAntagoniste=color,
                         behav=lgg, alert=alert,halo=0.1)
    svs+=svgg 

# contact law between different groups
lgg2 = pre.tact_behav(name='mpczm', law='MP3_CZM', cn=cn2, ct=ct2, smax=smax2, w=w2, stfr=stfr, dyfr=dyfr)
tacts+= lgg2

for i in range(len(colors)):
    for j in range(i+1,len(colors)):
        svgg = pre.see_table(CorpsCandidat='MAILx', candidat='CLxxx', colorCandidat=colors[i],
                             CorpsAntagoniste='MAILx', antagoniste='ALpxx', colorAntagoniste=colors[j],
                             behav=lgg2, alert=alert,halo=0.1)
        svs+=svgg


# wall-hexagon contact
lgp = pre.tact_behav(name='iqsc1', law='GAP_SGR_CLB', fric=dyfr)
tacts+= lgp

for color in colors:
    svgp = pre.see_table(CorpsCandidat='MAILx'   , candidat='CLxxx'   , colorCandidat=color,
                        CorpsAntagoniste='RBDY2', antagoniste='JONCx', colorAntagoniste='WALLx',
                        behav=lgp,  alert=alert)
    svs+=svgp

lptpt = pre.tact_behav(name='stiff', law='ELASTIC_ROD', stiffness=1e10, prestrain=0.)
tacts+= lptpt

svptpt = pre.see_table(CorpsCandidat='RBDY2', candidat='PT2Dx', colorCandidat='REDxx',
                       CorpsAntagoniste='RBDY2', antagoniste='PT2Dx', colorAntagoniste='REDxx',
                       behav=lptpt, alert=4*wall_width)
svs+=svptpt


post = pre.postpro_commands()


# writing files
pre.writeDatbox(dim, mats, mods, bodies, tacts, svs, post=post, gravy=[0., g, 0.])

print('Total groups:', len(colors))

try:
 pre.visuAvatars(bodies)
except:
 pass


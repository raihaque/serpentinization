import numpy as np 
import matplotlib.pyplot as plt


wall_pos = np.loadtxt('POSTPRO/top_wall_position.txt')
reaction_data = np.loadtxt('POSTPRO/reaction_data.txt')
reaction_data[:,1:] = reaction_data[:,1:]*100

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.plot(wall_pos[:,0], wall_pos[:,1], 'k-', label='wall position')

ax2.plot(reaction_data[:,0], reaction_data[:,1], '-o', label='nodes reacted')
ax2.plot(reaction_data[:,0], reaction_data[:,2], '-o', label='length reacted')
ax2.plot(reaction_data[:,0], reaction_data[:,3], '-o', label='porosity')


ax1.set_xlabel('time (s)')
ax1.set_ylabel('wall position (m)')
ax2.set_ylabel('porosity, nodes, length reacted (%)')

ax1.legend(loc='upper left')
ax2.legend(loc='lower right')
plt.savefig('wall_pos_with_time.png', bbox_inches='tight', dpi=300)
plt.show()
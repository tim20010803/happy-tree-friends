import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio
import matplotlib.font_manager as font_manager

#####overlay version of plot###########

font_size = 10

def set_overlay_figure_index( ):

    plt.axis('scaled')
    # plt.axis('equal')

    plt.xlabel( 'x axis', fontsize=font_size, labelpad=5 )
    plt.ylabel( 'y axis', fontsize=font_size, labelpad=5, rotation='vertical' )

    plt.tick_params(axis='x', labelsize=8)
    plt.tick_params(axis='y', labelsize=8)

    plt.xlim( -1 * x_range, x_range  )
    plt.ylim( -1 * x_range, x_range )

# input the data from files
data = pd.read_csv('./mainTree_data.csv')
data_non = pd.read_csv('./mainNonTree_data.csv')

times = data['Time'].unique()
times = data_non['Time'].unique()
E_total = data['SystemEnergy'][::len(data['Particle'].unique())]
Px_total = data['SystemMomentumX'][::len(data['Particle'].unique())]
Py_total = data['SystemMomentumY'][::len(data['Particle'].unique())]
Lx_total = data['SystemAngularMomentumX'][::len(data['Particle'].unique())]
Ly_total = data['SystemAngularMomentumY'][::len(data['Particle'].unique())]
P_total0 = (Px_total[0]**2 + Py_total[0]**2)**0.5
L_total0 = (Lx_total[0]**2 + Ly_total[0]**2)**0.5
E_total_non = data_non['SystemEnergy'][::len(data_non['Particle'].unique())]
Px_total_non = data_non['SystemMomentumX'][::len(data_non['Particle'].unique())]
Py_total_non = data_non['SystemMomentumY'][::len(data_non['Particle'].unique())]
Lx_total_non = data_non['SystemAngularMomentumX'][::len(data_non['Particle'].unique())]
Ly_total_non = data_non['SystemAngularMomentumY'][::len(data_non['Particle'].unique())]
P_total_non0 = (Px_total_non[0]**2 + Py_total_non[0]**2)**0.5
L_total_non0 = (Lx_total_non[0]**2 + Ly_total_non[0]**2)**0.5


# plotting parameters
nstep_per_image = 1           # plotting frequency
time_gap = 1 # the period per frame

def update_overlay( frame ):
    
    index = frame * len(data_non['Particle'].unique()) # the index of the momentum and the angular momentum

    frame_data = data[data['Time'] == times[frame]] # pick up the data at the same time
    positions_x = frame_data['PositionX'].tolist()
    positions_y = frame_data['PositionY'].tolist()

    frame_data_non = data_non[data_non['Time'] == times[frame]]
    positions_x_non = frame_data_non['PositionX'].tolist()
    positions_y_non = frame_data_non['PositionY'].tolist()

  
    # delete the old data on the figure
    plt.cla()
    plt.cla()

    set_overlay_figure_index()

    plt.scatter(positions_x, positions_y, c='red', s=3, marker='o', alpha=1, label='Tree', edgecolor='none')
    plt.scatter(positions_x_non, positions_y_non, c='#7FFF00', s=3, marker='o', alpha=0.7,label='NonTree', edgecolor='none')

    legend = plt.legend()

    handles = legend.legendHandles

    for handle in handles:
        handle.set_sizes([30])

    fig.suptitle('time = %.5f' % ( times[frame] ), c='blue')


# set the particle number
num_particles = len(data['Particle'].unique())

# set the range of the figure 
x_range = 3000.0

# create figure
fig = plt.figure( figsize=(6, 6), dpi=200 )
ax = plt.gca()
ax.set_facecolor('#D3D3D3') # #F0F0F0 (slightly lighter gray)
fig.subplots_adjust(left=0.15, right=0.9, bottom=0.1, top=0.95)

set_overlay_figure_index()

# create movie
nframe = len(data_non['Time'].unique()) # arbitrarily large
# print(nframe)
# print(nframe//time_gap)
# print("Tree data length:",len(data['Time']))
# print("NonTree data length:",len(data_non['Time']))
ani   = animation.FuncAnimation( fig, func=update_overlay, frames=nframe//time_gap, interval=1, repeat=False )

ani.save('ptc_pos_overlay-3000.gif', writer='pillow')
plt.show()

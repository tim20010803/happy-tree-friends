import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio
import matplotlib.font_manager as font_manager

font_size = 10

def set_figure_index():
    ax[0].set_aspect('equal')
    ax[1].set_aspect('equal')

    ax[0].set_xlabel( 'x axis', labelpad=5 )
    ax[0].set_ylabel( 'y axis', labelpad=5, rotation='vertical' )
    ax[1].set_xlabel( 'x axis', labelpad=5 )
    ax[1].set_ylabel( 'y axis', labelpad=5, rotation=-90 )

    ax[0].xaxis.label.set_fontsize(font_size)
    ax[0].yaxis.label.set_fontsize(font_size)
    ax[1].xaxis.label.set_fontsize(font_size)
    ax[1].yaxis.label.set_fontsize(font_size)
    ax[1].yaxis.set_label_position('right')

    ax[0].set_xlim( -1 * x_range, x_range )
    ax[0].set_ylim( -1 * x_range, x_range )
    ax[1].set_xlim( -1 * x_range, x_range )
    ax[1].set_ylim( -1 * x_range, x_range )
    ax[1].yaxis.tick_right()

    ax[0].tick_params(axis='both', labelsize=5)
    ax[1].tick_params(axis='both', labelsize=5)

    ax[0].set_title('Particle Animation\nwith Tree algorithm', y=1.2, fontweight='bold')
    ax[1].set_title('Particle Animation\nwithout Tree algorithm', y=1.2, fontweight='bold')



# input the data from files
data = pd.read_csv('mainTree_data.csv')
data_non = pd.read_csv('mainNonTree_data.csv')

#times = data['Time'].unique()
times = data_non['Time'].unique()
#E_total = data['SystemEnergy'][::len(data['Particle'].unique())]
#Px_total = data['SystemMomentumX'][::len(data['Particle'].unique())]
#Py_total = data['SystemMomentumY'][::len(data['Particle'].unique())]
#Lx_total = data['SystemAngularMomentumX'][::len(data['Particle'].unique())]
#Ly_total = data['SystemAngularMomentumY'][::len(data['Particle'].unique())]
#P_total0 = (Px_total[0]**2 + Py_total[0]**2)**0.5
#L_total0 = (Lx_total[0]**2 + Ly_total[0]**2)**0.5
E_total_non = data_non['SystemEnergy'][::len(data_non['Particle'].unique())]
Px_total_non = data_non['SystemMomentumX'][::len(data_non['Particle'].unique())]
Py_total_non = data_non['SystemMomentumY'][::len(data_non['Particle'].unique())]
Lx_total_non = data_non['SystemAngularMomentumX'][::len(data_non['Particle'].unique())]
Ly_total_non = data_non['SystemAngularMomentumY'][::len(data_non['Particle'].unique())]
P_total_non0 = (Px_total_non[0]**2 + Py_total_non[0]**2)**0.5
L_total_non0 = (Lx_total_non[0]**2 + Ly_total_non[0]**2)**0.5


# plotting parameters
nstep_per_image = 1           # plotting frequency
time_gap = 100 # the period per frame

def update( frame ):

    index = frame * len(data_non['Particle'].unique()) # the index of the momentum and the angular momentum

    frame_data = data[data['Time'] == times[frame]] # pick up the data at the same time
    positions_x = frame_data['PositionX'].tolist()
    positions_y = frame_data['PositionY'].tolist()

    frame_data_non = data_non[data_non['Time'] == times[frame]]
    positions_x_non = frame_data_non['PositionX'].tolist()
    positions_y_non = frame_data_non['PositionY'].tolist()

    
    #P_total = ((Px_total[index])**2 + (Py_total[index])**2)**0.5
    #L_total = ((Lx_total[index])**2 + (Ly_total[index])**2)**0.5
    #E_error = ( E_total[index] - E_total[0] ) / E_total[0]
    #P_error = ( P_total - P_total0 ) / P_total0
    #L_error = ( L_total - L_total0 ) / L_total0

    P_total_non = ((Px_total_non[index])**2 + (Py_total_non[index])**2)**0.5
    L_total_non = ((Lx_total_non[index])**2 + (Ly_total_non[index])**2)**0.5
    E_error_non = ( E_total_non[index] - E_total_non[0] ) / E_total_non[0]
    P_error_non = ( P_total_non - P_total_non0 ) / P_total_non0
    L_error_non = ( L_total_non - L_total_non0 ) / L_total_non0
    
    # delete the old data on the figure
    #ax[0].cla()
    ax[1].cla()

    #ax[0].scatter(positions_x, positions_y, c='green', s=0.5)
    ax[1].scatter(positions_x_non, positions_y_non, c='blue', s=0.5)

    set_figure_index()

    fig.suptitle('time = %.5f' % ( times[frame] ), c='blue')
    
    # (x)ax[0].text( 1000, 1200, 'the error of total energy= %10.3e\nthe error of total momentum = %10.3e\nthe error of total angular momentum = %10.3e' % (E_error, P_error, L_error), fontsize=8, color='black', ha='right',horizontalalignment='center', verticalalignment='center' )
    # (x)ax[1].text( 3750, 3700, 'the error of total energy= %10.3e\nthe error of total momentum = %10.3e\nthe error of total angular momentum = %10.3e' % (E_error_non, P_error_non, L_error_non), fontsize=8, color='black', ha='right', horizontalalignment='center', verticalalignment='center' )
    
    #fig.text( 0.47, 0.78, 'the error of total energy= %10.3e\nthe error of total momentum = %10.3e\nthe error of total angular momentum = %10.3e' % (E_error, P_error, L_error), fontsize=8, color='black', ha='right', va='center')
    fig.text( 0.95, 0.78, 'the error of total energy= %10.3e\nthe error of total momentum = %10.3e\nthe error of total angular momentum = %10.3e' % (E_error_non, P_error_non, L_error_non), fontsize=8, color='black', ha='right', va='center')

# set the particle number
num_particles = len(data['Particle'].unique())

# set the range of the figure 
x_range = 3000.0

# create figure
fig, ax = plt.subplots( 1, 2, sharex=False, sharey=False, dpi=200 )
#fig.subplots_adjust( hspace=0.0, wspace=0.5 )
fig.subplots_adjust(left=0.12, right=0.9, bottom=0.05, top=0.95, wspace=0.3)

set_figure_index()
'''
random_colors = np.random.randint(0, 256, size=(num_particles, 3))
random_colors = random_colors.astype(float)
random_colors /= 255.0
'''
# create movie
nframe = len(data_non['Time'].unique()) # arbitrarily large
ani   = animation.FuncAnimation( fig, func=update, frames=nframe//time_gap, interval=200, repeat=False )

plt.show()

#ani.save('ptc_pos.gif', writer='pillow')



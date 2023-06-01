import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio

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

def update( frame ):

    index = frame * len(data_non['Particle'].unique())
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

    ax[0].set_aspect('equal')
    ax[1].set_aspect('equal')

    ax[0].set_xlabel( 'x axis' )
    ax[0].set_ylabel( 'y axis' )
    ax[1].set_xlabel( 'x axis' )
    ax[1].set_ylabel( 'y axis' )

    ax[0].set_xlim( x_min, x_MAX )
    ax[0].set_ylim( x_min, x_MAX )
    ax[1].set_xlim( x_min, x_MAX )
    ax[1].set_ylim( x_min, x_MAX )

    fig.suptitle('time = %d' % ( times[frame] ))
    ax[0].set_title('Particle Animation\nwith Tree algorithm')
    #text0.set(text='the error of total energy= %10.3e\nthe error of total momentum = %10.3e\nthe error of total angular momentum = %10.3e' % (E_error, P_error, L_error))
    ax[1].set_title('Particle Animation\nwithout Tree algorithm')
    #text1.set(text='the error of total energy= %10.3e\nthe error of total momentum = %10.3e\nthe error of total angular momentum = %10.3e' % (E_error_non, P_error_non, L_error_non))
    #print(frame)
    

# set the particle number
num_particles = len(data['Particle'].unique())

# set the range of the figure 
x_min = -1000.0
x_MAX = 1000.0

# create figure
fig, ax = plt.subplots( 1, 2, sharex=False, sharey=False, dpi=200 )
fig.subplots_adjust( hspace=0.0, wspace=0.5 )

text0  = ax[1].text( 0.0, 1.3, '', fontsize=8, color='black',
                 ha='center', va='center' )
text1  = ax[1].text( 0.0, 1.3, '', fontsize=8, color='black',
                 ha='center', va='center' )

ax[0].set_aspect('equal')
ax[1].set_aspect('equal')

ax[0].set_xlabel( 'x axis' )
ax[0].set_ylabel( 'y axis' )
ax[1].set_xlabel( 'x axis' )
ax[1].set_ylabel( 'y axis' )

ax[0].set_xlim( x_min, x_MAX )
ax[0].set_ylim( x_min, x_MAX )
ax[1].set_xlim( x_min, x_MAX )
ax[1].set_ylim( x_min, x_MAX )

ax[0].set_title('Particle Animation\nwith Tree algorithm')
ax[1].set_title('Particle Animation\nwithout Tree algorithm')


random_colors = np.random.randint(0, 256, size=(num_particles, 3))
random_colors = random_colors.astype(float)
random_colors /= 255.0

# create movie
nframe = len(data_non['Time'].unique()) # arbitrarily large
ani   = animation.FuncAnimation( fig, func=update, frames=nframe, interval=200, repeat=False )

plt.show()

#ani.save('ptc_pos.gif', writer='pillow')



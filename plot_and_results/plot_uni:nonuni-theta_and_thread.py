import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

def input_data( mode ):
    if(mode == 0):
        data = pd.read_csv('./RuntimeTreeComplex_uni2.csv')
    else:
        data = pd.read_csv('./RuntimeTreeComplex_Nonuni2.csv')
    return data

def plot_thread_perf( mode ):
    data = input_data( mode )

    thread_num = np.array(data['NThread'].unique())
    one_thread = np.array(data[data['NThread'] == 1]['total_time(s)']) # the standard for thread compare
    particleNumber = data['particleNumber'].unique().tolist()

    ax        = plt.axes( xlim=(0,9), ylim=(0.9, 7.0) )
    ax.set_xlabel( 'thread numbers' )
    ax.set_ylabel( 'speed up (times)' )
    if (mode == 0):
        ax.set_title('Performance of multi-thread in uniform particle numbers')
    else:
        ax.set_title('Performance of multi-thread in nonuniform particle numbers')
    plt.xticks(np.arange(0, 10, 1))

    j = 0
    for i in particleNumber:

        time_origin = np.array(data[(data['particleNumber'] == i) & (data['theta'] == 1.0)]['total_time(s)'])
        time = np.delete(time_origin, -1)
        
        speed_up = one_thread[j] / time
        ax.plot( thread_num, speed_up, ls='-',  label='1e%d particles' % (np.log10(i)) , marker='o')
        j += 1

    ax.legend(loc='best')
    if(mode == 0):
        plt.savefig('uni-thread-perf2.png')
    else:
        plt.savefig('nonuni-thread-perf2.png')


def plot_theta_perf( mode ):
    data = input_data( mode )

    theta_num = np.sort(np.array(data['theta'].unique()))
    eight_thread = np.array(data[(data['NThread'] == 8) & (data['theta'] == 1.0)]['total_time(s)']) # the standard for theta compare
    particleNumber = data['particleNumber'].unique().tolist()

    ax        = plt.axes( xlim=(0.0 ,1.4), ylim=(0.0, 1.5) )
    ax.set_xlabel( 'theta numbers' )
    ax.set_ylabel( 'speed-up (times)' )
    if (mode == 0):
        ax.set_title('Performance of different theta in uniform particle numbers')
    else:
        ax.set_title('Performance of different theta in nonuniform particle numbers')
    fig.text(0.7, 0.15, '(thread = 8, theta = 1.0)', c='#444444', ha='center', va='center')
    # plt.xticks(np.arange(0, 10, 1))

    j = 1
    for i in particleNumber:

        time_origin = np.array(data[(data['particleNumber'] == i) & (data['NThread'] == 8.0)]['total_time(s)'])
        time = np.delete(time_origin, 0)
        
        speed_up = eight_thread[j] / time
        ax.plot( theta_num, speed_up, ls='-',  label='1e%d particles' % (np.log10(i)) , marker='o')
        j += 2

    ax.legend(loc='best')
    if(mode == 0):
        plt.savefig('uni-theta-perf2.png')
    else:
        plt.savefig('nonuni-theta-perf2.png')

def plot_thread_eff( mode ):
    data = input_data( mode )

    thread_num = np.array(data['NThread'].unique())
    theta_num = np.array(data['theta'].unique())

    one_thread = np.array(data[data['NThread'] == 1]['total_time(s)']) # the standard for thread compare
    eight_thread = np.array(data[(data['NThread'] == 8) & (data['theta'] == 1.0)]['total_time(s)']) # the standard for theta compare
    particleNumber = np.array(data['particleNumber'].unique())

    ax        = plt.axes( xlim=(2.5 ,6.5), ylim=(0.2, 1.2) )
    ax.set_xlabel( 'particle numbers (1og10)' )
    ax.set_ylabel( 'parallel efficiency' )
    if(mode == 0):
        ax.set_title('Parallel efficiency of multi-thread in uniform particle numbers')
    else:
        ax.set_title('Parallel efficiency of multi-thread in nonuniform particle numbers')
    plt.xticks(np.arange(3, 7, 1))
  
    for i in thread_num:
        time = np.array(data[(data['NThread'] == i) & (data['theta'] == 1.0)]['total_time(s)'])
        speed_up = np.empty(len(particleNumber))
        if (i == 8):
            time = time[::2]
        for j in range (0, len(particleNumber)):
            speed_up[j] = one_thread[j] / time[j]
            # print("%f / %f = %f" % (one_thread[j], time[j], speed_up[j]))
        eff = speed_up / i
        # print(eff)
        ax.plot( np.log10(particleNumber), eff, ls='-',  label='%d thread' % (i) , marker='o')

    ax.legend(loc='lower center', ncol=3)
    if(mode == 0):
        plt.savefig('uni-thread-eff.png')
    else:
        plt.savefig('nonuni-thread-eff.png')

def plot_theta_eff( mode ):
    data = input_data( mode )

    thread_num = np.array(data['NThread'].unique())
    theta_num = np.sort(np.array(data['theta'].unique()))

    one_thread = np.array(data[data['NThread'] == 1]['total_time(s)']) # the standard for thread compare
    eight_thread = np.array(data[(data['NThread'] == 8) & (data['theta'] == 1.0)]['total_time(s)']) # the standard for theta compare
    particleNumber = np.array(data['particleNumber'].unique())

    ax        = plt.axes( xlim=(2.5 ,6.5), ylim=(0.0, 1.0) )
    ax.set_xlabel( 'particle numbers (log10)' )
    ax.set_ylabel( 'parallel efficiency' )
    if(mode == 0):
        ax.set_title('Parallel efficiency of different theta in uniform particle numbers')
    else:
        ax.set_title('Parallel efficiency of different theta in nonuniform particle numbers')
    fig.text(0.8, 0.026, '(thread = 8, theta = 1.0)', c='#444444', ha='center', va='center')

    plt.xticks(np.arange(3, 7, 1))
  
    for i in theta_num:
        time = np.array(data[(data['NThread'] == 8) & (data['theta'] == i)]['total_time(s)'])
        if (i == 1):
            time = time[1::2]
        speed_up = np.empty(len(particleNumber))
        # print(time)
        if (i == 8):
            time = time[::2]
        for j in range (0, len(particleNumber)):
            speed_up[j] = one_thread[j] / time[j]
            # print("%f / %f = %f" % (one_thread[j], time[j], speed_up[j]))
        eff = speed_up / 8
        # print(theta_num)
        ax.plot(np.log10(particleNumber), eff, ls='-',  label='theta=%.1f' % (i) , marker='o')

    ax.legend(loc='upper right', ncol=2)
    if(mode == 0):
        plt.savefig('uni-theta-eff.png')
    else:
        plt.savefig('nonuni-theta-eff.png')

# create figure
fig       = plt.figure( figsize=(6,6), dpi=140 )

# plot_thread_perf(1)
# plot_theta_perf(1)
# plot_thread_eff(0)
plot_theta_eff(1)

# for i in range (0, 2):
#     plot_thread_perf(1)
#     plot_theta_perf(1)
#     plot_thread_eff(1)
#     plot_theta_eff(1)

plt.show()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('./RuntimeTree,theta=1,thread=8(2).csv')
data_non = pd.read_csv('./RuntimeNonTree.csv')

type = data['Type'].unique() # treeOMPNonUni(5) / treeOMPUni(5)

particleNumber = np.array(data['particleNumber'].unique()) # 5 points: 1e3, 1e4, 1e5, 1e6, 1e7
time_NonUni = np.array(data[data['Type'] == 'NonUni']['total_time(s)'])
time_Uni = np.array(data[data['Type'] == 'Uni']['total_time(s)'])

cuda_data = data_non[data_non['type'] == 'CUDA'] # particleNumber and time(s) of cuda

particleNumber_cuda = np.array(cuda_data['particleNumber'])# 5 points: 1e3, 1e4, 1e5, 1e6, 1e7
time_cuda = np.array(cuda_data['time(s)']) # time of cuda

nonOMP_data = data_non[data_non['type'] == 'NonTreeOMP'] # particleNumber and time(s) of nonOMP
particleNumber_nonOMP = np.array(nonOMP_data['particleNumber']) # 3 points: 1e5, 1e6, 1e7
time_nonOMP = np.array(nonOMP_data['time(s)']) # time of nonOMP

non_data = data_non[data_non['type'] == 'NonTree'] # particleNumber and time(s) of nonOMP
particleNumber_non = np.array(non_data['particleNumber']) # 3 points: 1e5, 1e6, 1e7
time_non = np.array(non_data['time(s)']) # time of nonOMP

# print(particleNumber_cuda)

prtN_log = np.empty(len(particleNumber))
order2 = np.empty(len(particleNumber))
x_data = np.logspace(3, 9, num=7, base=10.0)
NlogN = np.empty(len(x_data))
NlogN2 = np.empty(len(particleNumber))


for i in range ( 0, len(particleNumber) ):
    prtN_log[i] = np.log10(particleNumber[i])
    order2[i] = 2 * prtN_log[i] - 7.7
for i in range (0, len(x_data)):
    NlogN[i] = x_data[i] * np.log10(x_data[i])
    # NlogN2[i] = a * particleNumber[i] * np.log10(particleNumber[i]) + b
    # print("NlogN[%d] = %f; log(NlogN[%d]) = %f" % (i, NlogN[i], i, np.log10(NlogN[i])))

a = - 5.5

# create figure
fig       = plt.figure( figsize=(6,6), dpi=140 )
ax        = plt.axes( xlim=(3,9), ylim=(-3, 5) )
line_o2 = ax.plot( np.log10(particleNumber), order2, 'gray', ls='--',  label='$N^2$' )
line_NlogN = ax.plot( np.log10(x_data), np.log10(NlogN) + a, 'gray', linestyle='dotted', label='$NlogN$' )
ine_treeOMPNonUni, = ax.plot( prtN_log, np.log10(time_NonUni), 'r', ls='-', marker='o',  label='treeOMPNonUni' )
line_treeOMPUni, = ax.plot( prtN_log, np.log10(time_Uni), 'b', ls='-', marker='o',  label='treeOMPUni' )
line_CUDA, = ax.plot( np.log10(particleNumber_cuda), np.log10(time_cuda), 'green', ls='-', marker='o', label='CUDA' )
line_NonTreeOMP, = ax.plot( np.log10(particleNumber_nonOMP), np.log10(time_nonOMP), 'purple', ls='-', marker='o',  label='NonTreeOMP' )
line_NonTree, = ax.plot( np.log10(particleNumber_non), np.log10(time_non), 'orange', ls='-', marker='o',  label='NonTree' )
# line_NlogN2 = ax.plot( np.log10(particleNumber), np.log10(NlogN2) + a, 'brown', ls='--',  label='NlogN_type2' )
ax.set_xlabel( 'particle numbers (log)' )
ax.set_ylabel( 'time (log)' )
ax.set_title('Time complexity comparison')
ax.tick_params( top=False, right=True, labeltop=False, labelright=True )
ax.legend(loc='best')

plt.savefig('CodeType-Time.png')
plt.show()



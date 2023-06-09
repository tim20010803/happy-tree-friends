import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

data = pd.read_csv('./RuntimeTree,theta=1,thread=8(2).csv')
data_nonuni = pd.read_csv('./RuntimeTreeComplex_Nonuni2.csv')
data_uni = pd.read_csv('./RuntimeTreeComplex_uni2.csv')


particleNumber = np.array(data['particleNumber'].unique()) # 5 points: 1e3, 1e4, 1e5, 1e6, 1e7
Ctime_NonUni = np.array(data[data['Type'] == 'NonUni']['construct_time(s)'])
Ftime_NonUni = np.array(data[data['Type'] == 'NonUni']['force_time(s)'])
Ctime_Uni = np.array(data[data['Type'] == 'Uni']['construct_time(s)'])
Ftime_Uni = np.array(data[data['Type'] == 'Uni']['force_time(s)'])

#print(Ctime_NonUni)

Ctime_NonUni_pro = Ctime_NonUni / (Ctime_NonUni + Ftime_NonUni)
Ftime_NonUni_pro = Ftime_NonUni / (Ctime_NonUni + Ftime_NonUni)
Ctime_Uni_pro = Ctime_Uni / (Ctime_Uni + Ftime_Uni)
Ftime_Uni_pro = Ftime_Uni / (Ctime_Uni + Ftime_Uni)

NThread = np.array(data_nonuni['NThread'].unique())
Ctime_NonUni_thread = np.array(data_nonuni[(data_nonuni['particleNumber'] == 1000000) & (data_nonuni['theta'] == 1)]['construct_time(s)'])
Ftime_NonUni_thread = np.array(data_nonuni[(data_nonuni['particleNumber'] == 1000000) & (data_nonuni['theta'] == 1)]['force_time(s)'])
Ctime_Uni_thread = np.array(data_uni[(data_uni['particleNumber'] == 1000000) & (data_uni['theta'] == 1)]['construct_time(s)'])
Ftime_Uni_thread = np.array(data_uni[(data_uni['particleNumber'] == 1000000) & (data_uni['theta'] == 1)]['force_time(s)'])

Ctime_NonUni_thread = np.delete(Ctime_NonUni_thread, -1)
Ftime_NonUni_thread = np.delete(Ftime_NonUni_thread, -1)
Ctime_Uni_thread = np.delete(Ctime_Uni_thread, -1)
Ftime_Uni_thread = np.delete(Ftime_Uni_thread, -1)

Ctime_NonUni_thread_pro = Ctime_NonUni_thread / (Ctime_NonUni_thread + Ftime_NonUni_thread)
Ftime_NonUni_thread_pro = Ftime_NonUni_thread / (Ctime_NonUni_thread + Ftime_NonUni_thread)
Ctime_Uni_thread_pro = Ctime_Uni_thread / (Ctime_Uni_thread + Ftime_Uni_thread)
Ftime_Uni_thread_pro = Ftime_Uni_thread / (Ctime_Uni_thread + Ftime_Uni_thread)


bar_width = 0.5
index = np.arange(len(particleNumber))
delta = 0.01
diff = 0.235
fontsize = 6
font = FontProperties(weight='bold')

# create figure
fig, ax = plt.subplots( 1, 2, figsize=(7, 4), sharex=False, sharey=False, dpi=200 )
# fig.subplots_adjust(left=0.1, right=0.9, bottom=0.05, top=0.95, wspace=0.1)

ax[0].barh(index, Ctime_NonUni_pro, color='firebrick', label='construct time(nonuni)', height=bar_width, align='edge', edgecolor='black')
ax[0].barh(index, Ftime_NonUni_pro, color='lightcoral', label='force time(nonuni)', left=Ctime_NonUni_pro, height=bar_width, align='edge', edgecolor='black')

ax[1].barh(index, Ctime_Uni_pro, color='cornflowerblue', label='construct time(uni)', height=bar_width, align='edge', edgecolor='black')
ax[1].barh(index, Ftime_Uni_pro, color='skyblue', label='force time(uni)', left=Ctime_Uni_pro, height=bar_width, align='edge', edgecolor='black')

xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
xticklabels = ['0', '20', '40', '60', '80', '100']
ax[0].set_yticks(index + bar_width / 2)
ax[0].set_yticklabels(np.log10(particleNumber))
ax[0].set_xticks(xticks)
ax[0].set_xticklabels(xticklabels)

ax[1].set_yticks(index + bar_width / 2)
# ax[1].set_yticklabels(np.log10(particleNumber))
ax[1].set_xticks(xticks)
ax[1].set_xticklabels(xticklabels)
ax[1].set_yticklabels([])



for i, v in enumerate(Ctime_NonUni):
    string = "{:.2f}s\n{:.1f}%".format(v, Ctime_NonUni_pro[i]*100)
    ax[0].text(Ctime_NonUni_pro[i]/2, i+diff, string, ha='center', va='center', fontsize=fontsize, fontproperties=font)
    # print("i = ", i, "; v = ", v, "; index = ", index[i])

for i, v in enumerate(Ftime_NonUni):
    string = "{:.2f}s\n{:.1f}%".format(v, Ftime_NonUni_pro[i]*100)
    ax[0].text(Ftime_NonUni_pro[i]/2 + Ctime_NonUni_pro[i], i+diff, string, ha='center', va='center', fontsize=fontsize, fontproperties=font)

for i, v in enumerate(Ctime_Uni):
    string = "{:.2f}s\n{:.1f}%".format(v, Ctime_Uni_pro[i]*100)
    ax[1].text(Ctime_Uni_pro[i]/2, i+diff, string, ha='center', va='center', fontsize=fontsize, fontproperties=font)
    # print("i = ", i, "; v = ", v, "; index = ", index[i])

for i, v in enumerate(Ftime_Uni):
    string = "{:.2f}s\n{:.1f}%".format(v, Ftime_Uni_pro[i]*100)
    ax[1].text(Ftime_NonUni_pro[i]/2 + Ctime_Uni_pro[i], i+diff, string, ha='center', va='center', fontsize=fontsize, fontproperties=font)


ax[0].set_ylabel('particle numbers (log)')
ax[0].set_xlabel('Proportion of times (%)')
ax[1].set_xlabel('Proportion of times (%)')
fig.suptitle('Proportion of construct time and force time', fontsize=15, fontproperties=font )
ax[0].set_title('Nonuniform')
ax[1].set_title('Uniform')
ax[0].legend(bbox_to_anchor=(0.5, -0.2), loc='upper center', fontsize=fontsize)
ax[1].legend(bbox_to_anchor=(0.5, -0.2), loc='upper center', fontsize=fontsize)
for ax in ax.flatten():
    ax.spines['top'].set_visible(False) 
    ax.spines['right'].set_visible(False) 

plt.tight_layout()



# plt.savefig('bar_chart.png')
plt.savefig('bar_chart_pro_ptcnumber.png')
plt.show()



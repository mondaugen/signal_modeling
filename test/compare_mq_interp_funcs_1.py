import os
import numpy as np
import matplotlib.pyplot as plt

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

plotoutpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
plotoutpath+='mq_mod_err_comp_'

fpre=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
labels={
        'mq_mod_quintic':'DDM Quintic',
        'mq_mod_cubic':'DDM Cubic',
        'mq_cubic': 'MQ Cubic',
        'true_arg': 'True'
}
line_styles={
        'mq_mod_quintic':':',
        'mq_mod_cubic':'--',
        'mq_cubic':'-', 
        'true_arg':'-'
}
colors={
        'mq_mod_quintic':'k',
        'mq_mod_cubic':'k',
        'mq_cubic':'k',
        'true_arg':'grey'
}
d_ph={
    'mq_mod_quintic': None,
    'mq_mod_cubic': None,
    'mq_cubic': None
}
for k in d_ph.keys():
    with open(fpre+k+'_arg_ph.f64','r') as f:
        d_ph[k]=np.fromfile(f,dtype='double')
        N=len(d_ph[k])
with open(fpre+'mq_mod_quintic_arg_ph_x.f64','r') as f:
    d_ph['true_arg']=np.fromfile(f,dtype='double')
n=np.arange(N)
plt.figure(1)
for k in d_ph.keys():
    plt.plot(n,np.unwrap(d_ph[k]),color=colors[k],ls=line_styles[k],label=labels[k])
plt.legend()
plt.title('Phase function')
plt.xlim([0,8000])
plt.ylim([0,750])
plt.xlabel('Time (samples)')
plt.ylabel('Phase (radians)')
plt.savefig(plotoutpath+'phase_func.eps')

d_a={
    'mq_mod_quintic': None,
    'mq_mod_cubic': None,
    'mq_cubic': None
}
for k in d_a.keys():
    with open(fpre+k+'_arg_a.f64','r') as f:
        d_a[k]=np.fromfile(f,dtype='double')
        N=len(d_a[k])
with open(fpre+'mq_mod_quintic_arg_a_x.f64','r') as f:
    d_a['true_arg']=np.fromfile(f,dtype='double')
n=np.arange(N)
plt.figure(2)
for k in d_a.keys():
    plt.plot(n,np.unwrap(d_a[k]),
            color=colors[k],
            ls=line_styles[k],
            label=labels[k])
plt.legend()
plt.title('Log-amplitude function')
plt.xlim([2300,3900])
plt.ylim([0.22,0.32])
plt.xlabel('Time (samples)')
plt.ylabel('Log-amplitude')
plt.savefig(plotoutpath+'logamp_func.eps')

plt.figure(3)
for k in d_ph.keys():
    plt.plot(n,np.unwrap(d_ph[k])-d_ph['true_arg'],
            color=colors[k],
            ls=line_styles[k],
            label=labels[k])
plt.legend()
plt.title('Phase function error')
plt.xlim([0,7750])
plt.ylim([-2,7.75])
plt.xlabel('Time (samples)')
plt.ylabel('Error (radians)')
plt.savefig(plotoutpath+'phase_err.eps')

plt.figure(4)
for k in d_a.keys():
    plt.plot(n,np.unwrap(d_a[k])-d_a['true_arg'],
            color=colors[k],
            ls=line_styles[k],
            label=labels[k])
plt.legend()
plt.title('Log-amplitude function error')
plt.xlim([512,7750])
plt.ylim([-0.06,0.06])
plt.xlabel('Time (samples)')
plt.ylabel('Error (log-amplitude)')
plt.savefig(plotoutpath+'logamp_err.eps')

plt.show()

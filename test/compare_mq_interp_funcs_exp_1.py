import os
import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm

show_plots=False

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

plotoutpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
plotoutpath+='mq_exp_err_comp_'

phase_err_mod_2pi=True

fpre=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
labels={
        'mq_exp_mod_quintic':'DDM Quintic',
        'mq_exp_mod_cubic':'DDM Cubic',
        'mq_exp': 'MQ Cubic',
        'true_arg': 'True'
}
line_styles={
        'mq_exp_mod_quintic':':',
        'mq_exp_mod_cubic':'--',
        'mq_exp':'-', 
        'true_arg':'-'
}
colors={
        'mq_exp_mod_quintic':'k',
        'mq_exp_mod_cubic':'k',
        'mq_exp':'k',
        'true_arg':'grey'
}
d_ph={
    'mq_exp_mod_quintic': None,
    'mq_exp_mod_cubic': None,
    'mq_exp': None
}
d_x={
    'mq_exp_mod_quintic': None,
    'mq_exp_mod_cubic': None,
    'mq_exp': None,
}

for k in d_ph.keys():
    with open(fpre+k+'_arg_ph.f64','r') as f:
        d_ph[k]=np.fromfile(f,dtype='double')
with open(fpre+'mq_exp_mod_quintic_arg_ph_x.f64','r') as f:
    d_ph['true_arg']=np.fromfile(f,dtype='double')
for k in d_x.keys():
    d_x[k]=np.fromfile(fpre+k+'_est_x.dat',dtype='complex_')
d_x['true_arg']=np.fromfile(fpre+'mq_exp'+'_true_x.dat',dtype='complex_')

plt.figure(1)
for k in d_ph.keys():
    N=len(d_ph[k])
    n=np.arange(N)
    plt.plot(n,np.unwrap(d_ph[k]),color=colors[k],ls=line_styles[k],label=labels[k])
plt.legend()
plt.title('Phase function')
plt.xlim([0,8000])
plt.ylim([0,750])
plt.xlabel('Time (samples)')
plt.ylabel('Phase (radians)')
plt.savefig(plotoutpath+'phase_func.eps')

d_a={
    'mq_exp_mod_quintic': None,
    'mq_exp_mod_cubic': None,
    'mq_exp': None
}
for k in d_a.keys():
    with open(fpre+k+'_arg_a.f64','r') as f:
        d_a[k]=np.fromfile(f,dtype='double')
with open(fpre+'mq_exp_mod_quintic_arg_a_x.f64','r') as f:
    d_a['true_arg']=np.fromfile(f,dtype='double')
plt.figure(2)
for k in d_a.keys():
    N=len(d_a[k])
    n=np.arange(N)
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
    N=min(len(d_ph[k]),len(d_ph['true_arg']))
    n=np.arange(N)
    tmp=d_ph[k][:N]-d_ph['true_arg'][:N]
    if (phase_err_mod_2pi):
        tmp+=np.pi
        tmp%=(2.*np.pi)
        tmp-=np.pi
    plt.plot(n,tmp,
            color=colors[k],
            ls=line_styles[k],
            label=labels[k])
plt.title('Phase function error')
plt.xlim([0,7750])
plt.ylim([-np.pi*1.1,np.pi*1.1])
plt.xlabel('Time (samples)')
plt.ylabel('Error (radians)')
plt.legend(loc='best')
plt.savefig(plotoutpath+'phase_err.eps')

plt.figure(4)
for k in d_a.keys():
    N=min(len(d_a[k]),len(d_a['true_arg']))
    n=np.arange(N)
    tmp=d_a[k][:N]-d_a['true_arg'][:N]
    #plt.plot(n,np.unwrap(d_a[k])-d_a['true_arg'],
    plt.plot(n,tmp,
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

plt.figure(5)
for k in d_x.keys():
    N=min(len(d_x[k]),len(d_x['true_arg']))
    n=np.arange(N)
    tmp=20.*np.log10(np.abs(d_x[k][:N]-d_x['true_arg'][:N]))
    ma_,mai_=sm.lextrem(tmp,comp='max')
    #plt.plot(n,np.unwrap(d_a[k])-d_a['true_arg'],
    if (k != 'true_arg'):
        plt.plot(n[mai_],tmp[mai_],
                color=colors[k],
                ls=line_styles[k],
                label=labels[k])
plt.title('Original vs. estimated signals: upper error bound')
plt.xlabel('Time (samples)')
plt.ylabel('Error (dB power)')
plt.xlim([0,7750])
plt.legend(loc='lower right')
plt.savefig(plotoutpath+'true_vs_est_err.eps')

if (show_plots):
    plt.show()

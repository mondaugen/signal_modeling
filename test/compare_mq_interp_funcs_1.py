import os
import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm
import neplot as nep

# Color contrast config
# values further from 1, more contrast
clr_gamma=4.
clr_mapper=nep.PowerNormalize(clr_gamma)

show_plots=False

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

plotoutpath=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
plotoutpath+='mq_mod_err_comp_'

phase_err_mod_2pi=True

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

d_x={
    'mq_mod_quintic': None,
    'mq_mod_cubic': None,
    'mq_cubic': None,
}

d_xy={
    'true_arg':(0,0),
    'mq_cubic':(0,1),
    'mq_mod_cubic':(1,0),
    'mq_mod_quintic':(1,1)
}

for k in d_x.keys():
    d_x[k]=np.fromfile(fpre+k+'_est_x.dat',dtype='complex_')
d_x['true_arg']=np.fromfile(fpre+'mq_cubic'+'_true_x.dat',dtype='complex_')

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
    tmp=d_ph[k]-d_ph['true_arg']
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
    #plt.plot(n,np.unwrap(d_a[k])-d_a['true_arg'],
    plt.plot(n,d_a[k]-d_a['true_arg'],
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
plt.legend(loc='best')
plt.savefig(plotoutpath+'true_vs_est_err.eps')

# Hop size
H=128
# Analysis window / FFT size
N=512
# Sample rate
Fs=16000.

plt.figure(6)
fig6,asx5=plt.subplots(2,2,sharey=True,sharex=True,num=6)
for k in d_xy.keys():
    x_,y_=d_xy[k]
    asx5[x_,y_].specgram(d_x[k],NFFT=N,noverlap=(N-H),Fs=Fs,norm=clr_mapper,cmap="Greys")
    asx5[x_,y_].set_xlim([0,7750./Fs])
    asx5[x_,y_].set_ylim([50,300])
    asx5[x_,y_].set_title(labels[k])
    if (x_==1):
        asx5[x_,y_].set_xlabel('Time (seconds)')
    if (y_==0):
        asx5[x_,y_].set_ylabel('Frequency (Hz)')
fig6.suptitle('Spectrogram of original and resynthesized signals')
plt.savefig(plotoutpath+'all_spect.eps')

if (show_plots):
    plt.show()

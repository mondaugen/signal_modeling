import os
import numpy as np
import matplotlib.pyplot as plt
fpre=os.environ['HOME']+'/Documents/development/masters_thesis/reports/plots/'
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
    plt.plot(n,np.unwrap(d_ph[k]),label=k)
plt.legend()
plt.title('Phase function')

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
    plt.plot(n,np.unwrap(d_a[k]),label=k)
plt.legend()
plt.title('Log-amplitude function')
plt.show()

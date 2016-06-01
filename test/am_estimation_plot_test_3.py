# Here we estimate the amplitude and frequency of a real sound
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm

N=4096
M=32000*2+N
n=np.arange(N)
Fs=16000
H=512
m=np.arange(M)
# If the estimated values exceed these thresholds, they are considered erroneous
# Cannot change more than a1_r_th dB / hop
a1_r_th=20
# Absolute minimum value of acknowledged partial
th_db=-80 
a_r_e=[]
a_i_e=[]
b_ddm_hz=70.
o_ddm_hz=35.
b_ddm=np.round(b_ddm_hz/Fs*N)
o_ddm=np.round(o_ddm_hz/Fs*N)
print 'b_ddm, o_ddm = %d %d' %(b_ddm,o_ddm)
th_ddm=0.01
M_ddm=N/2
i_ddm=4

in_fname='tmp/ac_gtr_b3_gs4_sr16k.f64'
out_fname='/tmp/ac_gtr_b3_gs4_sr16k_aprx.f64'
res_fname='/tmp/ac_gtr_b3_gs4_sr16k_res.f64'
#in_fname='tmp/ac_gtr_a3_op_sr16k.f64'
#out_fname='/tmp/ac_gtr_a3_op_sr16k_aprx.f64'
#res_fname='/tmp/ac_gtr_a3_op_sr16k_res.f64'
#in_fname='tmp/close_guits_sr16k.f64'
#out_fname='/tmp/close_guits_aprx.f64'
#res_fname='/tmp/close_guits_res.f64'

with open(in_fname,'r') as f:
    x=np.r_[np.zeros(N),np.fromfile(f,dtype=np.float64)[:M-N]]

wc=np.r_[0.358735,0.488305,0.141265,0.011695]
#wc=np.r_[0.5,0.5]
w_=((wc*np.power(-1,np.arange(len(wc))))[:,np.newaxis]
        *np.cos(np.pi*2./N*np.outer(np.arange(len(wc)),n)))
w=np.sum(w_,0)
dw_=((2.*np.pi/N*wc[1:]*np.arange(1,len(wc))
        *np.power(-1,1+np.arange(1,len(wc))))[:,np.newaxis]
    *np.sin(np.pi*2./N*np.outer(np.arange(1,len(wc)),n)))
dw=np.sum(dw_,0)
# Coefficient to make normalize overlapped windows
c_ol=(np.r_[w,np.zeros(H)]+np.r_[np.zeros(H),w]).max()

a=[]
x_=np.zeros(M).astype('complex_')
h=0
while h <= (M-N):
    a.append(sm.ddm_p2_1_3_b(x[h:h+N],w,dw,
        b_ddm,o_ddm,th_ddm,M_ddm,i_ddm))
    a_r_e.append([])
    a_i_e.append([])
    for a_ in a[-1]:
        # Here find values of extremum and end of parabola
        # Sort them from highest to lowest, then subtract highest from lowest to
        # get total change in amplitude
        # check to see if this change is below a threshold
        a1_r_end=20.*np.real(a_[1]*N+a_[2]*(N**2.))/np.log(10.)
        n_ma=np.real(a_[1])/np.real((-2.*a_[2]))
        if (n_ma > 0.) and (n_ma < N):
            a1_r_ma=20.*np.real(a_[1]*n_ma+a_[2]*n_ma**2.)/np.log(10.)
            a1_extrm=np.sort(np.r_[0.,a1_r_ma,a1_r_end])
        else:
            a1_extrm=np.sort(np.r_[0.,a1_r_end])
        a1_r_ch=a1_extrm[-1]-a1_extrm[0]
#        print a1_r_ch  
        # if
        if ((a1_r_ch < a1_r_th) 
                and a1_extrm[-1]+20.*np.real(a_[0])/np.log(10.) > th_db):
#            print a1_extrm
            a_r_e[-1].append(np.real(a_[0]+a_[1]*n[:H]+a_[2]*n[:H]**2.))
            #print '\t',(a_r_e[-1][-1].min()-np.real(a_[0]))*20/np.log(10.)
            #print '\t',(a_r_e[-1][-1].max()-np.real(a_[0]))*20/np.log(10.)
            a_i_e[-1].append(np.imag(a_[1]+2.*a_[2]*n[:H]))
            x_[h:h+N]+=np.exp(a_[0]+a_[1]*n)*w/c_ol
    h+=H

mpl.rcParams['legend.fontsize'] = 10

fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')

h=0
for k in xrange(len(a_i_e)):
    for l in xrange(len(a_i_e[k])):
        ax.plot(np.arange(h,h+H).astype('d'),a_i_e[k][l]/(2*np.pi)*Fs,a_r_e[k][l]*20/np.log(10.),
            label='estimated',c='g')
    h+=H

#ax.plot(m, a_i_e, a_r_e*20/np.log(10.), label='estimated', c='g')
ax.legend()

fig2=plt.figure(2)
ax=fig2.gca()
ax.plot(m,np.real(x),label='true')
ax.plot(m,2.*np.real(x_),label='estimated')
ax.legend()

with open(out_fname,'w') as f:
    (2.*np.real(x_)).tofile(f)
with open(res_fname,'w') as f:
    (np.real(x)-2.*np.real(x_)).tofile(f)

plt.show()

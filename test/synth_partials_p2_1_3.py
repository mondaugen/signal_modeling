# Load in partials classified by plot_pep_3_3e.py
# Choose partial set.
# Synthesize partial set.

N_ARGS=3

import sys
import pickle
import numpy as np
import sigmod as sm

if (len(sys.argv) != N_ARGS):
    sys.stderr.write("""
    Arguments:
        input file
        partial set number

    Synthesizes to a file with same name as input file but with .f64 appended.
    The file is 1 channel of float64 samples. The sample rate of the output file
    is 16 KHz and it is assumed that the partials were also measured at the same
    rate.
    """)
    sys.exit(-1)

with open(sys.argv[1],'r') as f:
    ptls_all=pickle.load(f)

Fs=16000.

# Length of taper on first and last frame as percentage of frame size, to avoid
# clicking
a_tpr=0.01
grp=int(sys.argv[2])

ptls=ptls_all[grp]

# Find longest partial and make this the length of the sound
N=max(max(x_,key=lambda x__: x__[0])[0] for x_ in ptls)

y=np.zeros(N,dtype=np.float64)

for th in ptls:
    th.sort(key=lambda x_: x_[0])
    h=th[0][0]
    for i in xrange(len(th)-1):
        H=th[i+1][0]-th[i][0]
        N_tpr=int(H*a_tpr)
        w_tpr=np.hanning(2*N_tpr)
        y_=np.zeros(H,dtype=np.float64)
        phi_i0=np.imag(th[i][1][0])
        phi_i1=np.imag(th[i+1][1][0])
        w_i0=np.imag(th[i][1][1])
        w_i1=np.imag(th[i+1][1][1])
        psi_i0=np.imag(th[i][1][2])
        psi_i1=np.imag(th[i+1][1][2])
        # Compute M*
        M=np.round(((w_i1+w_i0)*H-2.*(phi_i1-phi_i0))/(4.*np.pi))
        # Compute polynomial coefficients
        q_=phi_i1-phi_i0+2.*np.pi*M
        r_=0.5*(w_i1+w_i0)
        s_=0.5*(psi_i1+psi_i0)
        # Done this way for more numerical stability, maybe?
        b=np.asmatrix(np.c_[q_,r_,s_].T)
        V=np.asmatrix(np.array([
            [4./H,-4.,0],
            [-6./H,6.,0.5*H],
            [3./H,-2.,-0.5*H]
        ]))
        D=np.diag([1./(H**2.),1./H,1.])
        c_=D*V*b
        c=np.array(np.vstack((c_,phi_i0))).flatten()
        y_=np.cos(np.polyval(c,np.arange(H)))
        a0_i0=np.real(th[i][1][0])
        a0_i1=np.real(th[i+1][1][0])
        a1_i0=np.real(th[i][1][1])
        a1_i1=np.real(th[i+1][1][1])
        a2_i0=np.real(th[i][1][2])
        a2_i1=np.real(th[i+1][1][2])
        # Find coefficients of cubic amplitude function:
        # mu(t) = d3*t^3 + d2*t^2 + d1*t + d0
        # at t=0, mu(0)= a0_i0 so
        d0=a0_i0
        # We find the rest of the coefficients using least squares
        # constrained so that
        # mu(H) = a0_i1
        # The least squares minimizes the error fit of the following:
        # mu'(0) = a1_i0
        # mu'(H) = a1_i1
        # mu''(0) = a2_i0
        # mu''(H) = a2_i1
        # Therefore the observation matrix is
        A_o=np.array([
            [0.,0.,1.],
            [3.*H*H,2*H,1.],
            [0.,2.,0.],
            [6.*H,2.,0.]
        ])
        # The constraint matrix is
        A_c=np.array([[H**3.,H**2.,H]])
        # Constraint vector
        b_c=np.array([[a0_i1-a0_i0]])
        # observation vector
        b_o=np.c_[a1_i0,a1_i1,a2_i0,a1_i1].T
        d_=sm.lstsq_c(A_o,b_o,A_c,b_c)
        d=np.array(np.vstack((d_,d0))).flatten()
        # Multiply by amplitude function
        y_*=np.exp(np.polyval(d,np.arange(H)))
        if (i==0):
            # Taper the beginning of the frame
            y_[:N_tpr]*=w_tpr[:N_tpr]
        elif (i==(len(th)-2)):
            # Taper the end of the frame
            y_[-N_tpr:]*=w_tpr[-N_tpr:]

        y[h:h+H]+=y_
        h+=H

fout_name=sys.argv[1]
fsuffidx=fout_name.rfind('.')
if (fsuffidx==-1):
    fout_name=fout_name+'.f64'
else:
    fout_name=fout_name[:fsuffidx]+'_p'+str(grp)+'.f64'

with open(fout_name,'w') as f:
    y.tofile(f)

print 'Wrote to: '+fout_name

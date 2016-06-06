# Perform STFT and estimate partial slopes of recording of
# flute

with open('tmp/flt_b4_sr16k.f64','r') as f:
    x=np.fromfile(f)
# Sample rate of file
Fs=16000
# Length of the file
N=len(x)
# sample indices of the file
n=np.arange(N)
# space to store the estimate ddm parameters
a=[]
# current hop
h=0
M=4096
Y=np.zeros((M,(N-M)/H+1)).astype('complex_')
m=np.arange(M)
# Windows for DDM
#w=0.5+0.5*np.cos(2.*np.pi*(m-M/2)/M)
#dw=-np.pi/M*np.sin(2.*np.pi*(m-M/2)/M)
wc=np.r_[0.358735,0.488305,0.141265,0.011695]
#wc=np.r_[0.5,0.5]
w_=((wc*np.power(-1,np.arange(len(wc))))[:,np.newaxis]
        *np.cos(np.pi*2./M*np.outer(np.arange(len(wc)),m)))
w=np.sum(w_,0)
dw_=((2.*np.pi/M*wc[1:]*np.arange(1,len(wc))
        *np.power(-1,1+np.arange(1,len(wc))))[:,np.newaxis]
    *np.sin(np.pi*2./M*np.outer(np.arange(1,len(wc)),m)))
dw=np.sum(dw_,0)
# Coefficient to make normalized overlapped windows
c_ol=(np.r_[w,np.zeros(H)]+np.r_[np.zeros(H),w]).max()
k=0
w_s=np.sum(w)
k_S=0

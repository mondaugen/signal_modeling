import numpy as np
import matplotlib.pyplot as plt

Fs=16000
T=1.
N_fft=1024
N_overlap=256
N=T*Fs
x=np.zeros((N,),dtype='complex_')
y=np.zeros((N,),dtype='complex_')
n=np.arange(N)
f0=100.
f1=150.
K_h=10
for k in xrange(K_h):
    x+=np.exp(1j*2.*np.pi*k*f0/Fs*n)*np.r_[n[:N/2]*2./N,n[N/2-1::-1]*2./N]
    y+=np.exp(1j*2.*np.pi*k*f1/Fs*n)*np.r_[n[:N/2]*2./N,n[N/2-1::-1]*2./N]
z=np.concatenate((x,y,x+y))
plt.figure(1)
S,F,T,IM=plt.specgram(z,NFFT=N_fft,Fs=Fs,noverlap=N_overlap)
K=2
V=S[N_fft/2:((0.5)+K_h*max(f0,f1)/Fs)*N_fft]
W=np.random.rand(V.shape[0],K)
H=np.random.rand(K,V.shape[1])
L=100
bet=1.
for l in xrange(L):
    H=(H*(np.array(np.mat(W.T)*np.mat((np.array(np.mat(W)*np.mat(H))**(bet-2.))*V))
            /np.array(np.mat(W.T)*np.mat(np.array(np.mat(W)*np.mat(H))**(bet-1.)))))
    W=(W*(np.array(np.mat((np.array(np.mat(W)*np.mat(H))**(bet-2.))*V)*np.mat(H.T))
            /np.array(np.mat(np.array(np.mat(W)*np.mat(H))**(bet-1.))*np.mat(H.T))))

plt.figure(2)
plt.imshow(np.multiply.outer(W[:,0],H[0,:]))
plt.figure(3)
plt.imshow(np.multiply.outer(W[:,1],H[1,:]))
plt.figure(4)
plt.imshow(np.array(np.mat(W)*np.mat(H)))
plt.show()

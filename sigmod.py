import numpy as np

def lextrem(a,comp='max'):
    """
    Output local extrema. Regions of equal values with different values
    surrounding are not considered extrema! Values on the edges of the array are
    not considered extrema either.
    """
    cmpfunc={'max':np.greater,'min':np.less}[comp]
    mab=np.r_[False,cmpfunc(a[1:],a[:-1])]&np.r_[cmpfunc(a[:-1],a[1:]),False]
    mai=np.transpose(np.nonzero(mab))
    ma=a[mai]
    return (ma,mai)

#def ddm_p2_3_3(x,H,w,dw,k_de):
#    """
#    Compute parameters of 2nd order polynomial using 8 bins surrounding the
#    maximum of the STFT at the centre of signal x.
#
#    x:  
#        the signal to analyse, must have length of at least 2*H+N_w where N_w is
#        the length of the window.
#    H:  
#        the hop size between analysis frames.
#    w:  
#        the analysis window.
#    dw: 
#        the derivative of the analysis window
#    k_de: 
#        the number of bins over which the maximum amplitude in the following
#        frame is searched
#
#    Returns 
#
#    a:
#        a vector containing the estimated parameters
#    """
#    N_x=len(x)
#    N_w=len(w)
#    if (N_w%2) == 0:
#        raise Exception('Window length must be odd.')
#    Ncx=(N_x-1)/2
#    Ncw=(N_w-1)/2
#    nx0=np.arange(-Ncw,Ncw+1)+Ncx
#    nx_H=nx0-H
#    nxH=nx0+H
#    x0=x[nx0]
#    xH=x[nxH]
#    x_H=x[nx_H]
#    # order is[-H;0:H]
#    # FFT are zero padded by one assuming N_w has odd length.
#    x_stk=np.vstack((x_H,x0,xH))
#    Xp1w=np.fft.fft(x_stk*w,N_w+1)
#    Xp2w=np.fft.fft(np.vstack((
#        2.*(nx_H-Ncx)*x_H,
#        2.*(nx0-Ncx)*x0,
#        2.*(nxH-Ncx)*xH))*w,
#        N_w+1)
#    Xdw_=np.fft.fft(x_stk*dw,N_w+1)
#    Xdw=Xp1w*(-2.*np.pi*1j*np.arange(0,N_w+1)/(N_w+1))+Xdw_
#    kma0=np.abs(Xp1w[1,:]).argmax()
#    kma_H=np.abs(Xp1w[0,(kma0-k_de):(kma0+k_de+1)]).argmax()+kma0-k_de
#    kmaH=np.abs(Xp1w[2,(kma0-k_de):(kma0+k_de+1)]).argmax()+kma0-k_de
#    A=np.c_[
#            np.r_[
#                Xp1w[0,(kma_H)-1:(kma_H+2)],
#                Xp1w[1,(kma0)-1:(kma0+2)],
#                Xp1w[2,(kmaH)-1:(kmaH+2)]
#                ],
#            np.r_[
#                Xp2w[0,(kma_H)-1:(kma_H+2)],
#                Xp2w[1,(kma0)-1:(kma0+2)],
#                Xp2w[2,(kmaH)-1:(kmaH+2)]
#                ]
#            ]
#    b=np.c_[
#            np.r_[
#                Xdw[0,(kma_H)-1:(kma_H+2)],
#                Xdw[1,(kma0)-1:(kma0+2)],
#                Xdw[2,(kmaH)-1:(kmaH+2)]
#                ]
#            ]
#    a=np.linalg.lstsq(A,-b)[0]
#    return a

def ddm_p2_3_3(x,H,w,dw,k_de):
    """
    Compute parameters of 2nd order polynomial using 8 bins surrounding the
    maximum of the STFT at the centre of signal x.

    x:  
        the signal to analyse, must have length of at least 2*H+N_w where N_w is
        the length of the window.
    H:  
        the hop size between analysis frames.
    w:  
        the analysis window.
    dw: 
        the derivative of the analysis window
    k_de: 
        the number of bins over which the maximum amplitude in the following
        frame is searched

    Returns 

    a:
        a vector containing the estimated parameters
    """
#    N_x=len(x)
    N_w=len(w)
    if (N_w%2) == 1:
        raise Exception('Window length must be even.')
    n=np.arange(N_w)
    nx_H=n
    nx0=n+H
    nxH=n+2*H
    x0=x[nx0]
    xH=x[nxH]
    x_H=x[nx_H]
    # order is[-H;0:H]
    # FFT are zero padded by one assuming N_w has odd length.
    x_stk=np.vstack((x_H,x0,xH))
    Xp1w=np.fft.fft(x_stk*w)
    Xp2w=np.fft.fft(np.vstack((
        2.*(n-H)*x_H,
        2.*(n)*x0,
        2.*(n+H)*xH))*w)
    Xdw_=np.fft.fft(x_stk*dw)
    Xdw=Xp1w*(-2.*np.pi*1j*n/N_w)+Xdw_
    kma0=np.abs(Xp1w[1,:]).argmax()
    kma_H=np.abs(Xp1w[0,(kma0-k_de):(kma0+k_de+1)]).argmax()+kma0-k_de
    kmaH=np.abs(Xp1w[2,(kma0-k_de):(kma0+k_de+1)]).argmax()+kma0-k_de
    A=np.c_[
            np.r_[
                Xp1w[0,(kma_H)-1:(kma_H+2)],
                Xp1w[1,(kma0)-1:(kma0+2)],
                Xp1w[2,(kmaH)-1:(kmaH+2)]
                ],
            np.r_[
                Xp2w[0,(kma_H)-1:(kma_H+2)],
                Xp2w[1,(kma0)-1:(kma0+2)],
                Xp2w[2,(kmaH)-1:(kmaH+2)]
                ]
            ]
    b=np.c_[
            np.r_[
                Xdw[0,(kma_H)-1:(kma_H+2)],
                Xdw[1,(kma0)-1:(kma0+2)],
                Xdw[2,(kmaH)-1:(kmaH+2)]
                ]
            ]
    a=np.linalg.lstsq(A,-b)[0]
    return a

def ddm_p2_1_3(x,w,dw):
    """
    Compute parameters of 2nd order polynomial using 2 bins surrounding the
    maximum of the STFT at the centre of signal x.

    x:  
        the signal to analyse, must have length of at least N_w where N_w is
        the length of the window.
    w:  
        the analysis window.
    dw: 
        the derivative of the analysis window

    Returns 

    a:
        a vector containing the estimated parameters
    """
    N_x=len(x)
    N_w=len(w)
    Ncx=(N_x-1)/2
    Ncw=(N_w-1)/2
    nx0=np.arange(N_x)
    x0=x[nx0]
    # FFT are zero padded by one assuming N_w has odd length.
    Xp1w=np.fft.fft(x0*w)
    Xp2w=np.fft.fft(2.*nx0*x0*w)
    Xdw_=np.fft.fft(x0*dw)
    Xdw=Xp1w*(-2.*np.pi*1j*nx0/N_x)+Xdw_
    kma0=np.abs(Xp1w).argmax()
    A=np.c_[
            np.r_[
                Xp1w[(kma0)-1:(kma0+2)],
                ],
            np.r_[
                Xp2w[(kma0)-1:(kma0+2)],
                ]
            ]
    b=np.c_[
            np.r_[
                Xdw[(kma0)-1:(kma0+2)],
                ]
            ]
    a=np.linalg.lstsq(A,-b)[0]
    gam=np.exp(a[0]*nx0+a[1]*nx0**2.)
    a0=(np.log(np.inner(x,np.conj(gam)))
        -np.log(np.inner(gam,np.conj(gam))))
    return np.vstack((a0,a))

def ddm_p2_1_3_b(x,w,dw,b,o,th,M):
    """
    Compute parameters of 2nd order polynomial using 2 bins surrounding the
    maximum of the STFT at the centre of signal x.

    x:  
        the signal to analyse, must have length of at least N_w where N_w is
        the length of the window.
    w:  
        the analysis window.
    dw: 
        the derivative of the analysis window
    b:
        the size of the bands (in samples) in which to search for maxima
    o:
        the overlap of the bands (in samples)
    th:
        the minimum amplitude of a true peak (peaks lower than this are not
        considered)
    M:
        the maximum index to consider. Can be used to only consider half the
        spectrum for example.

    Returns 

    a:
        a vector containing the estimated parameters
    """
    N_x=len(x)
    N_w=len(w)
    Ncx=(N_x-1)/2
    Ncw=(N_w-1)/2
    nx0=np.arange(N_x)
    x0=x[nx0]
    # FFT are zero padded by one assuming N_w has odd length.
    Xp1w=np.fft.fft(x0*w)
    Xp2w=np.fft.fft(2.*nx0*x0*w)
    Xdw_=np.fft.fft(x0*dw)
    Xdw=Xp1w*(-2.*np.pi*1j*nx0/N_x)+Xdw_
    b_=0
#    wb=np.hanning(b)
    result=[]
    while (b_ < M):
        kma0=(np.abs(Xp1w)[b_:b_+b]).argmax()+b_
        if (np.abs(Xp1w)[kma0] < th):
            b_+=o
            continue
        b_=kma0+o
        A=np.c_[
                np.r_[
                    Xp1w[(kma0)-1:(kma0+2)],
                    ],
                np.r_[
                    Xp2w[(kma0)-1:(kma0+2)],
                    ]
                ]
        c=np.c_[
                np.r_[
                    Xdw[(kma0)-1:(kma0+2)],
                    ]
                ]
        try:
            a=np.linalg.lstsq(A,-c)[0]
            gam=np.exp(a[0]*nx0+a[1]*nx0**2.)
            a0=(np.log(np.inner(x,np.conj(gam)))
                -np.log(np.inner(gam,np.conj(gam))))
            result.append(np.vstack((a0,a)))
        except ValueError:
            pass
    return result

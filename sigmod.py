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

def lextrem_win(x,b,o,i,th,ar,M,bi=0):
    """
    Search in windows of size b every hop o for values exceeding th and that are
    whose dB ratio with the average of the minimum values in the frame is
    greater than ar.

    x:
        the signal to search for maxima in.

    b:
        the size of the bands (in samples) in which to search for maxima.
    o:
        the hopsize between bands (in samples).
    th:
        the minimum amplitude of a true peak (peaks lower than this are not
        considered).
    i:
        If a maximum is found, the band in which maxima are searched is not
        shifted by o but rather set to begin at the index of the last maximum +
        some number of ignored bins, given by i (default 1).
    ar:
        The threshold ratio 20*log10(max_peak/avg_min_peak). If over this value,
        the peak is considered, otherwise it is not. If the value is less than 0
        the peak is always considered (default -1). For example, if ar=40, the
        max_peak has to be 100 times that of the avg_min_peak.
    M:
        the maximum index to consider. Can be used to only consider half the
        spectrum for example.
    bi:
        The initial index to start on (default 0).


    Returns:

    ks:
        Indices of local maxima.

    """

    b_=bi
    ks=[]
    if (o == 0):
        raise Exception('Hopsize cannot be 0.')
    while (b_ < M):
        tmp_=x[b_:b_+b]
        # Indices of values greater than their neighbours
        kma=np.where(np.r_[False,tmp_[1:]>tmp_[:-1]]
                & np.r_[tmp_[:-1]>tmp_[1:],False])[0]
        if len(kma) <= 0:
            b_+=o
            continue
        kma0=kma[tmp_[kma].argmax()]+b_
        if (x[kma0] < th):
            b_+=o
            continue
        if ar > 0.:
            l_min=x[b_:kma0].min()
            r_min=x[kma0+1:b_+b].min()
            alr_min=0.5*(l_min+r_min)
            if (20.*np.log10(x[kma0]/alr_min)<ar):
                b_+=o
                continue
        ks.append(kma0)
        b_=kma0+i
    return ks

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
    gam=np.exp(a[0]*nx0+a[1]*nx0**2.)
    a0=(np.log(np.inner(x0,np.conj(gam)))
        -np.log(np.inner(gam,np.conj(gam))))
    return np.vstack((a0,a))

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
    N_w=len(w)
    nx0=np.arange(N_w)
    x0=x[nx0]
    Xp1w=np.fft.fft(x0*w)
    Xp2w=np.fft.fft(2.*nx0*x0*w)
    Xdw_=np.fft.fft(x0*dw)
    Xdw=Xp1w*(-2.*np.pi*1j*nx0/N_w)+Xdw_
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
    try:
        a=np.linalg.lstsq(A,-b)[0]
    except ValueError:
        return np.zeros((3,1)).astype('complex_')
    gam=np.exp(a[0]*nx0+a[1]*nx0**2.)
    a0=(np.log(np.inner(x0,np.conj(gam)))
        -np.log(np.inner(gam,np.conj(gam))))
    return np.vstack((a0,a))

def ddm_p2_1_3_b(x,w,dw,b,o,th,M,i=1,ar=-1.):
    """
    Compute parameters of 2nd order polynomial using 2 bins surrounding the
    maximum of the STFT at the centre of signal x.

    x:  
        the signal to analyse, must have length of at least N_w where N_w is
        the length of the window.
    w:  
        the analysis window.
    dw: 
        the derivative of the analysis window.
    b:
        the size of the bands (in samples) in which to search for maxima.
    o:
        the hopsize between bands (in samples).
    th:
        the minimum amplitude of a true peak (peaks lower than this are not
        considered).
    M:
        the maximum index to consider. Can be used to only consider half the
        spectrum for example.
    i:
        If a maximum is found, the band in which maxima are searched is not
        shifted by o but rather set to begin at the index of the last maximum +
        some number of ignored bins, given by i (default 1).
    ar:
        The threshold ratio 20*log10(max_peak/avg_min_peak). If over this value,
        the peak is considered, otherwise it is not. If the value is less than 0
        the peak is always considered (default -1). For example, if ar=40, the
        max_peak has to be 100 times that of the avg_min_peak.

    Returns 

    a:
        a vector containing the estimated parameters.
    """
    N_w=len(w)
    nx0=np.arange(N_w)
    x0=x[nx0]
    # FFT are zero padded by one assuming N_w has odd length.
    Xp1w=np.fft.fft(x0*w)
    Xp2w=np.fft.fft(2.*nx0*x0*w)
    Xdw_=np.fft.fft(x0*dw)
    Xdw=Xp1w*(-2.*np.pi*1j*nx0/N_w)+Xdw_
    result=[]
    ks = lextrem_win(np.abs(Xp1w),b,o,i,th,ar,M)
    #b_=0
    #while (b_ < M):
    #    tmp_=np.abs(Xp1w)[b_:b_+b]
    #    # Indices of values greater than their neighbours
    #    kma=np.where(np.r_[False,tmp_[1:]>tmp_[:-1]]
    #            & np.r_[tmp_[:-1]>tmp_[1:],False])[0]
#   #     kma0=(np.abs(Xp1w)[b_:b_+b]).argmax()+b_
    #    if len(kma) <= 0:
    #        b_+=o
    #        continue
    #    kma0=kma[tmp_[kma].argmax()]+b_
    #    if (np.abs(Xp1w)[kma0] < th):
    #        b_+=o
    #        continue
    #    if ar > 0.:
    #        l_min=np.abs(Xp1w)[b_:kma0].min()
    #        r_min=np.abs(Xp1w)[kma0+1:b_+b].min()
    #        alr_min=0.5*(l_min+r_min)
    #        if (20.*np.log10(np.abs(Xp1w)[kma0]/alr_min)<ar):
    #            b_+=o
    #            continue
    #    b_=kma0+i
    for kma0 in ks:
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
            a0=(np.log(np.inner(x0,np.conj(gam)))
                -np.log(np.inner(gam,np.conj(gam))))
            result.append(np.vstack((a0,a)))
        except ValueError:
            pass
    return result

def pca_ne(X,mode='cov'):
    """
    Calculate the principal components of X.

    X:
        A (PxN) matrix of N observations on P variables (this is the transpose
        of the X in pca_ne.m).
    mode:
        Can be
        'cov':
            Uses the covariance matrix to compute PCs.
        'corr':
            Uses the correlation matrix to compute PCs.

    Returns:

    A:
        The principal components of X. This is a matrix of size (PxN) and is N
        realizations of P principal components. The first PCs are in the first
        row and the last in the last row.

    """
    if mode == 'cov':
        S=np.cov(X)
    elif mode == 'corr':
        S=np.corrcoef(X)
    else:
        raise Exception('Bad mode %s.\n' % (mode,))
    (l,V)=np.linalg.eig(S)
    li=np.flipud(l.argsort())
    V=V[:,li]
    A=np.inner(V.T,X.T)
    return A

def ddm_p2_3_1(x,H,w,dw,k_de):
    """
    Compute parameters of 2nd order polynomial using 2 bins in frames adjacent
    to the maximum of the STFT in time at the centre of signal x.

    x:  
        the signal to analyse, must have length of at least 2*H+N_w where N_w is
        the length of the window.
    H:  
        the hop size between analysis frames.
    w:  
        the analysis window.
        It is recommended the window be such that it is an even function of time
        with the 0th index of the window corresponding to time 0.
    dw: 
        the derivative of the analysis window
    k_de: 
        the number of bins over which the maximum amplitude in the adjacent 
        frames is searched

    Returns 

    a:
        a vector containing the estimated parameters
    """
#    N_x=len(x)
    N_w=len(w)
    #if (N_w%2) == 2:
    #    raise Exception('Window length must be even.')
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
                Xp1w[0,kma_H],
                Xp1w[1,kma0],
                Xp1w[2,kmaH]
                ],
            np.r_[
                Xp2w[0,kma_H],
                Xp2w[1,kma0],
                Xp2w[2,kmaH]
                ]
            ]
    b=np.c_[
            np.r_[
                Xdw[0,kma_H],
                Xdw[1,kma0],
                Xdw[2,kmaH]
                ]
            ]
    a=np.linalg.lstsq(A,-b)[0]
    return a

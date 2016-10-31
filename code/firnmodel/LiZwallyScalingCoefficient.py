import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
#matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

spot = os.path.dirname(sys.argv[0])
os.chdir(spot)


def LZbeta(A,T):
    beta1=-9.788+8.996*A-0.6165*T
    beta2=beta1/(-2.0178+8.4043*A-0.0932*T)
    return beta1#,beta2
    
    
if __name__ == '__main__': 
    As=np.arange(0.01,1.0,0.005)
    Ts=np.arange(-60,0,0.1)
    b1out=np.zeros([len(As),len(Ts)])
    b2out=b1out
    
    for ii in range(len(As)):
        for jj in range(len(Ts)):
            A=As[ii]
            T=Ts[jj]
            beta1=LZbeta(A,T)
            if beta1>100:
                print A
                print T
            b1out[ii,jj]=beta1
            #b2out[ii,jj]=beta2
            
    TT,AA=np.meshgrid(Ts,As)
    
    plt.figure(1)
    plt.clf()
    plt.contourf(TT,AA,b1out,256,cmap=plt.cm.seismic, vmin=-10, vmax=10)
    plt.colorbar()
    plt.xlabel('Temperature (C)')
    plt.ylabel('Accumulation Rate (m yr$^{-1}$)')
    plt.title('Value of Scaling coefficient for Li & Zwally (2011)')
    plt.savefig('LZcoeff.eps')
    plt.show()
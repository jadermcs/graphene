import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import linalg as la
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from cycler import cycler


#FUNCAO CONTINUA
def potv(xa,multi,lw):
    x=abs(xa)/lw
    return multi*700.*(1+np.tanh(5000*(-1+x)))
################
#FUNCAO DESCONTINUA
def potdv(xa):
    if abs(xa)<35.:
        return 0
    else:
        return 700.
###################
##FUNCAO MASSA
def potmm(xa):
    if abs(xa)<35.:
        return .043
    else:
        return .120
        

def mtxg(multi,m,lw):
    #m=potmm(t)
    inter = lw
    n=91
    t=-inter
    h=2*inter/n
    hc=6.5821*10e-16
    vp = np.zeros((n+1,n+1),float)
    for w in range(n):
        vp[w][w] = multi*(((potv(t,multi,lw)*2.*m)/(hc*hc))+2./(h*h))
        vp[w+1][w] = multi*(-1./(h*h))
        vp[w][w+1] = multi*(-1./(h*h))
        t+=h
    vp[n][n] = multi*(((potv(t,multi,lw)*2.*m)/(hc*hc))+2./(h*h))
    Av, Aw = la.eig(vp)
    return np.sort(Av/(2.*m/(hc*hc)))



def main():
#********CONSTANTES*****
    hc=6.582*10e-16
    mo = 1.66*10e-27
    mgaas = .067*mo
    mhh = .45 *mo
    mso = .154*mo
    mhh = .45*mo
    mlh = .082*mo
    eg = 1.424
    esp = .34
    x = np.linspace(-10,10,num=40)

    line = [eg/2 for xk in x]

    largurapc = np.linspace(5.,60.,num=100)


    dif1 = []
    dif2 = []
    dif3 = []
    dif4 = []

    for ax in largurapc:
        mtxgaas = mtxg(1,mgaas,ax)+eg/2
        bc = mtxgaas[0].real
        mtxbv = mtxg(-1,mhh,ax)-eg/2
        bvhh = mtxbv[0].real
        mtxbv2 = mtxg(-1,(mso*mhh)/(mso+mhh),ax)-eg/2
        bvhh2 = mtxbv2[0].real
        mtxspin = mtxg(-1,mlh,ax)-eg/2-esp
        spin = mtxspin[0].real

        dif1.append(bc - bvhh)
        dif2.append(bvhh - bvhh2)
        dif3.append(bvhh2 - spin)
        dif4.append(bvhh - spin)
        
    plt.plot(largurapc,dif1,'r-',label='Gap I')
    plt.plot(largurapc,dif2,'k-',label='Gap II')
    plt.plot(largurapc,dif3,'g-',label='Gap III')
    plt.plot(largurapc,dif4,'b-',label='Gap IV')
    plt.legend()
    plt.show()
   

if __name__ == "__main__":
    main()

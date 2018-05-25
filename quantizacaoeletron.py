import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from scipy import linalg as la

pto = 5.

#FUNCAO CONTINUA
def potv(xa,lw):
    x=abs(xa)/lw
    return 7000000.*(1+np.tanh(5000*(-1+x)))
################
#FUNCAO DESCONTINUA
def potdv(t):
    if abs(t)<35.:
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
        

def mtxg(lw):
    #m=potmm(t)
    inter = lw
    n=91
    t=-inter
    h=2*inter/n
    hc=6.582*10e-16
    mo = 1.66*10e-27
    mgaas = .067
    m = mgaas *mo
    eg = 1.424
    vp = np.zeros((n+1,n+1),float)
    for w in range(n):
        vp[w][w] = ((potv(t,lw)*2.*m)/(hc*hc))+2./(h*h)
        vp[w+1][w] = -1./(h*h)
        vp[w][w+1] = -1./(h*h)
        t+=h
    vp[n][n] = ((potv(t,lw)*2.*m)/(hc*hc))+2./(h*h)
    Av, Aw = la.eig(vp)
    return np.sort(Av/(2.*m/(hc*hc)))# + x*x*hc*hc/(2*m) + eg/2



def main():
#********CONSTANTES*****
    hc=6.582*10e-16
    mo = 1.66*10e-27
    mgaas = .067
    m = mgaas *mo
    eg = 1.424
    x = np.linspace(-2,2,num=40)
    print(x)
    fx,ax = plt.subplots(3,sharex=True)
    lista = [1.,5.,10.]
    count = 0
    for lw in lista:
        delz = mtxg(lw)

        y =  np.array([( delz[0].real + xk*xk*hc*hc/(2*m) + eg/2) for xk in x])
        y2 = np.array([(xk*xk*hc*hc/(2*m) + eg/2) for xk in x])
        n = np.array([( delz[1].real + xk*xk*hc*hc/(2*m) + eg/2) for xk in x])
        j = np.array([( delz[2].real + xk*xk*hc*hc/(2*m) + eg/2) for xk in x])
        ax[count].plot(x,y2,linestyle='--',color='#63f97c')
        ax[count].plot(x,y,'g-',x,n,'k-',x,j,'b-')
        count+=1
    
    g1label = mpatches.Patch(color='g', label='Spin-Up',linestyle='-')
    g2label = mpatches.Patch(color='#63f97c', label='Spin-Down',linestyle='-')
    blabel = mpatches.Patch(color='b', label='Heavy-Hole2',linestyle='-')
    klabel = mpatches.Patch(color='k', label='Heavy-Hole2',linestyle='-')

    #plt.plot(x,y,'g--',x,y2,'r--',x,n,'k--',x,j,'b--')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=4,handles=[g1label,g2label,klabel,blabel])
    plt.show()
   

if __name__ == "__main__":
    main()

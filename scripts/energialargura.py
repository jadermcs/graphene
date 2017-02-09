import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as la

#FUNCAO CONTINUA
def potv(xa,lw):
    x=abs(xa)/lw
    return 700.*(1+np.tanh(5000*(-1+x)))
################
#FUNCAO DESCONTINUA
def potdv(xa,lw):
    if abs(xa)<lw:
        return 0
    else:
        return 700.
###################
##FUNCAO MASSA
def potmm(xa,lw):
    if abs(xa)<lw:
        return .043
    else:
        return .120
        

def mtxg(n,t,h,m,hc,lw):
    #m=potmm(t)
    vp = np.zeros((n+1,n+1),float)
    for w in range(n):
        vp[w][w] = (potv(t,lw)*2.*m/(hc*hc))+2./(h*h)
        vp[w+1][w] = -1./(h*h)
        vp[w][w+1] = -1./(h*h)
        t+=h
    vp[n][n]=(potv(t,lw)*2.*m/(hc*hc))+2./(h*h)
    return vp

def analy(num,h,m,inter):
    return (num**2*h**2*np.pi**2)/(2*m*inter**2)

def main(M=250):
#********CONSTANTES*****
    file1 = open('dados.dat','w')
    hcort = 1
    mo = 1
    m = .067 * mo
    inter = 100
    passo = 50
#***********************
    N = 2*M+1
    h = (2.*inter)/N
    ebarra = 1./(2.*h*h)
    niv1=[]
    niv2=[]
    niv3=[]
    niv4=[]
    for lw in np.linspace(20,inter,num=passo):
        Ha = mtxg(N,-lw,h,m,hcort,lw)
        Av,Aw = la.eig(Ha)
        Av = np.sort(Av/(2.*m/(hcort*hcort)))#sorted
        e2 = Av[0]
        eanaly = np.array([analy(x+1,hcort,m,lw*2) for x in range(25)])
        niv1.append(Av[0].real)
        niv2.append(Av[1].real)
        niv3.append(Av[2].real)
        niv4.append(Av[3].real)
        print("%.1f\t%.7f\t%.7f\t%.7f" % (lw,Av.real[0],Av.real[1],Av.real[2]))
        file1.write("%.1f\t%.7f\t%.7f\t%.7f\n" % (lw,Av.real[0],Av.real[1],Av.real[2]))
    file1.close()
    x = np.linspace(20,inter,num=passo)
    #arq = open("energ_num.dat", "w")
    #for a in x:
    #        arq.write(str(y[a])+"\n")
    plt.plot(x, niv1, 'r--', x, niv2, 'b--', x, niv3, 'g--', x, niv4, 'k--')
    #ax[4].plot(x,potv(x),'c--')
    plt.show()

if __name__ == "__main__":
    main()

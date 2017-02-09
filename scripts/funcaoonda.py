import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as la

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
    return Aw

def analy(num,h,m,inter):
    return (num**2*h**2*np.pi**2)/(2*m*inter**2)

def main():
    hc=6.582*10e-16
    mo = 1.66*10e-27
    mgaas = .067*mo
    mhh = .45 *mo
    mso = .154*mo
    mhh = .45*mo
    mlh = .082*mo
    eg = 1.424
    esp = .34
    hole = 35.

    f, (ax) = plt.subplots(1, 1, sharex=True)

    #lista = [ax1,ax2,ax3,ax4]

    matrizes = [mtxg(1,mgaas,hole),mtxg(-1,mhh,hole),mtxg(-1,(mso*mhh)/(mso+mhh),hole),mtxg(-1,mlh,hole)]
    
    for j in range(1):
	    valor1 = 0
	    valor2 = 0
	    valor3 = 0
	    y = matrizes[j][:,1]**2
	    y2 = matrizes[j][:,2]**2
	    y3 = matrizes[j][:,3]**2
	    x = np.linspace(-hole,hole,num=len(y))
	    buraco = potv(x,1,hole)
	    h = hole*2/91
	    for i in y:
	    	valor1 += i*h
	    for i in y2:
	    	valor2+=i*h
	    for i in y3:
	    	valor3+=i*h
	    print(valor1,valor2,valor3)
	    ax.plot(x,y/valor1,'b:',x,y2/valor2,'g:',x,y3/valor3,'r:')
	    #lista[j].fill_between(x,buraco,0,facecolor='c',alpha=1)
	    ax.fill_between(x,y/valor1,0,facecolor='b',alpha=.75)
	    ax.fill_between(x,y2/valor2,0,facecolor='g',alpha=.75)
	    ax.fill_between(x,y3/valor3,0,facecolor='r',alpha=.75)

    plt.show()

if __name__ == "__main__":
    main()

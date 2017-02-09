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

    red_patch = mpatches.Patch(color='#cd0000', label='Heavy-Hole2')
    b_patch = mpatches.Patch(color='#0909d2' , label='Conduction Band')
    b2_patch = mpatches.Patch(color='#1c6bf2' , label='Light-Hole')
    g_patch = mpatches.Patch(color='#228b22' , label='Heavy-Hole1')

    bc = np.array([(+ (xk*xk*hc*hc)/(2*mgaas) + eg/2) for xk in x])
    bvhh = np.array([-hc*hc*xk*xk/(2*mhh)-eg/2 for xk in x])
    bvhh2 = np.array([ -hc*hc*xk*xk/(2*mhh)-eg/2-hc*hc*xk*xk/(2*mso) for xk in x])
    spin = np.array([ -hc*hc*xk*xk/(2*mlh)-eg/2-esp for xk in x])

    f, (ax1, ax2,ax3,ax4) = plt.subplots(1, 4, sharey=True)

    lista = [ax2,ax3,ax4]
    largurapc = [15,10,5]
    file1 = open("qb.dat","w")
    ax1.plot(x,0*x,'k--')
    ax1.fill_between(x,bc,3,facecolor='#0909d2',interpolate=True)
    ax1.fill_between(x,bvhh,bvhh2,where=bvhh>=bvhh2,facecolor='#228b22', interpolate=True)
    ax1.fill_between(x,bvhh2,spin,where=bvhh2>=spin,facecolor='#cd0000', interpolate=True)
    ax1.fill_between(x,spin,-5,facecolor='#1c6bf2',interpolate=True)
    for i in range(len(x)):
        file1.write("%s\t%s\t%s\t%s\t%s\t\n"%(x[i],bc[i],bvhh[i],bvhh2[i],spin[i]))
    file1.write("\n\n")
    ax1.grid(True)
    ax1.set_ylabel('E')
    ax1.set_title('Livre')
    for ax in range(3):
        mtxgaas = mtxg(1,mgaas,largurapc[ax])
        bc = np.array([(mtxgaas[0].real+ (xk*xk*hc*hc)/(2*mgaas) + eg/2) for xk in x])
        bc_2 = np.array([(mtxgaas[10].real+ (xk*xk*hc*hc)/(2*mgaas) + eg/2) for xk in x])
        bc_3 = np.array([(mtxgaas[20].real+ (xk*xk*hc*hc)/(2*mgaas) + eg/2) for xk in x])
        mtxbv = mtxg(-1,mhh,largurapc[ax])
        bvhh = np.array([mtxbv[0].real-hc*hc*xk*xk/(2*mhh)-eg/2 for xk in x])
        bvhh_2 = np.array([mtxbv[10].real-hc*hc*xk*xk/(2*mhh)-eg/2 for xk in x])
        bvhh_3 = np.array([mtxbv[20].real-hc*hc*xk*xk/(2*mhh)-eg/2 for xk in x])
        mtxbv2 = mtxg(-1,(mso*mhh)/(mso+mhh),largurapc[ax])
        bvhh2 = np.array([mtxbv2[0].real -hc*hc*xk*xk/(2*mhh)-eg/2-hc*hc*xk*xk/(2*mso) for xk in x])
        bvhh2_2 = np.array([mtxbv2[10].real -hc*hc*xk*xk/(2*mhh)-eg/2-hc*hc*xk*xk/(2*mso) for xk in x])
        bvhh2_3 = np.array([mtxbv2[20].real -hc*hc*xk*xk/(2*mhh)-eg/2-hc*hc*xk*xk/(2*mso) for xk in x])
        mtxspin = mtxg(-1,mlh,largurapc[ax])
        spin = np.array([mtxspin[0].real -hc*hc*xk*xk/(2*mlh)-eg/2-esp for xk in x])
        spin_2 = np.array([mtxspin[10].real -hc*hc*xk*xk/(2*mlh)-eg/2-esp for xk in x])
        spin_3 = np.array([mtxspin[20].real -hc*hc*xk*xk/(2*mlh)-eg/2-esp for xk in x])
        for i in range(len(x)):
            file1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n"%(x[i],bc[i],bc_2[i],bc_3[i],bvhh[i],bvhh_2[i],bvhh_3[i],bvhh2[i],bvhh2_2[i],bvhh2_3[i],spin[i],spin_2[i],spin_3[i]))
        file1.write("\n\n")
        lista[ax].plot(x,0*x,'k--')
        lista[ax].plot(x,bc,color='#0909d2')
        lista[ax].plot(x,bc_2,color='#0909d2',linestyle='--')
        lista[ax].plot(x,bc_3,color='#0909d2',linestyle=':')
        lista[ax].plot(x,bvhh,color='#228b22',linestyle=':')
        lista[ax].plot(x,bvhh_2,color='#228b22',linestyle='--')
        lista[ax].plot(x,bvhh_3,color='#228b22')
        lista[ax].plot(x,bvhh2,color='#cd0000',linestyle=':')
        lista[ax].plot(x,bvhh2_2,color='#cd0000',linestyle='--')
        lista[ax].plot(x,bvhh2_3,color='#cd0000')
        lista[ax].plot(x,spin,color='#1c6bf2',linestyle=':')
        lista[ax].plot(x,spin_2,color='#1c6bf2',linestyle='--')
        lista[ax].plot(x,spin_3,color='#1c6bf2')
        lista[ax].grid(True)
        lista[ax].set_title('LP='+str(largurapc[ax]))
    
    file1.close()
    ax4.set_ylim([-5,3])

    plt.legend(loc='lower center', bbox_to_anchor=(.0, 1.05),
          ncol=4, fancybox=True, shadow=True,handles=[b_patch,g_patch,red_patch,b2_patch])
    plt.show()
   

if __name__ == "__main__":
    main()

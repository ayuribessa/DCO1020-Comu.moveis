import numpy as np
import matplotlib.pyplot as plt
import math

#constantes

RAIO = 5e3
XGRID = 5*RAIO
YGRID = 6*np.sqrt(3/4)*RAIO
OFFSET = np.pi/6
HBSS = 30
PTDBM = 57
HMOB = 5
FC = 800
AHM = 3.2*(math.log10(11.75*HMOB))**2 - 4.97

def main():

    #Chamando a função que calcula os centros dos hexagonos
    centros = calcCentros(XGRID,YGRID)
    #Chamando função que plota os 7 hexagonos num grid
    MakeGrid(centros)
    plt.scatter(np.real(centros),np.imag(centros))
    #Mostrando o plot
    plt.show()

    #Criando o grid de pontos de medição para cada Bs
    xx,yy,posEachBs = calcPontosMedicao(XGRID,YGRID,centros)
    plt.scatter(np.real(centros),np.imag(centros))
     
    #Mostrando os pontos de medição em cima do grid para as 7 Erbs
    printPontosMedicao(posEachBs,centros)
    mtPotEachBsdBm = calcMatrizPotenciasBsdbm(posEachBs)
    printREMs(xx,yy,mtPotEachBsdBm,centros)


def calcPontosMedicao(XGRID,YGRID,centros):
    passo = math.ceil(RAIO/10)
    XGRID += (XGRID % passo)
    YGRID += (YGRID % passo)
    xx, yy = np.meshgrid( np.arange(0,XGRID,passo), np.arange(0,YGRID,passo) )
    posEachBs = []
    for i in range(7):
        posBs = (xx + yy*1j) - centros[i]
        posEachBs.append(posBs)
    posEachBs = np.array(posEachBs)
    return xx,yy,posEachBs

def calcMatrizPotenciasBsdbm(posEachBs):
    mtDisEachBsNorm = np.abs(posEachBs)
    np.where(mtDisEachBsNorm < RAIO, RAIO, mtDisEachBsNorm)
    mtPldb = 69.55+26.6*math.log10(FC)+(44.9-6.55*math.log10(HBSS))*np.log10(mtDisEachBsNorm/1e3) - 13.82*math.log10(HBSS) - AHM
    mtPotEachBsdBm = PTDBM - mtPldb
    return mtPotEachBsdBm

def printREMs(xx,yy,mtPotEachBsdBm,centros):
    for i in range(7):
        plt.pcolor(xx,yy, mtPotEachBsdBm[i],vmax = -115, vmin = -50)
        plt.colorbar()
        MakeGrid(centros) 
        plt.title('ERB {}'.format(i+1) )
        plt.autoscale()
        plt.show()

def printPontosMedicao(posEachBs,centros):
    for i in range(7):
        plt.scatter(posEachBs[i].real, posEachBs[i].imag)
        plt.title('ERB {}'.format(i+1))
        MakeGrid(centros - centros[i])
        plt.autoscale()
        plt.show()


def makeHex(centro):
    PontosHex = []
    for i in range(7):
        ponto = RAIO*np.exp(np.pi/3*(i-1)*1j)
        PontosHex.append(ponto)

    pontosHex = np.asarray(PontosHex)
    pontosHex += centro
    X = np.real(pontosHex)
    Y = np.imag(pontosHex)
    plt.plot(X,Y)
    #plt.scatter(centro.real,centro.imag, marker="s")
    #plt.show()

def calcCentros(XGRID,YGRID):
    Centros = [0] 
    for a in range(6):
        centro = RAIO*np.sqrt(3)*np.exp(((a-2)*np.pi/3 + OFFSET)*1j)
        Centros.append(centro)
    

    centros = np.asarray(Centros)
    centros += complex(XGRID/2,YGRID/2)
    return centros

#Função para desenhar o grid de hexagonos
def MakeGrid(centros):
    for point in centros:
        makeHex(point)
    
    
    
main()


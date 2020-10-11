
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
SIGMA_SHADOW = 8

def main():

    #1)Chamando a função que calcula os centros dos hexagonos
    centros = calcCentros(XGRID,YGRID)
    #2)Chamando função que plota os 7 hexagonos num grid
    #MakeGrid(centros)
    #Plotando os centros em cima do Grid
    #plt.scatter(centros.real , centros.imag)
    #Mostrando o plot
    #plt.show()

    #3)Criando o grid de pontos de medição para cada Bs
    xx,yy,posEachBs = calcPontosMedicao(XGRID,YGRID,centros)
    #plt.scatter(np.real(centros),np.imag(centros))
     
    #Mostrando os pontos de medição em cima do grid para as 7 Erbs
    #printPontosMedicao(posEachBs,centros)
    #4)Criando a matriz de potencias de cada ERB
    mtPotEachBsdBm = calcMatrizPotenciasBsdbm(posEachBs)
    #5)Mostra a REMs das 7 ERBS
    #printREMs(xx,yy,mtPotEachBsdBm,centros)
    calcShowMatrizPotenciasAll(mtPotEachBsdBm,xx,yy,centros)
    calcMtShadow(mtPotEachBsdBm, xx,yy,centros)
#função que calcula a malha de pontos de medição em cima do grid de hexagonos
def calcPontosMedicao(XGRID,YGRID,centros):
    passo = math.ceil(RAIO/20)
    XGRID += (XGRID % passo)
    YGRID += (YGRID % passo)
    #xx, yy = np.meshgrid( np.arange(0,XGRID,passo), np.arange(0,YGRID,passo) )
    xx, yy = np.meshgrid( np.linspace(0,XGRID,passo), np.linspace(0,YGRID,passo) )
    posEachBs = []
    for i in range(7):
        posBs = (xx + yy*1j) - centros[i]
        posEachBs.append(posBs)
    posEachBs = np.array(posEachBs)
    return xx,yy,posEachBs

def calcMtShadow(mtPotEachBsdBm, xx,yy,centros):
    MtPtShadowDbmFinal = np.NINF*np.ones((np.shape(yy)))
    dim = np.shape(yy)
    mtShadowing = SIGMA_SHADOW*np.random.randn(dim[0],dim[1])
    mtPtEachBsShadow = mtPotEachBsdBm + mtShadowing
    for i in range(7):
           MtPtShadowDbmFinal = np.maximum(MtPtShadowDbmFinal,mtPtEachBsShadow[i])
    plt.pcolor(xx,yy,MtPtShadowDbmFinal,cmap='hsv')
    plt.colorbar()
    MakeGrid(centros)
    plt.title('Todas as 7 Erbs com Shadowing')  #plt.axis('equal')
    plt.show()

def calcShowMatrizPotenciasAll(mtPotEachBsdBm,xx,yy,centros):
    MtPotenciasDbmFinal = np.NINF*np.ones(np.shape(YGRID))
    for i in range(7):
        MtPotenciasDbmFinal = np.maximum(MtPotenciasDbmFinal,mtPotEachBsdBm[i])
    plt.pcolor(xx,yy,MtPotenciasDbmFinal,cmap='hsv')
    plt.colorbar()
    MakeGrid(centros)
    #plt.axis('equal')
    plt.title('Todas as 7 Erbs sem shadowing')
    plt.show()

def calcMatrizPotenciasBsdbm(posEachBs):
    mtDisEachBsNorm = np.abs(posEachBs)
    np.where(mtDisEachBsNorm < RAIO, RAIO, mtDisEachBsNorm)
    mtPldb = 69.55+26.16*math.log10(FC)+(44.9-6.55*math.log10(HBSS))*np.log10(mtDisEachBsNorm/1e3) - 13.82*math.log10(HBSS) - AHM
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
        plt.scatter((posEachBs[i].real), (posEachBs[i].imag))
        plt.title('ERB {}'.format(i+1))
        MakeGrid(centros - centros[i])
        #plt.axis('equal')
        plt.show()



#Essa é a função que desenha o(s) hexagonos
def makeHex(centro):
    PontosHex = []
    for i in range(7):
        ponto = RAIO*np.exp(np.pi/3*i*1j)
        PontosHex.append(ponto)
    pontosHex = np.asarray(PontosHex)
    pontosHex += centro
    X = np.real(pontosHex)
    Y = np.imag(pontosHex)
    plt.plot(X,Y)
    #plt.scatter(centro.real,centro.imag, marker="s")
    #plt.show()

#Função que calcula os centros do 7 hexagonos do Grid
#Tem como entrada o tamanho do eixo x e y do Grid desejado
def calcCentros(XGRID,YGRID):
    #Centros inicia com um elemento sendo zero porque ao somar com o tamanho do grid/2 posterirmente
    #o primento elemento já vai ser o centro do hexagono central
    Centros = [0] 
    for a in range(6):
        centro = RAIO*np.sqrt(3)*np.exp((a*np.pi/3 + OFFSET)*1j)
        Centros.append(centro)
    centros = np.asarray(Centros)
    centros += complex(XGRID/2,YGRID/2)
    return centros

#Função para desenhar o grid de hexagonos
#Chama a função que desenha o hexagono o n° de vezes igual ao n° de elementos do array centros
#O array centros é um array de números complexos que são os centros dos 7 hexagonos
def MakeGrid(centros):
    for point in centros:
        makeHex(point)
    
    
    
main()


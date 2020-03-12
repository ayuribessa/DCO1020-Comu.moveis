import numpy as np
import matplotlib.pyplot as plt
import math

#constantes

RAIO = 5e3
XGRID = 5*RAIO
YGRID = 6*np.sqrt(3/4)*RAIO
OFFSET = np.pi/6

def main():


    #centro = float(input("digite o centro da do hexagono: "))   (Caso queira pedir um centro qualquer)
    #RAIO = float(input("Digite o RAIO do hexagono: "))          (Caso queria usar um centro qualquer)
    #makeHex(centro,RAIO)
    #distEntreSites= 2*np.sqrt(3/4)*RAIO
    #Definindo constantes:
   
    #Chamando a função que calcula os centros dos hexagonos
    centros = calcCentros(RAIO,OFFSET,XGRID,YGRID)
    #Chamando função que plota os 7 hexagonos num grid
    MakeGrid(RAIO,centros)
    #Mostrando o plot
    plt.show()

    #Criando o grid de pontos de medição para cada Bs
    posEachBs = calcPontosMedicao(RAIO,XGRID,YGRID,centros)
     
    #Mostrando os pontos de medição em cima do grid para as 7 Erbs
    printPontosMedicao(posEachBs,centros)


def calcPontosMedicao(RAIO,XGRID,YGRID,centros):
    passo = math.ceil(RAIO/10)
    XGRID += (XGRID % passo)
    YGRID += (YGRID % passo)
    xx, yy = np.meshgrid( np.arange(0,XGRID,passo), np.arange(0,YGRID,passo) )
    posEachBs = []
    for i in range(7):
        posBs = (xx + yy*1j) - centros[i]
        posEachBs.append(posBs)
    posEachBs = np.array(posEachBs)
    return posEachBs

def printPontosMedicao(posEachBs,centros):
    for i in range(7):
        plt.scatter(posEachBs[i].real, posEachBs[i].imag)
        plt.title('ERB {}'.format(i+1))
        MakeGrid(RAIO, centros - centros[i])
        plt.show()


def makeHex(centro,RAIO):
    PontosHex = []
    for i in range(7):
        #calcula um número complexo é o vértice do hexágono
        ponto = RAIO*np.exp(np.pi/3*(i-1)*1j)
        #ponto = RAIO*(np.cos((i-1)*np.pi/3) + np.sin((i-1)*np.pi/3)*1j)
        #Povoa a lista com esses pontos
        PontosHex.append(ponto)

    #Escalona os pontos pelo valor do centro( ajuste )
    #pontosHex = [x + centro for x in pontosHex]
    #X = [x.real for x in pontosHex]
    #Y = [y.imag for y in pontosHex]
    pontosHex = np.asarray(PontosHex)
    pontosHex += centro
    X = np.real(pontosHex)
    Y = np.imag(pontosHex)
    plt.plot(X,Y)
    #plota o centro também
    plt.scatter(centro.real,centro.imag, marker="s")
    #plt.show()

def calcCentros(RAIO, OFFSET,XGRID,YGRID):
    Centros = [0] #primeiro elemento como zero para não precisar somar o centro do hexagono do meio posteriormente
    for a in range(6):
        centro = RAIO*np.sqrt(3)*np.exp(((a-2)*np.pi/3 + OFFSET)*1j)
        Centros.append(centro)
    
    #centros = [x + complex(XGRID/2,YGRID/2) for x in centros]
    #print(centros)
    centros = np.asarray(Centros)
    centros += complex(XGRID/2,YGRID/2)
    return centros

#Função para desenhar o grid de hexagonos
def MakeGrid(RAIO,centros):
    for point in centros:
        makeHex(point,RAIO)
    
main()


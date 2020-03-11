import numpy as np
import matplotlib.pyplot as plt
import math

def main():


    #centro = float(input("digite o centro da do hexagono: "))   (Caso queira pedir um centro qualquer)
    #raio = float(input("Digite o raio do hexagono: "))          (Caso queria usar um centro qualquer)
    #makeHex(centro,raio)
    #distEntreSites= 2*np.sqrt(3/4)*raio
    #Definindo constantes:
    raio = 5e3
    xGrid = 5*raio
    yGrid = 6*np.sqrt(3/4)*raio
    offset = np.pi/6
    #Chamando a função que calcula os centros dos hexagonos
    centros = calcCentros(raio,offset,xGrid,yGrid)
    #Chamando função que plota os 7 hexagonos num grid
    MakeGrid(raio,centros)
    #Mostrando o plot
    plt.show()

    #Criando o grid de pontos de medição
    desenhaPontos(raio,xGrid,yGrid,centros)

def desenhaPontos(raio,xGrid,yGrid,centros):
    passo = math.ceil(raio/10)
    xGrid += (xGrid % passo)
    yGrid += (yGrid % passo)
    xx, yy = np.meshgrid( np.arange(0,xGrid,passo), np.arange(0,yGrid,passo) )
    posEachBs = []
    for i in range(7):
        posBs = (xx + yy*1j) - centros[i]
        posEachBs.append(posBs)
        plt.scatter(posEachBs[i].real, posEachBs[i].imag)
        plt.title('ERB {}'.format(i+1))
        MakeGrid(raio, centros - centros[i])
        plt.show()


def makeHex(centro,raio):
    PontosHex = []
    for i in range(7):
        #calcula um número complexo é o vértice do hexágono
        ponto = raio*np.exp(np.pi/3*(i-1)*1j)
        #ponto = raio*(np.cos((i-1)*np.pi/3) + np.sin((i-1)*np.pi/3)*1j)
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

def calcCentros(raio, offset,xGrid,yGrid):
    Centros = [0] #primeiro elemento como zero para não precisar somar o centro do hexagono do meio posteriormente
    for a in range(6):
        centro = raio*np.sqrt(3)*np.exp(((a-2)*np.pi/3 + offset)*1j)
        Centros.append(centro)
    
    #centros = [x + complex(xGrid/2,yGrid/2) for x in centros]
    #print(centros)
    centros = np.asarray(Centros)
    centros += complex(xGrid/2,yGrid/2)
    return centros

#Função para desenhar o grid de hexagonos
def MakeGrid(raio,centros):
    for point in centros:
        makeHex(point,raio)
    
main()
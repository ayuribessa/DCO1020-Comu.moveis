
import numpy as np
import math
import matplotlib.pyplot as plt

#constantes

PTDBM = 57
SENSIBILIDADE = -104
HMOB = 5
HBSS = 30
OFFSET = math.pi/6


def main():
    #Escolhe o modelo
    choice = getChoice()
    print('Aguarde um pouco...')
    vFreqs = np.array([800, 900, 1800, 1900, 2100])
    for freq in vFreqs:
        raio = 3000 #valor definido para ser o inicial do raio
        while True:
            flag = 0
            xgrid = 5*raio
            ygrid = 6*np.sqrt(3/4)*raio
            #calcula os centros das Bss
            centros = calcCentros(raio, xgrid, ygrid)
            #calcula o grid de pontos de medição
            xx, yy, ygrid = calcPontosMedicao(raio, centros, xgrid, ygrid)
            #calcula a posicao de cada Bss
            posEachBs = calcPosicoes(centros, xx, yy)
            #calcula a taxa de outage
            flag = calcOutage(choice,raio, freq, posEachBs, ygrid)
            if(flag != 0):
                break
            raio = raio + 5


def calcCentros(raio, xgrid, ygrid):

    Centros = [0]
    for a in range(6):
        centro = raio*np.sqrt(3)*np.exp((a*np.pi/3 + OFFSET)*1j)
        Centros.append(centro)
    centros = np.asarray(Centros)
    centros += complex(xgrid/2, ygrid/2)
    return centros


def calcPontosMedicao(raio, centros, xgrid, ygrid):

    passo = math.ceil(raio/50)
    xgrid += (xgrid % passo)
    ygrid = math.ceil(ygrid + (ygrid % passo))
    xx, yy = np.meshgrid(np.linspace(0, xgrid, passo),
                         np.linspace(0, ygrid, passo))
    return xx, yy, ygrid


def calcPosicoes(centros, xx, yy):

    posEachBs = []
    for i in range(7):
        posBs = (xx + yy*1j) - centros[i]
        posEachBs.append(posBs)
    posEachBs = np.array(posEachBs)
    return posEachBs


def calcMatrizPotenciasAll(mtPotEachBsdBm, ygrid):

    MtPotenciasDbmFinal = np.NINF*np.ones(np.shape(ygrid))
    for i in range(7):
        MtPotenciasDbmFinal = np.maximum(MtPotenciasDbmFinal, mtPotEachBsdBm[i])
    return MtPotenciasDbmFinal



def calcOutage(choice,raio, freq, posEachBs, ygrid):
   
    #ahm = ( 1.1*math.log10(fc) - 0.7) * HMOB - (1.56*math.log10(fc)-0.8) ahm para  cost na literatura
    ahm = 3.2*np.log10(11.75*HMOB)**2 - 4.97
    if  choice == 1: 
        mtPotEachBsdBm = calcMtPtBsdbmOkHata(raio, posEachBs, freq, ahm)
    else:
        mtPotEachBsdBm = calcMtPtBsdbmCost231(raio, posEachBs, freq, ahm)
    MtPotenciasDbmFinal = calcMatrizPotenciasAll(mtPotEachBsdBm, ygrid)
    x = ((MtPotenciasDbmFinal < SENSIBILIDADE).sum())
    y = np.size(MtPotenciasDbmFinal)
    outage = x/y

    if outage >= 0.0999998:
        printValues(freq,raio,outage)
        return 1
    return 0

def printValues(freq,raio,outage):
    print('------------------------------------------------------')
    print('Com a frequencia da portadora: {:d} MHz'.format(freq))
    print('O maior raio que deu taxa de outage com até 10% foi: {:.2f} m '.format(raio))
    print('Com outage calculado de: {:.4%}'.format(outage))

def calcMtPtBsdbmOkHata(raio, posEachBs, fc, ahm):
        
    mtDisEachBsNorm = np.abs(posEachBs)
    np.where(mtDisEachBsNorm < raio, raio, mtDisEachBsNorm)
    mtPldb = 69.55 + 26.16 * math.log10(fc) + (44.9-6.55*math.log10(HBSS))*np.log10(mtDisEachBsNorm/1e3) - \
        13.82*math.log10(HBSS) - ahm
    mtPotEachBsdBm = PTDBM - mtPldb
    return mtPotEachBsdBm

def calcMtPtBsdbmCost231(raio, posEachBs, fc,ahm):

    mtDisEachBsNorm = np.abs(posEachBs)
    np.where(mtDisEachBsNorm < raio, raio, mtDisEachBsNorm)
    mtPldb = 46.3 + 33.9*math.log10(fc) - 13.82*math.log10(HBSS) - ahm + \
        (44.9 - 6.55*math.log10(HBSS))*np.log10(mtDisEachBsNorm/1e3) + 3
    mtPotEachBsdBm = PTDBM - mtPldb  
    return mtPotEachBsdBm


def getChoice():
    while True:
        print('Deseja saber outage para Modelo de Okomura Hata (1) ou ' +
        'Cost 231 (2) ? ')
        choice = int(input('Digite <1> ou <2> '))
        if choice == 1 or choice == 2:
            return choice
        print('Só são válidos <1> ou <2>')
main()

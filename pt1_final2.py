import numpy as np
import math
import matplotlib.pyplot as plt


PTDBM = 57
SENSIBILIDADE = -104
HMOB = 5
HBSS = 30
AHM = 3.2*np.log10(11.75*HMOB)**2 - 4.97
OFFSET = np.pi/6

def main():
    vFreqs = np.array([800, 900, 1800, 1900, 2100])
    for freq in vFreqs:
        calcoutage2(freq)
  

    
def calcoutage2(freq):

    raio = 1000
    while True:
        flag = 0
        xgrid = 5*raio
        ygrid = 6*np.sqrt(3/4)*raio
        centros = calcCentros(raio,xgrid, ygrid)
        xx, yy, ygrid = calcPontosMedicao(raio,centros, xgrid, ygrid)
        posEachBs = calcPosicoes(centros, xx, yy)
        flag = calcOutage(raio,freq, posEachBs, ygrid)
        if(flag != 0):
            break
        raio = raio + 1




def calcCentros(raio,xgrid,ygrid):
   
    Centros = [0] 
    for a in range(6):
        centro = raio*np.sqrt(3)*np.exp((a*np.pi/3 + OFFSET)*1j)
        Centros.append(centro)
    

    centros = np.asarray(Centros)
    centros += complex(xgrid/2,ygrid/2)
    return centros

def calcPontosMedicao(raio,centros,xgrid,ygrid):
   
    passo = math.ceil(raio/50)
    xgrid +=  (xgrid % passo)
    ygrid = math.ceil( ygrid + (ygrid % passo))
    xx, yy = np.meshgrid( np.linspace(0,xgrid,passo), np.linspace(0,ygrid,passo) )
    return xx,yy,ygrid

def calcPosicoes(centros,xx,yy):
    posEachBs = []
    for i in range(7):
        posBs = (xx + yy*1j) - centros[i]
        posEachBs.append(posBs)
    posEachBs = np.array(posEachBs)
    return posEachBs

def calcMatrizPotenciasAll(mtPotEachBsdBm,ygrid):
    MtPotenciasDbmFinal = np.NINF*np.ones(np.shape(ygrid))
    for i in range(7):
        MtPotenciasDbmFinal = np.maximum(MtPotenciasDbmFinal,mtPotEachBsdBm[i])
    return MtPotenciasDbmFinal

def calcOutage(raio,freq,posEachBs,ygrid):
    mtPotEachBsdBm = calcMatrizPotenciasBsdbm(raio,posEachBs,freq)
    MtPotenciasDbmFinal = calcMatrizPotenciasAll(mtPotEachBsdBm,ygrid)
    x = ((MtPotenciasDbmFinal < SENSIBILIDADE).sum())
    y = np.size(MtPotenciasDbmFinal)
    outage = x/y

    if outage >= 0.099999:
        print('com a frequencia da portadora: {:.2f} mHz'.format(freq))
        print('O maior raio que deu taxa de outage com até 10%% foi: {:.2f} '.format(raio))
        print('Com outage de: {:.2%}'.format(outage))
        return 1 
    
    return 0
    

def calcMatrizPotenciasBsdbm(raio,posEachBs,fc):
    mtDisEachBsNorm = np.abs(posEachBs)
    np.where(mtDisEachBsNorm < raio, raio, mtDisEachBsNorm)
    mtPldb = 69.55+26.6*math.log10(fc)+(44.9-6.55*math.log10(HBSS))*np.log10(mtDisEachBsNorm/1e3) - 13.82*math.log10(HBSS) - AHM
    mtPotEachBsdBm = PTDBM - mtPldb
    return mtPotEachBsdBm

main()

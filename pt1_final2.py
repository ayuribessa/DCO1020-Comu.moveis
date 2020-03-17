import numpy as np
import math
import matplotlib.pyplot as plt


#constantes
RAIO = 10e3
PTDBM = 57
SENSIBILIDADE = -104
HMOB = 5
HBSS = 30
AHM = 3.2*np.log10(11.75*HMOB)**2 - 4.97
OFFSET = np.pi/6
YGRID = 6*np.sqrt(3/4)*RAIO 

def main():
    

    xgrid = 5*RAIO
    ygrid = 6*np.sqrt(3/4)*RAIO 
    
    centros = calcCentros(xgrid,ygrid)
    xx,yy,ygrid = calcPontosMedicao(centros,xgrid,ygrid)
    posEachBs =  calcPosicoes(centros,xx,yy)
    #fc = float(input("Digite fc: "))
    vFreqs = np.array([800,900,1800,1900,2100])
    calcOutage(vFreqs,posEachBs,ygrid)
    




def calcCentros(xgrid,ygrid):
    #Centros inicia com um elemento sendo zero porque ao somar com o tamanho do grid/2 posterirmente
    #o primento elemento já vai ser o centro do hexagono central
    Centros = [0] 
    for a in range(6):
        centro = RAIO*np.sqrt(3)*np.exp((a*np.pi/3 + OFFSET)*1j)
        Centros.append(centro)
    

    centros = np.asarray(Centros)
    centros += complex(xgrid/2,ygrid/2)
    return centros

def calcPontosMedicao(centros,xgrid,ygrid):
   
    passo = math.ceil(RAIO/50)
    xgrid +=  (xgrid % passo)
    ygrid = math.ceil( ygrid + (ygrid % passo))
    #xx, yy = np.meshgrid( np.arange(0,xgrid,passo), np.arange(0,ygrid,passo) )
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
    MtPotenciasDbmFinal = np.NINF*np.ones(np.shape(YGRID))
    for i in range(7):
        MtPotenciasDbmFinal = np.maximum(MtPotenciasDbmFinal,mtPotEachBsdBm[i])
    return MtPotenciasDbmFinal

def calcOutage(vFreqs,posEachBs,ygrid):
    for freq in vFreqs:
        mtPotEachBsdBm = calcMatrizPotenciasBsdbm(posEachBs,freq)
        MtPotenciasDbmFinal = calcMatrizPotenciasAll(mtPotEachBsdBm,ygrid)
        #z = np.size(x)
        #y = np.size(MtPotenciasDbmFinal)
        outage = 100*np.size(MtPotenciasDbmFinal[MtPotenciasDbmFinal < SENSIBILIDADE])/np.size(MtPotenciasDbmFinal)
        #outage = 100*np.size(MtPotenciasDbmFinal[MtPotenciasDbmFinal < SENSIBILIDADE])/np.size(MtPotenciasDbmFinal)
        #outage = 100 * (z/y)
        print("Freq da portadora = ",freq)
        print("outage = " ,outage)
        

def calcMatrizPotenciasBsdbm(posEachBs,fc):
    mtDisEachBsNorm = np.abs(posEachBs)
    np.where(mtDisEachBsNorm < RAIO, RAIO, mtDisEachBsNorm)
    mtPldb = 69.55+26.6*math.log10(fc)+(44.9-6.55*math.log10(HBSS))*np.log10(mtDisEachBsNorm/1e3) - 13.82*math.log10(HBSS) - AHM
    mtPotEachBsdBm = PTDBM - mtPldb
    return mtPotEachBsdBm

main()
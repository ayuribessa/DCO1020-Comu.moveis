import numpy as np
import matplotlib.pyplot as plt
import math



#constantes
RAIO = 200
SHAD = 50
PASSO = 10
XGRID0RI = 5*RAIO
YGRID0RI = 6*np.sqrt(3/4)*RAIO
OFFSET = np.pi/6
HBSS = 30
PTDBM = 57
HMOB = 5
FC = 800
AHM = 3.2*(math.log10(11.75*HMOB))**2 - 4.97
SIGMA_SHADOW = 8
# alphas_corr = 0


def MakeGrid(centros):
    for point in centros:
        makeHex(point)

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

def fCorrShadowing(alphas_corr, mtPontosMedicao):

    dimXS = math.ceil(XGRID0RI + (XGRID0RI % SHAD))
    dimYS = math.ceil(YGRID0RI + (YGRID0RI % SHAD))
    mtPosxShad, mtPosyShad = np.meshgrid(np.arange(0,dimXS + SHAD,SHAD), np.arange(0,dimYS + SHAD ,SHAD))
    
    mtShadowingSamples = []
    for erb in range(8):
        ShadowingSamples = SIGMA_SHADOW*np.random.randn(np.shape(mtPosyShad)[0],np.shape(mtPosyShad)[1])
        mtShadowingSamples.append(ShadowingSamples)
    mtShadowingSamples = np.asarray(mtShadowingSamples)

    sizeL, sizeC = np.shape(mtPontosMedicao) #pt4

    mtShadowingCorr = np.empty([7,sizeL,sizeC])
    for linha in range(sizeL ): #pt4
        for coluna in range(sizeC  ): #pt4
            dShadPoint = mtPontosMedicao[linha,coluna] #pt4

            dxIndexP1 = np.real(dShadPoint)/SHAD
            dyIndexP1 = np.imag(dShadPoint)/SHAD

            if dxIndexP1 % 1 == 0 and dyIndexP1 % 1 == 0:
                dxIndexP1 = math.floor(dxIndexP1) + 1
                dyIndexP1 = math.floor(dyIndexP1) + 1
              
                shadowingC = mtShadowingSamples[7][dyIndexP1 -1][dxIndexP1 -1 ] #pt3
                # amostra do sombreamento de cada ERB
                for i in range(7):
                    mtShadowingERB = mtShadowingSamples[i][dyIndexP1 -1][dxIndexP1 -1]
                    mtShadowingCorr[i][linha ][coluna ] = np.sqrt(alphas_corr)*shadowingC + np.sqrt(1-alphas_corr)*mtShadowingERB

            else:
                dxIndexP1 = math.floor(dxIndexP1) + 1
                dyIndexP1 = math.floor(dyIndexP1) + 1
                if dxIndexP1 == np.shape(mtPosyShad)[0] and dyIndexP1 == np.shape(mtPosyShad)[1]:
                    dxIndexP2 = dxIndexP1-1
                    dyIndexP2 = dyIndexP1
                    dxIndexP4 = dxIndexP1-1
                    dyIndexP4 = dyIndexP1-1
                    dxIndexP3 = dxIndexP1
                    dyIndexP3 = dyIndexP1-1        
    
                elif dyIndexP1 == np.shape(mtPosyShad)[0]:
                    dxIndexP2 = dxIndexP1+1
                    dyIndexP2 = dyIndexP1
                    dxIndexP4 = dxIndexP1+1
                    dyIndexP4 = dyIndexP1-1
                    dxIndexP3 = dxIndexP1
                    dyIndexP3 = dyIndexP1-1

                elif dxIndexP1 == np.shape(mtPosyShad)[1]:
                    dxIndexP2 = dxIndexP1-1
                    dyIndexP2 = dyIndexP1
                    dxIndexP4 = dxIndexP1-1
                    dyIndexP4 = dyIndexP1+1
                    dxIndexP3 = dxIndexP1
                    dyIndexP3 = dyIndexP1+1

                else:
                    dxIndexP2 = dxIndexP1+1
                    dyIndexP2 = dyIndexP1
                    dxIndexP4 = dxIndexP1+1
                    dyIndexP4 = dyIndexP1+1
                    dxIndexP3 = dxIndexP1
                    dyIndexP3 = dyIndexP1+1


                distX = (dShadPoint.real%SHAD)/SHAD
                distY = (dShadPoint.imag%SHAD)/SHAD
    
                
                    #----------------------------------------END PT.2--------------------------------------------------
                #ajuste do desvio padrão devido a regressão linear
                stdNormalFactor = np.sqrt( (1 - 2*distY + 2*(distY**2))*(1 - 2*distX + 2*(distX**2))) #pt3
                #amostras do sombreamento para os 4 pontos de gradedXIndexP1);
                Sample1 = mtShadowingSamples[7][dyIndexP1 -1,dxIndexP1 -1] #pt3
                Sample2 = mtShadowingSamples[7][dyIndexP2 -1,dxIndexP2 -1 ] #pt3
                Sample3 = mtShadowingSamples[7][dyIndexP3 -1,dxIndexP3 -1] #pt3
                Sample4 = mtShadowingSamples[7][dyIndexP4 -1,dxIndexP4 -1] #pt3
                shadowingC = ((1-distY)*(Sample1*(1-distX) + Sample2*distX) +  #pt4 #pt4
                               distY*(Sample3*(1 - distX) + Sample4*distX))/stdNormalFactor #pt4 #pt4
                for i in range(7):
                    Sample1 = mtShadowingSamples[i][dyIndexP1 -1,dxIndexP1 -1] #pt3
                    Sample2 = mtShadowingSamples[i][dyIndexP2 -1,dxIndexP2 -1] #pt3
                    Sample3 = mtShadowingSamples[i][dyIndexP3 -1,dxIndexP3 -1] #pt3
                    Sample4 = mtShadowingSamples[i][dyIndexP4 -1 ,dxIndexP4-1 ] #pt3
                    shadowingERB = ((1-distY)*(Sample1*(1-distX) + Sample2*distX) +  #pt4 #pt4
                                distY*(Sample3*(1 - distX) + Sample4*distX))/stdNormalFactor #pot4 #pt4
                    mtShadowingCorr[i][linha ][coluna ] = np.sqrt(alphas_corr)*shadowingC + np.sqrt(1-alphas_corr)*mtShadowingERB
    


    return mtShadowingCorr

def main():
    #calcula posição dos centros
    Centros = [0] 
    for a in range(6):
        centro = RAIO*np.sqrt(3)*np.exp((a*np.pi/3 + OFFSET)*1j)
        Centros.append(centro)
    centros = np.asarray(Centros)
    centros += complex(XGRID0RI/2, YGRID0RI /2)
    #calcula posição dos grids de medição
    dimX = math.ceil(XGRID0RI + (XGRID0RI % PASSO))
    dimY = math.ceil(YGRID0RI + (YGRID0RI % PASSO))
    mtPosX, mtPosY = np.meshgrid(np.arange(0,dimX + PASSO ,PASSO), np.arange(0,dimY + PASSO,PASSO) )
    mtPontosMedicao = mtPosX + mtPosY*1j
    #cria as matrizes 
    mtPowerFinalDbm  = np.NINF*np.ones(np.shape(mtPosY))
    mtPowerFinalShadDbm  = np.NINF*np.ones(np.shape(mtPosY))
    mtPowerFinalShadCorrDbm  = np.NINF*np.ones(np.shape(mtPosY))
    print("Sorteando 20 valores do coeficiente de correlação de sombreamento" +
           "entre 0 e 1: \n")
    alphas_corr = np.round(np.random.rand(20),decimals=1)

    for j in range(20):
        mtShadowingCorr = fCorrShadowing(alphas_corr[j], mtPontosMedicao)
        print(f"Para coef.correlação de sombreamento '{alphas_corr[j]:.1f}', o desvio-padrão das amostras" +
                f" de sombreamentos é: '{np.std(mtShadowingCorr):.1f}'")
    
    
main()
    
    
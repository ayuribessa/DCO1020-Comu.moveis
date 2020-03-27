import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

#constantes
#tirei 0 -1 na indexação da matrizshadowblahblah

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
ALPHA_CORR = 0.5


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

def fCorrShadowing(mtPontosMedicao):

    dimXS = math.ceil(XGRID0RI + (XGRID0RI % SHAD))
    dimYS = math.ceil(YGRID0RI + (YGRID0RI % SHAD))
    mtPosxShad, mtPosyShad = np.meshgrid(np.arange(0,dimXS + SHAD,SHAD), np.arange(0,dimYS + SHAD ,SHAD))
    mtPosShad = mtPosxShad + mtPosyShad*1j

    mtShadowingSamples = []
    for erb in range(8):
        ShadowingSamples = SIGMA_SHADOW*np.random.randn(np.shape(mtPosyShad)[0],np.shape(mtPosyShad)[1])
        mtShadowingSamples.append(ShadowingSamples)
    mtShadowingSamples = np.asarray(mtShadowingSamples)

    sizeL, sizeC = np.shape(mtPontosMedicao) #pt4

    mtShadowingCorr = np.empty([7,sizeL,sizeC])
    for linha in range(sizeL -1): #pt4
        for coluna in range(sizeC -1 ): #pt4
            dShadPoint = mtPontosMedicao[linha,coluna] #pt4

            dxIndexP1 = np.real(dShadPoint)/SHAD
            dyIndexP1 = np.imag(dShadPoint)/SHAD

            if dxIndexP1 % 1 == 0 and dyIndexP1 % 1 == 0:
                dxIndexP1 = math.floor(dxIndexP1) + 1
                dyIndexP1 = math.floor(dyIndexP1) + 1
              
                shadowingC = mtShadowingSamples[7][dyIndexP1 - 1][dxIndexP1 - 1] #pt3
                for i in range(7):
                    mtShadowingERB = mtShadowingSamples[i][dyIndexP1-1][dxIndexP1-1]
                    #linha -1 e col -1 em mtshadowingcorr?
                    mtShadowingCorr[i][linha ][coluna ] = np.sqrt(ALPHA_CORR)*shadowingC + np.sqrt(1-ALPHA_CORR)*mtShadowingERB

            else:
                dxIndexP1 = math.floor(dxIndexP1) + 1
                dyIndexP1 = math.floor(dyIndexP1) + 1
                if dxIndexP1 == np.shape(mtPosyShad)[1] and dyIndexP1 == np.shape(mtPosyShad)[0]:
                    dxIndexP2 = dxIndexP1-1
                    dyIndexP2 = dyIndexP1
                    dxIndexP4 = dxIndexP1-1
                    dyIndexP4 = dyIndexP1-1
                    dxIndexP3 = dxIndexP1
                    dyIndexP3 = dyIndexP1-1        
    
                elif dyIndexP1 == np.shape(mtPosyShad)[1]:
                    dxIndexP2 = dxIndexP1+1
                    dyIndexP2 = dyIndexP1
                    dxIndexP4 = dxIndexP1+1
                    dyIndexP4 = dyIndexP1-1
                    dxIndexP3 = dxIndexP1
                    dyIndexP3 = dyIndexP1-1

                elif dxIndexP1 == np.shape(mtPosyShad)[0]:
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
    
                # print('X = {:.2f} e Y = {:.2f}'.format(distX,distY))
                    #----------------------------------------END PT.2--------------------------------------------------
                #ajuste do desvio padrão devido a regressão linear
                stdNormalFactor = np.sqrt( (1 - 2*distY + 2*(distY**2))*(1 - 2*distX + 2*(distX**2))) #pt3
                #amostras do sombreamento para os 4 pontos de gradedXIndexP1);
                Sample1 = mtShadowingSamples[7][dyIndexP1 - 1,dxIndexP1 - 1] #pt3
                Sample2 = mtShadowingSamples[7][dyIndexP2 - 1,dxIndexP2 - 1] #pt3
                Sample3 = mtShadowingSamples[7][dyIndexP3 - 1,dxIndexP3 - 1] #pt3
                Sample4 = mtShadowingSamples[7][dyIndexP4 - 1,dxIndexP4 - 1] #pt3
                shadowingC = ((1-distY)*(Sample1*(1-distX) + Sample2*distX) +  #pt4 #pt4
                               distY*(Sample3*(1 - distX) + Sample4*distX))/stdNormalFactor #pt4 #pt4
                for i in range(7):
                    Sample1 = mtShadowingSamples[i][dyIndexP1 - 1,dxIndexP1 - 1] #pt3
                    Sample2 = mtShadowingSamples[i][dyIndexP2 - 1,dxIndexP2 - 1] #pt3
                    Sample3 = mtShadowingSamples[i][dyIndexP3 - 1,dxIndexP3 - 1] #pt3
                    Sample4 = mtShadowingSamples[i][dyIndexP4 - 1,dxIndexP4 - 1] #pt3
                    shadowingERB = ((1-distY)*(Sample1*(1-distX) + Sample2*distX) +  #pt4 #pt4
                                distY*(Sample3*(1 - distX) + Sample4*distX))/stdNormalFactor #pot4 #pt4
                    mtShadowingCorr[i][linha ][coluna ] = np.sqrt(ALPHA_CORR)*shadowingC + np.sqrt(1-ALPHA_CORR)*mtShadowingERB

    return mtShadowingCorr

Centros = [0] 
for a in range(6):
    centro = RAIO*np.sqrt(3)*np.exp((a*np.pi/3 + OFFSET)*1j)
    Centros.append(centro)
centros = np.asarray(Centros)
centros += complex(XGRID0RI/2, YGRID0RI /2)

dimX = math.ceil(XGRID0RI + (XGRID0RI % PASSO))
dimY = math.ceil(YGRID0RI + (YGRID0RI % PASSO))
mtPosX, mtPosY = np.meshgrid(np.arange(0,dimX + PASSO ,PASSO), np.arange(0,dimY + PASSO,PASSO) )
mtPontosMedicao = mtPosX + mtPosY*1j

mtPowerFinalDbm  = np.NINF*np.ones(np.shape(mtPosY))
mtPowerFinalShadDbm  = np.NINF*np.ones(np.shape(mtPosY))
mtPowerFinalShadCorrDbm  = np.NINF*np.ones(np.shape(mtPosY))
mtShadowingCorr = fCorrShadowing(mtPontosMedicao)
mtPowerEachBssShadowCorrDbm = np.empty((7,np.shape(mtPosY)[0],np.shape(mtPosY)[1]))
dim = np.shape(mtPosY)
for i in range(7):
    mtPosEachBs = mtPontosMedicao - centros[i]
    mtDistEachBs = np.abs(mtPosEachBs)
    np.where(mtDistEachBs < PASSO, PASSO, mtDistEachBs)
    mtPldb = 69.55+26.16*math.log10(FC)+(44.9-6.55*math.log10(HBSS))*np.log10(mtDistEachBs/1e3) - 13.82*math.log10(HBSS) - AHM
    mtShadowing = SIGMA_SHADOW*np.random.randn(dim[0],dim[1])
    mtPowerEachBsDbm = PTDBM - mtPldb
    mtPowerEachBsShadowDbm =  PTDBM - mtPldb + mtShadowing
    mtPowerFinalDbm = np.maximum(mtPowerFinalDbm,mtPowerEachBsDbm)
    mtPowerFinalShadDbm = np.maximum(mtPowerFinalShadDbm,mtPowerEachBsShadowDbm)
    mtPowerEachBssShadowCorrDbm[i] = PTDBM - mtPldb + mtShadowingCorr[i]
    mtPowerFinalShadCorrDbm = np.maximum(mtPowerFinalShadCorrDbm,mtPowerEachBssShadowCorrDbm[i])

# for i in range(7):
#     mtPowerFinalShadCorrDbm = np.maximum(mtPowerFinalShadCorrDbm,mtPowerEachBssShadowCorrDbm[i])


# plt.pcolor(mtPosX,mtPosY,mtPowerFinalDbm,cmap='hsv',vmax=10, vmin=-50)
# plt.colorbar()
# MakeGrid(centros) 
# plt.show()

# plt.pcolor(mtPosX,mtPosY,mtPowerFinalShadDbm,cmap='hsv',vmax=30,vmin=-60)
# plt.colorbar()
# MakeGrid(centros) 
# plt.show()

plt.pcolor(mtPosX,mtPosY,mtPowerFinalShadCorrDbm,cmap='hsv',vmax=20, vmin=-60)
plt.colorbar()
MakeGrid(centros) 
plt.show()



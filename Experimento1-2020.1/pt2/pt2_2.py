

import numpy as np
import matplotlib.pyplot as plt
import math

#constantes

RAIO = 200 
SHAD = 50
PASSO = 7
XGRID0RI = 5*RAIO
YGRID0RI = 6*np.sqrt(3/4)*RAIO
OFFSET = np.pi/6
HBSS = 30
PTDBM = 57
HMOB = 5
FC = 800
AHM = 3.2*(math.log10(11.75*HMOB))**2 - 4.97
SIGMA_SHADOW = 8

dimX = math.ceil(XGRID0RI + (XGRID0RI % PASSO))
dimY = math.ceil(YGRID0RI + (YGRID0RI % PASSO))
mtPosX, mtPosY = np.meshgrid(np.arange(0,dimX + PASSO ,PASSO), np.arange(0,dimY + PASSO,PASSO) )
mtPontosMedicao = mtPosX + mtPosY*1j

# dShadPoint = mtPontosMedicao[11,11]
#faltou testar esses outros :
dShadPoint = mtPontosMedicao[1,1] #retorna 0,0 nas D. de correlação 
# dShadPoint = mtPontosMedicao[0,(np.shape(mtPontosMedicao)[1])-1]
# dShadPoint = mtPontosMedicao[np.shape(mtPontosMedicao)[0]-1, 0]
# dShadPoint = mtPontosMedicao[ np.shape(mtPontosMedicao)[0] -1,  np.shape(mtPontosMedicao)[1] -1]

dimXS = math.ceil(XGRID0RI + (XGRID0RI % SHAD))
dimYS = math.ceil(YGRID0RI + (YGRID0RI % SHAD))
mtPosxShad, mtPosyShad = np.meshgrid(np.arange(0,dimXS + SHAD,SHAD), np.arange(0,dimYS + SHAD,SHAD))
mtPosShad = mtPosxShad + mtPosyShad*1j

mtShadowingSamples = SIGMA_SHADOW*np.random.randn(np.shape(mtPosyShad)[0],np.shape(mtPosyShad)[1])

dxIndexP1 = np.real(dShadPoint)/SHAD
dyIndexP1 = np.imag(dShadPoint)/SHAD

# x = np.shape(mtPosyShad)[0] 
# y = np.shape(mtPosyShad)[1] 
# z = np.shape(mtPosyShad) 
if dxIndexP1 % 1 == 0 and dyIndexP1 % 1 == 0:
    dxIndexP1 = math.floor(dxIndexP1) + 1
    dyIndexP1 = math.floor(dyIndexP1) + 1
    plt.plot(mtPosShad[dxIndexP1-1,dyIndexP1-1].imag, mtPosShad[dxIndexP1-1,dyIndexP1-1].real, 'g*')
    plt.plot(np.imag(mtPosShad),np.real(mtPosShad),'k.',fillstyle = 'none')
    plt.show()
    print("O ponto é um ponto de grade")
else:
    dxIndexP1 = math.floor(dxIndexP1) + 1
    dyIndexP1 = math.floor(dyIndexP1) + 1
    if dxIndexP1 == np.shape(mtPosyShad)[1] and dyIndexP1 == np.shape(mtPosyShad)[0]:
    # if dxIndexP1 == mtPosyShad[1][0] and dyIndexP1 == mtPosyShad[1][0]:
        dxIndexP2 = dxIndexP1-1
        dyIndexP2 = dyIndexP1
        dxIndexP4 = dxIndexP1-1
        dyIndexP4 = dyIndexP1-1
        dxIndexP3 = dxIndexP1
        dyIndexP3 = dyIndexP1-1        
    
    # elif dyIndexP1 == mtPosyShad[0][0]:
    elif dyIndexP1 == np.shape(mtPosyShad)[0]:
        dxIndexP2 = dxIndexP1+1
        dyIndexP2 = dyIndexP1
        dxIndexP4 = dxIndexP1+1
        dyIndexP4 = dyIndexP1-1
        dxIndexP3 = dxIndexP1
        dyIndexP3 = dyIndexP1-1

    # elif dxIndexP1 == mtPosyShad[1][0]:
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

    
    plt.plot(np.imag(mtPosShad),np.real(mtPosShad),'k.',fillstyle = 'none')
    plt.plot(np.imag(dShadPoint),np.real(dShadPoint),'sr')
    mt4points = np.array([ mtPosShad[dyIndexP1-1, dxIndexP1-1], mtPosShad[dyIndexP2-1, dxIndexP2-1],
                  mtPosShad[dyIndexP3-1, dxIndexP3-1], mtPosShad[dyIndexP4-1, dxIndexP4-1]])

    plt.plot((mt4points.imag),(mt4points.real),'b*')
    #Customizando eixos para dar zoom no ponto:
    # plt.axis([-2* SHAD + mtPosShad[dyIndexP3-1][dxIndexP3-1].real, 2*SHAD+ mtPosShad[dyIndexP4-1][dxIndexP4-1].real,
    #         -2 * SHAD + mtPosShad[dyIndexP3-1][dxIndexP3-1].imag, 2* SHAD + mtPosShad[dyIndexP1-1][dxIndexP1-1].imag])
    #plt.axis('equal')
    
    plt.show()

distX = (dShadPoint.real%SHAD)/SHAD
distY = (dShadPoint.imag%SHAD)/SHAD
    
print('X = {:.2f} e Y = {:.2f}'.format(distX,distY))
    #----------------------------------------END PT.2--------------------------------------------------
import canal
import numpy as np
import matplotlib.pyplot as plt

# d0 = d0
# p0 = p0
# nPoints = nPoints 
# totalLenght = totalLenght
# n = n
# sigma = sigma
# shadowingWindow = shadowingWindow 
# m = m  
# txPower = txPower
# nCDF = nCDF
# chFileName = chFileName


myChannel = canal.Channel(d0 = 5, p0 = 0, nPoints = 50000, totalLenght = 100, n = 4, sigma = 6,
                        shadowingWindow = 200, m = 4, txPower = 0, nCDF = 40, chFileName = 'prx_sintetico')

vtDist, vt_pathLoss, vtShadowCorr, vtFading, vtPrx = myChannel.GeraVetores()
vtPrxmW = 10**(vtPrx/10)

JanelaEstim = 100
dMeiaJanela = int(JanelaEstim/2)
nSamples = len(vtPrxmW)
#FIXME find a pythonic way to fill the following arrays and build the loop:
vtDesLarga = np.empty(nSamples-JanelaEstim)
vtDesPequenaEst = np.empty(nSamples-JanelaEstim)
ij = 0
for ik in range(dMeiaJanela,nSamples - dMeiaJanela):
    vtDesLarga[ij] = 10*np.log10(np.mean(vtPrxmW[ik-dMeiaJanela:ik+dMeiaJanela]))
    vtDesPequenaEst[ij] = vtPrx[ik] - vtDesLarga[ij]
    ij += 1

indexes = np.arange(dMeiaJanela,nSamples - dMeiaJanela)
#Valoes de vTpRx começando de 50
vtPtrxmWNew = 10 ** (vtPrx[indexes]/10)
#Todos valores de vtDesLarga
desLarga_lin = 10 ** (vtDesLarga[:len(indexes)]/10)
envNormal = np.sqrt(vtPtrxmWNew) / np.sqrt(desLarga_lin)

#ajuste do tamanho dos vetores devido a filtragem
vtDistEst = vtDist[dMeiaJanela : nSamples - dMeiaJanela]
vtPrxdBm = vtPrx[dMeiaJanela : nSamples - dMeiaJanela]

#cálculo da reta de perda de percurso
vtDistLog = np.log10(vtDist)
vtDistLogEst = np.log10(vtDistEst)
#Coeficientes da reta que melhor caracteriza P.L
CoefReta = np.polyfit(vtDistLogEst,vtPrxdBm,1)
#Expoente da perda de percurso
dNEst = -CoefReta[0]/10

print(dNEst)

#perda de percurso estimada
vtPathLossEst = np.polyval(CoefReta, vtDistLogEst)
#sombreamento
vtShadowCorrEst = vtDesLarga - vtPathLossEst
stdShad  = np.std(vtShadowCorrEst)
meanShad = np.mean(vtShadowCorrEst)

print(stdShad)
print(meanShad)

vtPathLossEst = -vtPathLossEst
vtPrxEst = myChannel.txPower - vtPathLossEst + vtShadowCorrEst + vtDesPequenaEst

#estimação da CDF do descanecimento de pequena escala   
vtn = np.arange(myChannel.nCDF)
xCDF = 1.2 ** (vtn) * 0.01

#cálculo da CDF
den = 0
cdffn = np.zeros(myChannel.nCDF)
for ik in range(myChannel.nCDF):
    for ij in range(len(envNormal)):
        if envNormal[ij] <= xCDF[ik]:
            den += 1
        cdffn[ik] = cdffn[ik] + den
        den = 0

#estrutura do histograma
xccdfEst = 20*np.log10(xCDF)
yccdfEst = cdffn/cdffn[-1]

# potencia recebida canal completo
plt.plot(vtDistLogEst, vtPrxEst)
plt.plot(vtDistLogEst, myChannel.txPower -vtPathLossEst)
plt.plot(vtDistLogEst, myChannel.txPower - vtPathLossEst + vtShadowCorrEst)
plt.show()

#path loss estimado x original
plt.plot(vtDistLogEst, -vtPathLossEst)
plt.plot(vtDistLog, - vt_pathLoss)
plt.show()

#sombreamento estimado x original
plt.plot(vtDistLog,vtShadowCorr)
plt.plot(vtDistLogEst,vtShadowCorrEst)
plt.show()

#fading original x estimado
plt.plot(vtDistLogEst, vtDesPequenaEst)
plt.plot(vtDistLog, vtFading)
plt.show()

    


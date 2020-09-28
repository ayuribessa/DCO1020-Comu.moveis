import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import nakagami

class Channel:
    def __init__(self, d0, p0,nPoints, totalLenght, n, sigma, shadowingWindow, m, txPower,nCDF, chFileName):
        self.d0 = d0
        self.p0 = p0
        self.nPoints = nPoints 
        self.totalLenght = totalLenght
        self.n = n
        self.sigma = sigma
        self.shadowingWindow = shadowingWindow 
        self.m = m  
        self.txPower = txPower
        self.nCDF = nCDF
        self.chFileName = chFileName
        self.dMed = totalLenght/nPoints


#Cria o objeto canal com os parâmetros do canal
channel = Channel(5, 0, 50000, 100, 4, 6, 200, 4, 0, 40,'prx_sintetico')
#Vetor de distancias do tx 
vtDist = np.arange(channel.d0 , channel.totalLenght, channel.dMed) # 1 47501
#Tamanho do vetor de amostras geradas no vetor
nSamples = np.size(vtDist)
# Geração do vetor com valores do Path Loss
vt_pathLoss = channel.p0 + 10*channel.n*np.log10(vtDist/channel.d0)
##Geração da V.A gaussiana com média zero e std sigma para o sombreamento
#Quant de pontos de amostras V.a
nShadowSamples = np.floor(nSamples/channel.shadowingWindow)
#Gera amostras, todos os argumentos devem ser int por isso o casting
vtShadowing = channel.sigma*np.random.randn(int(nShadowSamples)) # 1 237
#Amostras de última janela?? FIXME
restShadowing = channel.sigma*np.random.randn(1)*np.ones((int(nSamples%channel.shadowingWindow)))# 1 101
# restShadowing = channel.sigma*np.random.randn(1,1)*np.ones((1,int(nSamples%channel.shadowingWindow)))# 1 101
#Repetição do mesmo valor de sombreamento durante a janela de correlação?? FIXME
vtShadowing = np.ones((channel.shadowingWindow,1))*vtShadowing #200, 237
#Redimensiona o array multi-d vtShadowing para unidimensional, ordenado por cada coluna 'Fortran style'.
vtShadowing = np.ravel(vtShadowing,order='F')
#Acrescenta restShadowing
vtShadowing = np.hstack((vtShadowing,restShadowing))
#filtragem(filtro médio móvel)
jan = int(channel.shadowingWindow/2)
icont = 0
#Criando o array especificando o tamanho
vtShadowCorr = np.empty((nSamples - 2*jan))
#aplicando o filtro da média móvel
for i in range(jan,nSamples-jan):
    vtShadowCorr[icont] = np.mean(vtShadowing[i-jan:i+jan])
    icont += 1

vtShadowCorr = vtShadowCorr*np.std(vtShadowing)/np.std(vtShadowCorr);

vtShadowCorr = vtShadowCorr - np.mean(vtShadowCorr)+ np.mean(vtShadowing);

vtNakagamiNormEnvelope = nakagami.rvs(channel.m, size = len(vtShadowing))
vtNakagamiSampledB = 20*np.log10(vtNakagamiNormEnvelope)
vTxPower = channel.txPower*np.ones(nSamples)

vTxPower = vTxPower[jan:nSamples-jan]
vt_pathLoss = vt_pathLoss[jan:nSamples-jan]
vtFading = vtNakagamiSampledB[jan:nSamples-jan]
vtDist = vtDist[jan:nSamples-jan]
vtPrx = vTxPower - vt_pathLoss + vtShadowCorr + vtFading

log_distancia = np.log10(vtDist)
plt.plot(log_distancia, vtPrx, label = 'Prx canal completo')
plt.plot(log_distancia,channel.txPower - vt_pathLoss,linewidth = 1, label = 'Prx somente perda de percurso')
plt.plot(log_distancia,channel.txPower - vt_pathLoss + vtShadowCorr,linewidth = 2, label = 'Prx perda de percurso + sombreamento')
plt.legend()
plt.xlabel('$log_{10}(d)$')
plt.ylabel('Potência (dBm)')
plt.show()
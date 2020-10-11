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

def geraPLeShadowing(sigma, nSamples, shadowingWindow):
    nShadowSamples = np.floor(nSamples/shadowingWindow)
    #Gera amostras, todos os argumentos devem ser int por isso o casting
    vtShadowing = sigma*np.random.randn(int(nShadowSamples)) # 1 237
    restShadowing = sigma*np.random.randn(1)*np.ones((int(nSamples%shadowingWindow)))# 1 101
    vtShadowing = np.ones((shadowingWindow,1))*vtShadowing #200, 237
    #Redimensiona o array multi-d vtShadowing para unidimensional, ordenado por cada coluna 'Fortran style'.
    vtShadowing = np.ravel(vtShadowing,order='F')
    #Acrescenta restShadowing
    vtShadowing = np.hstack((vtShadowing,restShadowing))

    return vtShadowing

def geraShadowingCorr(shadowingWindow, vtShadowing,nSamples):
    jan = int(shadowingWindow/2)
    icont = 0
    #Criando o array especificando o tamanho
    vtShadowCorr = np.empty((nSamples - 2*jan))
    #aplicando o filtro da média móvel
    for i in range(jan,nSamples-jan):
        vtShadowCorr[icont] = np.mean(vtShadowing[i-jan:i+jan])
        icont += 1

    vtShadowCorr = vtShadowCorr*np.std(vtShadowing)/np.std(vtShadowCorr);
    vtShadowCorr = vtShadowCorr - np.mean(vtShadowCorr)+ np.mean(vtShadowing);
    return vtShadowCorr

def printSinal(vtDist,vtPrx,txPower, vt_pathLoss, vtShadowCorr):
    log_distancia = np.log10(vtDist)
    plt.plot(log_distancia, vtPrx, label = '$P_{rx}$ canal completo')
    plt.plot(log_distancia,txPower - vt_pathLoss,linewidth = 1, label = '$P_{rx}$ somente perda de percurso')
    plt.plot(log_distancia,txPower - vt_pathLoss + vtShadowCorr,linewidth = 1, label = '$P_{rx}$ perda de percurso + sombreamento')
    plt.legend()
    plt.xlabel('$log_{10}(d)$')
    plt.ylabel('Potência (dBm)')
    plt.show()

def printHistPdf(vtNakagamiNormEnvelope,m, vtShadowing):
    n_bins = 100
    plt.hist(vtNakagamiNormEnvelope,bins=n_bins,density=True,label='Histograma normalizado das amostras')
    x = np.linspace(nakagami.ppf(0.001,m),nakagami.ppf(0.999,m), len(vtShadowing))
    rv = nakagami(m)
    plt.plot(x,rv.pdf(x),'r',label='PDF')
    plt.legend()
    plt.show()  
    
def main():
        #Cria o objeto canal com os parâmetros do canal
        channel = Channel(d0 = 5, p0 = 0, nPoints = 50000, totalLenght = 100, n = 4, sigma = 6,
                        shadowingWindow = 200, m = 4, txPower = 0, nCDF = 40, chFileName = 'prx_sintetico')
        
        #Vetor de distancias do tx 
        vtDist = np.arange(channel.d0 , channel.totalLenght, channel.dMed) 
        #Tamanho do vetor de amostras geradas no vetor
        nSamples = np.size(vtDist)
        #Gera vetor de pathloss em dB
        vt_pathLoss = channel.p0 + 10*channel.n*np.log10(vtDist/channel.d0)
        #Gera vetor de shadowing
        vtShadowing = geraPLeShadowing(channel.sigma, nSamples, channel.shadowingWindow)
        #Gera vetor de shadowing correlacionado
        vtShadowCorr = geraShadowingCorr(channel.shadowingWindow, vtShadowing,nSamples)
        #Geração desv. pequena escala
        vtNakagamiNormEnvelope = nakagami.rvs(channel.m, size = len(vtShadowing))
        vtNakagamiSampledB = 20*np.log10(vtNakagamiNormEnvelope)
        vTxPower = channel.txPower*np.ones(nSamples)
        #ajustes 
        jan = int(channel.shadowingWindow/2)
        vTxPower = vTxPower[jan:nSamples-jan]
        vt_pathLoss = vt_pathLoss[jan:nSamples-jan]
        vtFading = vtNakagamiSampledB[jan:nSamples-jan]
        vtDist = vtDist[jan:nSamples-jan]
        #potencia resultante contabilizando todos os efeitos
        vtPrx = vTxPower - vt_pathLoss + vtShadowCorr + vtFading
        printSinal(vtDist,vtPrx,channel.txPower, vt_pathLoss, vtShadowCorr)
        printHistPdf(vtNakagamiNormEnvelope,channel.m,vtShadowing)

main()






#------------------------------------------------------------------ Sem funções -----------
# d0 = 5
# p0 = 0
# nPoints = 5000 
# totalLenght = 100
# n = 4
# sigma = 6
# shadowingWindow = 200 
# m = 4  
# txPower = 0
# nCDF = 40
# chFileName = 'prx_sintetico'
# #Cria o objeto canal com os parâmetros do canal
# channel = Channel(d0, p0, nPoints, totalLenght, n, sigma, shadowingWindow, m, txPower, nCDF, chFileName)

# vtDist = np.arange(channel.d0 , channel.totalLenght, channel.dMed) 
# nSamples = np.size(vtDist)
# vt_pathLoss = channel.p0 + 10*channel.n*np.log10(vtDist/channel.d0)
# ##Geração da V.A gaussiana com média zero e std sigma para o sombreamento
# #Quant de pontos de amostras V.a
# nShadowSamples = np.floor(nSamples/channel.shadowingWindow)
# #Gera amostras, todos os argumentos devem ser int por isso o casting
# vtShadowing = channel.sigma*np.random.randn(int(nShadowSamples)) # 1 237
# restShadowing = channel.sigma*np.random.randn(1)*np.ones((int(nSamples%channel.shadowingWindow)))# 1 101
# vtShadowing = np.ones((channel.shadowingWindow,1))*vtShadowing #200, 237
# #Redimensiona o array multi-d vtShadowing para unidimensional, ordenado por cada coluna 'Fortran style'.
# vtShadowing = np.ravel(vtShadowing,order='F')
# #Acrescenta restShadowing
# vtShadowing = np.hstack((vtShadowing,restShadowing))
# #filtragem(filtro médio móvel)
# jan = int(channel.shadowingWindow/2)
# icont = 0
# #Criando o array especificando o tamanho
# vtShadowCorr = np.empty((nSamples - 2*jan))
# #aplicando o filtro da média móvel
# for i in range(jan,nSamples-jan):
#     vtShadowCorr[icont] = np.mean(vtShadowing[i-jan:i+jan])
#     icont += 1

# vtShadowCorr = vtShadowCorr*np.std(vtShadowing)/np.std(vtShadowCorr);

# vtShadowCorr = vtShadowCorr - np.mean(vtShadowCorr)+ np.mean(vtShadowing);

# #Geração desv. pequena escala
# vtNakagamiNormEnvelope = nakagami.rvs(channel.m, size = len(vtShadowing))
# vtNakagamiSampledB = 20*np.log10(vtNakagamiNormEnvelope)
# vTxPower = channel.txPower*np.ones(nSamples)
# #ajuste 
# vTxPower = vTxPower[jan:nSamples-jan]
# vt_pathLoss = vt_pathLoss[jan:nSamples-jan]
# vtFading = vtNakagamiSampledB[jan:nSamples-jan]
# vtDist = vtDist[jan:nSamples-jan]
# #potencia resultante contabilizando todos os efeitos
# vtPrx = vTxPower - vt_pathLoss + vtShadowCorr + vtFading

# #######prints
# #print sinal completo no tempo
# log_distancia = np.log10(vtDist)
# plt.plot(log_distancia, vtPrx, label = '$P_{rx}$ canal completo')
# plt.plot(log_distancia,channel.txPower - vt_pathLoss,linewidth = 2, label = '$P_{rx}$ somente perda de percurso')
# plt.plot(log_distancia,channel.txPower - vt_pathLoss + vtShadowCorr,linewidth = 2, label = '$P_{rx}$ perda de percurso + sombreamento')
# plt.legend()
# plt.xlabel('$log_{10}(d)$')
# plt.ylabel('Potência (dBm)')
# plt.show()

# #Print histograma + pdf
# plt.hist(vtNakagamiNormEnvelope,bins='auto',density=True,label='Histograma normalizado das amostras')
# n_bins = 100
# x = np.linspace(nakagami.ppf(0.001,channel.m),nakagami.ppf(0.999,channel.m), len(vtShadowing))
# rv = nakagami(channel.m)
# plt.plot(x,rv.pdf(x),'r',label='PDF')
# plt.legend()
# plt.show()  


## plot com subplots(fica pequeno)::
# fig, (ax1,ax2) = plt.subplots(2,1)
# log_distancia = np.log10(vtDist)
# ax1.plot(log_distancia, vtPrx, label = '$P_{rx}$ canal completo')
# ax1.plot(log_distancia,channel.txPower - vt_pathLoss,linewidth = 2, label = '$P_{rx}$ somente perda de percurso') # plt.plot(log_distancia,channel.txPower - vt_pathLoss + vtShadowCorr,linewidth = 1, label = '$P_{rx}$ perda de percurso + sombreamento')
# ax1.legend()
# ax1.set_xlabel('$log_{10}(d)$')
# ax1.set_ylabel('Potência (dBm)')

# ax2.hist(vtNakagamiNormEnvelope,bins='auto',density=True,label='Histograma normalizado das amostras')
# n_bins = 100
# x = np.linspace(nakagami.ppf(0.001,channel.m),nakagami.ppf(0.999,channel.m), len(vtShadowing))
# rv = nakagami(channel.m)
# ax2.plot(x,rv.pdf(x),'r',label='PDF')
# ax2.legend()
# plt.show()
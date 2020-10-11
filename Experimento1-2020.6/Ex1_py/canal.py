#ho1.1.py em forma de função para ser usada em outro programa
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
    
    def GeraVetores(self):
    #Vetor de distancias do tx 
        vtDist = np.arange(self.d0 , self.totalLenght, self.dMed) 
        #Tamanho do vetor de amostras geradas no vetor
        nSamples = np.size(vtDist)
        #Gera vetor de pathloss em dB
        vt_pathLoss = self.p0 + 10*self.n*np.log10(vtDist/self.d0)
        #Gera vetor de shadowing
        vtShadowing = self.geraPLeShadowing(self.sigma, nSamples, self.shadowingWindow)
        #Gera vetor de shadowing correlacionado
        vtShadowCorr = self.geraShadowingCorr(self.shadowingWindow, vtShadowing,nSamples)
        #Geração desv. pequena escala
        vtNakagamiNormEnvelope = nakagami.rvs(self.m, size = len(vtShadowing))
        vtNakagamiSampledB = 20*np.log10(vtNakagamiNormEnvelope)
        vTxPower = self.txPower*np.ones(nSamples)
        #ajustes 
        jan = int(self.shadowingWindow/2)
        vTxPower = vTxPower[jan:nSamples-jan]
        vt_pathLoss = vt_pathLoss[jan:nSamples-jan]
        vtFading = vtNakagamiSampledB[jan:nSamples-jan]
        vtDist = vtDist[jan:nSamples-jan]
        #potencia resultante contabilizando todos os efeitos
        vtPrx = vTxPower - vt_pathLoss + vtShadowCorr + vtFading
        # printSinal(vtDist,vtPrx,self.txPower, vt_pathLoss, vtShadowCorr)
        # printHistPdf(vtNakagamiNormEnvelope,self.m,vtShadowing)
        return vtDist, vt_pathLoss, vtShadowCorr, vtFading, vtPrx

    def geraPLeShadowing(self,sigma, nSamples, shadowingWindow):
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

    def geraShadowingCorr(self,shadowingWindow, vtShadowing,nSamples):
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

# def printSinal(vtDist,vtPrx,txPower, vt_pathLoss, vtShadowCorr):
#     log_distancia = np.log10(vtDist)
#     plt.plot(log_distancia, vtPrx, label = '$P_{rx}$ canal completo')
#     plt.plot(log_distancia,txPower - vt_pathLoss,linewidth = 2, label = '$P_{rx}$ somente perda de percurso')
#     plt.plot(log_distancia,txPower - vt_pathLoss + vtShadowCorr,linewidth = 2, label = '$P_{rx}$ perda de percurso + sombreamento')
#     plt.legend()
#     plt.xlabel('$log_{10}(d)$')
#     plt.ylabel('Potência (dBm)')
#     plt.show()

# def printHistPdf(vtNakagamiNormEnvelope,m, vtShadowing):
#     n_bins = 100
#     plt.hist(vtNakagamiNormEnvelope,bins=n_bins,density=True,label='Histograma normalizado das amostras')
#     x = np.linspace(nakagami.ppf(0.001,m),nakagami.ppf(0.999,m), len(vtShadowing))
#     rv = nakagami(m)
#     plt.plot(x,rv.pdf(x),'r',label='PDF')
#     plt.legend()
#     plt.show()  

# #Cria o objeto canal com os parâmetros do canal
# channel = Channel(d0 = 5, p0 = 0, nPoints = 5000, totalLenght = 100, n = 4, sigma = 6,
#                 shadowingWindow = 200, m = 4, txPower = 0, nCDF = 40, chFileName = 'prx_sintetico')






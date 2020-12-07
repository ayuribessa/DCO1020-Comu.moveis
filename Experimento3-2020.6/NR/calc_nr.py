#Constantes para este projeto
SCALING_FACTOR = 1
DIRECTION = 'Downlink'
MODE = 'FDD'
RMAX = 948/1024 #Max Coding Rate da tabela MCS


def calcTputNR(numComponetCarriers, ModOrder, LayersMIMO, numerology, nRBs, overhead):

    if (nRBs < 1):
        nRBs = 1
    elif (nRBs > 275):
        nRBs = 275

    ts = 10e-3 / (14 * (2 ** numerology))
    maxTput = 10e-6 * numComponetCarriers * ModOrder * LayersMIMO * RMAX * (nRBs * 12 / ts) * (1 - overhead)
    return maxTput



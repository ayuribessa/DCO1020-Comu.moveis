#Constantes para este projeto
SCALING_FACTOR = 1
DIRECTION = 'Downlink'
MODE = 'FDD'
RMAX = 948/1024 #Max Coding Rate da tabela MCS


# carrierAgregation = {
#     '1': 1,                    
#     '2': 2,                  
#     '3': 3,                
#     '4': 4,                  
#     '5': 5,                  
#     '6': 6,                  
#     '7': 7,                  
#     '8': 8,                  
#     '9': 9,                  
#     '10': 10,                 
#     '11': 11,                 
#     '12': 12,                 
#     '13': 13,                 
#     '14': 14,                 
#     '15': 15,                   
#     '16': 16,                 
# }                           

# mimo = {
#     '1': 1,
#     '2': 2,
#     '4': 4,
#     '8': 8,
#  }

# modulationOrder = {
#     '2': 2, #2QPSK
#     '4': 4, #16QAM
#     '6': 6, #64QAM
#     '8': 8, #256QAM
# }

# numerology = { 
#     #Espaçamento entre subportadoras.
#     '0': (0,15),
#     '1': (1,30),
#     '2': (2,60),
#     '3': (3,120),
#     '4': (4,240),
# }

# overhead = {
#     '0.14': 0.14,
#     '0.18': 0.18,
# }


def calcTputNR(numComponetCarriers, ModOrder, LayersMIMO, numerology, nRBs, overhead):


    if (nRBs < 1):
        nRBs = 1
    elif (nRBs > 275):
        nRBs = 275

    ts = 10e-3 / (14 * (2 ** numerology))
    maxTput = 10e-6 * numComponetCarriers * ModOrder * LayersMIMO * RMAX * (nRBs * 12 / ts) * (1 - overhead)
    return maxTput



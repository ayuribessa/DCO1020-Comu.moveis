

import numpy as np
import matplotlib.pyplot as plt
import begin
from colorama import init, Fore, Style 
init(autoreset=True)

#constantes
PTDBM = 57
HMOB = 5
AHM = 3.2*np.log10(11.75*HMOB)**2 - 4.97
HBS = 30
SENSITIVITY = -104 

def calc_centers_bss(r):
    ''' Calcula o centro de cada hexágono e retorna um array de centros'''
    dim_x = 5*r
    dim_y = 6*np.sqrt(3/4) * r
    offset = np.pi/6
    center_hexes = [0]
    for i in range(6):
        center_hexes.append(r*np.sqrt(3)*np.exp(1j*(i*np.pi/3 + offset)))
    center_hexes = np.asarray(center_hexes) + complex(dim_x/2, dim_y/2)
    return center_hexes, dim_x, dim_y


def calc_pontos_medicao(center_hexes, dim_x, dim_y, r):
    '''Calcula a matriz de pontos de medição de referência, depois cada Bss tem seus pontos de medição 
    ajustados de acordo com seu centro '''
    passo = int(np.ceil(r/50))
    dim_y = dim_y + dim_y % passo
    dim_x = dim_x + dim_x % passo
    posx_mat, posy_mat = np.meshgrid(np.linspace(0, dim_x + passo, passo), np.linspace(0, dim_y + passo, passo))
    mt_pts_medicao = posx_mat + 1j*posy_mat
    mt_pts_medicao_cada_erb = ajuste_pontos_medicao_cada_bss(mt_pts_medicao, center_hexes, r)
    return mt_pts_medicao, mt_pts_medicao_cada_erb


def ajuste_pontos_medicao_cada_bss(mt_pts_medicao, centros, r):
    ''' Ajusta os pontos de medição em relação a cada ERB, criando um array com os pontos
    de medição de cada ERB '''
    mt_pts_medicao_cada_erb = np.empty([7, np.shape(mt_pts_medicao)[0], np.shape(mt_pts_medicao)[1]],dtype=complex)
    for i in range(7):
        mt_pts_medicao_cada_erb[i] = mt_pts_medicao - centros[i]
       
    return mt_pts_medicao_cada_erb


def calc_mts_pt(mt_pts_medicao_cada_erb, mt_pts_medicao, center_hexes, freq, r, x):
    ''' Cria  as matrizes de potencias finais e popula com os valores de potencia para as ERBs, e
     cria uma matriz para as 7 ERBs em conjunto.'''
    mt_pt_dbm_final = np.ones(np.shape(mt_pts_medicao))*np.NINF
    mt_pt_dbm_bss = np.empty_like(mt_pts_medicao_cada_erb,dtype=float)
    for i in range(7):
        mt_dist_bss = np.abs(mt_pts_medicao_cada_erb[i])
        np.where(mt_dist_bss < r, r, mt_dist_bss)
        mt_pt_dbm = x + (44.9-6.55*np.log10(HBS))*np.log10(mt_dist_bss/1e3) - 13.82*np.log10(HBS) - AHM
        mt_pt_dbm_bss[i] = PTDBM - mt_pt_dbm
        mt_pt_dbm_final = np.maximum(mt_pt_dbm_final, mt_pt_dbm_bss[i])

    return mt_pt_dbm_bss , mt_pt_dbm_final
    
@begin.start
def main(modelo:'MODELO = Modelo a ser usado. Opções: okomura ou cost' = 'okomura'):
    ''' Calcula raio para outage máximo de 10% usando modelo okomura-hata ou cost231'''
    print(f'Para modelo: {modelo} ')
    freqs = [800, 900, 1800, 1900, 2100]
    for freq in freqs:
        #iniciando o raio em 1e3 por convenção
        raio = 1e3
        if modelo != 'okomura' and modelo != 'cost':
            print('Opção não disponível, usando <okomura> ')
            modelo = 'okomura'
        if modelo == 'okomura':
            x = 69.55 + 26.16*np.log10(freq) 
        elif modelo == 'cost':
            x = 46.3 + 33.9*np.log10(freq) + 3

        while True:
            # calcula os centros dos hexagonos(posição das BSs)
            centers_bss, dim_x, dim_y = calc_centers_bss(raio)
            # Cria os pontos de medição de referência, e os pontos de medição ajustados de cada ERB
            mt_pts_medicao, mt_pts_medicao_cada_erb = calc_pontos_medicao(centers_bss, dim_x, dim_y, raio)
            # calculas as matrizes de potencia
            mt_pt_dbm_bss, mt_pt_dbm_final = calc_mts_pt(mt_pts_medicao_cada_erb, mt_pts_medicao, centers_bss, freq,raio, x)
            outage_rate = np.count_nonzero(mt_pt_dbm_final < SENSITIVITY) / np.size(mt_pt_dbm_final) 
            if(outage_rate >= 0.09988):
                break
            raio += 20
        
        print('Com a frequencia'+ Fore.BLUE + Style.BRIGHT +  f' {freq} MHz',end='')
        print(', a taxa de outage de'+Fore.GREEN + Style.BRIGHT + f' {outage_rate:.2%} ', end='')
        print(', foi atingida com raio'+ Fore.LIGHTYELLOW_EX + Style.BRIGHT +  f' {raio} m')


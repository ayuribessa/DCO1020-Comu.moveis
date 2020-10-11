
import numpy as np
import matplotlib.pyplot as plt

#constantes
PTDBM = 57
HMOB = 5
AHM = 3.2*np.log10(11.75*HMOB)**2 - 4.97
HBS = 30
# FC = 800
RAIO_HEX = 10e3
DIM_X = 5*RAIO_HEX
DIM_Y = 6*np.sqrt(3/4) * RAIO_HEX
SENSITIVITY = -104 
# sem casting para int no ceil, np.linspace levanta exceção no numpy 1.18.1
PASSO = int(np.ceil(RAIO_HEX/50))


def main():

    freqs = [800, 900, 1800, 1900, 2100]
    for freq in freqs:
        # calcula os centros dos hexagonos(posição das BSs)
        centers_bss = calc_centers_bss()
        # Cria os pontos de medição de referência, e os pontos de medição ajustados de cada ERB
        mt_pts_medicao, mt_pts_medicao_cada_erb = calc_pontos_medicao(centers_bss)
        # calculas as matrizes de potencia
        mt_pt_dbm_bss, mt_pt_dbm_final = calc_mts_pt(mt_pts_medicao_cada_erb, mt_pts_medicao, centers_bss, freq)
        # Exibe os RME's
        # print_RME(mt_pt_dbm_bss, mt_pt_dbm_final, mt_pts_medicao, centers_bss)
        outage_rate = np.count_nonzero(mt_pt_dbm_final < SENSITIVITY) / np.size(mt_pt_dbm_final) 
        print(f'Para frequencia {freq} Mhz, a taxa de outage foi de {outage_rate:.2%}')


def draw_hex(center, r=RAIO_HEX):
    ''' Desenha um hexagono com um centro e um raio informados, calculando
    seus vértices como números complexos e ligando os pontos'''
    vertices_hex = [r*np.exp(np.pi/3*i*1j) for i in range(7)]
    vertices_hex = [vertices_hex[i] + center for i in range(7)]
    plt.plot(np.real(vertices_hex), np.imag(vertices_hex))


def draw_grid(center_hexes, r=RAIO_HEX):
    ''' Desenha uma grade com 7 hexagonos. Recebe o centro de cada um e seu raio '''
    for i in range(len(center_hexes)):
        draw_hex(center_hexes[i])
    # plt.scatter(np.real(center_hexes),np.imag(center_hexes))
    plt.axis('equal')
    plt.show()


def calc_centers_bss(r=RAIO_HEX):
    ''' Calcula o centro de cada hexágono e retorna um array de centros'''
    offset = np.pi/6
    center_hexes = [0]
    for i in range(6):
        center_hexes.append(r*np.sqrt(3)*np.exp(1j*(i*np.pi/3 + offset)))
    center_hexes = np.asarray(center_hexes) + complex(DIM_X/2, DIM_Y/2)
    # draw_grid(center_hexes)
    return center_hexes


def calc_pontos_medicao(center_hexes, r=RAIO_HEX):
    '''Calcula a matriz de pontos de medição de referência, depois cada Bss tem seus pontos de medição 
    ajustados de acordo com seu centro '''
    dim_y = DIM_Y + DIM_Y % PASSO
    dim_x = DIM_X + DIM_X % PASSO
    posx_mat, posy_mat = np.meshgrid(np.linspace(0, dim_x + PASSO, PASSO), np.linspace(0, dim_y + PASSO, PASSO))
    # posx_mat, posy_mat = np.meshgrid(np.linspace(0,dim_x ,PASSO),np.linspace(0,dim_y ,PASSO))
    #FIXME mt_pts_medicao deveria ser float ????
    mt_pts_medicao = posx_mat + 1j*posy_mat
    mt_pts_medicao_cada_erb = ajuste_pontos_medicao_cada_bss(mt_pts_medicao, center_hexes, r=RAIO_HEX)
    return mt_pts_medicao, mt_pts_medicao_cada_erb


def ajuste_pontos_medicao_cada_bss(mt_pts_medicao, centros, r=RAIO_HEX):
    ''' Ajusta os pontos de medição em relação a cada ERB, criando um array com os pontos
    de medição de cada ERB '''
    mt_pts_medicao_cada_erb = np.empty([7, np.shape(mt_pts_medicao)[0], np.shape(mt_pts_medicao)[1]],dtype=complex)
    for i in range(7):
        mt_pts_medicao_cada_erb[i] = mt_pts_medicao - centros[i]
       
    return mt_pts_medicao_cada_erb


def calc_mts_pt(mt_pts_medicao_cada_erb, mt_pts_medicao, center_hexes, freq, r=RAIO_HEX):
    ''' Cria  as matrizes de potenciasi finais e popula com os valores de potencia para as ERBs
     e uma para as 7 ERBs em conjunto.'''
    mt_pt_dbm_final = np.ones(np.shape(mt_pts_medicao))*np.NINF
    mt_pt_dbm_bss = np.empty_like(mt_pts_medicao_cada_erb,dtype=float)
    for i in range(7):
        mt_dist_bss = np.abs(mt_pts_medicao_cada_erb[i])
        np.where(mt_dist_bss < r, r, mt_dist_bss)
        x = (44.9 - 6.55*np.log10(HBS))*np.log10(mt_dist_bss/1e3)
        mt_pl_dbm = 69.55 + 26.16*np.log10(freq) + x - 13.82*np.log10(HBS) - AHM
        mt_pt_dbm_bss[i] = PTDBM - mt_pl_dbm
        mt_pt_dbm_final = np.maximum(mt_pt_dbm_final, mt_pt_dbm_bss[i])

    return mt_pt_dbm_bss , mt_pt_dbm_final


def print_RME(mt_pt_dbm_bss, mt_pt_dbm_final, mt_pts_medicao, center_hexes):
    for i in range(7):
        plt.pcolor(np.real(mt_pts_medicao), np.imag(mt_pts_medicao), mt_pt_dbm_bss[i],vmax = -115, vmin = -50)
        plt.colorbar()
        plt.title(f"Erb {i + 1}")
        draw_grid(center_hexes)
    #     plt.show()
    plt.pcolor(np.real(mt_pts_medicao), np.imag(mt_pts_medicao), mt_pt_dbm_final,vmax = -115, vmin = -50)
    plt.colorbar()
    plt.title("Todas as 7 Erbs")
    draw_grid(center_hexes)
    plt.show()
    

main()

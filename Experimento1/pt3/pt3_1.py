import numpy as np 
import matplotlib.pyplot as plt 

PTDBM = 57
HMOB = 5
AHM = 3.2*np.log10(11.75*HMOB)**2 - 4.97
HBS = 30
FC = 800
RAIO_HEX = 5E3
DIM_X = 5*RAIO_HEX
DIM_Y = 6*np.sqrt(3/4)* RAIO_HEX
PASSO = np.ceil(RAIO_HEX/10)


def main():
    centers_bss = calc_centers_bss()
    mt_pts_medicao, mt_pts_medicao_cada_erb = calc_pontos_medicao(centers_bss)
    # mt_pt_dbm_bss, mt_pt_dbm_final = calc_mts_pt(mt_pts_medicao_cada_erb, mt_pts_medicao,centers_bss)
    mt_pt_dbm_bss, mt_pt_dbm_final = calc_mts_pt(mt_pts_medicao_cada_erb, mt_pts_medicao,centers_bss)
    print_RME(mt_pt_dbm_bss, mt_pt_dbm_final,mt_pts_medicao, centers_bss)

def draw_hex(center,r=RAIO_HEX):
    vertices_hex = [r*np.exp(np.pi/3*i*1j) for i in range(7)]
    vertices_hex = [vertices_hex[i] + center for i in range(7)]
    plt.plot(np.real(vertices_hex), np.imag(vertices_hex))

def draw_grid(center_hexes,r=RAIO_HEX):
    for i in range(len(center_hexes)):
        draw_hex(center_hexes[i])
    # plt.scatter(np.real(center_hexes),np.imag(center_hexes))
    plt.axis('equal')
    plt.show()

def calc_centers_bss(r=RAIO_HEX):
    offset = np.pi/6
    center_hexes = [0]
    for i in range(6):
        center_hexes.append(r*np.sqrt(3)*np.exp(1j*(i*np.pi/3 + offset)))
    center_hexes = np.asarray(center_hexes) +  complex(DIM_X/2,DIM_Y/2)
    # draw_grid(center_hexes)
    return center_hexes

#matriz de referência, depois cada Bss tem seus pontos de medição ajustados de acordo
#com seu centro
def calc_pontos_medicao(center_hexes, r=RAIO_HEX):
    dim_y =  DIM_Y + DIM_Y % PASSO
    dim_x =  DIM_X + DIM_X % PASSO
    posx_mat, posy_mat = np.meshgrid(np.linspace(0,dim_x ,PASSO),np.linspace(0,dim_y ,PASSO))
    mt_pts_medicao = posx_mat + 1j*posy_mat
    mt_pts_medicao_cada_erb = ajuste_pontos_medicao_cada_bss(mt_pts_medicao,center_hexes,r=RAIO_HEX)
    return mt_pts_medicao, mt_pts_medicao_cada_erb

def ajuste_pontos_medicao_cada_bss(mt_pts_medicao,centros,r=RAIO_HEX):
    mt_pts_medicao_cada_erb = []
    for i in range(7):
        pontos_bs = mt_pts_medicao - centros[i]
        mt_pts_medicao_cada_erb.append(pontos_bs)
        # plt.scatter(np.real(mt_pts_medicao_cada_erb[i]),np.imag(mt_pts_medicao_cada_erb[i]))
        # plt.title(f'ERB {i}')
        # draw_grid(centros - centros[i])
    mt_pts_medicao_cada_erb = np.asarray(mt_pts_medicao_cada_erb)
    return mt_pts_medicao_cada_erb

def calc_mts_pt(mt_pts_medicao_cada_erb, mt_pts_medicao, center_hexes,r=RAIO_HEX):
    #TODO verificar o os calculos dessa função.
    mt_pt_dbm_final = np.ones(np.shape(mt_pts_medicao))*np.NINF
    mt_pt_dbm_bss = np.empty(np.shape(mt_pts_medicao_cada_erb))
    for i in range(7):
        mt_dist_bss = np.abs(mt_pts_medicao_cada_erb[i])
        np.where(mt_dist_bss < r, r, mt_dist_bss)
        x = 44.9 - 6.55*np.log10(HBS)*np.log10(mt_dist_bss/1e3)
        mt_pl_dbm = 69.55 - 26.16*np.log10(FC) + x - 13.82*np.log10(HBS) - AHM
        mt_pt_dbm_bss[i] = PTDBM - mt_pl_dbm
        mt_pt_dbm_final = np.maximum(mt_pt_dbm_final,mt_pt_dbm_bss[i])
        # plt.pcolor(mt_pts_medicao.real, mt_pts_medicao.imag, mt_pt_dbm_bss[i])
        # plt.colorbar()
        # plt.title(f"Erb {i + 1}")
        # draw_grid(center_hexes)
        plt.pcolor(np.real(mt_pts_medicao), np.imag(mt_pts_medicao), mt_pt_dbm_final)
        plt.colorbar()
        plt.title("Todas as 7 Erbs")
        draw_grid(center_hexes)
        plt.show()
    return mt_pt_dbm_bss , mt_pt_dbm_final

def print_RME(mt_pt_dbm_bss, mt_pt_dbm_final, mt_pts_medicao,center_hexes):
    # for i in range(7):
    #     plt.pcolor(mt_pts_medicao.real, mt_pts_medicao.imag, mt_pt_dbm_bss[i])
    #     plt.colorbar()
    #     plt.title(f"Erb {i + 1}")
    #     draw_grid(center_hexes)
    #     plt.show()
    plt.pcolor(mt_pts_medicao.real, mt_pts_medicao.imag, mt_pt_dbm_final)
    plt.colorbar()
    plt.title("Todas as 7 Erbs")
    draw_grid(center_hexes)
    plt.show()
main()

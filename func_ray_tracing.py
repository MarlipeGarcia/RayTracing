# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 15:09:17 2022

Conjuntos de funções para o programa de raios acústicos

@author: Marlipe
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Normal do plano
# =============================================================================

def normal_plano(plano):
    
    normal_plano = np.zeros((1,3))
    
    x_max, y_max, z_max = np.max(np.max(plano, axis=0), axis=0)
    ponto_centro = np.array([x_max, y_max, z_max])/2
    
    for pos in plano:
        v1 = pos[1,:] - pos[0,:]
        v2 = pos[2,:] - pos[0,:]
        n = np.cross(v1, v2)
        n = n / np.linalg.norm(n)
        
        n_dir = np.linalg.norm((n*1 + pos[0,:]) - ponto_centro)
        n_inv = np.linalg.norm((n*-1 + pos[0,:]) - ponto_centro)
        
        if n_inv < n_dir:
            
            n = n*-1
        
        normal_plano = np.append(normal_plano, [n], axis=0)
    
    normal_plano = normal_plano[-len(plano):]
    return normal_plano

# =============================================================================
# Distancia ponto e plano
# =============================================================================

def dist_ponto_plano(ponto, vetor_direcao, plano):
       
    nl = vetor_direcao / np.linalg.norm(vetor_direcao)
    
    normal = normal_plano(plano)
    
    C_ii = np.inf
    plano_inter = 0
    ponto_inter = np.zeros(3)
    dist_perc = 0
        
    for ii in range (0, len(plano)):

        a, b, c = normal[ii]
        d = np.dot(normal[ii],plano[ii][0])
    
        A = np.array([[a, b, c, 0],
                      [1, 0, 0, -nl[0]],
                      [0, 1, 0, -nl[1]],
                      [0, 0, 1, -nl[2]]])
    
        B = np.array([d, ponto[0], ponto[1], ponto[2]])
    
        try:
            C = np.linalg.solve(A,B)
        except:
            # print('Singular matrix para o plano %d' %ii)
            C = [np.inf, np.inf, np.inf, np.inf]
        
        # ax.plot(C[0],C[1],C[2],'*g')
        
        PlanoValido = C[3] > 0 and C[3] < C_ii
        
        PlanoPontoDentro = C[0] <= np.max(plano[ii][:,0]) and C[0] >= np.min(plano[ii][:,0]) \
                        and C[1] <= np.max(plano[ii][:,1]) and C[1] >= np.min(plano[ii][:,1]) \
                        and C[2] <= np.max(plano[ii][:,2]) and C[2] >= np.min(plano[ii][:,2])
                
        if PlanoValido and PlanoPontoDentro:
                                    
            ponto_inter = C[0:3]
            plano_inter = ii
            dist_perc = np.linalg.norm(ponto-ponto_inter)
            C_ii = C[3]
            # print('O ponto de intersecção é: ', ponto_inter)
            # print('O plano é: ', plano_inter)
            # print('A distância é: ', dist_perc)
    
    plot_raio = False
    
    if plot_raio:
        ax = plt.gca()
        ax.plot([ponto[0], ponto_inter[0]],[ponto[1], ponto_inter[1]],[ponto[2], ponto_inter[2]],'--g')
        plt.pause(1)
    
    return ponto_inter, plano_inter, dist_perc

# =============================================================================
# Distancia ponto e reta
# =============================================================================

def dist_ponto_reta(ponto,vetor_direcao,receptor):
    
    nl = vetor_direcao / np.linalg.norm(vetor_direcao)
       
    distRecep = np.cross((receptor - ponto),nl)
    
    distRecep = np.linalg.norm(distRecep)
    
    # print('A distância com receptor é: ', distRecep)
   
    return distRecep

# =============================================================================
# Parâmetros objetivos
# =============================================================================

def calc_parametros_objetivos(tempo, energia):
    
    energia_dB = 10*np.log10(energia)
    
    T0_ponto = np.where(energia_dB < 0)
    T0_ponto = T0_ponto[0]
    
    T5_ponto = np.where(energia_dB >= -5)
    T5_ponto = T5_ponto[0]
    
    T10_ponto = np.where(energia_dB >= -10)
    T10_ponto = T10_ponto[0]
    
    T25_ponto = np.where(energia_dB >= -25)
    T25_ponto = T25_ponto[0]
    
    T35_ponto = np.where(energia_dB >= -35)
    T35_ponto = T35_ponto[0]
    
    PolyCoef = np.polyfit(energia_dB[T5_ponto[-1]:T25_ponto[-1]], tempo[T5_ponto[-1]:T25_ponto[-1]], 1)
    T20 = (np.polyval(PolyCoef, -25) - np.polyval(PolyCoef, -5)) * 3
    
    PolyCoef = np.polyfit(energia_dB[T5_ponto[-1]:T35_ponto[-1]], tempo[T5_ponto[-1]:T35_ponto[-1]], 1)
    T30 = (np.polyval(PolyCoef, -35) - np.polyval(PolyCoef, -5)) * 2
    
    PolyCoef = np.polyfit(energia_dB[T0_ponto[0]:T10_ponto[-1]], tempo[T0_ponto[0]:T10_ponto[-1]], 1)
    EDT = (np.polyval(PolyCoef, -10) - np.polyval(PolyCoef, 0)) * 6
    
    T50ms_ponto = np.where(tempo <= 50e-3)
    T50ms_ponto = T50ms_ponto[0]
    
    C50 = 10*np.log10(np.sum(energia[:T50ms_ponto[-1]]) / np.sum(energia[T50ms_ponto[-1]:-1]))
    
    D50 = np.sum(energia[:T50ms_ponto[-1]]) / np.sum(energia)
    
    return T20, T30, EDT, C50, D50

# =============================================================================
# Leitura do arquivo DXF
# =============================================================================

def leitura_dxf(sala):
    
    file = open(sala, 'r')
    
    plano_count = 0
    plano_pos = []
    plano_layer = []
    
    receptor_count = 0
    receptor_pos = []
    receptor_layer = []
    
    fonte_count = 0
    fonte_pos = []
    fonte_layer = []
    
    while True:
        
        line = file.readline()
        
        if line.strip() == '3DFACE':
            
            pos = np.zeros((4,3))
            
            for ii in range(7):
                line = file.readline()
            
            line = file.readline()
            layer = line
            
            for ii in range(2):
                line = file.readline()
        
            line = file.readline()
            line = file.readline()
            pos[0][0] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[0][1] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[0][2] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[1][0] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[1][1] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[1][2] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[2][0] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[2][1] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[2][2] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[3][0] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[3][1] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[3][2] = line.strip()
            
            if layer[0:2] == 'M_':
            
                plano_pos.append(pos)
                plano_layer.append(layer)
                plano_count += 1
            
        if line.strip() == 'POINT':
            
            pos = np.zeros((3))
            
            for ii in range(7):
                line = file.readline()
            
            line = file.readline()
            layer = line
            
            for ii in range(2):
                line = file.readline()
        
            line = file.readline()
            line = file.readline()
            pos[0] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[1] = line.strip()
            
            line = file.readline()
            line = file.readline()
            pos[2] = line.strip()                   
            
            if layer[0:2] == 'R_':
            
                receptor_pos.append(pos)
                receptor_layer.append(layer)
                receptor_count += 1
                
            if layer[0:2] == 'F_':
            
                fonte_pos.append(pos)
                fonte_layer.append(layer)
                fonte_count += 1
        
        if not line:
            break    
    
    file.close()
    
    return plano_pos, plano_layer, receptor_pos, receptor_layer, fonte_pos, fonte_layer
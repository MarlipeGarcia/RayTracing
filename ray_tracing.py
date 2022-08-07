# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 11:38:38 2022

@author: Marlipe
"""

import numpy as np
import matplotlib.pyplot as plt
import func_ray_tracing as func

# =============================================================================
# Dados de entrada
# =============================================================================

# Leitura da sala

plano_pos, plano_layer, receptor_pos, receptor_layer, fonte_pos, fonte_layer = func.leitura_dxf('Sala1.dxf')

# Dados da receptor

# ReceptorPosicao = np.array([5,1,2])
NumReceptor = len(receptor_pos)
ReceptorPosicao = receptor_pos[0]
ReceptorRaioInicial = .1

# Dados da fonte

# FontePosicao = np.array([1,1,1])
FontePosicao = fonte_pos[2]

FonteNumDiv = 36

FonteTheta = np.linspace(0, 2*np.pi, FonteNumDiv, endpoint=False)
FontePhi = np.linspace(-np.pi/2, np.pi/2, FonteNumDiv)
FonteTheta, FontePhi = np.meshgrid(FonteTheta, FontePhi)

FonteVetorDirecao = np.array([np.cos(FonteTheta), np.sin(FonteTheta), np.sin(FontePhi)])
FonteVetorDirecao = np.reshape(FonteVetorDirecao.T, (FonteNumDiv ** 2,3))

numRaio = len(FonteVetorDirecao)

potFonte = 1
distFonteRecep = np.sqrt(np.sum((FontePosicao-ReceptorPosicao)**2))  # Distância entre fonte e receptor
angAber = np.arcsin(ReceptorRaioInicial/distFonteRecep)  # Ângulo de abertura entre fonte e esfera receptor

EnergiaInicial = potFonte/(2*np.pi*numRaio*(distFonteRecep**2)*(1-np.cos(angAber)))
EnergiaInicial = potFonte/numRaio

# Dados do plano

# Lx = 5.0
# Ly = 7.3
# Lz = 3.0

# plano = np.array([
#                   [[0, 0, 0], [Lx, 0, 0], [Lx, Ly, 0], [0, Ly, 0]],
#                   [[0, 0, 0], [0, 0, Lz], [Lx, 0, Lz], [Lx, 0, 0]],
#                   [[0, 0, 0], [0, Ly, 0], [0, Ly, Lz], [0, 0, Lz]],
#                   [[0, 0, Lz], [0, Ly, Lz], [Lx, Ly, Lz], [Lx, 0, Lz]],
#                   [[0, Ly, 0], [Lx, Ly, 0], [Lx, Ly, Lz], [0, Ly, Lz]],
#                   [[Lx, 0, 0], [Lx, 0, Lz], [Lx, Ly, Lz], [Lx, Ly, 0]],
#                 ])


plano = plano_pos
Lx, Ly, Lz = np.max(np.max(plano, axis=0), axis=0)

# plano_absorcao = np.array([.2, .3, .4, .4, .3, .2])
# plano_espalhamento = np.array([.2, .3, .4, .4, .3, .2])

plano_absorcao = np.array([.02, .02, .02, .02, .02, .02, .02, .1])
plano_espalhamento = np.array([.2, .2, .2, .2, .2, .2, .2, .5])

normal = func.normal_plano(plano)

# Parâmetros de parada

maxRefl = round(-5/np.log10(1-np.min(plano_absorcao)))

maxDist = np.inf

minEnergia = 1e-6

# Parâmetros energia decaimento

vel_propag = 341
absorcao_ar = 0.4e-3
FreqSimul = 16000

# =============================================================================
# Plot da sala
# =============================================================================

plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Acústica de Raios - Sala de Aula')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')

# Plot fonte

ax.plot(FontePosicao[0],FontePosicao[1],FontePosicao[2],'*r')
ax.text(FontePosicao[0],FontePosicao[1],FontePosicao[2], 'Fonte')

# Plot receptor

for ii in range(NumReceptor):
    ReceptorPosicao = receptor_pos[ii]
    ax.plot(ReceptorPosicao[0],ReceptorPosicao[1],ReceptorPosicao[2],'og')
    ax.text(ReceptorPosicao[0],ReceptorPosicao[1],ReceptorPosicao[2],'Receptor %d' %ii)

# Plot plano

num_plano = 0
for pos_plano in plano:
    pos_plano = np.append(pos_plano, [pos_plano[0,:]], axis=0)
    ax.plot(pos_plano[:,0],pos_plano[:,1],pos_plano[:,2],'-b')
    
    pos_text = (pos_plano[2,:] + pos_plano[0,:])/2
    label_plano = 'Plano %d' %num_plano
    ax.text(pos_text[0], pos_text[1], pos_text[2], label_plano)
    num_plano += 1

plt.pause(.1)
plt.show()

# =============================================================================
# Raios acústicos
# =============================================================================

DistanciaPercorrida = np.zeros((NumReceptor,numRaio))
EnergiaInicial = np.ones((NumReceptor,numRaio)) * EnergiaInicial
EnergiaDissipada = np.ones((NumReceptor,numRaio))

for ir in range(NumReceptor):
    
    print('\n')
    print('Receptor número: %d' %(ir)) 
    
    ReceptorPosicao = receptor_pos[ir]
    

    for ii in range(numRaio):
    
        ponto = FontePosicao 
        RaioDirecao = FonteVetorDirecao[ii]
        # print('\n')
        # print('Raio número: %d' %(ii))     

        ReceptorRaio = ReceptorRaioInicial      

        for num_refl in range(maxRefl):
            
            dist_recep = func.dist_ponto_reta(ponto, RaioDirecao, ReceptorPosicao)
    
            ponto_inter, plano_inter, dist_aux = func.dist_ponto_plano(ponto, RaioDirecao, plano)
    
            ponto = ponto_inter
            DistanciaPercorrida[ir,ii] += dist_aux
       
            # print('Distância percorrida: ', DistanciaPercorrida[ii])
            # print('Número de reflexão: ', num_refl)
            # print('Plano interceptado: ', plano_inter)
        
            if DistanciaPercorrida[ir,ii] > maxDist or dist_recep < ReceptorRaio or EnergiaDissipada[ir,ii] < minEnergia:
            
                break
    
            ReceptorRaio *= 1.01 
    
            RaioDirecao = RaioDirecao - 2 * np.dot(RaioDirecao, normal[plano_inter]) * normal[plano_inter]
        
            g = np.random.random(1)
        
            if g > plano_espalhamento[plano_inter]:
            
                RaioDirecao = RaioDirecao
            
            elif g <= plano_espalhamento[plano_inter]:
            
                g1, g2 = np.random.random(2)
                EspalhamentoTheta = 2*np.pi*g1
                EspalhamentoPhi = np.arccos(np.sqrt(g2))
                EspalhamentoDirecao = np.array([np.cos(EspalhamentoTheta), np.sin(EspalhamentoTheta), np.sin(EspalhamentoPhi)])
                EspalhamentoDirecao = EspalhamentoDirecao / np.linalg.norm(EspalhamentoDirecao)    
        
                RaioDirecao = EspalhamentoDirecao
        
            EnergiaDissipada[ir,ii] *= 1 - plano_absorcao[plano_inter] 
    
        # print('Número de reflexões: ', num_refl)
        # print('Distância percorrida: ', DistanciaPercorrida[ii])
        # print('Energia raio: ', EnergiaDissipada[ii])     
        
print('Acabou')

# =============================================================================
# Energia Decaimento
# =============================================================================

EnergiaRaio = EnergiaDissipada * EnergiaInicial * np.exp(- absorcao_ar * DistanciaPercorrida)

TempoRaio = DistanciaPercorrida / vel_propag

TempoFreqSimul = np.around(FreqSimul * TempoRaio)

from scipy.sparse import csr_matrix

EnergiaDec = {}
TempoEnergia = {}
EnergiaTotal = {}

for ir in range(NumReceptor):

    EnergiaCol = np.zeros(len(TempoFreqSimul[ir]), dtype = int)
    EnergiaRow = np.array(TempoFreqSimul[ir], dtype = int)
    EnergiaTotal[ir] = csr_matrix((EnergiaRaio[ir], (EnergiaRow, EnergiaCol))).toarray()

    TempoEnergia[ir] = np.arange(len(EnergiaTotal[ir])) / FreqSimul

    aux1 = np.sum(EnergiaTotal[ir]) * TempoEnergia[ir][1]
    aux2 = np.cumsum(EnergiaTotal[ir]) * TempoEnergia[ir][1]

    EnergiaDec[ir] = ((aux1) - (aux2))/ aux1

# =============================================================================
# T60 Teórico
# =============================================================================

T60Sabine = 0.161*Lx*Ly*Lz/(np.mean(plano_absorcao)*2*(Lx*Ly + Lx*Lz + Ly*Lz))
print('\n')
print('T60 Sabine igual a %.2f s' %T60Sabine)

# =============================================================================
# Parâmetro objetivos
# =============================================================================

T20 = np.zeros(NumReceptor)
T30 = np.zeros(NumReceptor)
EDT = np.zeros(NumReceptor)
C50 = np.zeros(NumReceptor)
D50 = np.zeros(NumReceptor)

for ir in range(NumReceptor):

    T20[ir], T30[ir], EDT[ir], C50[ir], D50[ir] = func.calc_parametros_objetivos(TempoEnergia[ir], EnergiaDec[ir])

    print('\n')
    print('Parâmetro objetivos - Ray Tracing')
    print('T20 igual a %.2f s' %T20[ir])
    print('T30 igual a %.2f s' %T30[ir])
    print('EDT igual a %.2f s' %EDT[ir])
    print('C50 igual a %.2f dB' %C50[ir])
    print('D50 igual a %.2f %%' %(D50[ir]*100))

dados = np.array([T20, T30, EDT, C50, D50])

# =============================================================================
# Plot Energia
# =============================================================================

for ir in range(NumReceptor):
    
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.set_title('Acústica de Raios - Sala de Aula')
    ax.set_xlabel('Tempo [s]')
    ax.set_ylabel('Energia Total [W]')
    ax.set_xlim([0, T60Sabine*1.1])
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax.plot(TempoEnergia[ir], EnergiaTotal[ir],'-*')
    ax.grid()
    plt.show()

    ax = fig.add_subplot(212)
    ax.set_title('Acústica de Raios - Sala de Aula')
    ax.set_xlabel('Tempo [s]')
    ax.set_ylabel('Energia Decaimento [dB]')
    ax.set_xlim([0, T60Sabine*1.1])
    ax.plot(TempoEnergia[ir], 10*np.log10(EnergiaDec[ir]))
    ax.grid()
    plt.show()
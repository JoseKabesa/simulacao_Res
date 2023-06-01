from EDH_explicita import EDH_explicita
from EDH_analitica import EDH_analitica
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math as mt
import time
'''DADOS DE ENTRADA'''
k = 20  # permeabilidade md
mi = 1.2  # viscosidade cp
phi = 0.3  # porosidade
ct = 150e-6  # compressibilidade total (kgf/cm2)-1
p0 = 3000  # psi
N = 10
L = 200  # m
anos = 20  # anos
tf = 365 * anos
dt = 0.5


'''CALCULO DA SOLUÇÃO NÚMERICA'''
A, qw, P1 = EDH_explicita(k, mi, phi, ct, p0, N, L, tf, dt)

'''CALCULO DA SOLUÇÃO ANALÍTICA'''
P2 = EDH_analitica(k, mi, phi, ct, p0, N, L, tf, dt, A, -qw)

plt.xlabel("Distância (m)")
plt.ylabel("Pressão (psi)")
plt.title("Solução Explícita e Analítica")
plt.show()

erro = abs(np.delete(P2, tf, 1) - np.delete(P1, tf, 1)) / np.delete(P2, tf, 1) * 100

plt.plot(erro)
plt.xlabel("Blocos")
plt.ylabel("Erro em %")
plt.title("Erro entre as soluções")
plt.show()

print('Modificação')
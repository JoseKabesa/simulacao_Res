def EDH_analitica(k,mi,phi,ct,p0,N,L,tf,dt,A,qw):

    import numpy as np
    import math as mt
    import matplotlib.pyplot as plt
    """
        Função que cálcula EDH de 1 Ordem para fluxo linear usando solução analítica
        ---------------------------------------------------------
        Parâmetros de entrada:
        k = permeabilidade em mD
        mi = viscosidade em cp
        phi = porosidade
        ct =  compressibilidade total em (kgf/cm2)^-1
        p0 = pressão inicial do sistema em psi
        N = número de blocos da malha
        L = comprimento total em m
        tf = tempo final de simulação
        dt = passo de tempo
        OBS: O tempo está para ser calculado em dias, mas uma conversão pode ser feita
        Por exemplo: 1h = 1/24
        A =  área em m2
        qw = vazão de produção poço em bbl/d
        ---------------------------------------------------------
        """
    c1 = 270.4325  # constante de conversão de unidade para o termo de fluxo da condição de contorno
    c2 = 9.6777778e-8  # constante de conversão para o eta
    dx = L / N
    tempo = np.arange(1, tf + dt, dt, dtype=int)
    distancia = np.arange(dx / 2, L, dx, dtype=int)
    eta = c2 * k / (mi * phi * ct)
    p = np.zeros((N, int(tf / dt) + 1))
    p[:, 0] = p0

    for i in tempo:
        for j in range(N):
            x = distancia[j]
            p[j,i] = p0 - c1 * qw * mi * L/(k * A) * (mt.sqrt(4 * eta * i / (mt.pi * L ** 2))*mt.exp(-x ** 2/(4 * eta *i))-(x / L)*mt.erfc(x / mt.sqrt((4 * eta* i))))
        if i % 365 == 0:
            plt.plot(distancia,p[:,i],'g--')

    return(p)
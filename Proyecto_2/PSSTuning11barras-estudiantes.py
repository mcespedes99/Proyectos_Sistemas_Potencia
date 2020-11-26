#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
import cmath
from numpy import linalg as LA
from scipy import signal
import control
control.use_matlab_defaults()
import matplotlib.pyplot as plt


# Cargar matrices A, B, C y D
files = ['Proyecto_2/A.txt', 'Proyecto_2/B.txt', 'Proyecto_2/C.txt']
matrices = []
for file in files:
    matrix = []
    with open(file, 'r') as f:
        for row in f:
            data = [float(x) for x in row.split(',')]
            matrix.append(data)
    f.close()
    matrix = np.array(matrix)
    matrices.append(matrix)

A = matrices[0]
B = matrices[1]
C = matrices[2]
D = np.array(0)

# Obtener valores propios de A
lamb, v = LA.eig(A)

# Definir cantidad de decimales
places = 6

# Mostrar valores propios inestables
for eigval in lamb:
	if np.real(eigval) > 0 and np.imag(eigval) != 0:
		print('Valor propio inestable de A: ' + str(np.around(eigval, places)))

# Calcular función de transferencia G(s) a partir de A, B, C y D
num, den = signal.ss2tf(A, B, C, D)

# Obtener residuos r y polos p de la función de transferencia
num = num.reshape((np.size(num),))
den = den.reshape((np.size(den),))
r, p, k = signal.residue(num, den)

# Determinar índices de polos inestables
ind_unstable = []
for (i, pole) in enumerate(p):
	if np.real(pole) > 0 and np.imag(pole) != 0:
		print('Polo inestable de G(s): ' + str(np.around(pole, places)))
		ind_unstable.append(i)

# Leer residuo asociado a polo inestable con parte real positiva
pos = 0
for i in ind_unstable:
	if np.imag(p[i]) > 0:
		pos = i
ri = r[pos]

# Desplegar polos y residuos
msg = 'Polo con parte imaginaria positiva: ' + str(np.around(p[pos], places))
print(msg)
# Desplegar residuo
msg = 'Residuo asociado a polo con parte imaginaria positiva: ' \
      + str(np.around(ri, places))
print(msg)

# cálculo de ángulo del residuo de interes

thetaRi=cmath.phase(ri)

# Diseño del PSS
#Para Caso C:
if thetaRi < 0:
    thetaRi += 2*np.pi
# Caso A y B
if 0 <= thetaRi and thetaRi <= np.pi:
    thetaPSSpos = np.pi - thetaRi
else:
    thetaPSSpos = thetaRi - np.pi

# Calculo de numero de bloques de compensacion requeridos
Np=math.ceil(abs(thetaPSSpos)/(math.pi/3))

# calculo de valor a en bloques de compensacion, para realm. positiva
ap=(1-math.sin(thetaPSSpos/Np))/(1+math.sin(thetaPSSpos/Np))

# frecuencia de oscilacion a amortiguar
wi=lmbda[pos].imag

# valor de T del bloque con realim. positivo
Tp=1/(wi*math.sqrt(ap))

# calculo de constantes de tiempo en bloques de compensación realm. positiva
T1p=Tp
T2p=ap*Tp

Tw=1.5 #Revisar si debe ser 1.5 o 10

# creacion de FT de los bloques para realm pos.
bloquepos=control.tf([T1p, 1],[T2p, 1])


# creacion de bloque Washout (se usa igual a modelo generic1 de RAMSES)
WS=control.tf([Tw, 0],[Tw, 1])

# se grafica un ejemplo de diagrama de bode para bloque de comp. en realm. positiva
bodeplot = control.bode_plot(bloquepos)
plt.savefig('bode.pdf')

# funcion de transferencia de lazo cerrado (realm. positiva)
TFpos=sys/(1-bloquepos**Np*sys*WS)

# Se grafica el Lugar geometrico de las raices para el sistema de lazo abierto
rlocus=control.root_locus(-bloquepos**Np*sys*WS,xlim=(-40,40),ylim=(-40,40))
plt.savefig('root_locus.pdf')

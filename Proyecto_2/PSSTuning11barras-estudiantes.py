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

sys = control.tf(num,den)

r, p, k = signal.residue(num, den)
print("\n")
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
msg = '\nPolo con parte imaginaria positiva: ' + str(np.around(p[pos], places))
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
wi=p[pos].imag
print("\nFrecuencia de oscilación a amortiguar: "+str(np.around(wi/(2*math.pi), places))+" Hz")

# valor de T del bloque con realim. positivo
Tp=1/(wi*math.sqrt(ap))

# calculo de constantes de tiempo en bloques de compensación realm. positiva
print("\nConstantes de tiempo de bloques de compensación:")
T1p=Tp
print("T_1= "+str(np.around(T1p, places))+" s")
T2p=ap*Tp
print("T_2= "+str(np.around(T1p, places))+" s")

Tw=10 #Revisar si debe ser 1.5 o 10
print("\nConstante de tiempo de filtro Washout:")
print("Tw= "+str(np.around(Tw, places))+" s")

# creacion de FT de los bloques para realm pos.
bloquepos=control.tf([T1p, 1],[T2p, 1])

# creacion de bloque Washout (se usa igual a modelo generic1 de RAMSES)
WS=control.tf([1, 0],[Tw, 1])

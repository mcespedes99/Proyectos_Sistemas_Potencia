#Daniel: es necesario que las admitancias se ingresen en la forma G+j*B para la definición de las ecuaciones de potencia. Puede ver el siguiente link: https://docs.julialang.org/en/v1/manual/complex-and-rational-numbers/
using LinearAlgebra
println("Inserte el número de barras")
barras = parse(Int8, readline()) #cantidad de barras del problema
println("Inserte el número de barras oscilantes")
oscilante = parse(Int8, readline()) #cantidad de barras oscilantes del problema
println("Inserte el número de barras generadoras")
generadoras = parse(Int8, readline()) #cantidad de barras generadoras del problema
println("Inserte el número de barras de carga")
cargas = parse(Int8, readline()) #cantidad de barras cargas del problema
admitancia = zeros(Int8,barras,barras) #crea una matriz de dimensión de la cantidad de barras y la llena de ceros.

for i in 1:barras
  for n in 1:barras
    if n == i
      println("Inserte la admitancia entre $i y tierra")
      admitancia[i,i] = admitancia[i,i] + parse(Int8, readline())
    else
      println("Inserte la admitancia entre $i y $n")
      admitancia[i,n] = -parse(Int8, readline())
      admitancia[i,i] = admitancia[i,i] -admitancia[i,n]
    end
  end 
end# Este código arma la matriz de admitancias.

#Cada barra tiene una P, Q, V y cita asociada. Solo se conocen 2 de estas

#Definición de potencias:
#Barra 1: tipo Swing (N.A.)
#V1 = 1.05
#cita_1 = 0 #Referencia

#Barra 2: Tipo Gen
#V2 = 1.05 #p.u.
P2 = 0.5 #MW
#cita_2 = 0; #Valor inicial asumido

#Barra 3: tipo Gen
#V3 = 1.05 #p.u.
P3 = 600 #MW (Aquí la potencia es positiva porque entra al sistema)
#cita_3 = 0; #Valor inicial asumido

#Barra 4: tipo Load
P4 = -700 #MW (Aquí la potencia es negativa porque está siendo consumida por la barra)
Q4 = -400 #MW 
#V4 = 1 #Valor inicial asumido
#cita_4 = 0 #Valor inicial asumido


#Barra 5: tipo Load
P5 = -700 #MW
Q5 = -400 #MW 
#V5 = 1 #Valor inicial asumido
#cita_5 = 0 #Valor inicial asumido

#Barra 6: tipo Load
P6 = -700 #MW
Q6 = -400 #MW 
#V6 = 1 #Valor inicial asumido
#cita_6 = 0 #Valor inicial asumido

#Definición de los valores de las tensiones y ángulos de las barras en una lista que va desde la barra 1 hasta la barra 6 
#Los valores desconocidos de V se asumen en 1:
V = [1.05 1.05 1.07 1 1 1]

#Lista ecuaciones de potencias activas de las barras
ecuaciones_potencia = zeros(Int8,barras,barras)

for i in 1:barras
  for n in 1:barras
    if n == i
      Vi
    else
      println("Inserte la admitancia entre $i y $n")
      admitancia[i,n] = -parse(Int8, readline())
      admitancia[i,i] = admitancia[i,i] -admitancia[i,n]
    end
  end 
end
# real(1 + 2im) 
# imag(1 + 2im)
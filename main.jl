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

# real(1 + 2im) 
# imag(1 + 2im)

# 1. Definición del vector X: (vectores de tensiones y ángulos de las barras)
# X = [V1 V2 V3 V4 V5 V6 cita_1 cita_2 cita_3 cita_4 cita_5 cita_6]
# Aquellas V desconocidas se asumen como 1 en la primera iteración y las cita 0.
X = [1.05 1.05 1.07 1 1 1; 0 0 0 0 0 0] 

# 2. Cálculo del primer missmatch: se utilizan las ecuaciones implícitas

vector_missmatch = zeros(8) #Definicón del vector de missmatch

for j = 2:6
  if j < 4
    for i in 1:barras 
      vector_missmatch[j-1] += X[i,1]*X[j,1]*(real(matriz[i,j])*cos(X[i,2]-X[i,2])+imag(matriz[i,j])*sin(X[i,2]-X[i,2]))
    end
  elseif j == 4
    vector_missmatch[j-1] += X[i,1]*X[j,1]*(real(matriz[i,j])*cos(X[i,2]-X[i,2])+imag(matriz[i,j])*sin(X[i,2]-X[i,2]))
    vector_missmatch[j+1] += X[i,1]*X[j,1]*(real(matriz[i,j])*cos(X[i,2]-X[i,2])+imag(matriz[i,j])*sin(X[i,2]-X[i,2]))
  elseif j == 5
    vector_missmatch[j] += X[i,1]*X[j,1]*(real(matriz[i,j])*cos(X[i,2]-X[i,2])+imag(matriz[i,j])*sin(X[i,2]-X[i,2]))
    vector_missmatch[j+1] += X[i,1]*X[j,1]*(real(matriz[i,j])*cos(X[i,2]-X[i,2])+imag(matriz[i,j])*sin(X[i,2]-X[i,2]))
  else
    vector_missmatch[j+1] += X[i,1]*X[j,1]*(real(matriz[i,j])*cos(X[i,2]-X[i,2])+imag(matriz[i,j])*sin(X[i,2]-X[i,2]))
    vector_missmatch[j+2] += X[i,1]*X[j,1]*(real(matriz[i,j])*cos(X[i,2]-X[i,2])+imag(matriz[i,j])*sin(X[i,2]-X[i,2]))
  end
end

# 3. Cálculo del Jacobiano para la iteración correspondiente:
J = zeros(Int8,8,8)
# Primero se agrega la submatriz de dP/d_cita que es 5x5
for i = 2:6
  for j = 2:6
    if i==j
      for k = 1:6
        J[i-1,j-1] += -X[1,i]*(X[1,k]*(real(matriz[i,k])*sin(X[2,i]-X[2,k])-imag(matriz[i,k])*cos(X[2,i]-X[2,k]))-X[1,i]^2 * imag(matriz[i,i]))
      end
    else
      J[i-1,j-1] = X[1,i]*X[1,j]*(real(matriz[i,j])*sin(X[2,i]-X[2,j]) - imag(matriz[i,j])*cos(X[2,i]-X[2,j]))
    end
  end
end

# A continuación se añade la submatriz dP/dV que es de 5x3
for i = 2:6
  for j = 4:6
    if i==j
      for k = 1:6
        J[i-1,j-1+5] += X[1,k]*(real(matriz[i,k])*cos(X[2,i]-X[2,k])+imag(matriz[i,k])*sin(X[2,i]-X[2,k])) + X[1,i]*real(matriz[i,i])
      end
    else
      J[i-1,j-1+5] = X[1,i]*(real(matriz[i,j])*cos(X[2,i]-X[2,j]) + imag(matriz[i,j])*sin(X[2,i]-X[2,j]))
    end
  end
end

# Como tercer punto, de añade la submatriz dQ/d_cita que es 3x5
for i = 4:6
  for j = 2:6
    if i==j
      for k = 1:6
        J[i-1+5,j-1] += X[1,i]*(X[1,k]*(real(matriz[i,k])*cos(X[2,i]-X[2,k])+imag(matriz[i,k])*sin(X[2,i]-X[2,k])) - X[1,i]^2 *real(matriz[i,i]))
      end
    else
      J[i-1+5,j-1] = -X[1,i]*X[1,j]*(real(matriz[i,j])*cos(X[2,i]-X[2,j]) + imag(matriz[i,j])*sin(X[2,i]-X[2,j]))
    end
  end
end

#Por último, se añade la submatriz dQ/dV que es 3x3
for i = 4:6
  for j = 4:6
    if i==j
      for k = 1:6
        J[i-1+5,j-1+5] += X[1,k]*(real(matriz[i,k])*sin(X[2,i]-X[2,k])-imag(matriz[i,k])*cos(X[2,i]-X[2,k])) - X[1,i]*imag(matriz[i,i])
      end
    else
      J[i-1+5,j-1+5] = X[1,i]*(real(matriz[i,j])*sin(X[2,i]-X[2,j]) - imag(matriz[i,j])*cos(X[2,i]-X[2,j]))
    end
  end
end
#Fin del Jacobiano
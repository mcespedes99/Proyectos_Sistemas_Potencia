using LinearAlgebra
barras = 6 #se define la cantidad de barras
r = zeros(barras,barras) #es una matriz con las reistencias, donde r[x,y] es la resitencia entre la barra x y la y.
x = zeros(barras,barras) #es una matriz con las impedancias, donde x[x,y] es la impedancia entre la barra x y la y.
bcap = zeros(barras,barras) #es una matriz con los bcap, donde bcap[x,y] es el bcap entre la barra x y la y.
r[1,2]=0.1
x[1,2]=0.2
bcap[1,2]=0.02
r[1,4]=0.05
x[1,4]=0.2
bcap[1,4]=0.02
r[1,5]=0.08
x[1,5]=0.3
bcap[1,5]=0.03
r[2,3]=0.05
x[2,3]=0.25
bcap[2,3]=0.03
r[2,4]=0.05
x[2,4]=0.1
bcap[2,4]=0.01
r[2,5]=0.1
x[2,5]=0.3
bcap[2,5]=0.02
r[2,6]=0.07
x[2,6]=0.2
bcap[2,6]=0.025
r[3,5]=0.12
x[3,5]=0.26
bcap[3,5]=0.025
r[3,6]=0.02
x[3,6]=0.1
bcap[3,6]=0.01
r[4,5]=0.2
x[4,5]=0.4
bcap[4,5]=0.04
r[5,6]=0.1
x[5,6]=0.3
bcap[5,6]=0.03
# debido a que estas tres matrices son simétricas se añade esta simetría, sin la necesidad de otras 32 líneas de código.
for i in 1:barras
  for j in 1:barras
    if i > j
      r[i,j]= r[j,i]
      x[i,j]= x[j,i]
      bcap[i,j]= bcap[j,i]
    end
  end
end
# ahora se forma la matriz de admitancias utilizando las fórmulas conocidas, donde Y=r/(r^2+x^2)+jx/(r^2+x^2)
#Y utilizando el modelo pi entre las líneas, por lo que bcap es una admitancia que va de la línea a tierra.
matriz = zeros(Complex{Float64},barras,barras)
for i in 1:barras
  for j in 1:barras
    if j == i
      matriz[i,i] = matriz[i,i] -(bcap[i,i])im
    else
      if r[i,j] != 0
        matriz[i,j]= -r[i,j]/(r[i,j]^2+x[i,j]^2) + (x[i,j]/(r[i,j]^2+x[i,j]^2))im
        matriz[i,i]=matriz[i,i]-(bcap[i,j])im -matriz[i,j]
      end    
    end
  end
end
# 1. Definición del vector X: (vectores de tensiones y ángulos de las barras)
# X = [V1 V2 V3 V4 V5 V6 cita_1 cita_2 cita_3 cita_4 cita_5 cita_6]
# Aquellas V desconocidas se asumen como 1 en la primera iteración y las cita 0.
X = [1.05 1.05 1.07 1 1 1; 0 0 0 0 0 0] 

#Potencias de las barras. Las que no se conocen, se asignan como cero; de manera que la dimensión del vector sea la correcta:
Pbarras= [0 500 600 -700 -700 -700]
Qbarras = [0 0 0 -700 -700 -700]

# 2. Cálculo del primer missmatch: se utilizan las ecuaciones implícitas

vector_missmatch = zeros(8) #Definicón del vector de missmatch

for j = 2:6
  if j < 4
    for i in 1:barras 
      vector_missmatch[j-1] += X[1,i]*X[1,j]*(real(matriz[j,i])*cos(X[2,j]-X[2,i])+imag(matriz[j,1])*sin(X[2,j]-X[2,i])) - Pbarras[j]
    end
  elseif j == 4
    for i in 1:barras 
      vector_missmatch[j-1] += X[1,i]*X[1,j]*(real(matriz[j,i])*cos(X[2,j]-X[2,i])+imag(matriz[j,1])*sin(X[2,j]-X[2,i])) - Pbarras[j]
    end

    for i in 1:barras 
      vector_missmatch[j] += X[1,i]*X[1,j]*(real(matriz[j,i])*sin(X[2,j]-X[2,i])+imag(matriz[j,1])*cos(X[2,j]-X[2,i])) - Qbarras[j]
    end
  elseif j == 5
    for i in 1:barras 
      vector_missmatch[j] += X[1,i]*X[1,j]*(real(matriz[j,i])*cos(X[2,j]-X[2,i])+imag(matriz[j,1])*sin(X[2,j]-X[2,i])) - Pbarras[j]
    end
    for i in 1:barras 
      vector_missmatch[j+1] += X[1,i]*X[1,j]*(real(matriz[j,i])*sin(X[2,j]-X[2,i])+imag(matriz[j,1])*cos(X[2,j]-X[2,i])) - Qbarras[j]
    end
  else
    for i in 1:barras 
      vector_missmatch[j+1] += X[1,i]*X[1,j]*(real(matriz[j,i])*cos(X[2,j]-X[2,i])+imag(matriz[j,1])*sin(X[2,j]-X[2,i])) - Pbarras[j]
    end
    for i in 1:barras 
      vector_missmatch[j+2] += X[1,i]*X[1,j]*(real(matriz[j,i])*sin(X[2,j]-X[2,i])+imag(matriz[j,1])*cos(X[2,j]-X[2,i])) - Qbarras[j]
    end
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

